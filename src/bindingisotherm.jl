
####################################################
"""
$(TYPEDSIGNATURES)

Least squares linear regression.
"""
function lssqregr(points)
    
    ## G = [protc ones(length(protc))]
    ## dobs = sdsc
    G = [points[:,1] ones(size(points,1))]
    dobs = points[:,2]

    ## overdetermined least squares
    ## mpost = (G'*G)*G'*dobs
    ## (G'*G) * h = G'
    lstmp = (G'*G) \ G'
    mpost = lstmp * dobs

    angcoe = mpost[1]
    intercept = mpost[2]

    return angcoe,intercept
end

####################################################

function calcfreeSDSNbound(lstlines::Array{<:Array{<:Real,2},1}) #Vector{Array{<:Real,2}})
    nlines = length(lstlines)
    intercept = Vector{Float64}(undef,nlines)
    angcoe = Vector{Float64}(undef,nlines)
    for i=1:nlines
        angcoe[i],intercept[i] = lssqregr(lstlines[i])
    end
    
    ##
    ## totSDS = freeSDS + Nbound * protc
    ##
    # y = angcoe.*x .+ intercept
    # x = (-intercept .+ 0)./angcoe
    Nbound = 1.0./angcoe
    freeSDS = -intercept./angcoe

    ## Sort in ascending order of SDS concentration
    idxsort = sortperm(freeSDS)
    freeSDS = freeSDS[idxsort] 
    Nbound  = Nbound[idxsort] 

    return freeSDS,Nbound
end

####################################################

"""
$(TYPEDSIGNATURES)

Define a set of features on the Beta mix to subsequently compute the binding isotherm.
It uses the stationary and inflection points at given protein concentrations (protcon). 
"""
function findcurvefeatures(betamix::BetaMix2D,protcon::Vector{<:Real})
    
    function fwd1ptxy(xpt::Real,ypt::Real)
        xy = hcat(xpt,ypt)
        out = ANISPROU.forwmod2D(betamix.betpar,xy,betamix.modkonamp)
        return out[1] 
    end

    xpt,ypt = 0.0,0.0
    deriv1f(xpt::Real) = ForwardDiff.derivative(xpt->fwd1ptxy(xpt,ypt),xpt)
    deriv2f(xpt::Real) = ForwardDiff.derivative(xpt->deriv1f(xpt),xpt)

    ##==================================================
    ## get stationary and inflection points using
    ##   first and second derivatives    
    lowlim = betamix.betpar.a 
    uplim  = betamix.betpar.b 
    ny = length(protcon)
    statpts = Array{Vector{Float64}}(undef,ny)
    inflpts = Array{Vector{Float64}}(undef,ny)
    for i=1:ny
        ypt = protcon[i] ## define ypt for closure derivf(xpt)
        statpts[i] = find_zeros(deriv1f,lowlim,uplim)
        inflpts[i] = find_zeros(deriv2f,lowlim,uplim)
    end

    return statpts,inflpts
end

##################################################################
##################################################################

"""
$(TYPEDSIGNATURES)

Define a set of features on a single individual Beta function to subsequently compute the binding isotherm.
It uses a ratio of amplitude compared to maximum amplitude (ampratios).
"""
function freeboundsds_singlebeta(betpar::ScaledBeta2DParams,mcur::Vector{<:Real} ;
                                 ampratios::Union{Vector{<:Real},Nothing}=nothing,
                                 nptsprot::Integer=10)

    # init
    sdsfreeNbound = Array{Float64,2}(undef,0,2) #npoints,2)

    ## selectd values of amplitudes to trace points
    if ampratios==nothing
        ampratios = collect(LinRange(0.1,1.0,6))
    end
    
    ## sample points for protein concentration (vector)
    ypos = collect(LinRange(betpar.ymin,betpar.ymax,nptsprot))

    namppoints = length(ampratios)
    for p=1:namppoints
        ## find the straight line(s) for a given amplitude 
        fsds,nbou = linefreesdsNbound(betpar,mcur,ypos,ampratios[p])
        # store results
        sdsfreeNbound = vcat(sdsfreeNbound,[fsds nbou])
    end

    ## sort data
    sdsfreeNbound = sortslices(sdsfreeNbound,dims=1)

    return sdsfreeNbound #,sdsfNb_mode
end

###########################################################

function linefreesdsNbound(betpar::ScaledBeta2DParams,mcur::Vector{<:Real},
                           ypos::Vector{<:Real},ampratio::Real)

    ## find the point where ampl/max(ampl)==selectedampl
    nypts = length(ypos)
    ## swapping for least squares, y first and x last, [protein] and [SDS]...
    sdsc = Dict()
    protc = Dict()
    ## left, right with respect to the mode!!
    lssides = ["left","right"]
 
    ##------------------------------------------------
    ## find the x,y values for given amplitudes
    ## loop on all y values  
    for i=1:nypts
        ## the following returns 1 or 2-element array
        xpos,ytmp,sidevec = findxbeta(betpar,mcur,ampratio,ypos[i])
        ## y first and x last, [protein] and [SDS]...
        for s=1:length(xpos)
            try
                push!(sdsc[sidevec[s]],xpos[s])
                push!(protc[sidevec[s]],ytmp[s])
            catch
                sdsc[sidevec[s]]  = Float64[]
                protc[sidevec[s]] = Float64[]
                push!(sdsc[sidevec[s]],xpos[s])
                push!(protc[sidevec[s]],ytmp[s])
            end
        end       
    end

    nkeys = length(keys(sdsc))

    ## find the straight line by least squares
    intercept = Vector{Float64}(undef,nkeys)
    angcoe = Vector{Float64}(undef,nkeys)
    for (i,sd) in enumerate(keys(sdsc))
        ## least squares on swapped x and y, i.e., [SDS] and [protein]
        G = [protc[sd] ones(length(protc[sd]))]
        dobs = sdsc[sd]
        ## overdetermined least squares
        ## mpost = (G'*G)*G'*dobs
        ## (G'*G) * h = G'
        lstmp = (G'*G) \ G'
        mpost = lstmp * dobs

        angcoe[i] = mpost[1]
        intercept[i] = mpost[2]
    end    

    return intercept,angcoe
end


###########################################################

function findxbeta(betpar::ScaledBeta2DParams,mcur::Vector{<:Real},
                   ampratio::Real,yval::Real)

    ## yval -> protein concentration

    ##------------------------------------------------------
    function misfroot(x::Real)
        bval = scaledbeta(mode,kon,betpar.a,betpar.b,amp,x)
        ## abs(f(x)) because we have also negative amplitudes,
        ##  however, the beta is always either all positive or
        ##  all negative
        msr = bval-ampval 
        return msr 
    end
    ##------------------------------------------------------

    @assert 0.0<ampratio<=1.0

    ## mode,kon,amp
    ## get the values of model parameters at y=ycur
    mode,kon,amp = getmodparbeta(betpar,mcur,yval)

    ## left, right with respect to the mode!!
    maxamp = scaledbeta(mode,kon,betpar.a,betpar.b,amp,mode)
    ## actual amplitude value
    ampval = ampratio*maxamp

    lssides = ["left","right"]
    if ampval==maxamp
        xpos,ypos,sides = [mode],[yval],["left"]
    else
        xpos = Vector{Float64}(undef,2)
        ypos = Vector{Float64}(undef,2)
        sides = Vector{String}(undef,2)
        for (i,sd) in enumerate(lssides)
            sd=="left" ? (abracket=betpar.a) : (abracket=mode)
            sd=="left" ? (bbracket=mode) : (bbracket=betpar.b)
            ## https://github.com/JuliaMath/Roots.jl
            xpoint = fzero(misfroot,abracket,bbracket)
 #@show i,sd,abracket,bbracket,xpoint,sd
            xpos[i],ypos[i],sides[i] = xpoint,yval,sd
        end
    end
#@show xpos,ypos,sides
    return xpos,ypos,sides
end
 
###########################################################
