

##################################################################

function freeboundsds_betamix(betmix::BetaMix2D ;
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

##################################################################

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