

###############################################
"""
$(TYPEDEF)

Structure containing the parameters of the 2D Beta functions.

# Fields 

$(TYPEDFIELDS)

"""
struct ScaledBeta2DParams
    "number of model parameters"
    nummodpar::Integer
    "type of function on y (protein concentration) for the mode"
    modefuny::String 
    "type of function on y (protein concentration) for the confidence parameter"
    konfuny::String 
    "type of function on y (protein concentration) for the amplitude"
    ampfuny::String 
    "lower bound for Beta function (along x, i.e. SDS concentration)"
    a::Real
    "upper bound for Beta function (along y, i.e. protein concentration)"
    b::Real
    "user defined minimum y (protein concentration)"
    ymin::Real
    "user defined maximum y (protein concentration)"
    ymax::Real
    "starting index in the vector of model parameters for mode"
    idxmode::Integer
    "starting index in the vector of model parameters for confidence parameter"
    idxkon::Integer
    "starting index in the vector of model parameters for amplitude parameter"
    idxamp::Integer

    ##----------------------------------------------------------------------
    function ScaledBeta2DParams(; modefuny,konfuny,ampfuny,a,b,ymin,ymax)
        nummodpar::Integer=0
        if modefuny=="linear"
            nummodpar+=2
            idxmode=1
        elseif modefuny=="quadratic"
            nummodpar+=3
            idxmode=1
        end
        if konfuny=="linear"
            idxkon=nummodpar+1
            nummodpar+=2
        elseif konfuny=="quadratic"
            idxkon=nummodpar+1
            nummodpar+=3
        end
        if ampfuny=="linear"
            idxamp=nummodpar+1
            nummodpar+=2
        elseif ampfuny=="quadratic"
            idxamp=nummodpar+1
            nummodpar+=3
        end
        return new(nummodpar,modefuny,konfuny,ampfuny,a,b,ymin,ymax,idxmode,idxkon,idxamp)
    end
    ##----------------------------------------------------------------------
end

###########################################################

"""
$(TYPEDEF)

Structure containing the observed (measured) data, i.e. the enthalpy, the protein and SDS concentration.

# Fields 

$(TYPEDFIELDS)

"""
Base.@kwdef struct ITCObsData
    "Enthalpy values"
    enthalpy::Vector{<:Real}
    "Indices pointing to single experiments data"
    idxdata::Vector{UnitRange{Int64}}
    "Concentration of SDS (first column) and protein (second column)"
    sdsprotcon::Array{<:Real,2}
    "Name of the protein"
    protein::String
end

###########################################################

"""
$(TYPEDEF)

Structure containing the observed (measured) data, i.e. the enthalpy, the protein and SDS concentration.

# Fields 

$(TYPEDFIELDS)

"""
struct BetaMix2D
    "structure containing the parameters of the 2D Beta functions"
    betpar::ScaledBeta2DParams
    "2D array where each column is a set of mode, confidence and amplitude parameters"
    modkonamp::Matrix{<:Real}
    "name of the protein"
    protein::String
end
  
##############################################################

"""
$(TYPEDSIGNATURES)

Set the constraints for the Newton optimization.
"""
function setconstraints(betpar::ScaledBeta2DParams,mstart::Matrix{<:Real},
                        lowc::Vector{<:Real},upc::Vector{<:Real})

    lenm = size(mstart,1)
    ncomp = size(mstart,2) #div(lenm,betpar.nummodpar)

    lowconstr = zeros(betpar.nummodpar,ncomp)
    upconstr = zeros(betpar.nummodpar,ncomp)

    #@show size(lowc),size(upc)
    for c=1:ncomp
        #c1 = (c-1)*betpar.nummodpar
        for i=1:betpar.nummodpar
            lowconstr[i,c] = lowc[i]
            upconstr[i,c] = upc[i]
        end
    end
    
    # @show mstart.<=lowconstr
    # @show mstart.>=upconstr

    ## check mstart...
    @assert all(mstart.>=lowconstr)
    @assert all(mstart.<=upconstr)

    return vec(lowconstr),vec(upconstr)
end

######################################################################

"""
$(TYPEDSIGNATURES)

Calculate misfit between observed/measured and calculated data.
"""
function misfitbeta2D(betpar::ScaledBeta2DParams,dobs::ITCObsData,
                      invCd::Matrix{<:Real},mcur::Matrix{<:Real})
    
    # compute synthetic data
    # xy = dobs.sdsprotcon
    dcalc = forwmod2D(betpar,dobs.sdsprotcon,mcur)

    # calculate likelihood value
    difcalobs = dcalc.-dobs.enthalpy
    tmp1 = invCd * difcalobs
    misflkl = 0.5 * dot(difcalobs,tmp1)

    # total misfit
    misf = misflkl 

    return misf
end

###########################################################

"""
$(TYPEDSIGNATURES)

Compute the forward response for given input parameters, i.e. the enthalpy 
 values (2D) for given Beta functions (2D) parameters.
"""
function forwmod2D(betpar::ScaledBeta2DParams,xy::Array{<:Real,2},mcur::Matrix{<:Real})
                   
    npts = size(xy,1)
    nummodpar = size(mcur,1)
    ncomp = size(mcur,2)

    sscbeta = zeros(Real,npts) #zeros(eltype(mcur),npts)
    # sscbeta .= zero(eltype(mcur))   #0.0 ## make sure it's zeroed

    for n=1:ncomp
        # sum all the components
        sscbeta .+= singlescaledbeta2D(betpar,xy,mcur[:,n]) 
    end

    return sscbeta
end

###########################################################

"""
 $(TYPEDSIGNATURES)

  2D-fied scaled Beta function for given x (SDS concentration) and y (protein concentration) values.
"""
function singlescaledbeta2D(betpar::ScaledBeta2DParams,xy::Array{<:Real,2},mcur::Vector{<:Real})

    # pre-allocate array of arrays (one for each y value)
    npts = size(xy,1)
    scbeta2d = zeros(Real,npts) #zeros(eltype(mcur),npts)
    for l=1:npts
        ## xy[l,1] -> SDS concentration
        ## xy[l,2] -> protein concentration
        ycur = xy[l,2] # protein concentration
        xcur = xy[l,1]
        
        ## get the values of model parameters at y=ycur
        mode,kon,amp = getmodparbeta(betpar,mcur,ycur)

        # beta pdf along x
        ## xy[l,1] -> SDS concentration
        scbeta2d[l] = scaledbeta(mode,kon,betpar.a,betpar.b,amp,xcur)
    end

    return scbeta2d
end

###########################################################

"""
 $(TYPEDSIGNATURES)

Get the parameters (mode, confidence, amplitude) of the Beta functions for given y (protein concentration).
"""
function getmodparbeta(betpar::ScaledBeta2DParams,mcur::Vector{<:Real},ycur::Real)

    ##-------------------------------------------
    ## Linear: Straight line through 2 points:
    ##
    ## y = (y2-y1)/(x2-x1) * (x-x1) + y1
    ##

    ##-------------------------------------------
    ## Quadratic: to be implemented...
    ##
    ## 
    ##

    ## limits of protein concentration
    x1 = betpar.ymin
    x2 = betpar.ymax
    @assert (x2-x1)!=0.0

    ## mode
    if betpar.modefuny=="linear"
        mode1,mode2 = mcur[betpar.idxmode:betpar.idxmode+1]
        mode = ((mode2-mode1)/(x2-x1))*(ycur-x1) + mode1
    elseif betpar.modefuny=="quadratic"
        error("Quadratic representation of mode along y not yet implemented.")
    end
    
    ## confidence parameter
    if betpar.konfuny=="linear"
        kon1,kon2 = mcur[betpar.idxkon:betpar.idxkon+1]
        kon  = ((kon2-kon1)/(x2-x1))*(ycur-x1) + kon1
    elseif betpar.konfuny=="quadratic"
        error("Quadratic representation of confidence parameter along y not yet implemented.")
    end
    
    ## amplitude
    if betpar.ampfuny=="linear"
        amp1,amp2 = mcur[betpar.idxamp:betpar.idxamp+1]
        amp = ((amp2-amp1)/(x2-x1))*(ycur-x1) + amp1
    elseif betpar.ampfuny=="quadratic"
        error("Quadratic representation of amplitude parameter along y not yet implemented.")
    end

    return mode,kon,amp
end

###########################################################

"""
 $(TYPEDSIGNATURES)

1D modified scaled Beta function.
"""
@inline function scaledbeta(mo::Real,kon::Real,a::Real,b::Real,
                            amplscale::Real,x::Real)

    if a > mo || mo > b
        @show mo,kon,amplscale,a,b
    end
    
    @assert b>a
    @assert a <= mo <= b
    ## keep kon at >2.0 instead of >=2.0 to avoid weird shapes...
    @assert kon>2.0 

    rmo = (mo-a)/(b-a)
    α = rmo * (kon-2.0) + 1.0
    β = (1.0-rmo) * (kon-2.0) + 1.0
    @assert α>0
    @assert β>0

    betatyp = Beta(α,β)

    xresc = (x-a)/(b-a)
    tmpscbeval = pdf(betatyp,xresc)/(b-a)

    # btype = "absampl"

    # if btype == "pdfkind"
    #     ## pdf version with scaling, the value of amplscale does NOT correspond
    #     ##   to the actual value of the peak
    #     scbeval = amplscale * tmpscbeval

    # elseif btype == "absampl"
        ## here the value of amplscale corresponds to the actual value of the peak
        xrescmode = (mo-a)/(b-a)
        ## value at the mode
        peakval = pdf(betatyp,xrescmode)/(b-a)
        ## normalize first, then rescale
        scbeval = amplscale * (tmpscbeval/peakval)
    # end

    return scbeval
end

###############################################################################

"""
$(TYPEDSIGNATURES)

Solve the inverse problem, i.e., fit the measured enthalpy data,
   using an Interior Point Newton method from the Optim.jl package.

# Arguments

- `betpar`: a struct containing the parameters for the Beta functions
            See [`ScaledBeta2DParams`](@ref)
- `dobs`: a struct containing the observed (measured) data and concentrations
            See [`ITCObsData`](@ref)
- `invCd`: inverse of the covariance matrix on observations (precision matrix)
- `mstart`: the starting model
- `lowconstr`: lower constraints
- `upconstr`: upper constraints
- `outdir`: output directory to save results

"""
function solveinvprob(betpar::ScaledBeta2DParams,dobs::ITCObsData,
                      invCd::Matrix{<:Real},mstart::Matrix{<:Real},
                      lowconstr::Vector{<:Real},upconstr::Vector{<:Real},
                      outdir::String)  

    # References

    # The algorithm was originally written by Tim Holy (@timholy, tim.holy@gmail.com).
    # J Nocedal, SJ Wright (2006), Numerical optimization, second edition. Springer.
    # A Wächter, LT Biegler (2006), On the implementation of an interior-point filter line-search algorithm for large-scale nonlinear programming. Mathematical Programming 106 (1), 25-57.

    ### https://github.com/JuliaNLSolvers/NLSolversBase.jl/blob/0f07421dbef63e8f53ca658e6be2981c09ca9c5b/src/objective_types/constraints.jl 
    ### Constraints
    #
    # Constraints are specified by the user as
    #    lx_i ≤   x[i]  ≤ ux_i  # variable (box) constraints
    #    lc_i ≤ c(x)[i] ≤ uc_i  # linear/nonlinear constraints
    # and become equality constraints with l_i = u_i. ±∞ are allowed for l
    # and u, in which case the relevant side(s) are unbounded.
    #
    # The user supplies functions to calculate c(x) and its derivatives.
    #
    # Of course we could unify the box-constraints into the
    # linear/nonlinear constraints, but that would force the user to
    # provide the variable-derivatives manually, which would be silly.

    # struct TwiceDifferentiableConstraints{F,J,H,T} <: Constraints
    #     c!::F # c!(storage, x) stores the value of the constraint-functions at x
    #     jacobian!::J # jacobian!(storage, x) stores the Jacobian of the constraint-functions
    #     h!::H   # h!(storage, x) stores the hessian of the constraint functions
    #     bounds::ConstraintBounds{T}
    # end
    
    ## Newton from Optim.jl
    #------------------------------------------------
    
    function closmisfitOPTIM(clmcur::Vector{<:Real})
        clmcur2d=reshape(clmcur,betpar.nummodpar,ncomp)
        msf = misfitbeta2D(betpar,dobs,invCd,clmcur2d) #invCm,mprior,mcur)
        return msf
    end

    #------------------------------------------------

    function fun_grad!(gr,mcur) 
        ## gradient of the misfit function
        gr .= ForwardDiff.gradient(closmisfitOPTIM,mcur)
        return nothing
    end

    #------------------------------------------------
    
    function fun_hess!(hess,mcur) 
        ## Hessian of the misfit function
        hess .= ForwardDiff.hessian(closmisfitOPTIM,mcur)
        return nothing
    end
 
    ##===================================================
    ## Some checks
    @assert issymmetric(invCd)
    @assert isposdef(invCd)

    
    ##===================================================
    ## https://julianlsolvers.github.io/Optim.jl/stable/#examples/generated/ipnewton_basics/
    
    ## 2d -> 1D
    mstartvec = vec(mstart)
    ncomp = size(mstart,2)
    
    ## Objective function without constraints
    df = TwiceDifferentiable(closmisfitOPTIM, fun_grad!, fun_hess!, mstartvec)

    ## Add constraints
    dfc = TwiceDifferentiableConstraints(lowconstr, upconstr)

    ## Run the Newton inversion with box minimization
    println("\nRunning optimization with IPNewton...")
    result = optimize(df, dfc, mstartvec, IPNewton())
    mpostvec = Optim.minimizer(result)

    ## reshape mpost
    mpost = reshape(mpostvec,size(mstart))

    ## structure holding Beta parameters and the solution in terms of
    ##   mode, confidence parameter and amplitude
    betmix = BetaMix2D(betpar,mpost,dobs.protein)

    ## show some info
    println(result)

    ##================================
    ## save data
    try
        mkdir(outdir)
    catch
        nothing
    end

    outfile = joinpath(outdir,dobs.protein*"_ITCresults.h5")
    println("Saving results to $outfile\n")

    hf = h5open(outfile,"w")
    hf["protein"] = dobs.protein
    hf["mpost"] = mpost
    hf["enthalpy"] = dobs.enthalpy
    hf["sdscon"] = dobs.sdsprotcon[:,1] # sdscon
    hf["procon"] = dobs.sdsprotcon[:,2] # procon
    hf["beta_a"] = betpar.a
    hf["beta_b"] = betpar.b
    hf["beta_ymin"] = betpar.ymax
    hf["beta_ymax"] = betpar.ymin
    hf["nummodpar"] = betpar.nummodpar
    hf["mstart"] = mstart
    hf["lowconstr"] = lowconstr
    hf["upconstr"] = upconstr
    hf["invCd"] = invCd
    hf["modefuncy"] = betpar.modefuny
    hf["konfuncy"]  = betpar.konfuny
    hf["ampfuncy"] = betpar.ampfuny
    close(hf)

    return betmix
end

########################################################################
########################################################################
