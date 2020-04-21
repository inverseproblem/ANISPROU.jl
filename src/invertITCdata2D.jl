


using LinearAlgebra
using Distributions
using ForwardDiff
#using ReverseDiff
#using PositiveFactorizations
using Optim


# https://en.wikipedia.org/wiki/Protein_folding

#
# Derivatives incomplete gamma function
# https://en.wikipedia.org/wiki/Incomplete_gamma_function

###############################################

Base.@kwdef  struct ScaledBeta2DParams
    "number of model parameters"
    nummodpar::Integer
    "sample points at (x,y)"
    xy::AbstractArray{<:Real,2}
    "lower bound for Beta function (along x)"
    a::Real
    "upper bound for Beta function (along y)"
    b::Real
    "user defined minimum y"
    ymin::Real
    "user defined maximum y"
    ymax::Real
end


###########################################################

function setconstraints(betpar,ncomp,mstart,lowc,upc)

    lowconstr = zeros(betpar.nummodpar*ncomp)
    upconstr = zeros(betpar.nummodpar*ncomp)

    @show size(lowc),size(upc)
    for c=1:ncomp
        c1 = (c-1)*betpar.nummodpar
        for i=1:betpar.nummodpar
            lowconstr[c1+i] = lowc[i]
            upconstr[c1+i]  = upc[i]
        end
    end
    
    @show mstart.<=lowconstr
    @show mstart.>=upconstr

    ## check mstart...
    @assert all(mstart.>=lowconstr)
    @assert all(mstart.<=upconstr)

    return lowconstr,upconstr
end

######################################################################

function misfitbeta2D(betpar::ScaledBeta2DParams,dobs::AbstractVector{<:Real},
                      invCd::AbstractMatrix{<:Real},mcur::AbstractVector{<:Real})
    
    # compute synthetic data
    dcalc = forwmod2D(betpar,mcur)

    # calculate likelihood value
    difcalobs = dcalc.-dobs
    tmp1 = invCd * difcalobs
    misflkl = 0.5 * dot(difcalobs,tmp1)

    # total misfit
    misf = misflkl 

    return misf
end

###########################################################

function forwmod2D(betpar::ScaledBeta2DParams,mcur::AbstractVector{<:Real})
                   
    nummodparam = betpar.nummodpar
    npts = size(betpar.xy,1)
    lenm = length(mcur)
    # make sure that mcur has a multiple of N number of elements
    @assert mod(lenm,nummodparam)==0

    ncomp = div(lenm,nummodparam)

    sscbeta = zeros(eltype(mcur),npts)
    sscbeta .= zero(eltype(mcur))   #0.0 ## make sure it's zeroed

    for n=1:ncomp
        c1 = (n-1)*betpar.nummodpar
        # sum all the components
        sscbeta .+= singlescaledbeta2D(betpar,mcur[c1+1],mcur[c1+2],mcur[c1+3],
                                       mcur[c1+4],mcur[c1+5],mcur[c1+6])                                      
    end

    return sscbeta
end

###########################################################

"""
  2-Dfied scaled Beta pdf for given x and y locations.
"""
function singlescaledbeta2D(betpar::ScaledBeta2DParams,
                            mode1::Real,mode2::Real,kon1::Real,
                            kon2::Real,ampl1::Real,ampl2::Real)

    ## Straight line through 2 points:
    ##
    ## y = (y2-y1)/(x2-x1) * (x-x1) + y1
    ##

    ## limits of protein concentration
    x1 = betpar.ymin
    x2 = betpar.ymax
    @assert (x2-x1)!=0.0

    # pre-allocate array of arrays (one for each y value)
    npts = size(betpar.xy,1)
    scbeta2d = zeros(typeof(mode1),npts)
    for l=1:npts
        ## xy[l,1] -> SDS concentration
        ## xy[l,2] -> protein concentration
        xcur = betpar.xy[l,2] # protein concentration

        mo   = ((mode2-mode1)/(x2-x1))*(xcur-x1) + mode1
        kon  = ((kon2-kon1)/(x2-x1))*(xcur-x1) + kon1
        ampl = ((ampl2-ampl1)/(x2-x1))*(xcur-x1) + ampl1
        
        # beta pdf along x
        ## xy[l,1] -> SDS concentration
        scbeta2d[l] = scaledbeta(mo,kon,betpar.a,betpar.b,
                                 ampl,betpar.xy[l,1])
    end

    return scbeta2d
end

###########################################################

@inline function scaledbeta(mo::Real,kon::Real,a::Real,b::Real,
                            amplscale::Real,x::Real)
    
    @assert b>a
    @assert a <= mo <= b
    @assert kon>=2.0

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

###########################################################
###########################################################

function solveNewton2D(xy::AbstractArray{<:Real,2},a::Real,b::Real,
                       minprotcon::Real,maxprotcon::Real,
                       dobs::AbstractVector{<:Real},
                       invCd::AbstractMatrix{<:Real},
                       #invCm::AbstractMatrix{<:Real},
                       mstart::AbstractVector{<:Real},niter::Integer,newtonstep::Real)

    #------------------------------------------------
    function closmisfit(clmcur::AbstractVector{<:Real})
        msf = misfitbeta2D(xy,a,b,minprotcon,maxprotcon,dobs,invCd,clmcur) #invCm,mprior,mcur)
        return msf
    end
    #------------------------------------------------

    @assert newtonstep>0.0

    doplot =false

    mcur = copy(mstart)
    mall = zeros(length(mcur),niter)
    misfall = zeros(niter+1)
    misfall[1] = misfitbeta2D(xy,a,b,minprotcon,maxprotcon,dobs,invCd,mcur)
    println("Initial misfit: $(misfall[1]) ")
    print("Newton step factor: $newtonstep ")


    if doplot
        figure()
        subplot(223)
        title("model parameters")
        plot(mcur,".-",label="mcur")
        legend()
    end


    for it=1:niter
        print("\nNewton iteration #$it     ")

        ## gradient of the misfit function
        grad = ForwardDiff.gradient(closmisfit,mcur)

        ## Hessian matrix
        hess = ForwardDiff.hessian(closmisfit,mcur)

        ##---------------------------------
        ## 1. solve the linear system
        # H p = -g
        # Ax = y
        # x = A\y
        
        ## check positive definitess of Hessian
        # pdhess = isposdef(hess)
        # print("isposdef(Hessian): $(pdhess)")

        # if !pdhess
            ## Hessian is not positive definite,
            ##   so find a pos. def. approximation using 
            ##  PositiveFactorizations.jl
            F = cholesky(Positive, hess )
            # H p = -g
            ## Julia should use the appropriate algorithm
            ##  given F as a lower triangular Cholesky decomposition
            pvec =  -(F \ grad)

            # pdhess = isposdef(F.L*F.U)
            # print(" $(pdhess) ")

        # else
            
        #     # H p = -g
        #     pvec = -(hess \ grad)
            
        # end


        ##---------------------------------
        ## 2. update the current solution
        mcur = mcur .+ newtonstep .* pvec
        
        mall[:,it] .= mcur


        if doplot
            figure()
            subplot(221)
            title("Gradient")
            plot(grad,".-")
            subplot(222)
            title("Hessian")
            imshow(hess)
            colorbar()
            subplot(223)
            title("model parameters")
            plot(mcur,".-",label="mcur")
            legend()
            subplot(224)
            title("inv(Hessian)")
            imshow(inv(hess))
            colorbar()
        end
        

        ## misfit
        misfall[it+1] = misfitbeta2D(xy,a,b,minprotcon,maxprotcon,dobs,invCd,mcur)
        print("misfit: $(misfall[it+1]) ")

    end
    println()

    return mall,misfall
end

###########################################################


function graddescent2D(betpar::ScaledBeta2DParams,dobs::AbstractVector{<:Real},
                       invCd::AbstractMatrix{<:Real},mstart::AbstractVector{<:Real},
                       niter::Integer,graddescstep::Real)

    #------------------------------------------------
    function closmisfit(clmcur::AbstractVector{<:Real})
        msf = misfitbeta2D(betapar,dobs,invCd,clmcur) #invCm,mprior,mcur)
        return msf
    end
    #------------------------------------------------

    @assert graddescstep>0.0

    doplot =false

    mcur = copy(mstart)
    mall = zeros(length(mcur),niter)
    misfall = zeros(niter+1)
    misfall[1] = misfitbeta2D(betapar,dobs,invCd,mcur)
    println("Initial misfit: $(misfall[1]) ")
    print("Gradient descent, step factor: $graddescstep ")

    if doplot
        figure()
        subplot(223)
        title("model parameters")
        plot(mcur,".-",label="mcur")
        legend()
    end


    for it=1:niter
        print("\nGradient descent iteration #$it     ")

        ## gradient of the misfit function
        grad = ForwardDiff.gradient(closmisfit,mcur)

        ## update current model
        mcur = mcur .- graddescstep .* grad 
        
        mall[:,it] .= mcur


        if doplot 
            figure()
            subplot(121)
            title("Gradient")
            plot(grad,".-")
            subplot(122)
            title("model parameters")
            plot(mcur,".-",label="mcur")
            legend()
        end

        ## misfit
        misfall[it+1] = misfitbeta2D(betapar,dobs,invCd,mcur)
        print("misfit: $(misfall[it+1]) ")

    end
    println()

    return mall,misfall
end


########################################################################

function IPNewtonOptim(protein::String,betpar::ScaledBeta2DParams,dobs::AbstractVector{<:Real},
                       invCd::AbstractMatrix{<:Real},mstart::AbstractVector{<:Real},
                       lowconstr::AbstractVector{<:Real},upconstr::AbstractVector{<:Real})  

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

    # struct TwiceDifferentiableConstraints{F,J,H,T} <: AbstractConstraints
    #     c!::F # c!(storage, x) stores the value of the constraint-functions at x
    #     jacobian!::J # jacobian!(storage, x) stores the Jacobian of the constraint-functions
    #     h!::H   # h!(storage, x) stores the hessian of the constraint functions
    #     bounds::ConstraintBounds{T}
    # end

    
    ## Newton from Optim.jl
    #------------------------------------------------
    
    function closmisfitOPTIM(clmcur::AbstractVector{<:Real})
        msf = misfitbeta2D(betpar,dobs,invCd,clmcur) #invCm,mprior,mcur)
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
    ## https://julianlsolvers.github.io/Optim.jl/stable/#examples/generated/ipnewton_basics/

    ## Objective function without constraints
    df = TwiceDifferentiable(closmisfitOPTIM, fun_grad!, fun_hess!, mstart)

    ## Add constraints
    dfc = TwiceDifferentiableConstraints(lowconstr, upconstr)

    ## Run the Newton inversion with box minimization
    res = optimize(df, dfc, mstart, IPNewton())

    println(res)
    ##===================================================

    mpost = Optim.minimizer(result)
    @show mstart
    @show mpost

    ##================================
    ## save data
    hf = h5open("../output/"*protein*"_ITCresults.h5","w")
    hf["protein"] = protein
    hf["mpost"] = mpost
    hf["dobs"] = dobs
    hf["sdscon"] = sdscon
    hf["procon"] = procon
    hf["beta_a"] = betpar.a
    hf["beta_b"] = betpar.b
    hf["beta_ymin"] = betpar.ymax
    hf["beta_ymax"] = betpar.ymin
    hf["nummodpar"] = betpar.nummodpar
    hf["mstart"] = mstart
    hf["lowconstr"] = lowconstr
    hf["upconstr"] = upconstr
    hf["invCd"] = invCd
    close(hf)

    return mpost
end

########################################################################
########################################################################

function SimulatingAnnealingDrea(betaparxy::ScaledBeta2DParams,dobs::AbstractVector{<:Real},
                                 invCd::AbstractMatrix{<:Real},mstart::AbstractVector{<:Real},
                                 lowconstr::AbstractVector{<:Real},upconstr::AbstractVector{<:Real},
                                 niter::Integer)  

    #------------------------------------------------
      function closmisfitOPTIM(clmcur::AbstractVector{<:Real})
        msf = misfitbeta2D(betapar,dobs,invCd,clmcur) #invCm,mprior,mcur)
        return msf
    end
    #------------------------------------------------

    res = optimize(closmisfitOPTIM, lowconstr, upconstr, mstart, SAMIN(), Optim.Options(iterations=niter))

    return res
end

# ########################################################################

