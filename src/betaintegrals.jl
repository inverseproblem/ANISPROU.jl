

###########################################################

function area_betamix(betamix::BetaMix2D,protcon::Real)

    ##-------------------------------------------------
    ## betaval(x::Real) = scaledbeta(mode,kon,abeta,bbeta,amp,x)

    ##------------------------------------
    ncomp=size(betamix.modkonamp,2)
    integr = zeros(ncomp)
    err = zeros(ncomp)
    for c=1:ncomp

        betpar = betamix.betpar
        mcur = betamix.modkonamp[:,c]

        @assert length(mcur)==betpar.nummodpar
        
        ## get the values of model parameters at y=ycur
        mode,kon,amp = getmodparbeta(betpar,mcur,protcon)

        ##---------------------------------------------------------
        # (val,err) = hquadrature(f::Function, xmin::Real, xmax::Real;
        #                     reltol=1e-8, abstol=0, maxevals=0)
        abeta = betpar.a
        bbeta = betpar.b
        yval = protcon
        ## get the values of model parameters at y=ycur
        mode,kon,amp = getmodparbeta(betpar,mcur,yval)

        integr[c],err[c] = hquadrature(x->scaledbeta(mode,kon,abeta,bbeta,amp,x),abeta,bbeta,
                                       reltol=1e-8,abstol=0,maxevals=0)

    end

    return integr,err
end

###########################################################

function volume_betamix(betamix::BetaMix2D,minprotcon::Real,maxprotcon::Real)

    ##---------------------------------------------------
    function betaval(xy::Vector{<:Real})::Real
        ## get the values of model parameters at y=ycur
        x=xy[1]
        y=xy[2]
        mode,kon,amp = getmodparbeta(betpar,mcur,y)
        ## calculate value of beta function
        bval = scaledbeta(mode,kon,betpar.a,betpar.b,amp,x)
        return bval 
    end
    ##---------------------------------------------------
  
    ##---------------------------------------------------------
    # (val,err) = hquadrature(f::Function, xmin::Real, xmax::Real;
    #                     reltol=1e-8, abstol=0, maxevals=0)
    betpar = betamix.betpar
    xmin = [betpar.a,minprotcon]
    xmax = [betpar.b,maxprotcon]

    ncomp = size(betamix.modkonamp,2)
    integr = zeros(ncomp)
    err = zeros(ncomp)
    mcur=nothing # given the closure of betaval() and the private scope of for loop
    for c=1:ncomp
        mcur = betamix.modkonamp[:,c]
        integr[c],err[c] = hcubature(betaval,xmin,xmax,
                                     reltol=1e-8,abstol=0,maxevals=0)
    end

    return integr,err
end

###########################################################
