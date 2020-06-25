
###########################################################

"""
$(TYPEDSIGNATURES)

Calculate the area of each individual Beta function for a given protein concentration. 
The SDS concentration axis needs a unit conversion which is handled by the `volumescal` 
arguments.

# Arguments

- `betamix`: structure of type `BetaMix2D` containing the parameters of the Beta functions, 
            the modes, confidence parameters and amplitudes and protein name
- `protcon`: protein concentration value (y axis) at which to perform the calculation of the area
- `volumescal`: scaling factor in μl to convert from mM to mole, instrument dependent. The default 
               is 203.0μl.
"""
function area_enthalpy(betamix::BetaMix2D,protcon::Real ; volumescal::Real=203.0)

    # @assert betamix.betpar.ymin<=protcon<=betamix.betpar.ymax
        
    # scale the SDS concentration from mM to mole 10^(-3)*203*10^(-6) =2.03*10^(-7).
    # the 203 µl is specific for certain instruments
    scalfact = volumescal*10^(-9)

    ##------------------------------------
    ## Define a new ScaledBeta2DParams with the new units
    betpar_inp = deepcopy(betamix.betpar)
    betpar_scal = ScaledBeta2DParams(modefuny=betpar_inp.modefuny,konfuny=betpar_inp.konfuny,
                                     ampfuny=betpar_inp.ampfuny,a=betpar_inp.a*scalfact,b=betpar_inp.b*scalfact,
                                     ymin=betpar_inp.ymin,ymax=betpar_inp.ymax)
    
    ## Define a new BetaMix2D with the new units
    modkonamp_scal = copy(betamix.modkonamp)
    # scale the mode in terms of SDS concentration:
    modkonamp_scal[1:2,:] .*= scalfact
    betamix_scal = BetaMix2D(betpar_scal,modkonamp_scal,betamix.protein)
    ##------------------------------------

    ncomp = size(betamix_scal.modkonamp,2)
    integr = zeros(ncomp)
    err = zeros(ncomp)
    for c=1:ncomp

        ## betamix_scal  _scal!!!
        bp = betamix_scal.betpar
        mcur = betamix_scal.modkonamp[:,c]

        @assert length(mcur)==bp.nummodpar
    
        ##---------------------------------------------------------
        # (val,err) = hquadrature(f::Function, xmin::Real, xmax::Real;
        #                     reltol=1e-8, abstol=0, maxevals=0)
        abeta = bp.a
        bbeta = bp.b
        ## get the values of model parameters at y=ycur
        mode,kon,amp = getmodparbeta(bp,mcur,protcon)

        integr[c],err[c] = hquadrature(x->scaledbeta(mode,kon,abeta,bbeta,amp,x),abeta,bbeta,
                                       reltol=1e-8,abstol=0,maxevals=0)

    end

    return integr,err
end

###########################################################

"""
$(TYPEDSIGNATURES)

Calculate the area of each individual Beta function for a set of  protein concentrations.

# Arguments

- `betamix`: structure of type `BetaMix2D` containing the parameters of the Beta functions, 
            the modes, confidence parameters and amplitudes and protein name
- `protcosn`: array of protein concentrations at which to calculate the areas
- `volumescal`: scaling factor in μl to convert from mM to mole, instrument dependent. The default 
               is 203.0μl.
"""
function areasvsprotcon(betamix::BetaMix2D,protcons::Vector{Float64}; volumescal::Real=203.0,
                        doplot::Bool=true)

    N=length(protcons)
    ncomp = size(betamix.modkonamp,2)
    areas = zeros(N,ncomp)
    errar = zeros(N,ncomp)
    for i=1:N
        areas[i,:],errar[i,:] = area_enthalpy(betamix,protcons[i],volumescal=volumescal)
    end
    
    if doplot
        plotareavsprotcon(protcons,areas)
    end

    return areas,errar
end

###########################################################

"""
$(TYPEDSIGNATURES)

Calculate the volume of each single Beta function within given bounds of protein concentration.
The SDS concentration axis needs a unit conversion which is handled by the `volumescal` 
arguments.

# Arguments

- `betamix`: structure of type `BetaMix2D` containing the parameters of the Beta functions, 
            the modes, confidence parameters and amplitudes and protein name
- `minprotcon`: lower bound of protein concentration value (y axis) to perform the integral
- `maxprotcon`: upper bound of protein concentration value (y axis) to perform the integral
- `volumescal`=203.0: scaling factor in μl to convert from mM to mole, instrument dependent. The default 
               is 203.0μl.
"""
function volume_enthalpy(betamix::BetaMix2D,minprotcon::Real,maxprotcon::Real ; volumescal::Real=203.0)

    @assert betamix.betpar.ymin<=minprotcon<=betamix.betpar.ymax
    @assert betamix.betpar.ymin<=maxprotcon<=betamix.betpar.ymax
    @assert minprotcon<maxprotcon

    # scale the SDS concentration from mM to mole 10^(-3)*203*10^(-6) =2.03*10^(-7).
    # the 203 µl is specific for certain instruments
    scalfact = volumescal*10^(-9)

    ##------------------------------------
    ## Define a new ScaledBeta2DParams with the new units
    betpar_inp = deepcopy(betamix.betpar)
    betpar_scal = ScaledBeta2DParams(modefuny=betpar_inp.modefuny,konfuny=betpar_inp.konfuny,
                                     ampfuny=betpar_inp.ampfuny,a=betpar_inp.a*scalfact,b=betpar_inp.b*scalfact,
                                     ymin=betpar_inp.ymin,ymax=betpar_inp.ymax)
    
    ## Define a new BetaMix2D with the new units
    modkonamp_scal = copy(betamix.modkonamp)
    # scale the mode in terms of SDS concentration:
    modkonamp_scal[1:2,:] .*= scalfact
    betamix_scal = BetaMix2D(betpar_scal,modkonamp_scal,betamix.protein)
    ##------------------------------------


    ##---------------------------------------------------
    function betaval(xy::Vector{<:Real})::Real
        ## get the values of model parameters at y=ycur
        x=xy[1]
        y=xy[2]
        mode,kon,amp = getmodparbeta(bp,mcur,y)
        ## calculate value of beta function
        bval = scaledbeta(mode,kon,bp.a,bp.b,amp,x)
        return bval 
    end
    ##---------------------------------------------------
    
    ##---------------------------------------------------------
    # (val,err) = hquadrature(f::Function, xmin::Real, xmax::Real;
    #                     reltol=1e-8, abstol=0, maxevals=0)
    bp = betamix_scal.betpar
    xmin = [bp.a,minprotcon]
    xmax = [bp.b,maxprotcon]

    ncomp = size(betamix_scal.modkonamp,2)
    integr = zeros(ncomp)
    err = zeros(ncomp)
    mcur=nothing # given the closure of betaval() and the private scope of for loop
    for c=1:ncomp
        mcur = betamix_scal.modkonamp[:,c]
        integr[c],err[c] = hcubature(betaval,xmin,xmax,reltol=1e-8,abstol=0,maxevals=0)
    end

    return integr,err
end


############################################################

