

using Revise
using ANISPROU
using LinearAlgebra

###########################################################

function launchall()

    ## do the fitting of data with Beta functions
    inpdir="inputdata/"
    protein = "IM7"
    betamix,dobs = invertITCdata(inpdir,protein)

    ## construct binding isotherm using selected amplitude points
    sdsfNbou = bindingisotherm_singlebetas(protein,betamix)
    
    ## compute area for all Beta components at
    ##   requested protein concentration
    protcon = 0.08
    area,errarea = area_betamix(betamix,protcon)
    println("\nArea for all components: $area\n")

    ## compute volume for all Beta components within
    ##   requested bounds of protein concentration
    minprotcon = betamix.betpar.ymin
    maxprotcon = betamix.betpar.ymax
    volume,errvol = volume_betamix(betamix,minprotcon,maxprotcon)
    println("\nVolume for all component: $volume\n")

    return #mpost,betpar,dobs,area,sdsfNbou,volume
end

###########################################################

function bindingisotherm_singlebetas(protein::String,betamix::BetaMix2D)

    betpar = betamix.betpar
    mcur = betamix.modkonamp

    ncomp = size(mcur,2)
    sdsfNbou = Array{Matrix{<:Real}}(undef,ncomp)

    ## number of sampling points along x, i.e., [SDS]
    nxpts = 5

    ## loop on all Beta functions (components)
    for ico=1:ncomp

        ## define custom amplitude ratios values for each component
        if ico in [1,3,4]
            # ampratios = collect(LinRange(0.3,1.0,6))
            ampratios = (collect(LinRange(0.4,1.0,nxpts))).^(1/2)
        elseif ico==2
            ampratios = (collect(LinRange(0.8,1.0,nxpts))).^(1/2)
        end

        ## calculate intercept and angular coefficient, i.e., [SDS]_free and N_bound
        sdsfNbou[ico] = freeboundsds_singlebeta(betpar,mcur[:,ico],ampratios=ampratios)

    end

    ## plot results
    plotbindingisotherm(protein,betpar,mcur,sdsfNbou,"figs",Npts=100)

    return sdsfNbou
end

######################################################################

function invertITCdata(inpdir::String,protein::String)

    ##===========================================
    ## read data
    ## data = readallexperiments("inputdata/",["IM7","FN3","sFN3"])
    data = readallexperiments(inpdir,[protein])
    
    #protein = "FN3" # "IM7" "FN3" "sFN3"
    sdscon = data[protein]["sdscon"]
    procon = data[protein]["procon"]
    enout  = data[protein]["enout"] 

    dobs = enout

    ##=========================
    ## forward test
    a = minimum(sdscon) # lower bound for beta domain
    b = 1.2*maximum(sdscon) # upper bound for beta domain

    minprotcon = minimum(procon) 
    maxprotcon = maximum(procon)
    # x axis is SDS concentration, y axis is protein concentration
    xy = [sdscon procon] 

    ## define the parameters of Beta 2D functions
    modefuny = "linear"
    konfuny  = "linear"
    ampfuny  = "linear"
    betpar = ScaledBeta2DParams(modefuny=modefuny,konfuny=konfuny,
                                ampfuny=ampfuny,xy=xy,a=a,b=b,
                                ymin=minprotcon,ymax=maxprotcon)

    ##=======================================
    ## Only testing initial guess?
    ## Used to visualise the calculated data from the initial guess
    ##   compared with the actual observed data, in order to manually
    ##   find a good initial model
    onlytestinitialguess = false

    ##===========================================
    ## Starting model

    if protein=="IM7"
        ## starting model IM7
        ## Each row is one component (one 2D Beta function)
        ##  Just add (remove) rows to increase (decrease) the number of components
        ## Elements are: 2 for mode, 2 for the confidence parameter and
        ##   2 for the amplitude parameter
        comp1 = [0.37, 2.5,  50.0, 30.0, -500.0, -1250.0]
        comp2 = [1.4,  4.8,  60.0, 40.0, -500.0,  -800.0]
        comp3 = [4.5,  10.0, 100.0, 80.0,   30.0,    40.0]
        comp4 = [6.0,  15.0,  80.0, 50.0, -400.0,  -500.0]

    elseif protein=="FN3"
        ## starting model FN3
        ## Each row is one component (one 2D Beta function)
        ##  Just add (remove) rows to increase (decrease) the number of components
        ## Elements are: 2 for mode, 2 for the confidence parameter and
        ##   2 for the amplitude parameter
        comp1 = [0.4,  1.4,  50.0, 30.0, -500.0, -1250.0]
        comp2 = [1.8,  4.0,  60.0, 40.0, -500.0,  -800.0]
        comp3 = [3.2,  6.4, 100.0, 80.0,   30.0,    40.0]
        comp4 = [4.5,  8.0,  80.0, 50.0, -400.0,  -500.0]

    elseif protein=="sFN3"
        ## starting model sFN3
        ## Each row is one component (one 2D Beta function)
        ##  Just add (remove) rows to increase (decrease) the number of components
        ## Elements are: 2 for mode, 2 for the confidence parameter and
        ##   2 for the amplitude parameter
        comp1 = [0.4,  1.4,  50.0, 30.0, -500.0, -1250.0]
        comp2 = [1.8,  4.0,  60.0, 40.0, -500.0,  -800.0]
        comp3 = [3.2,  6.4, 100.0, 80.0,   30.0,    40.0]
        comp4 = [4.5,  8.0,  80.0, 50.0, -400.0,  -500.0]

    end

    ## mstart is a 2D array where each column represents one component
    mstart = [comp1 comp2 comp3 comp4]
    @show size(mstart)

    ############################################
    if onlytestinitialguess==true
        # Visualize initial guess
        plotinitialguess(protein,betpar,dobs,mstart)

        ## stop here when testing initial models
        return nothing
    end

    ###############################################
    ## Set constraints
    ## Elements are: 2 for mode, 2 for the confidence parameter and
    ##   2 for the amplitude parameter
    ## lower constraints 
    lowc = [betpar.a,  betpar.a, 2.0, 2.0, -2000.0, -2000.0]
    ## upper constraints
    upc  = [betpar.b, betpar.b, 500.0, 500.0, 2000.0, 2000.0]

    lowconstr,upconstr = setconstraints(betpar,mstart,lowc,upc)

    ###############################################
    ## run the Newton optimization
    nobs = length(dobs)
    stdobs = 2.0 .* ones(nobs)
    invCd = inv(diagm(stdobs.^2))

    ## IPNewton from Optim.jl, box constraints
    outdir = "output"
    betamix = solveinvprob(protein,betpar,dobs,invCd,mstart,lowconstr,upconstr,outdir)
 
    ##================================
    ## plot results
    outdir = "figs"
    plotresults(protein,betamix.betpar,dobs,mstart,betamix.modkonamp,outdir)

    return betamix,dobs
end

######################################################################


