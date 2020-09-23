

using Revise
using ANISPROU
using LinearAlgebra

###########################################################

function launchall()

    ##===========================================================
    ## do the fitting of data with Beta functions
    inpdir="inputdata/"
    protein = "IM7"
    betamix,dobs = invertITCdata(inpdir,protein)

    # ## plot straight lines defined by parameters as a function of concentration
    #plotparamlines(betamix)

    # ## plot single components for fixed protein concentration
    #plotbetacomp1D(betamix,0.148)
   
    # # plot 3D surface from results
    # if isdefined(@__MODULE__,:Makie)
    #     plotsurface3D(dobs,betamix,yscal=70,ymin=0.0)
    # end

    # ## plot fit to single experiments
    # outdir="figs"
    # plotsingleexperiments(outdir,dobs,betamix)


    # ##===========================================================
    # ## compute area for all Beta components at
    # ##   requested protein concentration
    # protcon = 0.08
    # area,errarea = area_enthalpy(betamix,protcon)
    # println("\nArea at protein concentration $protcon for each component: $area, error $errarea\n")

    ## Area for a set of protein concentrations
    N = 15
    protcons = collect(LinRange(0.0,0.14,N))
    outdir = "output"
    areas,erras,linfitres = areasvsprotcon(betamix,protcons,outdir,protein)

    # ## plot areas only
    # plotareavsprotcon(proteinname,protcons,areas,linfitres,outdir)
    # ## plot parameters lines and areas
    # plotparamlines(betamix,protcons,areas)


    # ## compute volume for all Beta components within
    # ##   requested bounds of protein concentration
    # minprotcon = betamix.betpar.ymin
    # maxprotcon = betamix.betpar.ymax
    # volume,errvol = volume_enthalpy(betamix,minprotcon,maxprotcon)
    # println("\nVolume within bounds for each component: $volume, error $errvol\n")


    ##===========================================================
    ## construct binding isotherm for Beta mix using stationary and inflection points
    freeSDS,Nbound = bindingisotherm_betamix(betamix,dobs)


    # ## construct binding isotherm for each Beta function using selected amplitude points
    # # sdsfNbou_singlebetas = bindingisotherm_singlebetas(dobs,betamix)

    return betamix,dobs #,freeSDS,Nbound
end

###########################################################

function invertITCdata(inpdir::String,protein::String)

    ##===========================================
    ## read data
    ## data = readallexperiments("inputdata/",["IM7","FN3","sFN3"])
    data = readallexperiments(inpdir,[protein],discninitrows=0)
     

    #protein = "FN3" # "IM7" "FN3" "sFN3"
    sdscon = data[protein]["sdscon"]
    procon = data[protein]["procon"]
    enthalpy = data[protein]["enout"] 
    idxdata = data[protein]["idxdata"]


    ## x axis is SDS concentration, y axis is protein concentration
    ## xy = [sdscon procon]
    dobs = ITCObsData(protein=protein,enthalpy=enthalpy,idxdata=idxdata,
                      sdsprotcon=[sdscon procon])
    # plot only the observed data
    #plotobsdata(dobs)

    ##=========================
    ## forward test
    a = minimum(dobs.sdsprotcon[:,1]) # lower bound for beta domain
    b = 1.05*maximum(dobs.sdsprotcon[:,1]) # upper bound for beta domain

    minprotcon = minimum(dobs.sdsprotcon[:,2]) 
    maxprotcon = 1.05*maximum(dobs.sdsprotcon[:,2]) # 1.1* 
    
    @show a,b,minprotcon,maxprotcon

    ## define the parameters of Beta 2D functions
    modefuny = "linear"
    konfuny  = "linear" #"constant"
    ampfuny  = "linear"
    betpar = ScaledBeta2DParams(modefuny=modefuny,konfuny=konfuny,
                                ampfuny=ampfuny,a=a,b=b,
                                ymin=minprotcon,ymax=maxprotcon)

    ##=======================================
    ## Only testing initial guess?
    ## Used to visualise the calculated data from the initial guess
    ##   compared with the actual observed data, in order to manually
    ##   find a good initial model
    onlytestinitialguess = false #true

    ##===========================================
    ## Starting model

    if protein=="IM7"
        ## starting model IM7
        ## Each row is one component (one 2D Beta function)
        ##  Just add (remove) rows to increase (decrease) the number of components
        ## Elements are: 2 for mode, 2 for the confidence parameter and
        ##   2 for the amplitude parameter
   
        if betpar.konfuny=="constant"
            # constant kon
            comp1 = [0.6,  1.5,  30.0, -2.5,  -5.0 ]  
            comp2 = [1.7,  4.8,  50.0, -1.6,  -4.0 ] 
            comp3 = [4.5,  10.0, 35.0, 0.12, 0.16 ] 
            comp4 = [6.2,  15.3, 50.0, -1.6, -2.0 ]

        else
            comp1 = [0.6,  1.5,  35.0, 30.0, -2.5,  -5.0 ]  
            comp2 = [1.7,  4.8,  60.0, 40.0, -1.6,  -4.0 ] 
            comp3 = [4.5,  10.0, 40.0, 30.0, 0.12, 0.16 ] 
            comp4 = [6.2,  15.3, 60.0, 40.0, -1.6, -2.0 ]
        end

     end

    ## mstart is a 2D array where each column represents one component
    mstart = [comp1 comp2 comp3 comp4] 
    @show size(mstart)

    ## plotparamlines(BetaMix2D(betpar,mstart,dobs.protein))

    ############################################
    if onlytestinitialguess==true
        # Visualize initial guess
        plotinitialguess(betpar,dobs,mstart)

        ## stop here when testing initial models
        return nothing,nothing
    end

    ###############################################
    ## Set constraints
    ## Elements are: 2 for mode, 2 for the confidence parameter and
    ##   2 for the amplitude parameter
    # ## lower constraints 
    # lowc = [betpar.a,  betpar.a, 2.1, 2.1, -20.0, -20.0]
    # ## upper constraints
    # upc  = [betpar.b, betpar.b, 500.0, 500.0, 10.0, 10.0]
    # lowconstr,upconstr = setconstraints(betpar,mstart,lowc,upc)

    if betpar.konfuny=="constant"
        ## constant kon
        # lower constraints [confidence must be >2.0]
        lcs1 = [betpar.a, betpar.a, 2.1, -20.0, -20.0]
        lcs2 = [betpar.a, betpar.a, 2.1, -20.0, -20.0]
        lcs3 = [betpar.a, betpar.a, 2.1,   0.0,   0.0]
        lcs4 = [betpar.a, betpar.a, 2.1, -20.0, -20.0]
        lowconstr = [lcs1 lcs2 lcs3 lcs4]
        # upper constraints
        ucs1 = [betpar.b, betpar.b, 1000.0,  0.0,  0.0]
        ucs2 = [betpar.b, betpar.b, 1000.0,  0.0,  0.0]
        ucs3 = [betpar.b, betpar.b, 1000.0, 10.0, 10.0]
        ucs4 = [betpar.b, betpar.b, 1000.0,  0.0,  0.0]
        upconstr = [ucs1 ucs2 ucs3 ucs4]

    else
        # lower constraints [confidence must be >2.0]
        lcs1 = [betpar.a, betpar.a, 2.1, 2.1, -20.0, -20.0]
        lcs2 = [betpar.a, betpar.a, 2.1, 2.1, -20.0, -20.0]
        lcs3 = [betpar.a, betpar.a, 2.1, 2.1,   0.0,   0.0]
        lcs4 = [betpar.a, betpar.a, 2.1, 2.1, -20.0, -20.0]
        lowconstr = [lcs1 lcs2 lcs3 lcs4]
        # upper constraints
        ucs1 = [betpar.b, betpar.b, 1000.0, 1000.0,  0.0,  0.0]
        ucs2 = [betpar.b, betpar.b, 1000.0, 1000.0,  0.0,  0.0]
        ucs3 = [betpar.b, betpar.b, 1000.0, 1000.0, 10.0, 10.0]
        ucs4 = [betpar.b, betpar.b, 1000.0, 1000.0,  0.0,  0.0]
        upconstr = [ucs1 ucs2 ucs3 ucs4]

    end

    ## check mstart...
    @show mstart.<=lowconstr
    @show mstart.>=upconstr
    @assert all(mstart.>=lowconstr)
    @assert all(mstart.<=upconstr)

    
    ###############################################
    ## run the Newton optimization
    nobs = length(dobs.enthalpy)
    stdobs = 0.5 .* ones(nobs)  
    @show stdobs
    Cd = diagm(stdobs.^2)
    invCd = inv(Cd)

    ## IPNewton from Optim.jl, box constraints
    outdir = "output"
    betamix = solveinvprob(betpar,dobs,invCd,mstart,lowconstr,upconstr,
                           outdir,applynonlinconstr=true,constrarea=6.0) ## 2.0

    ##================================
    ## plot results
    outdir = "figs"
    plotresults(betamix,dobs,mstart,outdir)

    return betamix,dobs
end

#######################################################################

function bindingisotherm_betamix(betamix::BetaMix2D,dobs::ITCObsData)

    ##===========================================
    ## Find stationary and inflection points
    ny = 4 
    protcon = collect(LinRange(1.1*betamix.betpar.ymin,0.8*betamix.betpar.ymax,ny))
    statpts,inflpts = findcurvefeatures(betamix,protcon)

    # plot found points/features
    outdir = "figs"
    plotfoundfeatures(betamix,protcon,statpts,inflpts,outdir)

    ##===========================================
    ## Selection of local minima and maxima
    selectstatpts = Vector{Array{<:Real,2}}(undef,0)
    
    push!(selectstatpts, [ statpts[1][2] protcon[1];
                           statpts[2][2] protcon[2];
                           statpts[3][2] protcon[3];
                           statpts[4][2] protcon[4] ] )
    
    push!(selectstatpts, [ statpts[3][3] protcon[3];
                           statpts[4][3] protcon[4] ] )

    push!(selectstatpts, [ statpts[2][4] protcon[2];
                           statpts[3][4] protcon[3];
                           statpts[4][4] protcon[4] ] )

    push!(selectstatpts, [ statpts[1][5] protcon[1];
                           statpts[2][5] protcon[2];
                           statpts[3][5] protcon[3];
                           statpts[4][5] protcon[4] ] )

    push!(selectstatpts, [ statpts[1][6] protcon[1];
                           statpts[2][6] protcon[2];
                           statpts[3][6] protcon[3];
                           statpts[4][6] protcon[4] ] )
    

    ##===========================================
    ## Selection of inflection points
    selectinflpts = Vector{Array{<:Real,2}}(undef,0)
    
    push!(selectinflpts, [ inflpts[1][2] protcon[1];
                           inflpts[2][2] protcon[2];
                           inflpts[3][2] protcon[3];
                           inflpts[4][2] protcon[4] ] )
    
    push!(selectinflpts, [ inflpts[1][3] protcon[1];
                           inflpts[2][3] protcon[2];
                           inflpts[3][3] protcon[3];
                           inflpts[4][3] protcon[4] ] )

    push!(selectinflpts, [ inflpts[1][4] protcon[1];
                           inflpts[2][4] protcon[2];
                           inflpts[3][4] protcon[3];
                           inflpts[4][4] protcon[4] ] )
    
    push!(selectinflpts, [ inflpts[1][5] protcon[1];
                           inflpts[2][5] protcon[2];
                           inflpts[3][5] protcon[3];
                           inflpts[4][5] protcon[4] ] )

    push!(selectinflpts, [ inflpts[1][6] protcon[1];
                           inflpts[2][6] protcon[2];
                           inflpts[3][6] protcon[3];
                           inflpts[4][6] protcon[4] ] )

    push!(selectinflpts, [ inflpts[1][7] protcon[1];
                           inflpts[2][7] protcon[2];
                           inflpts[3][7] protcon[3];
                           inflpts[4][7] protcon[4] ] )


    ##===========================================
    ## Linear regression to get the parameters of straight lines

    ## find the straight line by least squares
    ## lstlines = [selectstatpts..., selectinflpts...]
    outdir = "output"
    freeSDS,Nbound = calcfreeSDSNbound(protcon,statpts,inflpts,
                                       selectstatpts,selectinflpts,
                                       outdir,dobs.protein)

    ##===========================================
    ## Plot stuff
    outdir = "figs"
    plotbindisotherm(betamix,protcon,dobs,statpts,inflpts,freeSDS,Nbound,outdir)

    return freeSDS,Nbound
end

##########################################################

# ######################################################################

# function bindingisotherm_singlebetas(dobs::ITCObsData,betamix::BetaMix2D)

#     betpar = betamix.betpar
#     mcur = betamix.modkonamp
#     xy = dobs.sdsprotcon
#     protein = betamix.protein

#     ncomp = size(mcur,2)
#     sdsfNbou = Array{Matrix{<:Real}}(undef,ncomp)

#     ## number of sampling points along x, i.e., [SDS]
#     nxpts = 5

#     ## loop on all Beta functions (components)
#     for ico=1:ncomp

#         ## define custom amplitude ratios values for each component
#         if ico in [1,3,4]
#             # ampratios = collect(LinRange(0.3,1.0,6))
#             ampratios = (collect(LinRange(0.4,1.0,nxpts))).^(1/2)
#         elseif ico==2
#             ampratios = (collect(LinRange(0.8,1.0,nxpts))).^(1/2)
#         end

#         ## calculate intercept and angular coefficient, i.e., [SDS]_free and N_bound
#         sdsfNbou[ico] = freeboundsds_singlebeta(betpar,mcur[:,ico],ampratios=ampratios)

#     end

#     ## plot results
#     outdir="figs"
#     plotbindisotherm_singlebetas(protein,betpar,dobs.sdsprotcon,mcur,sdsfNbou,outdir,Npts=100)

#     return sdsfNbou
# end

# ######################################################################

