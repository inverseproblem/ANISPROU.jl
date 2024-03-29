

################################################################
"""
$(TYPEDSIGNATURES)

Plot the observed data as a scatter plot.
"""
function plotobsdata(dobs)
    
    sdscon = dobs.sdsprotcon[:,1]
    procon = dobs.sdsprotcon[:,2]
    protein = dobs.protein

    PyPlot.figure(figsize=(11,6))
    PyPlot.title("$protein observed data")
    #tricontour(sdscon,procon,dobs,cmap=PyPlot.get_cmap("rainbow"))
    PyPlot.scatter(sdscon,procon,c=dobs.enthalpy,cmap=PyPlot.get_cmap("rainbow"))
    PyPlot.colorbar(label="Enthalpy [kJ/mol]")
    PyPlot.xlabel("SDS concentration [mM]")
    PyPlot.ylabel("Protein concentration [mM]")    
   
    PyPlot.tight_layout()
    return
end

####################################################################

"""
$(TYPEDSIGNATURES)

Plot the fit to the enthalpy data as a result of the initial guess, i.e., the starting model parameters of the Beta functions.
"""
function plotinitialguess(betpar,dobs,mstart)
    
    dcalcstart = forwmod2D(betpar,dobs.sdsprotcon,mstart)
    sdscon = dobs.sdsprotcon[:,1]
    procon = dobs.sdsprotcon[:,2]
    protein = dobs.protein

    #lenm = length(mstart)
    nummodparam = size(mstart,1) #betpar.nummodpar
    #ncomp = size(mstart,2) #div(lenm,nummodparam)
    ym = LinRange(betpar.ymin,betpar.ymax,40)

    PyPlot.figure(figsize=(12,5))
    #suptitle("Observed data and initial guess")

    PyPlot.subplot(121)
    PyPlot.title("$protein observed data")
    #tricontour(sdscon,procon,dobs,cmap=PyPlot.get_cmap("rainbow"))
    PyPlot.scatter(sdscon,procon,c=dobs.enthalpy,cmap=PyPlot.get_cmap("rainbow"))
    PyPlot.colorbar(label="Enthalpy [kJ/mol]")
    PyPlot.xlabel("SDS concentration [mM]")
    PyPlot.ylabel("Protein concentration [mM]")    
    plotmodelines(betpar,mstart,"start. m.")

    PyPlot.subplot(122)
    PyPlot.title("Calculated data from starting model")
    #tricontour(sdscon,procon,dcalcstart,cmap=PyPlot.get_cmap("rainbow"))
    PyPlot.scatter(sdscon,procon,c=dcalcstart,cmap=PyPlot.get_cmap("rainbow"))
    PyPlot.colorbar(label="Enthalpy [kJ/mol]")
    PyPlot.xlabel("SDS concentration [mM]")
    PyPlot.ylabel("Protein concentration [mM]")
    plotmodelines(betpar,mstart,"start. m.")
   
    PyPlot.tight_layout()
    return
end

############################################
"""
$(TYPEDSIGNATURES)

Plot the results of inverting the ITC data to fit the enthalpy function in 2D.
"""
function plotresults(betamix,dobs,mstart,outdir)

    mpost = betamix.modkonamp
    betpar = betamix.betpar
    dcalcstart = forwmod2D(betpar,dobs.sdsprotcon,mstart)
    dcalccur = forwmod2D(betpar,dobs.sdsprotcon,mpost)
    sdscon = dobs.sdsprotcon[:,1]
    procon = dobs.sdsprotcon[:,2]
    protein = dobs.protein

    ##==============================
    ## show data
    lenm = length(mstart)
    nummodparam = betpar.nummodpar
    #ncomp = div(lenm,nummodparam)
    ym = LinRange(betpar.ymin,betpar.ymax,40)


    PyPlot.figure(figsize=(14,7))
    PyPlot.subplot(231)
    PyPlot.title("$protein observed data")
    #tricontour(sdscon,procon,dobs,cmap=PyPlot.get_cmap("rainbow"))
    PyPlot.scatter(sdscon,procon,c=dobs.enthalpy,cmap=PyPlot.get_cmap("rainbow"))
    PyPlot.colorbar(label="Enthalpy [kJ/mol]")
    PyPlot.xlabel("SDS concentration [mM]")
    PyPlot.ylabel("Protein concentration [mM]")
    plotmodelines(betpar,mstart,"start. m.")
    # xlim([-2,betpar.b])
    # ylim([0.9*betpar.ymin,1.05*betpar.ymax])
    #plotmodelines(betpar,mpost,"mpost")
    
    PyPlot.subplot(232)
    PyPlot.title("Calculated data from starting model")
    #tricontour(sdscon,procon,dcalcstart,cmap=PyPlot.get_cmap("rainbow"))
    PyPlot.scatter(sdscon,procon,c=dcalcstart,cmap=PyPlot.get_cmap("rainbow"))
    PyPlot.colorbar(label="Enthalpy [kJ/mol]")
    PyPlot.xlabel("SDS concentration [mM]")
    PyPlot.ylabel("Protein concentration [mM]")
    plotmodelines(betpar,mstart,"start. m.")

    PyPlot.subplot(233)
    PyPlot.title("Model parameters")
    PyPlot.plot(vec(mstart),".-",label="starting model")
    PyPlot.plot(vec(mpost),".-",label="posterior model")
    PyPlot.legend()

    PyPlot.subplot(234)
    PyPlot.title("Calculated data from posterior model")
    #tricontour(sdscon,procon,dcalccur,cmap=PyPlot.get_cmap("rainbow"))
    PyPlot.scatter(sdscon,procon,c=dcalccur,cmap=PyPlot.get_cmap("rainbow"))
    PyPlot.colorbar(label="Enthalpy [kJ/mol]")
    PyPlot.xlabel("SDS concentration [mM]")
    PyPlot.ylabel("Protein concentration [mM]")
    plotmodelines(betpar,mpost,"post. m.")    
    
    
    PyPlot.subplot(235)
    PyPlot.title("Calculated minus observed data")
    #tricontour(sdscon,procon,dcalccur,cmap=PyPlot.get_cmap("rainbow"))
    dcmo = dcalccur-dobs.enthalpy
    vmax = maximum(abs.(dcmo))
    PyPlot.scatter(sdscon,procon,c=dcmo,vmin=-vmax,vmax=vmax,cmap=PyPlot.get_cmap("RdBu"))
    PyPlot.colorbar(label="ΔEnthalpy [kJ/mol]")
    PyPlot.xlabel("SDS concentration [mM]")
    PyPlot.ylabel("Protein concentration [mM]")
    # xlim([-2,betpar.b])
    # ylim([0.9*betpar.ymin,1.05*betpar.ymax])
    
    PyPlot.subplot(236)
    PyPlot.title("Calculated vs. observed data")
    PyPlot.plot(dcalcstart,"--",label="calc. data from start. mod.",linewidth=0.8)
    PyPlot.plot(dobs.enthalpy,"-k",label="observed data")
    PyPlot.plot(dcalccur,"-r",label="calc. data from post. mod.")
    PyPlot.legend()

    # subplot(235)
    # title("Misfit")
    # plot(misfall,"o-")
    # xlabel("Iteration number")
    # ylabel("misfit")

    PyPlot.tight_layout()
    
    try
        mkdir(outdir)
    catch
        nothing
    end
    PyPlot.savefig(joinpath(outdir,protein*"_results.pdf"))
    return nothing
end

############################################

"""
$(TYPEDSIGNATURES)

Plot the lines defined by the modes (default) or, optionally, other 
  parameters (see `firstidpar`).
"""
function plotmodelines(betpar,mcur,modname; firstidpar=1)
    #Npoints = 40
    ncomp = size(mcur,2)
    ym = [betpar.ymin,betpar.ymax] #LinRange(betpar.ymin,betpar.ymax,Npoints)
    y1 = betpar.ymin
    y2 = betpar.ymax
    for i=1:ncomp
        i2 = firstidpar+1
        i1 = firstidpar
        xm = (mcur[i2,i].-mcur[i1,i])./(y2.-y1).*(ym.-y1) .+ mcur[i1,i]
        if i==ncomp
            PyPlot.plot(xm,ym,"-k",linewidth=0.6,label="mode $modname")
        else
            PyPlot.plot(xm,ym,"-k",linewidth=0.6)
        end
        PyPlot.plot(xm[[1,end]],ym[[1,end]],"ok",markersize=2.0)
    end
    PyPlot.legend(loc="lower right")
    # xlim([-2,betpar.b])
    # ylim([0.9*betpar.ymin,1.05*betpar.ymax])
    return
end

############################################

"""
$(TYPEDSIGNATURES)

Plot the lines defined by the model parameters as a function of protein concentration.
"""
function plotparamlines(betamix,protcons=nothing,areas=nothing)

    betpar = betamix.betpar
    mcur = betamix.modkonamp
    
    @assert betpar.modefuny=="linear"
    #@assert betpar.konfuny=="linear"
    if betpar.konfuny=="constant"
        
    elseif betpar.konfuny=="linear"

    end
    @assert betpar.ampfuny=="linear"

    Npoints = 100
    ncomp = size(mcur,2)
    xm = LinRange(0.0,betpar.ymax,Npoints)
    parname = ["Mode","Confidence par.","Amplitude"]


    PyPlot.figure(figsize=(12,9))
    ylab = ["SDS concentration [mM]",
            "Confidence parameter",
            "Enthalpy"]
    
    for ip=1:3

        PyPlot.subplot(2,2,ip)
        PyPlot.title("$(parname[ip])")
        PyPlot.xlabel("Protein concentration [mM]")
        PyPlot.axvline(betpar.ymin,color="black")
        PyPlot.axvline(betpar.ymax,color="black")

        if ip==1
            PyPlot.axhline(betpar.a,color="black")
            PyPlot.axhline(betpar.b,color="black")
        end

        for i=1:ncomp
            if ip==1
                firstidpar = betpar.idxmode
                i2 = firstidpar+1
                i1 = firstidpar
                x1 = betpar.ymin
                x2 = betpar.ymax
                ym = (mcur[i2,i].-mcur[i1,i])./(x2.-x1).*(xm.-x1) .+ mcur[i1,i] 
                
            elseif ip==2
                if betpar.konfuny=="constant"
                    firstidpar = betpar.idxkon
                    ym = mcur[firstidpar,i].*ones(length(xm))

                elseif betpar.konfuny=="linear"
                    firstidpar = betpar.idxkon
                    i2 = firstidpar+1
                    i1 = firstidpar
                    x1 = betpar.ymin
                    x2 = betpar.ymax
                    ym = (mcur[i2,i].-mcur[i1,i])./(x2.-x1).*(xm.-x1) .+ mcur[i1,i] 

                end
                    
            elseif ip==3
                firstidpar = betpar.idxamp
                i2 = firstidpar+1
                i1 = firstidpar
                x1 = betpar.ymin
                x2 = betpar.ymax
                ym = (mcur[i2,i].-mcur[i1,i])./(x2.-x1).*(xm.-x1) .+ mcur[i1,i] 

            end

            PyPlot.plot(xm,ym,"-",linewidth=1,label="comp. $i")
            PyPlot.plot(xm[[1,end]],ym[[1,end]],".k",markersize=5.0)
            PyPlot.ylabel(ylab[ip])

        end
        PyPlot.legend()
        
        # xlim([-2,betpar.b])
        # ylim([0.9*betpar.ymin,1.05*betpar.ymax])
    end

    if areas!=nothing
        PyPlot.subplot(224)
        PyPlot.title("Areas of single components")
        for i=1:ncomp
            PyPlot.plot(protcons,areas[:,i],"o-",label="Beta comp. #$i")
        end
        PyPlot.legend()
        PyPlot.xlabel("Protein concentration [mM]")
        PyPlot.ylabel("Area [kJ]")
    end

    PyPlot.tight_layout()

    return
end

#####################################################

# """
# $(TYPEDSIGNATURES)

# Plot the results of binding isotherm calculations using the individual Beta functions.
# """
# function plotbindisotherm_singlebetas(protein,betpar,xy,mpost,sdsfNb,outdir; Npts=100)

#     dcalc = forwmod2D(betpar,xy,mpost)    
#     sdscon = xy[:,1]
#     procon = xy[:,2]
    
#     figure(figsize=(13,5))
    
#     subplot(121)
#     title("Protein $protein fitted lines on given amplitude ratios")
#     scatter(sdscon,procon,c=dcalc,cmap=PyPlot.get_cmap("rainbow"))
#     colorbar()
#     xlabel("SDS concentration [mM]")
#     ylabel("Protein concentration [mM]")
#     ##
#     ##  y = m*x + y0
#     ##
#     ## One line per row
#     ##N = 100
#     firstplot=true
#     for ico=1:length(sdsfNb)
#         for r=1:size(sdsfNb[ico],1)
#             x = collect(LinRange(betpar.ymin,betpar.ymax,Npts))
#             y = sdsfNb[ico][r,2].*x .+ sdsfNb[ico][r,1]
#             if firstplot
#                 plot(y,x,"--r",linewidth=0.7,label="selected fitted amplitudes")
#                 firstplot = false
#             else
#                 plot(y,x,"--r",linewidth=0.7)
#             end
#         end
#     end

#     plotmodelines(betpar,mpost,"")

#     subplot(122)
#     title("Protein $protein: binding isotherm from individual Beta functions")
#     for ico=1:length(sdsfNb)
#         plot(sdsfNb[ico][:,1],sdsfNb[ico][:,2],"o-")
#     end
#     xlabel("free SDS conc. [mM]")
#     ylabel("Nbound")

#     tight_layout()

#     savefig(joinpath(outdir,protein*"_bind-isoth_singlebetas.pdf"))
#     return nothing
# end

######################################################################

"""
$(TYPEDSIGNATURES)

Save the fitting surface in the VTK format for Paraview.
"""
function saveresultVTK(protein,betpar,mpost)

    ## VTK stuff
    Nx=500
    Ny=500
    xvtk = collect(LinRange(betpar.a,betpar.b,Nx))
    yvtk = collect(LinRange(betpar.ymin,betpar.ymax,Ny))
    xy = zeros(Nx*Ny,2)
    for j=1:Ny, i=1:Nx
        xy[i+Nx*(j-1),:] .= (xvtk[i],yvtk[j])
    end
    vtkfile = vtk_grid("../output/"*protein*"_gridITCdata", xvtk, 100.0 .* yvtk)
    betparvtk = ScaledBeta2DParams(nummodpar=nummodpar,a=betpar.a,b=betpar.b,
                                   ymin=betpar.ymin,ymax=betpar.ymax)
    dcalcvtk = forwmod2D(betparvtk,xy,mpost)
    dcvtk2 = reshape(dcalcvtk,Nx,Ny)
    vtkfile["energyout",VTKPointData()] = dcvtk2
    outfiles = vtk_save(vtkfile)

    return 
end

####################################################
"""
$(TYPEDSIGNATURES)

Plot the results of binding isotherm calculations using the Beta mix, i.e., sum of all Beta fitting functions.

"""
function plotbindisotherm(betamix,protcon,dobs,statpts,inflpts,freeSDS,Nbound,outdir;
                          resstdev=nothing)

    protein = betamix.protein
    N = 1000
    x = collect(LinRange(betamix.betpar.a,betamix.betpar.b,N))

    ##-------------------------------------
    PyPlot.figure(figsize=(12,10))
    
    PyPlot.subplot(211)
    PyPlot.title("Protein $protein: binding isotherm from stationary points of Beta mix")
    PyPlot.scatter(dobs.sdsprotcon[:,1],dobs.sdsprotcon[:,2],c=dobs.enthalpy,cmap=PyPlot.get_cmap("rainbow"))
    PyPlot.colorbar(label="Enthalpy [kJ/mol]")
    PyPlot.xlabel("SDS concentration [mM]")
    PyPlot.ylabel("Protein concentration [mM]")    

    ny = length(protcon)
    for i in 1:ny
        for z in 1:length(statpts[i])
            PyPlot.plot(statpts[i][z],protcon[i],"or")
            PyPlot.text(statpts[i][z],protcon[i],string(z),color="red")
        end
        
        for z in 1:length(inflpts[i])
            PyPlot.plot(inflpts[i][z],protcon[i],"ob")
            PyPlot.text(inflpts[i][z],protcon[i],string(z),color="blue")
        end
    end

    nlines = length(freeSDS)
    for i=1:nlines
        ## totSDS = freeSDS + Nbound * protc
        # y = angcoe.*x .+ intercept
        # x = (-intercept .+ 0)./angcoe
        ypl = (1.0/Nbound[i]).*x .+ (-freeSDS[i]/Nbound[i])
        PyPlot.plot(x,ypl,"-k",linewidth=0.5)
    end

    # xlim([betamix.betpar.a,betamix.betpar.b])
    PyPlot.ylim([0.92*betamix.betpar.ymin,1.05betamix.betpar.ymax])

    ##
    ## totSDS = freeSDS + Nbound * protc
    ##
    PyPlot.subplot(212)
    PyPlot.title("Binding isotherm for protein $protein")
    xbi = [0.0, freeSDS...]
    ybi = [0.0, Nbound...]
    if resstdev!=nothing
        resstdevbi = [0.0,resstdev...]
        PyPlot.scatter(xbi,ybi,c=resstdevbi,cmap=PyPlot.cm_get_cmap("rainbow"))
        PyPlot.colorbar(label="residuals standard dev.")
    else
        PyPlot.scatter(xbi,ybi)
    end
    PyPlot.plot(xbi,ybi,"-")
    PyPlot.xlabel("free SDS concentration [mM]")
    PyPlot.ylabel("Nbound")

    PyPlot.tight_layout()
    try
        mkdir(outdir)
    catch
        nothing
    end
    PyPlot.savefig(joinpath(outdir,protein*"_bind-isoth_betamix.pdf"))
    return
end

#########################################################

"""
$(TYPEDSIGNATURES)

Plot the points/features defined on the enthalpy curves to find the binding isotherm.
"""
function plotfoundfeatures(betamix,protcon,statpts,inflpts,outdir)
    
    protein = betamix.protein
    N = 1000
    x = collect(LinRange(betamix.betpar.a,betamix.betpar.b,N))
    ny = length(protcon)

    function fwd1ptxy(xpt::Real,ypt::Real)
        xy = hcat(xpt,ypt)
        out = forwmod2D(betamix.betpar,xy,betamix.modkonamp)
        return out[1] 
    end

    PyPlot.figure(figsize=(12,6))
    PyPlot.title("Protein $protein: stationary and inflection points")
    for i in 1:ny
        y = [protcon[i] for _ in 1:N]
        slbet = fwd1ptxy.(x,y)
        PyPlot.plot(x,slbet,label="betamix at prot. conc. $(round(y[1],digits=4)) ")
        
        PyPlot.plot(statpts[i],fwd1ptxy.(statpts[i],protcon[i]),"or")
        for z in 1:length(statpts[i])
            PyPlot.text(statpts[i][z],fwd1ptxy.(statpts[i][z],protcon[i]),
                        string(z),color="red")
        end
        
        PyPlot.plot(inflpts[i],fwd1ptxy.(inflpts[i],protcon[i]),"ob")
        for z in 1:length(inflpts[i])
            PyPlot.text(inflpts[i][z],fwd1ptxy.(inflpts[i][z],protcon[i]),
                        string(z),color="blue")
        end
    end
    PyPlot.legend()
    PyPlot.xlabel("SDS concentration [mM]")
    PyPlot.ylabel("Enthalpy [kJ/mol]")
    PyPlot.tight_layout()
    
    try
        mkdir(outdir)
    catch
        nothing
    end
    PyPlot.savefig(joinpath(outdir,protein*"_foundfeatures_betamix.pdf"))
    return
end

#####################################################

"""
$(TYPEDSIGNATURES)

Plot each single experiment, i.e., enthalpy for an initial protein concentration and increasing SDS concentration.
If the third argument `betamix` is passed, shows a comparison of measured and calculate data (from results of inversion) including each single Beta component.
"""
function plotsingleexperiments(outdir,dobs,betamix=nothing) #; expernumber=nothing)

    nexper = length(dobs.idxdata)
    obsdata = dobs.enthalpy
    protein = dobs.protein

    for i in 1:nexper

        idx = dobs.idxdata[i]
        xy = dobs.sdsprotcon[idx,:]
        obsexper = dobs.enthalpy[idx]

        PyPlot.figure(figsize=(12,6.4))
        PyPlot.title("Protein: $protein, initial concentration: $(xy[1,2])")
        PyPlot.plot(xy[:,1],obsexper,".-k",linewidth=2,label="measured data")

        if betamix!=nothing
            dcalcall = forwmod2D(betamix.betpar,xy,betamix.modkonamp)
            PyPlot.plot(xy[:,1],dcalcall,".-r",linewidth=2,label="calculated data")
            ncomp = size(betamix.modkonamp,2)
            for c=1:ncomp
                dcalc1 = singlescaledbeta2D(betamix.betpar,xy,betamix.modkonamp[:,c])
                PyPlot.plot(xy[:,1],dcalc1,".--",linewidth=0.7,markersize=0.75,label="Beta comp. $c")
            end
        end
        PyPlot.legend()
        PyPlot.xlabel("SDS concentration [mM]")
        PyPlot.ylabel("Enthalpy [kJ/mol")
        PyPlot.tight_layout()

        try
            mkdir(outdir)
        catch
            nothing
        end
        if betamix!=nothing
            PyPlot.savefig(joinpath(outdir,protein*"_experiment$(i)_fit.pdf"))
        else
            PyPlot.savefig(joinpath(outdir,protein*"_experiment$(i).pdf"))
        end
    end

    return
end

#####################################################
"""
$(TYPEDSIGNATURES)

"""
function plotareavsprotcon(proteinname,protcons,areas,linfitres,resstdev,outdir )

    ncomp = size(areas,2)
    PyPlot.figure()
 
    PyPlot.println("### Protein $proteinname, area vs. protein conc. ###")
    ##------------------------
    # linear fits
    ncomp = size(linfitres,1)
    x1 = protcons[1]
    x2 = protcons[end]
    xm = [x1,x2]
    for i=1:ncomp
        acoe   = linfitres[i,1]
        interc = linfitres[i,2]
        println("Linear fit, Beta comp. #$i, ang.coef.=$acoe,\n                           intecept=$interc")
        ym = acoe.*xm .+ interc
        PyPlot.plot(xm,ym,"-k",linewidth=1.5)
    end

    ##------------------------------
    # areas 
    for i=1:ncomp
        rsd = round(resstdev[i],sigdigits=3)
        PyPlot.plot(protcons,areas[:,i],"o",
                    markersize=8,linewidth=2.5,label="Beta comp. #$i, res. std $rsd")
    end
    PyPlot.legend()
    PyPlot.xlabel("Protein concentration [mM]")
    PyPlot.ylabel("Area [kJ]")
    PyPlot.tight_layout()

    try
        mkdir(outdir)
    catch
        nothing
    end
    PyPlot.savefig(joinpath(outdir,proteinname*"_areavsprotcon.pdf"))
    return
end

#####################################################

"""
$(TYPEDSIGNATURES)

Plot the components and the sum of Beta functions for a given protein concentration.
"""
function plotbetacomp1D(betamix,protcon) #; expernumber=nothing)

    protein = betamix.protein

    # function fwd1ptxy(xpt::Real,ypt::Real)
    #     xy = hcat(xpt,ypt)
    #     out = forwmod2D(betamix.betpar,xy,betamix.modkonamp)
    #     return out[1] 
    # end

    PyPlot.figure(figsize=(12,6.4))
    PyPlot.title("Protein: $protein, concentration: $protcon")

    N = 100
    x = collect(LinRange(betamix.betpar.a,betamix.betpar.b,N))
    y = protcon.*ones(N)
    xy = [x y]

    ncomp = size(betamix.modkonamp,2)
    for c=1:ncomp
        #dcalc1 = singlescaledbeta2D(betamix.betpar,xy,betamix.modkonamp[:,c])
        xy = [x y]
        enth = forwmod2D(betamix.betpar,xy,betamix.modkonamp[:,c:c])
        bpars = round.(betamix.modkonamp[:,c:c],digits=2)
        PyPlot.plot(x,enth,"--",linewidth=1,label="Beta comp. $c, $(bpars)")
        #@show c,betamix.modkonamp[:,c:c]
    end

    enth = forwmod2D(betamix.betpar,xy,betamix.modkonamp)
    PyPlot.plot(x,enth,"-",linewidth=2,color="red",label="sum of Beta func.")

    PyPlot.legend()
    PyPlot.xlabel("SDS concentration [mM]")
    PyPlot.ylabel("Enthalpy [kJ/mol")
    PyPlot.tight_layout()

    return
end

#############################################################################3
