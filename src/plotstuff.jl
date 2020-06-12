

################################################################
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

    figure(figsize=(12,5))
    #suptitle("Observed data and initial guess")

    subplot(121)
    title("$protein observed data")
    #tricontour(sdscon,procon,dobs,cmap=PyPlot.get_cmap("rainbow"))
    scatter(sdscon,procon,c=dobs.enthalpy,cmap=PyPlot.get_cmap("rainbow"))
    colorbar(label="Enthalpy [kJ/mol]")
    xlabel("SDS concentration [mM]")
    ylabel("Protein concentration [mM]")    
    plotmodelines(betpar,mstart,"start. m.")

    subplot(122)
    title("Calculated data from starting model")
    #tricontour(sdscon,procon,dcalcstart,cmap=PyPlot.get_cmap("rainbow"))
    scatter(sdscon,procon,c=dcalcstart,cmap=PyPlot.get_cmap("rainbow"))
    colorbar(label="Enthalpy [kJ/mol]")
    xlabel("SDS concentration [mM]")
    ylabel("Protein concentration [mM]")
    plotmodelines(betpar,mstart,"start. m.")
   
    tight_layout()
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


    figure(figsize=(14,7))
    subplot(231)
    title("$protein observed data")
    #tricontour(sdscon,procon,dobs,cmap=PyPlot.get_cmap("rainbow"))
    scatter(sdscon,procon,c=dobs.enthalpy,cmap=PyPlot.get_cmap("rainbow"))
    colorbar(label="Enthalpy [kJ/mol]")
    xlabel("SDS concentration [mM]")
    ylabel("Protein concentration [mM]")
    # xlim([-2,betpar.b])
    # ylim([0.9*betpar.ymin,1.05*betpar.ymax])
    #plotmodelines(betpar,mpost,"mpost")
    
    subplot(232)
    title("Calculated data from starting model")
    #tricontour(sdscon,procon,dcalcstart,cmap=PyPlot.get_cmap("rainbow"))
    scatter(sdscon,procon,c=dcalcstart,cmap=PyPlot.get_cmap("rainbow"))
    colorbar(label="Enthalpy [kJ/mol]")
    xlabel("SDS concentration [mM]")
    ylabel("Protein concentration [mM]")
    plotmodelines(betpar,mstart,"start. m.")

    subplot(233)
    title("Model parameters")
    plot(vec(mstart),".-",label="starting model")
    plot(vec(mpost),".-",label="posterior model")
    legend()

    subplot(234)
    title("Calculated data from posterior model")
    #tricontour(sdscon,procon,dcalccur,cmap=PyPlot.get_cmap("rainbow"))
    scatter(sdscon,procon,c=dcalccur,cmap=PyPlot.get_cmap("rainbow"))
    colorbar(label="Enthalpy [kJ/mol]")
    xlabel("SDS concentration [mM]")
    ylabel("Protein concentration [mM]")
    plotmodelines(betpar,mpost,"post. m.")    
    
    
    subplot(235)
    title("Calculated minus observed data")
    #tricontour(sdscon,procon,dcalccur,cmap=PyPlot.get_cmap("rainbow"))
    dcmo = dcalccur-dobs.enthalpy
    vmax = maximum(abs.(dcmo))
    scatter(sdscon,procon,c=dcmo,vmin=-vmax,vmax=vmax,cmap=PyPlot.get_cmap("RdBu"))
    colorbar(label="ΔEnthalpy [kJ/mol]")
    xlabel("SDS concentration [mM]")
    ylabel("Protein concentration [mM]")
    # xlim([-2,betpar.b])
    # ylim([0.9*betpar.ymin,1.05*betpar.ymax])
    
    subplot(236)
    title("Calculated vs. observed data")
    plot(dcalcstart,"--",label="calc. data from start. mod.",linewidth=0.8)
    plot(dobs.enthalpy,"-k",label="observed data")
    plot(dcalccur,"-r",label="calc. data from post. mod.")
    legend()

    # subplot(235)
    # title("Misfit")
    # plot(misfall,"o-")
    # xlabel("Iteration number")
    # ylabel("misfit")

    tight_layout()
    
    try
        mkdir(outdir)
    catch
        nothing
    end
    savefig(joinpath(outdir,protein*"_results.pdf"))
    return nothing
end

############################################

"""
$(TYPEDSIGNATURES)

Plot the lines defined by the modes.
"""
function plotmodelines(betpar,mcur,modname)
    Npoints = 40
    ncomp = size(mcur,2)
    ym = LinRange(betpar.ymin,betpar.ymax,Npoints)
    for i=1:ncomp
        # s1 = (i-1)*betpar.nummodpar + 1 
        # s2 = (i-1)*betpar.nummodpar + 2 
        x1 = betpar.ymin
        x2 = betpar.ymax
        xm = (mcur[2,i].-mcur[1,i])./(x2.-x1).*(ym.-x1) .+ mcur[1,i] 
        if i==ncomp
            plot(xm,ym,"-k",linewidth=0.7,label="mode $modname")
        else
            plot(xm,ym,"-k",linewidth=0.7)
        end
    end
    legend()
    # xlim([-2,betpar.b])
    # ylim([0.9*betpar.ymin,1.05*betpar.ymax])
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
function plotbindisotherm(betamix,protcon,dobs,statpts,inflpts,
                          freeSDS,Nbound,outdir)

    protein = betamix.protein
    N = 1000
    x = collect(LinRange(betamix.betpar.a,betamix.betpar.b,N))

    ##-------------------------------------
    figure(figsize=(12,10))
    
    subplot(211)
    title("Protein $protein: binding isotherm from stationary points of Beta mix")
    scatter(dobs.sdsprotcon[:,1],dobs.sdsprotcon[:,2],c=dobs.enthalpy,cmap=PyPlot.get_cmap("rainbow"))
    colorbar(label="Enthalpy [kJ/mol]")
    xlabel("SDS concentration [mM]")
    ylabel("Protein concentration [mM]")    

    ny = length(protcon)
    for i in 1:ny
        for z in 1:length(statpts[i])
            plot(statpts[i][z],protcon[i],"or")
            text(statpts[i][z],protcon[i],string(z),color="red")
        end
        
        for z in 1:length(inflpts[i])
            plot(inflpts[i][z],protcon[i],"ob")
            text(inflpts[i][z],protcon[i],string(z),color="blue")
        end
    end

    nlines = length(freeSDS)
    for i=1:nlines
        ## totSDS = freeSDS + Nbound * protc
        # y = angcoe.*x .+ intercept
        # x = (-intercept .+ 0)./angcoe
        ypl = (1.0/Nbound[i]).*x .+ (-freeSDS[i]/Nbound[i])
        plot(x,ypl,"-k",linewidth=0.5)
    end

    # xlim([betamix.betpar.a,betamix.betpar.b])
    ylim([0.92*betamix.betpar.ymin,1.05betamix.betpar.ymax])

    ##
    ## totSDS = freeSDS + Nbound * protc
    ##
    subplot(212)
    title("Binding isotherm for protein $protein")
    plot(freeSDS,Nbound,"o-")
    xlabel("free SDS concentration [mM]")
    ylabel("Nbound")

    tight_layout()
    savefig(joinpath(outdir,protein*"_bind-isoth_betamix.pdf"))
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
        out = ANISPROU.forwmod2D(betamix.betpar,xy,betamix.modkonamp)
        return out[1] 
    end

    figure(figsize=(12,6))
    title("Protein $protein: stationary and inflection points")
    for i in 1:ny
        y = [protcon[i] for _ in 1:N]
        slbet = fwd1ptxy.(x,y)
        plot(x,slbet,label="betamix at prot. conc. $(round(y[1],digits=4)) ")
        
        plot(statpts[i],fwd1ptxy.(statpts[i],protcon[i]),"or")
        for z in 1:length(statpts[i])
            text(statpts[i][z],fwd1ptxy.(statpts[i][z],protcon[i]),
                 string(z),color="red")
        end
        
        plot(inflpts[i],fwd1ptxy.(inflpts[i],protcon[i]),"ob")
        for z in 1:length(inflpts[i])
            text(inflpts[i][z],fwd1ptxy.(inflpts[i][z],protcon[i]),
                 string(z),color="blue")
        end
    end
    legend()
    xlabel("SDS concentration [mM]")
    ylabel("Enthalpy [kJ/mol]")
    tight_layout()
    
    savefig(joinpath(outdir,protein*"_foundfeatures_betamix.pdf"))
    return
end

#####################################################

"""
$(TYPEDSIGNATURES)

Plot each single experiments, i.e., enthalpy for an initial protein concentration and increasing SDS concentration, comparing measured and calculate data (from results of inversion).
"""
function plotsingleexperiments(dobs,betamix,outdir) #; expernumber=nothing)

    nexper = length(dobs.idxdata)
    obsdata = dobs.enthalpy
    protein = betamix.protein
    ncomp = size(betamix.modkonamp,2)

    for i in 1:nexper

        idx = dobs.idxdata[i]
        xy = dobs.sdsprotcon[idx,:]
       
        dcalcall = forwmod2D(betamix.betpar,xy,betamix.modkonamp)
        obsexper = dobs.enthalpy[idx]

        figure(figsize=(12,6.4))
        title("Protein: $protein, initial concentration: $(xy[1,2])")
        plot(xy[:,1],obsexper,".-r",linewidth=2,label="measured data")
        plot(xy[:,1],dcalcall,".-k",linewidth=2,label="calculated data")
        for c=1:ncomp
            dcalc1 = singlescaledbeta2D(betamix.betpar,xy,betamix.modkonamp[:,c])
            plot(xy[:,1],dcalc1,".--",linewidth=0.5,markersize=0.5,label="Beta comp. $c")
        end
        legend()
        xlabel("SDS concentration [mM]")
        ylabel("Enthalpy [kJ/mol")
        tight_layout()

        savefig(joinpath(outdir,protein*"_experiment$(i).pdf"))
    end

    return
end

#####################################################
