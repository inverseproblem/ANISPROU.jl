
################################################################

function plotinitialguess(protein,betpar,sdscon,procon,dobs,mstart)
    
    dcalcstart = forwmod2D(betpar,mstart)


    lenm = length(mstart)
    nummodparam = betpar.nummodpar
    ncomp = div(lenm,nummodparam)
    ym = LinRange(betpar.ymin,betpar.ymax,40)

    figure(figsize=(12,5))
    suptitle("Observed data and initial guess")

    subplot(121)
    title("$protein observed data")
    #tricontour(sdscon,procon,dobs,cmap=PyPlot.get_cmap("rainbow"))
    scatter(sdscon,procon,c=dobs,cmap=PyPlot.get_cmap("rainbow"))
    colorbar()
    xlabel("SDS concentration [mM]")
    ylabel("Protein concentration [mM]")    
    plotmodelines(ncomp,betpar,mstart,"mstart")

    subplot(122)
    title("dcalc start")
    #tricontour(sdscon,procon,dcalcstart,cmap=PyPlot.get_cmap("rainbow"))
    scatter(sdscon,procon,c=dcalcstart,cmap=PyPlot.get_cmap("rainbow"))
    colorbar()
    xlabel("SDS concentration [mM]")
    ylabel("Protein concentration [mM]")
    plotmodelines(ncomp,betpar,mstart,"mstart")
   
    #tight_layout()
    return
end

############################################

function plotresults(protein,betpar,sdscon,procon,dobs,mstart,mpost,outdir)

    dcalcstart = forwmod2D(betpar,mstart)
    dcalccur = forwmod2D(betpar,mpost)

    ##==============================
    ## show data
    lenm = length(mstart)
    nummodparam = betpar.nummodpar
    ncomp = div(lenm,nummodparam)
    ym = LinRange(betpar.ymin,betpar.ymax,40)


    figure(figsize=(14,7))
    subplot(231)
    title("$protein observed data")
    #tricontour(sdscon,procon,dobs,cmap=PyPlot.get_cmap("rainbow"))
    scatter(sdscon,procon,c=dobs,cmap=PyPlot.get_cmap("rainbow"))
    colorbar()
    xlabel("SDS concentration [mM]")
    ylabel("Protein concentration [mM]")
    xlim([-2,betpar.b])
    ylim([0.9*betpar.ymin,1.05*betpar.ymax])
    #plotmodelines(ncomp,betpar,mpost,"mpost")
    
    subplot(232)
    title("dcalc start")
    #tricontour(sdscon,procon,dcalcstart,cmap=PyPlot.get_cmap("rainbow"))
    scatter(sdscon,procon,c=dcalcstart,cmap=PyPlot.get_cmap("rainbow"))
    colorbar()
    xlabel("SDS concentration [mM]")
    ylabel("Protein concentration [mM]")
    plotmodelines(ncomp,betpar,mstart,"mstart")

    subplot(233)
    title("model parameters")
    plot(mstart,".-",label="mstart")
    plot(mpost,".-",label="mpost")
    legend()

    subplot(234)
    title("dcalc mpost")
    #tricontour(sdscon,procon,dcalccur,cmap=PyPlot.get_cmap("rainbow"))
    scatter(sdscon,procon,c=dcalccur,cmap=PyPlot.get_cmap("rainbow"))
    colorbar()
    xlabel("SDS concentration [mM]")
    ylabel("Protein concentration [mM]")
    plotmodelines(ncomp,betpar,mpost,"mpost")    
    
    
    subplot(235)
    title("dcalc-dobs")
    #tricontour(sdscon,procon,dcalccur,cmap=PyPlot.get_cmap("rainbow"))
    dcmo = dcalccur-dobs
    vmax = maximum(abs.(dcmo))
    scatter(sdscon,procon,c=dcmo,vmin=-vmax,vmax=vmax,cmap=PyPlot.get_cmap("RdBu"))
    colorbar()
    xlabel("SDS concentration [mM]")
    ylabel("Protein concentration [mM]")
    xlim([-2,betpar.b])
    ylim([0.9*betpar.ymin,1.05*betpar.ymax])
    
    subplot(236)
    title("dcalc vs. dobs")
    plot(dcalcstart,"--",label="dcalc mstart",linewidth=0.8)
    plot(dobs,"-k",label="dobs")
    plot(dcalccur,"-r",label="dcalc mpost")
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

function plotmodelines(ncomp,betpar,mcur,modname)
    Npoints = 40
    ym = LinRange(betpar.ymin,betpar.ymax,Npoints)
    for i=1:ncomp
        s1 = (i-1)*betpar.nummodpar + 1 
        s2 = (i-1)*betpar.nummodpar + 2 
        x1 = betpar.ymin
        x2 = betpar.ymax
        xm = (mcur[s2].-mcur[s1])./(x2.-x1).*(ym.-x1) .+ mcur[s1] 
        if i==ncomp
            plot(xm,ym,"-k",linewidth=0.7,label="mode $modname")
        else
            plot(xm,ym,"-k",linewidth=0.7)
        end
    end
    legend()
    xlim([-2,betpar.b])
    ylim([0.9*betpar.ymin,1.05*betpar.ymax])
end

############################################

function saveresultVTK(  )

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
    betparvtk = ScaledBeta2DParams(nummodpar=nummodpar,xy=xy,a=betpar.a,b=betpar.b,
                                   ymin=betpar.ymin,ymax=betpar.ymax)
    dcalcvtk = forwmod2D(betparvtk,mpost)
    dcvtk2 = reshape(dcalcvtk,Nx,Ny)
    vtkfile["energyout",VTKPointData()] = dcvtk2
    outfiles = vtk_save(vtkfile)

    return 
end

############################################

## Makie plot    
#scene = Scene(resolution = (500, 500))
# scene=Makie.surface(xvtk,100*yvtk,0.01*dcvtk2,colormap = Reverse(:amp))
# wireframe!(scene,xvtk,100*yvtk,0.01*dcvtk2,overdraw = true, transparency = true, color = (:black, 0.1))
# display(scene)
