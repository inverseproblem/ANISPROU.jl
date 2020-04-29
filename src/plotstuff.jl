
################################################################

function plotinitialguess(protein,betpar,dobs,mstart)
    
    dcalcstart = forwmod2D(betpar,mstart)
    sdscon = betpar.xy[:,1]
    procon = betpar.xy[:,2]
    
    #lenm = length(mstart)
    nummodparam = size(mstart,1) #betpar.nummodpar
    #ncomp = size(mstart,2) #div(lenm,nummodparam)
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
    plotmodelines(betpar,mstart,"mstart")

    subplot(122)
    title("dcalc start")
    #tricontour(sdscon,procon,dcalcstart,cmap=PyPlot.get_cmap("rainbow"))
    scatter(sdscon,procon,c=dcalcstart,cmap=PyPlot.get_cmap("rainbow"))
    colorbar()
    xlabel("SDS concentration [mM]")
    ylabel("Protein concentration [mM]")
    plotmodelines(betpar,mstart,"mstart")
   
    #tight_layout()
    return
end

############################################

function plotresults(protein,betpar,dobs,mstart,mpost,outdir)

    dcalcstart = forwmod2D(betpar,mstart)
    dcalccur = forwmod2D(betpar,mpost)
   sdscon = betpar.xy[:,1]
    procon = betpar.xy[:,2]
    
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
    scatter(sdscon,procon,c=dobs,cmap=PyPlot.get_cmap("rainbow"))
    colorbar()
    xlabel("SDS concentration [mM]")
    ylabel("Protein concentration [mM]")
    xlim([-2,betpar.b])
    ylim([0.9*betpar.ymin,1.05*betpar.ymax])
    #plotmodelines(betpar,mpost,"mpost")
    
    subplot(232)
    title("dcalc start")
    #tricontour(sdscon,procon,dcalcstart,cmap=PyPlot.get_cmap("rainbow"))
    scatter(sdscon,procon,c=dcalcstart,cmap=PyPlot.get_cmap("rainbow"))
    colorbar()
    xlabel("SDS concentration [mM]")
    ylabel("Protein concentration [mM]")
    plotmodelines(betpar,mstart,"mstart")

    subplot(233)
    title("model parameters")
    plot(vec(mstart),".-",label="mstart")
    plot(vec(mpost),".-",label="mpost")
    legend()

    subplot(234)
    title("dcalc mpost")
    #tricontour(sdscon,procon,dcalccur,cmap=PyPlot.get_cmap("rainbow"))
    scatter(sdscon,procon,c=dcalccur,cmap=PyPlot.get_cmap("rainbow"))
    colorbar()
    xlabel("SDS concentration [mM]")
    ylabel("Protein concentration [mM]")
    plotmodelines(betpar,mpost,"mpost")    
    
    
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
    xlim([-2,betpar.b])
    ylim([0.9*betpar.ymin,1.05*betpar.ymax])
end

############################################

function plotbindingisotherm(protein,betpar,mpost,sdsfNb,outdir; Npts=100)

    dcalc = forwmod2D(betpar,mpost)    
    sdscon = betpar.xy[:,1]
    procon = betpar.xy[:,2]
    
    figure(figsize=(13,5))
    subplot(121)
    scatter(sdscon,procon,c=dcalc,cmap=PyPlot.get_cmap("rainbow"))
    colorbar()
    xlabel("SDS concentration [mM]")
    ylabel("Protein concentration [mM]")
    ##
    ##  y = m*x + y0
    ##
    ## One line per row
    ##N = 100
    firstplot=true
    for ico=1:length(sdsfNb)
        for r=1:size(sdsfNb[ico],1)
            x = collect(LinRange(betpar.ymin,betpar.ymax,Npts))
            y = sdsfNb[ico][r,2].*x .+ sdsfNb[ico][r,1]
            if firstplot
                plot(y,x,"--r",linewidth=0.7,label="selected fitted amplitudes")
                firstplot = false
            else
                plot(y,x,"--r",linewidth=0.7)
            end
        end
    end

    plotmodelines(betpar,mpost,"")

    subplot(122)
    for ico=1:length(sdsfNb)
        plot(sdsfNb[ico][:,1],sdsfNb[ico][:,2],"o-")
    end
    xlabel("free SDS conc. [mM]")
    ylabel("Nbound")

    tight_layout()

    savefig(joinpath(outdir,protein*"_binding-isotherm.pdf"))
    return nothing
end

######################################################################

function saveresultVTK(betpar,mpost)

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
