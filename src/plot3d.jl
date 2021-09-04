

"""
$(TYPEDSIGNATURES)

Plot the a 3D surface from the Beta mix together measured data as circles.
"""
function plotsurface3D(dobs,betamix ; yscal=1e2, zscal=2.0, markersize=3500,
                       displayfig=true, ymin=nothing,ymax=nothing,
                       savefig=false,outdir="")

    # markersize=350Makie.px

    println("\nPlotting 3D surface from Beta mix and measured data as circles.")
    println(" Scaling factors are $yscal for [$(dobs.protein)] and $zscal for enthalpy.")
    
    N = 300
    xmin = betamix.betpar.a
    xmax = betamix.betpar.b
    if ymin==nothing
        ymin = betamix.betpar.ymin
    end
    if ymax==nothing
        ymax = betamix.betpar.ymax
    end

    xv = collect(LinRange(xmin,xmax,N))
    yv = collect(LinRange(ymin,ymax,N))

    nx = length(xv)
    ny = length(yv)
    dcalc = zeros(nx,ny)
    for j=1:ny
        xy = [xv yv[j]*ones(nx)]
        dcalc[:,j] = forwmod2D(betamix.betpar,xy,betamix.modkonamp)
    end
    zmin = minimum(dcalc)
    zmax = maximum(dcalc)
       

    xobs = dobs.sdsprotcon[:,1]
    yobs = dobs.sdsprotcon[:,2]
    zobs = dobs.enthalpy


    ## Makie plot    
    fig = Makie.Figure(resolution = (1000, 800))
    ax1 = Makie.Axis3(fig[1,1]) # Axis3 -> 3D

    Makie.surface!(ax1,xv,yscal*yv,zscal*dcalc,alpha=1.0,colormap=:rainbow1)

    Makie.scatter!(ax1,xobs,yscal*yobs,zscal*zobs,color=:black,markersize=markersize)

    ax1.xlabel = "[SDS]"
    ax1.ylabel = "[$(dobs.protein)]x$yscal"
    ax1.zlabel = "Enthalpy x $zscal"

    if savefig
        flname = joinpath(outdir,"surfaceplot.png")
        GLMakie.save(flname,fig)
    end

    Makie.center!(fig.scene)
    if displayfig
        display(fig)
    end
    
    return fig
end

