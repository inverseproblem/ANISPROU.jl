
#using Makie,AbstractPlotting

"""
$(TYPEDSIGNATURES)

Plot the a 3D surface from the Beta mix together measured data as circles.
"""
function plotsurface3D(dobs,betamix ; yscal=1e2, zscal=2.0, markersize=350Makie.px,
                       displayscene=true, ymin=nothing,ymax=nothing)

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

    # yscal = 110
    # zscal = 0.01

    ## Makie plot    
    scene = Makie.Scene(resolution = (1000, 1000))
    Makie.surface!(scene,xv,yscal*yv,zscal*dcalc,alpha=1.0)
    
    # lim = FRect3D((xmin,ymin,zmin),(xmax,ymax,zmax))
    # Makie.surface!(scene,xv,yv,dcalc,limits=lim)

    # ,axis=(showgrid=false,)) #,colormap=Reverse(:amp))

    # julia> scatter(rand(10), axis = (showgrid = false,))
    # julia> scatter(rand(10), axis = (showgrid = (false,true),))
    # julia> scatter(rand(10), axis = (showticks = (false,true),))
    # julia> scatter(rand(10), axis = (showticks = false,))

    #Makie.wireframe!(scene,xv,yscal*yv,zscal*dcalc,overdraw = true, transparency = true, color = (:black, 0.1))
    Makie.scatter!(scene,xobs,yscal*yobs,zscal*zobs,color=:black,markersize=markersize)

    axis = scene[Makie.Axis] # get axis
    axis[:names, :axisnames] = ("[SDS]","[$(dobs.protein)]x$yscal","Enthalpy x $zscal")

    #Makie.scale!(scene,1,yscal,zscal)
    Makie.center!(scene)
    if displayscene
        Makie.display(scene)
    end
    #Makie.save("surfaceplot.png", scene)

    return scene
end

