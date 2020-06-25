
```@meta
Author = "Andrea Zunino"
```

# ANISPROU

ANalysis of ISothermal titration calorimetry data on sodium dodecyl sulphate mediated PRotein Unfolding.

Analysis of isothermal titration calorimetry (ITC) data on sodium dodecyl sulphate (SDS) mediated protein unfolding (ANISPROU) is a tool developed to globally fit an entire dataset and to extract thermodynamic values from this fit. ITC data on SDS mediated protein unfolding, at different protein concentrations, is used as an input for the fitting. The linearity of the features in the ITC data as a function of protein concentration allows the data to be fitted using a number of 3D beta functions, each representing a thermodynamic event. Besides the enthalpy of unfolding, the binding isotherm is also among the outputs. 


Contents:
```@contents
Pages = ["index.md"]
```


# Installation

To install the package simple enter into the package manager mode in Julia by typing "`]`" at the REPL prompt and then use `add`, i.e.,
```
(v1.4) pkg> add ANISPROU
```
The package will be automatically downloaded from the web and installed.




# Tutorial
In the following a step-by-step tutorial illustrating how to process the ITC data is shown. The complete code is also available in the folder `examples`.

## Read the data

First of all, import the package and read the data. The function [`readallexperiments`](@ref) takes care of reading the measred data from a ASCII/text file produced by some common instruments. 

However, the user can read/obtain the data in any other way, as long as it is possible to create a structure [`ITCObsData`](@ref), which it is what is used in all the calculations (e.g., data fitting).

If using the [`readallexperiments`](@ref) function, the directory containing the data set and the file names must have a certain structure in order for the function reading the data to work. Such structure should be like the following example, where the directory name is the protein name `IM7`:
```
100uMIM7highSDS_2.DAT
125uMIM7highSDS_1.DAT
25uMIM7highSDS_5.DAT
50uMIM7highSDS_4.DAT
75uMIM7highSDS_3.DAT
```
The file names must follow the following rules:
- they must end with ".DAT"
- they must contain the name of protein (e.g., "IM7") in the file name
- the file format must be the same than the provided example (a commonly used instrument).

Two optional parameters can be passed to [`readallexperiments`](@ref), namely `scalfactor` (defaulting to 0.004184 to convert Cal/mol to kJ/mol) which scales the enthalpy values, and `discninitrows` (defaulting to 2) which skips a certain  number of initial rows from the data set because usually initial data are affected by strong instrument noise which could bias the subsequent fitting process.  See [`readallexperiments`](@ref) to adapt it to a different case.

Next step is to extract the concentration of SDS and protein, the measured enthalpy and the indices for each experiments (`idxdata`) in the global data set. Then instantiate the `ITCObsData` structure containing the measured (observed) enthalpy values along with other information. See [`ITCObsData`](@ref).
```@example procITC
using PyPlot # hide
using ANISPROU

inpdir="../../examples/inputdata/" # directory containing input data
protein = "IM7" # protein name
data = readallexperiments(inpdir,[protein])

sdscon = data[protein]["sdscon"]
procon = data[protein]["procon"]
enthalpy = data[protein]["enout"] 
idxdata = data[protein]["idxdata"]

dobs = ITCObsData(protein=protein,enthalpy=enthalpy,idxdata=idxdata,sdsprotcon=[sdscon procon])
```
The `x` axis represents the SDS concentration, while the `y` axis the protein concentration.
To plot the observed data one can use the function [`plotobsdata`](@ref):
```@example procITC
plotobsdata(dobs)
savefig("plotobsdata.svg") # hide
nothing # hide
```
![](plotinitguess.svg)
Another way of visualising the observed data is plotting each single experiment in 1D using [`plotsingleexperiments`](@ref):
```@example procITC
outdir="figs"
plotsingleexperiments(outdir,dobs)
savefig("plotsingexpobs.svg") # hide
nothing # hide
```
![](plotsingexpobs.svg)


## [Setup the inverse problem: fitting the enthalpy data](@id setupip)

Now define the parameters of the 2D Beta functions: the type of function for the mode, concentration parameter and amplitude (e.g., linear), the `x` and `y` limits in terms of minimum and maximum of SDS and protein concentrations. See [`ScaledBeta2DParams`](@ref). Often enlarging by a little bit the domain on the `x` axis helps the Beta functions to better fit the data.
```@example procITC
a = 0.99*minimum(dobs.sdsprotcon[:,1]) # lower bound for beta domain
b = 1.05*maximum(dobs.sdsprotcon[:,1]) # upper bound for beta domain

minprotcon = minimum(dobs.sdsprotcon[:,2]) 
maxprotcon = maximum(dobs.sdsprotcon[:,2]) 

## define the parameters of Beta 2D functions
modefuny = "linear"
konfuny  = "linear"
ampfuny  = "linear"
betpar = ScaledBeta2DParams(modefuny=modefuny,konfuny=konfuny,ampfuny=ampfuny,a=a,b=b,ymin=minprotcon,ymax=maxprotcon)

```

Now we need to define a so-called starting model, i.e., a set of parameters for the Beta functions (mode, concentration and amplitude) which constitutes our initial guess in order to fit the measured enthalpy data.

 The starting model is represented by a 2D array where the number of rows is the number of Beta functions components. Each column contains the set of parameters necessary to define a  single Beta component, namely mode, confidence parameter and amplitude.
 If we set the all the functions for mode, concentration and amplitude to be *linear* (see above), then the parameters of the Beta functions (mode, confidence and amplitude) will vary along `y` (protein concentration) following the equation of a straight line passing through two poins:
```math
 y = \dfrac{(y_2-y_1)}{(x_2-x_1)} (x-x_1) + y_1
``` 
 where ``y_1``, ``y_2``, ``x_1`` and  ``x_2`` are given as explained in the following. ``x_1`` and  ``x_2`` represent the value of the Beta parameters at the minimum and maximum protein concentration ``y_1`` and ``y_2`` specified in the structure [`ScaledBeta2DParams`](@ref), called `ymin` and `ymax`.
 Then the elements of the column vector represent the following:
- elements 1 and 2: value of the mode at the two points where the protein concentration equals `betpar.ymin` `betpar.ymax`, part of the structure [`ScaledBeta2DParams`](@ref)
- elements 3 and 4: value of the confidence parameter at the two points where the protein concentration equals `betpar.ymin` `betpar.ymax`, part of the structure [`ScaledBeta2DParams`](@ref)
- elements 5 and 6: value of the mode at the two points where the protein concentration equals `betpar.ymin` `betpar.ymax`, part of the structure [`ScaledBeta2DParams`](@ref).

In the following we create a starting model use 4 Beta components. To add more (remove) components we can simple add more (remove) columns in mstart.
```@example procITC
# Elements are: 2 for mode, 2 for the confidence parameter and
#   2 for the amplitude parameter
comp1 = [0.6,  1.5,  30.0, 30.0, -2.5,  -5.0 ]
comp2 = [1.7,  4.8,  60.0, 40.0, -1.6,  -4.0 ]
comp3 = [4.0,  9.0,  40.0, 80.0, 0.12, 0.16 ]
comp4 = [5.0,  12.0, 20.0, 50.0, -1.6, -2.0 ]

# mstart is a 2D array where each column represents one component
mstart = [comp1 comp2 comp3 comp4]
	
nothing # hide
```
To visually check the goodness of our first guess we can plot it using
```@example procITC
plotinitialguess(betpar,dobs,mstart)
savefig("plotinitguess.svg") # hide
nothing # hide
```
![](plotinitguess.svg)


## [Solve the inverse problem: find optimal mix of Beta functions](@id invprobsect)

The solution of the inverse problem, that is, finding the set of model parameters which produces a "best" fit to the observed data is based on a constrained Newton method. The Newton method requires the computation of both the gradient of the misfit function with respect to model parameters and the Hessian matrix. Both gradient and Hessian matrix are calculated using automatic differentiation, specifically using the "forward mode" approach provided by the [`ForwardDiff.jl`](https://github.com/JuliaDiff/ForwardDiff.jl) package.

First, we need to set the constraints for the Newton optimization. That is done by specifying the lower and upper bounds for each parameter and and passing it to the function [`setconstraints`](@ref), as shown in the following:
```@example procITC
lowc = [betpar.a, betpar.a, 2.1, 2.1, -20.0, -20.0] # lower constraints [confidence must be >2.0]
upc  = [betpar.b, betpar.b, 500.0, 500.0, 10.0, 10.0] # upper constraints
lowconstr,upconstr = setconstraints(betpar,mstart,lowc,upc)
nothing # hide
```
In order to solve the inverse problem we need a covariance matrix (symmetric positive-definite) representing the uncertainty on the observed data. What is actually required by the software is the inverse of such covariance matrix `invCd`, i.e., if ``\mathbf{C}_D`` is the covariance matrix on the observations, representing the noise on the data, we need to input the code ``\mathbf{C}^{-1}_D``, sometimes called the precision matrix.
```@example procITC
using LinearAlgebra
nobs = length(dobs.enthalpy)
stdobs = 0.2 .* ones(nobs) # standard deviation of the error on measured data
invCd = inv(diagm(stdobs.^2)) # in this case a diagonal precision matrix
nothing # hide
```
Now we can run the Newton optimization algorithm to solve the inverse problem, provided also a starting model as explained in the previous [`section`](@ref setupip). The algorithm is from the package [`Optim.jl`](https://github.com/JuliaNLSolvers/Optim.jl), specifically it is an interior-point primal-dual Newton algorithm [](https://julianlsolvers.github.io/Optim.jl/stable/#algo/ipnewton/).
The optization is launched by the following code, where `outdir` is the directory where the output results will be written. 
```@example procITC
outdir = "output"
betamix = solveinvprob(betpar,dobs,invCd,mstart,lowconstr,upconstr,outdir)
nothing # hide
```
The output `betamix` is a structure of type [`BetaMix2D`](@ref) which holds the optimized parameters and other additional information.

Finally, it is possible to visualize the results as following:
```@example procITC
outdir = "figs"
plotresults(betamix,dobs,mstart,outdir)
savefig(outdir*"/IM7_results.svg") # hide
nothing # hide
```
![](figs/IM7_results.svg)

The function [`plotsingleexperiments`](@ref) provides a way to show the fit of the solution to the single experiments, i.e., enthalpy for different initial protein concentrations. 
```@example procITC
# plot fit to single experiments
outdir="figs"
plotsingleexperiments(outdir,dobs,betamix)
savefig(outdir*"/IM7_experiment1.svg") # hide
nothing # hide
```
![](figs/IM7_experiment1.svg)

### Plot the 3D surface vs. observed data

In case the package `Makie` is installed, it is possible to make a 3D plot showing the surface defined by the Beta mix and, in addition, the set of observed data as circles.
```@example procITC
using Makie
plotsurface3D(dobs,betamix)
scene=plotsurface3D(dobs,betamix,displayscene=false) # hide
Makie.save("figs/surfaceplot.png", scene) # hide
nothing # hide
```
![](figs/surfaceplot.png)


## Estimation of the binding isotherm

The binding isotherm can be estimated by defining a set of points or "features" on the 2D Beta functions such that given features for a certain protein concentration can be related to corresponding features at different protein concentrations. Such features could be, for instance, the peaks of the curves or their inflection points. A set of the same feature for different protein concentrations should produce a trend close to a straight line, which can be used to estimate the binding isotherm. The angular coefficient of such straight line will then represent the binding number ``N_{\rm bound}`` and the intercept the concentration of free SDS ``[SDS]_{\rm free}``. The two are in fact related by the following relation:
```math
 [\mathrm{SDS}]_{\rm total} =  [\mathrm{SDS}]_{\rm free} + N_{\rm bound} [\mathrm{Protein}]
```
where ``[...]`` represents the concentration.

Thus, to estimate the binding isotherm, we start by finding a set of features, in this case stationary and inflection points at given values of protein concentration (in this case 4 values):
```@example procITC
ny = 4 # number of protein concentrations to investigate
protcon = collect(LinRange(betamix.betpar.ymin,betamix.betpar.ymax,ny)) # set of protein concentrations
statpts,inflpts = findcurvefeatures(betamix,protcon) # find features
nothing # hide
```
To plot the found features one can do the following:
```@example procITC
# plot found points/features
outdir = "figs"
plotfoundfeatures(betamix,protcon,statpts,inflpts,outdir)
savefig("plotfeat.svg") # hide
nothing # hide
```
![](plotfeat.svg)

The next step involves selecting a subset of the found points to construct the binding isotherm: this can be done by looking at the previous plot and picking only desired points. The numbers on the plot nearby the found points correspond to the index in the array of found points.
```@example procITC
# ===========================================
# Selection of local minima and maxima
locminmax = [] ##Vector{Array{<:Real,2}}(undef,0)

push!(locminmax, [ statpts[1][2] protcon[1];
                   statpts[2][2] protcon[2];
                   statpts[3][2] protcon[3];
                   statpts[4][2] protcon[4] ] )

push!(locminmax, [ statpts[2][3] protcon[2];
                   statpts[3][3] protcon[3];
                   statpts[4][3] protcon[4] ] )

push!(locminmax, [ statpts[1][4] protcon[1];
                   statpts[2][6] protcon[2];
                   statpts[3][6] protcon[3];
                   statpts[4][6] protcon[4] ] )


# ===========================================
# Selection of inflection points
inflectionpts = [] ##Vector{Array{<:Real,2}}(undef,0)
    
push!(inflectionpts, [ inflpts[1][2] protcon[1];
                       inflpts[2][2] protcon[2];
                       inflpts[3][2] protcon[3];
                       inflpts[4][2] protcon[4] ] )

push!(inflectionpts, [ inflpts[1][3] protcon[1];
                       inflpts[2][3] protcon[2];
                       inflpts[3][3] protcon[3];
                       inflpts[4][3] protcon[4] ] )

push!(inflectionpts, [ inflpts[1][4] protcon[1];
                       inflpts[2][4] protcon[2];
                       inflpts[3][4] protcon[3];
                       inflpts[4][4] protcon[4] ] )

push!(inflectionpts, [ inflpts[1][5] protcon[1];
                       inflpts[2][5] protcon[2];
                       inflpts[3][5] protcon[3];
                       inflpts[4][5] protcon[4] ] )

push!(inflectionpts, [ inflpts[1][6] protcon[1];
                       inflpts[2][6] protcon[2];
                       inflpts[3][6] protcon[3];
                       inflpts[4][8] protcon[4] ] )

push!(inflectionpts, [ inflpts[1][7] protcon[1];
                       inflpts[2][7] protcon[2];
                       inflpts[3][7] protcon[3];
                       inflpts[4][9] protcon[4] ] )

nothing # hide
```

Once we have a set of set of features we can estimate the best fitting straight lines for each set by performing a least squares linear regression. 
```@example procITC
# find the straight line by least squares
lstlines = [locminmax..., inflectionpts...] # list of set of points
freeSDS,Nbound = calcfreeSDSNbound(lstlines) # do the linear regression
nothing # hide
```
Finally, the resulting binding isotherm is plotted by the following:
```@example procITC 
outdir = "figs"
plotbindisotherm(betamix,protcon,dobs,statpts,inflpts,freeSDS,Nbound,outdir)
savefig("plotbindisoth.svg") # hide
nothing # hide
```
![](plotbindisoth.svg)


## Calculating areas and volumes

Compute the area for each Beta component at requested protein concentration:
```@example procITC
protcon = 0.08 # requested protein concentration
area,errarea = area_enthalpy(betamix,protcon) # compute area
println("Area at protein concentration $protcon for each component: \n$area, \nintegration error\n $errarea\n")
nothing # hide
```
It is also possible to calculate the area for each Beta component for a set of different protein concentrations:
```@example procITC
N = 15
protcons = collect(LinRange(0.0,0.14,N)) # set of protein concentrations
areas,erras = areasvsprotcon(betamix,protcons,doplot=true)
savefig("plotsetareas.svg") # hide
nothing # hide
```	
![](plotsetareas.svg)
	
Compute the volume for all Beta components within requested bounds of protein concentration:
```@example procITC
minprotcon = betamix.betpar.ymin  # lower bound for integral
maxprotcon = betamix.betpar.ymax  # upper bound for integral
volume,errvol = volume_enthalpy(betamix,minprotcon,maxprotcon) # compute volume
println("Volume within bounds for each component: \n$volume, \nintegration error\n $errvol\n")
nothing # hide
```


# Public API

## General calculations
```@docs
ANISPROU
readallexperiments
ITCObsData
ScaledBeta2DParams
BetaMix2D
setconstraints
solveinvprob
findcurvefeatures
calcfreeSDSNbound
area_enthalpy
areasvsprotcon
volume_enthalpy
```
## Plotting
```@docs
plotobsdata
plotinitialguess
plotresults
plotsingleexperiments
plotfoundfeatures
plotbindisotherm
saveresultVTK
plotsurface3D
plotareavsprotcon
plotbetacomp1D
```

# Other non exported functions
```@docs
ANISPROU.lssqregr
ANISPROU.plotmodelines
ANISPROU.getmodparbeta
ANISPROU.forwmod2D
ANISPROU.readsingleexperiment
ANISPROU.scaledbeta
ANISPROU.misfitbeta2D
ANISPROU.singlescaledbeta2D 
```

