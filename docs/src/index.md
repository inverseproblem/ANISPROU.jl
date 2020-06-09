
```@meta
Author = "Andrea Zunino"
```

# ANISPROU

```@contents
Pages = ["index.md"]
```

# User guide


## Installation

To install the package simple enter into the package manager mode in Julia by typing "`]`" at the REPL prompt and then use `add`, i.e.,
```
(v1.4) pkg> add ANISPROU
```
The package will be automatically downloaded from the web and installed.



## Usage examples


# Public API

# General
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
area_betamix
volume_betamix
```
## Plotting
```@docs
plotresults
plotfoundfeatures
plotbindisotherm_betamix
plotbindisotherm_singlebetas
```

