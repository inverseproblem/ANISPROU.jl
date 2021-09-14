[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://inverseproblem.github.io/ANISPROU.jl/dev)
![CI](https://github.com/inverseproblem/ANISPROU.jl/workflows/CI/badge.svg)

# ANISPROU
Analysis of isothermal titration calorimetry (ITC) data on sodium dodecyl sulphate (SDS) mediated protein unfolding (ANISPROU) is a tool developed to globally fit an entire dataset and to extract thermodynamic values from this fit. ITC data on SDS mediated protein unfolding, at different protein concentrations, is used as an input for the fitting. The linearity of the features in the ITC data as a function of protein concentration allows the data to be fitted using a number of 3D beta functions, each representing a thermodynamic event. Besides the enthalpy of unfolding, the binding isotherm is also among the outputs. 

Documentation: https://inverseproblem.github.io/ANISPROU.jl/dev

![surf3d](surf3d1_paper.png)

If you use ANISPROU for a scientific publication or else, please cite our paper:
 
Tidemand F., Zunino A., Johansen N., Hansen A. F., Westh P., Mosegaard K., Arleth L. (2021), **A semi-empirical analysis of complex ITC data from protein-surfactant interactions**, *Analytical Chemistry*, [DOI:10.1021/acs.analchem.1c02558](http://dx.doi.org/10.1021/acs.analchem.1c02558)
