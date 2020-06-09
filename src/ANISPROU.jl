#
# This file is a part of ANISPROU. License is MIT
# Copyright (c) 2020 Andrea Zunino
#


"""
ANISPROU

Analysis of isothermal titration calorimetry (ITC) data on sodium dodecyl sulphate (SDS) mediated protein unfolding.

# Exports

$(EXPORTS)
"""
module ANISPROU

using DelimitedFiles
using Distributions 
using ForwardDiff
#using ReverseDiff
using HDF5
using LinearAlgebra
using Roots
using Cubature
using Optim
using PyPlot 
using WriteVTK
using DocStringExtensions

export ScaledBeta2DParams,BetaMix2D,ITCObsData
export readallexperiments
export setconstraints,solveinvprob
export plotinitialguess,plotresults
export saveresultVTK
export findcurvefeatures
export calcfreeSDSNbound
export freeboundsds_singlebeta,freeboundsds_betamix
export area_betamix,volume_betamix
export plotsingleexperiments,plotbindisotherm_singlebetas
export plotbindisotherm_betamix,plotfoundfeatures


include("readITCdata.jl")
include("invertITCdata2D.jl")
include("betaintegrals.jl")
include("bindingisotherm.jl")
include("plotstuff.jl")

end # module
