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
#using HDF5
using JLD2
using LinearAlgebra
using Roots
using Cubature,QuadGK
using Optim
using PyPlot
PyPlot.rc("font", family="serif", size=12)
using WriteVTK
using DocStringExtensions
using Requires
using GLMakie


export ScaledBeta2DParams,BetaMix2D,ITCObsData
export readallexperiments
export solveinvprob
export plotobsdata,plotinitialguess,plotresults
export plotsingleexperiments
export saveresultVTK
export findcurvefeatures
export calcfreeSDSNbound
export area_enthalpy,volume_enthalpy,areasvsprotcon
export plotbindisotherm,plotfoundfeatures
export plotareavsprotcon,plotbetacomp1D
export plotparamlines

include("utils.jl")
include("readITCdata.jl")
include("invertITCdata2D.jl")
include("betaintegrals.jl")
include("bindingisotherm.jl")
include("plotstuff.jl")

# 3D plotting
include("plot3d.jl")

# # 3D plotting
# function __init__()
#     @require GLMakie="e9467ef8-e4e7-5192-8a1a-b1aee30e663a"  begin
#         # @require AbstractPlotting="537997a7-5e4e-5d89-9595-2241ea00577e" begin
#         import .GLMakie
#         include("plot3d.jl")
#         export plotsurface3D
#        # end
#     end
# end



end # module
