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
using Cubature,QuadGK
using Optim
using PyPlot 
using WriteVTK
using DocStringExtensions
using Requires

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

include("readITCdata.jl")
include("invertITCdata2D.jl")
include("betaintegrals.jl")
include("bindingisotherm.jl")
include("plotstuff.jl")

# 3D plotting
function __init__()
    @require Makie="ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"  begin
        # @require AbstractPlotting="537997a7-5e4e-5d89-9595-2241ea00577e" begin
        import .Makie
        include("plot3d.jl")
        export plotsurface3D
       # end
    end
end



end # module
