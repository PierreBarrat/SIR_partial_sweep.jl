module SIR_partial_sweep

using OrdinaryDiffEq
using LinearAlgebra
using Parameters
using StatsBase

import Base: vec, length, size, getindex

include("objects.jl")
export SIRParameters, SIRRegion, SIRState, regions, parameters

include("dynamics.jl")
export SIRSolution, simulate

include("analytics.jl")
export equilibrium

end
