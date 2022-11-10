module SIR_partial_sweep

using Chain
using LinearAlgebra
using OrdinaryDiffEq
using Parameters
using StatsBase

import Base: vec, length, size, getindex

include("objects.jl")
export SIRParameters, SIRRegion, SIRState, regions, parameters, cross_immunity

include("dynamics.jl")
export SIRSolution, simulate

include("tools.jl")
export set_infected!, frequency

include("analytics.jl")
export equilibrium

end
