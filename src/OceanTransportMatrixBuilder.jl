module OceanTransportMatrixBuilder

# stdlibs
using SparseArrays
using LinearAlgebra

# deps
using Distances: haversine
using StaticArrays: FieldVector

# Main functions to be exposed
include("matrixbuilding.jl")

# Grid stuff
include("gridcellgeometry.jl")
include("gridtopology.jl")

# Derivatives
include("classicderivatives.jl")
include("dyads.jl")
include("triads.jl")

# Velocities and fluxes and Redi/GM stuff
include("velocities.jl")
include("RediGM.jl")

# Extras
include("extratools.jl") # <- I think this should be in a separate "base" repo

export velocity2fluxes
export facefluxesfrommasstransport
export facefluxesfromvelocities
export makegridmetrics
export makeindices
export transportmatrix

end
