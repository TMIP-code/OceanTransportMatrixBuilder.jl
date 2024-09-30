module OceanTransportMatrixBuilder

# stdlibs
using SparseArrays
using LinearAlgebra

# deps
using ProgressMeter: @showprogress
using Distances: haversine

include("preprocessing.jl")
include("matrixbuilding.jl")
include("extratools.jl") # <- I think this should be in a separate "base" repo

export facefluxesfrommasstransport
export facefluxesfromvelocities
export makemodelgrid
export makeindices
export transportmatrix

end
