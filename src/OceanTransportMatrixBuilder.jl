module OceanTransportMatrixBuilder

# stdlibs
using SparseArrays
using LinearAlgebra

# deps
using Distances: haversine

include("preprocessing.jl")
include("matrixbuilding.jl")
include("grids.jl")
include("topology.jl")
include("extratools.jl") # <- I think this should be in a separate "base" repo

export velocity2fluxes
export facefluxesfrommasstransport
export facefluxesfromvelocities
export makemodelgrid
export makeindices
export transportmatrix

end
