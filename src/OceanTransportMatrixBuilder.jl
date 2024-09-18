module OceanTransportMatrixBuilder

# stdlibs
using SparseArrays

# deps
using ProgressMeter: @showprogress
using Distances: haversine

include("preprocessing.jl")
include("matrixbuilding.jl")

export makeualldirections
export makemodelgrid
export makeindices
export transportmatrix

end
