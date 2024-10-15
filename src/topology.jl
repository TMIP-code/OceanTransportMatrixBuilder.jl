abstract type AbstractGrid end
struct BipolarGrid <: AbstractGrid
    nx::Int64
    ny::Int64
    nz::Int64
end
struct TripolarGrid <: AbstractGrid
    nx::Int64
    ny::Int64
    nz::Int64
end
struct UnknownGrid <: AbstractGrid
    nx::Int64
    ny::Int64
    nz::Int64
end

"""
    gridtopology(lon_vertices, lat_vertices, lev)

Returns the type of grid based on how the vertices connect at the north pole.
"""
function gridtopology(lon_vertices, lat_vertices, lev)
    nx = size(lon_vertices, 2)
    ny = size(lon_vertices, 3)
    nz = length(lev)
    # North pole vertices
    NPlon = @view lon_vertices[3:4,:,end]
    NPlat = @view lat_vertices[3:4,:,end]
    # If all the latitudes of the "northmost" vertices are 90, then it's a "regular" grid with 2 poles
    if all(NPlat .== 90)
        return BipolarGrid(nx,ny,nz)
    # Otherwise check if the north pole is split in two
    elseif NPlon == rot180(NPlon) && NPlat == rot180(NPlat)
        return TripolarGrid(nx,ny,nz)
    else
        return UnknownGrid(nx,ny,nz)
    end
end

# Default behavior for grids (assumed bipolar)
# wrap around in the i direction (longitude)
i₊₁(C, g::AbstractGrid) = C.I[1] < g.nx ? C + CartesianIndex(1, 0, 0) : CartesianIndex(1, C.I[2], C.I[3])
i₋₁(C, g::AbstractGrid) = C.I[1] > 1 ? C + CartesianIndex(-1, 0, 0) : CartesianIndex(g.nx, C.I[2], C.I[3])
i₊₁(C::CartesianIndex{2}, g::AbstractGrid) = C.I[1] < g.nx ? C + CartesianIndex(1, 0) : CartesianIndex(1, C.I[2])
i₋₁(C::CartesianIndex{2}, g::AbstractGrid) = C.I[1] > 1 ? C + CartesianIndex(-1, 0) : CartesianIndex(g.nx, C.I[2])
# No connection by default in the j direction (latitude)
j₊₁(C, g::AbstractGrid) = C.I[2] < g.ny ? C + CartesianIndex(0, 1, 0) : nothing
j₋₁(C, ::AbstractGrid) = C.I[2] > 1 ? C + CartesianIndex(0, -1, 0) : nothing
j₊₁(C::CartesianIndex{2}, g::AbstractGrid) = C.I[2] < g.ny ? C + CartesianIndex(0, 1) : nothing
j₋₁(C::CartesianIndex{2}, ::AbstractGrid) = C.I[2] > 1 ? C + CartesianIndex(0, -1) : nothing
# No connection by default in the k direction (depth)
k₊₁(C, g::AbstractGrid) = C.I[3] < g.nz ? C + CartesianIndex(0, 0, 1) : nothing
k₋₁(C, ::AbstractGrid) = C.I[3] > 1 ? C + CartesianIndex(0, 0, -1) : nothing



function ishift(C, g::AbstractGrid, n=0)
    i, j, k = C.I
    CartesianIndex(mod1(i + n, g.nx), j, k)
end
function jshift(C, g::AbstractGrid, n=0)
    i, j, k = C.I
    j2 = j + n
    (1 ≤ j2 ≤ g.ny) ? CartesianIndex(i, j2, k) : nothing
end
function kshift(C, g::AbstractGrid, n=0)
    i, j, k = C.I
    k2 = k + n
    (1 ≤ k2 ≤ g.nz) ? CartesianIndex(i, j, k2) : nothing
end

# Special behavior for tripolar grids where the seam connects the top and "folds" it around
# ┌─────┬─────┬─────┬─────┬─────┬─────┐
# │  6  │  5  │  4  │  3  │  2  │  1  │
# ├─────┼─────┼─────┼─────┼─────┼─────┤ ◄─ seam
# │  1  │  2  │  3  │  4  │  5  │  6  │
# └─────┴─────┴─────┴─────┴─────┴─────┘
#            ─► increasing i ("east")
j₊₁(C, g::TripolarGrid) = C.I[2] < g.ny ? C + CartesianIndex(0, 1, 0) : CartesianIndex(g.nx - C.I[1] + 1, g.ny, C.I[3])
j₊₁(C::CartesianIndex{2}, g::TripolarGrid) = C.I[2] < g.ny ? C + CartesianIndex(0, 1) : CartesianIndex(g.nx - C.I[1] + 1, g.ny)

function jshift(C, g::TripolarGrid, n=0)
    i, j, k = C.I
    j2 = j + n
    if !(1 ≤ j2)
        return nothing
    elseif j2 ≤ g.ny
        return CartesianIndex(i, j2, k)
    else
        i2 = g.nx - i + 1
        return CartesianIndex(i2, j, k)
    end
end

# error for unknown grids
i₊₁(C, ::UnknownGrid) = error("Unknown grid type")
i₋₁(C, ::UnknownGrid) = error("Unknown grid type")
j₊₁(C, ::UnknownGrid) = error("Unknown grid type")
j₋₁(C, ::UnknownGrid) = error("Unknown grid type")
k₊₁(C, ::UnknownGrid) = error("Unknown grid type")
k₋₁(C, ::UnknownGrid) = error("Unknown grid type")


