

abstract type ArakawaGridCell end

struct AGridCell <: ArakawaGridCell
    u_pos::Symbol
    v_pos::Symbol
end

struct BGridCell <: ArakawaGridCell
    u_pos::Symbol
    v_pos::Symbol
end

struct CGridCell <: ArakawaGridCell
    u_pos::Symbol
    v_pos::Symbol
end


"""
    arakawa = getarakawagrid(u_lon, u_lat, v_lon, v_lat, gridmetrics)

Returns the type of the grid (A, B, or C), the grid position of the velocity points,
and the error of that position relative to the perimeter of the cell.
The positions are relative to the cell center C, and are as follows:

```
                NW───N───NE
increasing j    │         │
      ▲         W    C    E
      │         │         │
                SW───S───SE

               ──► increasing i
```

The different Arakawa grids recongnized here are:

```
    ┌─────────┐    uv───────uv    ┌────v────┐
    │         │    │         │    │         │
    │   Cuv   │    │    C    │    u    C    u
    │         │    │         │    │         │
    └─────────┘    uv───────uv    └────v────┘
      A-grid         B-grid         C-grid
```

where each grid is centered in C.
Also returns the distances from the velocity points to the grid points.
"""
function getarakawagrid(u_lon, u_lat, v_lon, v_lat, gridmetrics)

    # Unpack gridmetrics
    (; lon, lat, lon_vertices, lat_vertices) = gridmetrics

    i = j = 1
    u_point = (u_lon[i, j], u_lat[i, j])
    v_point = (v_lon[i, j], v_lat[i, j])

    C = (lon[i, j], lat[i, j])
    SW = (lon_vertices[1, i, j], lat_vertices[1, i, j])
    SE = (lon_vertices[2, i, j], lat_vertices[2, i, j])
    NE = (lon_vertices[3, i, j], lat_vertices[3, i, j])
    NW = (lon_vertices[4, i, j], lat_vertices[4, i, j])
    S = midpointonsphere(SW, SE)
    N = midpointonsphere(NE, NW)
    W = midpointonsphere(SW, NW)
    E = midpointonsphere(SE, NE)

    cell = (; C, SW, SE, NE, NW, S, N, W, E)

    u_distances = (; (k => haversine(P, u_point) for (k,P) in pairs(cell))...)
    v_distances = (; (k => haversine(P, v_point) for (k,P) in pairs(cell))...)

    u_distance, u_pos = findmin(u_distances)
    v_distance, v_pos = findmin(v_distances)

    # Arakawa grid type
    if u_pos == v_pos == :C
        arakawa = AGridCell(u_pos, v_pos)
    elseif u_pos == v_pos && u_pos ∈ (:NE, :NW, :SE, :SW)
        arakawa = BGridCell(u_pos, v_pos)
    elseif u_pos ∈ (:E, :W) && v_pos ∈ (:N, :S)
        arakawa = CGridCell(u_pos, v_pos)
    else
        error("Unknown Arakawa grid type")
    end

    cellperimeter = haversine(SW, SE) + haversine(SE, NE) + haversine(NE, NW) + haversine(NW, SW)
    relerr = (u_distance + v_distance) / cellperimeter

    relerr > 0.01 && warn("Relative error in grid positions in $arakawa is $relerr")

    return arakawa

end

"""
    u, u_lon, u_lat, v, v_lon, v_lat = interpolateontodefaultCgrid(u, u_lon, u_lat, v, v_lon, v_lat, gridmetrics)

Interpolates the velocity fields `u` and `v` from B- or C-grid
onto the default C-grid (centered on the cell faces).
"""
interpolateontodefaultCgrid(u, u_lon, u_lat, v, v_lon, v_lat, gridmetrics) = interpolateontodefaultCgrid(u, u_lon, u_lat, v, v_lon, v_lat, gridmetrics, getarakawagrid(u_lon, u_lat, v_lon, v_lat, gridmetrics))
interpolateontodefaultCgrid(u, u_lon, u_lat, v, v_lon, v_lat, gridmetrics, ::CGridCell) = u, u_lon, u_lat, v, v_lon, v_lat
interpolateontodefaultCgrid(u, u_lon, u_lat, v, v_lon, v_lat, gridmetrics, ::AGridCell) = error("Interpolation not implemented for A-grid type")
function interpolateontodefaultCgrid(u, u_lon, u_lat, v, v_lon, v_lat, gridmetrics, arakawagrid::BGridCell)

    (; u_pos, v_pos) = arakawagrid
    u_pos == v_pos == :NE || error("Interpolation not implemented for this B-grid($u_pos,$v_pos) type")

    _FillValue = u.properties["_FillValue"]

    # Make sure u and v are in memory
    u = u |> Array
    v = v |> Array

    nx, ny, nz = size(u)

    # unpack gridmetrics
    (; lon_vertices, lat_vertices) = gridmetrics

    # It seems that u/v is NaN on boundaries (for ACCESS-ESM1-5)
    # and that umo/vmo were computed as if u/v were 0 on boundaries
    # so that's what we do here
    u2 = replace(u, _FillValue => 0.0)
    v2 = replace(v, _FillValue => 0.0)
    u2 = 0.5(u2 + [fill(0.0, nx, 1, nz);; u2[:, 1:end-1, :]])
    v2 = 0.5(v2 + [fill(0.0, 1, ny, nz); v2[1:end-1, :, :]])
    SE_points = [(lon, lat) for (lon, lat) in zip(lon_vertices[2, :, :], lat_vertices[2, :, :])]
    NE_points = [(lon, lat) for (lon, lat) in zip(lon_vertices[3, :, :], lat_vertices[3, :, :])]
    NW_points = [(lon, lat) for (lon, lat) in zip(lon_vertices[4, :, :], lat_vertices[4, :, :])]
    u2_points = [midpointonsphere(SE, NE) for (SE, NE) in zip(NE_points, SE_points)]
    v2_points = [midpointonsphere(NE, NW) for (NE, NW) in zip(NW_points, NE_points)]
    u2_lon = [P[1] for P in u2_points]
    u2_lat = [P[2] for P in u2_points]
    v2_lon = [P[1] for P in v2_points]
    v2_lat = [P[2] for P in v2_points]
    return u2, u2_lon, u2_lat, v2, v2_lon, v2_lat

end





"""
    vertexpermutation(lon_vertices, lat_vertices)

Returns the permutation that sorts the vertices in the default orientation,
where the default orientation is the following:

```
    4 ────────┐ 3
              │   ▲
              │   │ j ("North")
              │
    1 ────────┘ 2
    ─► i ("East")
```
"""
function vertexpermutation(lon_vertices, lat_vertices)
    # Make sure the vertices are in the right shape (4, nx, ny)
    @assert size(lon_vertices, 1) == size(lat_vertices, 1) == 4
    # Take the first grid cell
    i = j = 1
    # Turn the vertices into points
    points = collect(zip(lon_vertices[:, i, j], lat_vertices[:, i, j]))
    points_east = collect(zip(lon_vertices[:, i+1, j], lat_vertices[:, i+1, j]))
    points_north = collect(zip(lon_vertices[:, i, j+1], lat_vertices[:, i, j+1]))
    # Find the common points
    common_east = Set(points) ∩ Set(points_east)
    common_noth = Set(points) ∩ Set(points_north)
    # Find the indices of the common points
    idx_east = findall(in(common_east), points)
    idx_north = findall(in(common_noth), points)
    idx3 = only(idx_east ∩ idx_north) # common to all 3 cells
    idx2 = only(setdiff(idx_east, idx3)) # common to (i,j) and (i+1,j) only
    idx4 = only(setdiff(idx_north, idx3)) # common to (i,j) and (i,j+1) only
    idx1 = only(setdiff(1:4, idx2, idx3, idx4)) # only in (i,j)
    return [idx1, idx2, idx3, idx4]
end




# distances
function horizontaldistance(lon, lat, I, J)
    I = horizontalindex(I)
    J = horizontalindex(J)
    PI = (getindexornan(lon, I), getindexornan(lat, I))
    PJ = (getindexornan(lon, J), getindexornan(lat, J))
    return haversine(PI, PJ)
end
horizontaldistance(lon, lat, I, ::Nothing) = NaN
function verticaldistance(Z::Vector, I, J)
    I = verticalindex(I)
    J = verticalindex(J)
    return abs(getindexornan(Z, J) - getindexornan(Z, I))
end
verticaldistance(Z::Array, I, J) = abs(getindexornan(Z, J) - getindexornan(Z, I))

horizontalindex(I::CartesianIndex{2}) = I
horizontalindex(I::CartesianIndex{3}) = CartesianIndex(I.I[1], I.I[2])
verticalindex(I::CartesianIndex{1}) = I
verticalindex(I::CartesianIndex{3}) = CartesianIndex(I.I[3])




# The default orientation is the following:
#
#     4 ────┐ 3
#           │
#     1 ────┘ 2
#
function vertexindices(dir)
    dir == :south ? (1, 2) :
    dir == :east ? (2, 3) :
    dir == :north ? (3, 4) :
    dir == :west ? (1, 4) :
    error()
end
vertexpoint(vlon, vlat, i, j, vertexidx) = (vlon[vertexidx,i,j], vlat[vertexidx,i,j])
function verticalfacewidth(vlon, vlat, i, j, dir)
    a, b = vertexindices(dir)
    A = vertexpoint(vlon, vlat, i, j, a)
    B = vertexpoint(vlon, vlat, i, j, b)
    haversine(A, B)
end
verticalfacewidth(edge_length_2D, i, j, dir) = edge_length_2D[dir][i, j]

function verticalfacearea(vlon, vlat, lev_bnds_or_thkcello, i, j, k, dir)
    height = cellthickness(lev_bnds_or_thkcello, i, j, k)
    width = verticalfacewidth(vlon, vlat, i, j, dir)
    height * width
end
function verticalfacearea(edge_length_2D, lev_bnds_or_thkcello, i, j, k, dir)
    height = cellthickness(lev_bnds_or_thkcello, i, j, k)
    width = verticalfacewidth(edge_length_2D, i, j, dir)
    height * width
end

cellthickness(lev_bnds::Matrix, i, j, k) = abs(lev_bnds[2,k] - lev_bnds[1,k])
cellthickness(thkcello, i, j, k) = thkcello[i, j, k]


function centroid2edgedistance(lon, lat, vlon, vlat, i, j, dir)
    a, b = vertexindices(dir)
    C = (lon[i, j], lat[i, j])
    A = vertexpoint(vlon, vlat, i, j, a)
    B = vertexpoint(vlon, vlat, i, j, b)
    M = midpointonsphere(A, B)
    haversine(C, M)
end

function midpointonsphere(A, B)
    if abs(A[1] - B[1]) < 180
        (A .+ B) ./ 2
    else # if the edge crosses the longitudinal edge of the map
        (A .+ B) ./ 2 .+ (180, 0)
    end
end


function horizontalcentroiddistance(lon, lat, iA, jA, iB, jB)
    A = (lon[iA, jA], lat[iA, jA])
    B = (lon[iB, jB], lat[iB, jB])
    return haversine(A, B)
end




function makegridmetrics(; areacello, volcello, lon, lat, lev, lon_vertices, lat_vertices)

	# volume (3D)
    FillValue = volcello.properties["_FillValue"]
    v3D = volcello |> Array{Union{Missing, Float64}}
	v3D = replace(v3D, missing => NaN, 0 => NaN, FillValue => NaN)

    # area (2D)
    FillValue = areacello.properties["_FillValue"]
    area2D = areacello |> Array{Union{Missing, Float64}}
    area2D = replace(area2D, missing => NaN, 0 => NaN, FillValue => NaN)

	# depth and cell height (3D)
	thkcello = v3D ./ area2D
    ZBOT3D = cumsum(thkcello, dims = 3)
    Z3D = ZBOT3D - 0.5 * thkcello
	zt = lev |> Array

    lat = lat |> Array
    lon = lon |> Array

    # same with lon_vertices
    lon_vertices = lon_vertices |> Array{Float64}
    lat_vertices = lat_vertices |> Array{Float64}

    # sort the vertices to match the default orientation
    vertexidx = vertexpermutation(lon_vertices, lat_vertices)
    lon_vertices = lon_vertices[vertexidx,:,:]
    lat_vertices = lat_vertices[vertexidx,:,:]

    C = CartesianIndices(size(lon))

    gridtopology = getgridtopology(lon_vertices, lat_vertices, zt)

    dirs = (:south, :east, :north, :west)
    𝑗s = (j₋₁, i₊₁, j₊₁, i₋₁)
    edge_length_2D = Dict(d=>[verticalfacewidth(lon_vertices, lat_vertices, 𝑖.I[1], 𝑖.I[2], d) for 𝑖 in C] for d in dirs)
    distance_to_edge_2D = Dict(d=>[centroid2edgedistance(lon, lat, lon_vertices, lat_vertices, 𝑖.I[1], 𝑖.I[2], d) for 𝑖 in C] for d in dirs)
    distance_to_neighbour_2D = Dict(d=>[horizontaldistance(lon, lat, 𝑖, 𝑗(𝑖, gridtopology)) for 𝑖 in C] for (d,𝑗) in zip(dirs,𝑗s))

	return (; area2D, v3D, thkcello, lon_vertices, lat_vertices, lon, lat, Z3D, zt, edge_length_2D, distance_to_edge_2D, distance_to_neighbour_2D, gridtopology)
end

