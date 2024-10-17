

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
    arakawa = getarakawagrid(u_lon, u_lat, v_lon, v_lat, modelgrid)

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
function getarakawagrid(u_lon, u_lat, v_lon, v_lat, modelgrid)

    # Unpack modelgrid
    (; lon, lat, lon_vertices, lat_vertices) = modelgrid

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
    u, u_lon, u_lat, v, v_lon, v_lat = interpolateontodefaultCgrid(u, u_lon, u_lat, v, v_lon, v_lat, modelgrid)

Interpolates the velocity fields `u` and `v` from B- or C-grid
onto the default C-grid (centered on the cell faces).
"""
interpolateontodefaultCgrid(u, u_lon, u_lat, v, v_lon, v_lat, modelgrid) = interpolateontodefaultCgrid(u, u_lon, u_lat, v, v_lon, v_lat, modelgrid, getarakawagrid(u_lon, u_lat, v_lon, v_lat, modelgrid))
interpolateontodefaultCgrid(u, u_lon, u_lat, v, v_lon, v_lat, modelgrid, ::CGridCell) = u, u_lon, u_lat, v, v_lon, v_lat
interpolateontodefaultCgrid(u, u_lon, u_lat, v, v_lon, v_lat, modelgrid, ::AGridCell) = error("Interpolation not implemented for A-grid type")
function interpolateontodefaultCgrid(u, u_lon, u_lat, v, v_lon, v_lat, modelgrid, arakawa::BGridCell)

    (; u_pos, v_pos) = arakawa
    u_pos == v_pos == :NE || error("Interpolation not implemented for this B-grid($u_pos,$v_pos) type")

    _FillValue = u.properties["_FillValue"]

    # Make sure u and v are in memory
    u = u |> Array
    v = v |> Array

    nx, ny, nz = size(u)

    # unpack modelgrid
    (; lon_vertices, lat_vertices) = modelgrid

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
