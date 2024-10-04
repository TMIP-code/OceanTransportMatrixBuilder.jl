"""
    arakawa, uo_pos, vo_pos, uo_relerr, vo_relerr = gridtype(uo_lon, uo_lat, vo_lon, vo_lat, modelgrid)

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
function gridtype(uo_lon, uo_lat, vo_lon, vo_lat, modelgrid)

    # Unpack modelgrid
    (; lon, lat, lon_vertices, lat_vertices) = modelgrid

    i = j = 1
    uo_point = (uo_lon[i, j], uo_lat[i, j])
    vo_point = (vo_lon[i, j], vo_lat[i, j])

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

    uo_distances = (; (k => haversine(P, uo_point) for (k,P) in pairs(cell))...)
    vo_distances = (; (k => haversine(P, vo_point) for (k,P) in pairs(cell))...)

    uo_distance, uo_pos = findmin(uo_distances)
    vo_distance, vo_pos = findmin(vo_distances)

    # Arakawa grid type
    if uo_pos == vo_pos == :C
        arakawa = :A
    elseif uo_pos == vo_pos && uo_pos ∈ (:NE, :NW, :SE, :SW)
        arakawa = :B
    elseif uo_pos ∈ (:E, :W) && vo_pos ∈ (:N, :S)
        arakawa = :C
    else
        error("Unknown Arakawa grid type")
    end

    cellperimeter = haversine(SW, SE) + haversine(SE, NE) + haversine(NE, NW) + haversine(NW, SW)
    uo_relerr = uo_distance / cellperimeter
    vo_relerr = vo_distance / cellperimeter

    return arakawa, uo_pos, vo_pos, uo_relerr, vo_relerr

end

function interpolateontodefaultCgrid(uo, uo_lon, uo_lat, vo, vo_lon, vo_lat, modelgrid)

    _FillValue = uo.properties["_FillValue"]

    # Make sure uo and vo are in memory
    uo = uo |> Array
    vo = vo |> Array

    nx, ny, nz = size(uo)

    arakawa, uo_pos, vo_pos = gridtype(uo_lon, uo_lat, vo_lon, vo_lat, modelgrid)

    # unpack modelgrid
    (; lon_vertices, lat_vertices) = modelgrid

    if arakawa == :C && uo_pos == :E && vo_pos == :N
        return uo, uo_lon, uo_lat, vo, vo_lon, vo_lat
    elseif arakawa == :B && uo_pos == vo_pos == :NE
        # It seems that uo/vo is NaN on boundaries (for ACCESS-ESM1-5)
        # and that umo/vmo were computed as if uo/vo were 0 on boundaries
        # so that's what we do here
        uo2 = replace(uo, _FillValue => 0.0)
        vo2 = replace(vo, _FillValue => 0.0)
        uo2 = 0.5(uo2 + [fill(0.0, nx, 1, nz);; uo2[:, 1:end-1, :]])
        vo2 = 0.5(vo2 + [fill(0.0, 1, ny, nz); vo2[1:end-1, :, :]])
        SE_points = [(lon, lat) for (lon, lat) in zip(lon_vertices[2, :, :], lat_vertices[2, :, :])]
        NE_points = [(lon, lat) for (lon, lat) in zip(lon_vertices[3, :, :], lat_vertices[3, :, :])]
        NW_points = [(lon, lat) for (lon, lat) in zip(lon_vertices[4, :, :], lat_vertices[4, :, :])]
        uo2_points = [midpointonsphere(SE, NE) for (SE, NE) in zip(NE_points, SE_points)]
        vo2_points = [midpointonsphere(NE, NW) for (NE, NW) in zip(NW_points, NE_points)]
        uo2_lon = [P[1] for P in uo2_points]
        uo2_lat = [P[2] for P in uo2_points]
        vo2_lon = [P[1] for P in vo2_points]
        vo2_lat = [P[2] for P in vo2_points]
        return uo2, uo2_lon, uo2_lat, vo2, vo2_lon, vo2_lat
    else
        error("Interpolation not implemented for this grid type")
    end

end