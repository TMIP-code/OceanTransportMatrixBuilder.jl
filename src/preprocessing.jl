
"""
    velocity2fluxes(u, u_lon, u_lat, v, v_lon, v_lat, modelgrid, ρ)

Return the fluxes (integrated over cell faces) given veloticies `u` and `v` and their (lon,lat).

The fluxes are calculated from the density ρ (number or 3D array) and from the cell face areas between the cells.
The mean density from the two cells that share the face is used.
The minimum thickness of the two cells that share the face is used.
"""
function velocity2fluxes(u, u_lon, u_lat, v, v_lon, v_lat, modelgrid, ρ)

    # Unpack modelgrid
    (; thkcello, edge_length_2D, v3D, lon_vertices, lat_vertices, zt) = modelgrid

    # make indices
    indices = makeindices(v3D)

    # Interpolate to C-grid
    u, _, _, v, _, _ = interpolateontodefaultCgrid(u, u_lon, u_lat, v, v_lon, v_lat, modelgrid)

    # grid type
    gridtopology = getgridtopology(lon_vertices, lat_vertices, zt)

    ϕᵢ = zeros(size(u))
    ϕⱼ = zeros(size(v))

    # Calculate fluxes from velocities
    for 𝑖 in indices.C
        i, j = 𝑖.I # indices to access the 2D array
        𝑗 = i₊₁(𝑖, gridtopology) # grid cell to the "east"
        ϕᵢ[𝑖] = u[𝑖] * twocellnanmean(ρ, 𝑖, 𝑗) * twocellnanmin(thkcello, 𝑖, 𝑗) * edge_length_2D[:east][i,j]
        𝑗 = j₊₁(𝑖, gridtopology) # grid cell to the "north"
        ϕⱼ[𝑖] = v[𝑖] * twocellnanmean(ρ, 𝑖, 𝑗) * twocellnanmin(thkcello, 𝑖, 𝑗) * edge_length_2D[:north][i,j]
    end
    # TODO check that this orientation is correct for all models?
    # I think this only applies to Arakawa C-grids...

    return ϕᵢ, ϕⱼ
end

"""
    twocellnanmean(x, 𝑖, 𝑗)

Return the nanmean of `x` at the two cell indices `𝑖` and `𝑗`.
"""
twocellnanmean(x::Number, 𝑖, 𝑗) = x
twocellnanmean(x, 𝑖, 𝑗) = nanmean2(x[𝑖], x[𝑗])

"""
    nanmean2(a, b)

Return the nanmean of `a` and `b` (scalars).
"""
function nanmean2(a, b)
    wa = !isnan(a)
    wb = !isnan(b)
    (wa * a + wb * b) / (wa + wb)
end

"""
    twocellnanmin(x, 𝑖, 𝑗)

Return the nanmin of `x` at the two cell indices `𝑖` and `𝑗`.
"""
twocellnanmin(x::Number, 𝑖, 𝑗) = x
twocellnanmin(x, 𝑖, 𝑗) = nanmin2(x[𝑖], x[𝑗])

"""
    nanmin2(a, b)

Return the nanmin of `a` and `b` (scalars).
"""
nanmin2(a, b) = isnan(a) ? b : isnan(b) ? a : min(a, b)

"""
    facefluxesfrommasstransport(; umo, vmo)

Return the fluxes integrated over each cell face (east, west, north, south, top, bottom)
given mass transport `umo` and `vmo` (fluxes across faces).

See also `facefluxes`.
"""
function facefluxesfrommasstransport(; umo, vmo)

    FillValue = umo.properties["_FillValue"]
    @assert isequal(FillValue, vmo.properties["_FillValue"])

    # Convert to in-memory Array to avoid slow getindex
    # Convert to Float64 for double-precision mass conservation
    umo = umo |> Array{Float64}
    vmo = vmo |> Array{Float64}

    return facefluxes(umo, vmo; FillValue)

end

"""
    facefluxesfromvelocities(; uo, uo_lon, uo_lat, vo, vo_lon, vo_lat, modelgrid, ρ)

Return the fluxes integrated over each cell face (east, west, north, south, top, bottom)
given either `uo`, `vo`, their lon/alt locations, and `modelgrid` and `ρ`.

See also `facefluxes`.
"""
function facefluxesfromvelocities(; uo, uo_lon, uo_lat, vo, vo_lon, vo_lat, modelgrid, ρ)

    FillValue = uo.properties["_FillValue"]
    @assert isequal(FillValue, vo.properties["_FillValue"])

    # Convert to in-memory Array to avoid slow getindex
    # Convert to Float64 for double-precision mass conservation
    umo, vmo = velocity2fluxes(uo, uo_lon, uo_lat, vo, vo_lon, vo_lat, modelgrid, ρ)

    return facefluxes(umo, vmo; FillValue)

end

"""
    facefluxes(umo, vmo; FillValue)

Return the fluxes integrated over each cell face (east, west, north, south, top, bottom)
given the east and north fluxes, `umo` and `vmo`.

The west and south fluxes are computed by simple shift in coordinates.
The top and bottom fluxes are computed by mass conservation from the seafloor up.
"""
function facefluxes(umo, vmo; FillValue)

	iseastborder = isnan.(umo[[2:end;1],:,:]) .& .!isnan.(umo)
	isnorthborder = isnan.([vmo[:,2:end,:] vmo[end:-1:1,end:end,:]]) .& .!isnan.(vmo)
	umo[iseastborder] .= 0
	vmo[isnorthborder] .= 0

	@info "Making ϕeast and ϕwest"
	ϕeast = replace(umo, NaN=>0.0, FillValue=>0.0) .|> Float64
	ϕwest = ϕeast[[end;1:end-1],:,:] .|> Float64

	@info "Making ϕnorth and ϕsouth"
	# Check that south pole ϕnorth is zero (so that it causes no issues with circular shift)
	ϕnorth = replace(vmo, NaN=>0.0, FillValue=>0.0) .|> Float64
	ϕsouth = ϕnorth[:,[end;1:end-1],:] .|> Float64

	@info "Making ϕtop and ϕbottom"
	# Then build ϕtop and ϕbottom from bottom up
	# Mass conservation implies that
	# 	ϕwest + ϕsouth + ϕbottom - ϕeast - ϕnorth - ϕtop = 0
	# except at the top.
	# We could build w from the bottom up, by looking at each
	# water column's sea floor, but it's simpler to go from the top down,
	# and then remove ϕbottom
	ϕbottom = similar(ϕeast)
	ϕtop = similar(ϕeast)
	for k in reverse(eachindex(axes(ϕeast, 3)))
		if k == lastindex(axes(ϕeast, 3))
			@views @. ϕbottom[:,:,k] = 0 # seafloor ϕbottom is zero
		else
			@views @. ϕbottom[:,:,k] = ϕtop[:,:,k+1] # otherwise it's ϕtop from below
		end
		@views @. ϕtop[:,:,k] = ϕbottom[:,:,k] + ϕwest[:,:,k] + ϕsouth[:,:,k] - ϕeast[:,:,k] - ϕnorth[:,:,k]
	end

    ϕ = (
        east = ϕeast,
        west = ϕwest,
        north = ϕnorth,
        south = ϕsouth,
        top = ϕtop,
        bottom = ϕbottom
    )

	return ϕ
end


function makemodelgrid(; areacello, volcello, lon, lat, lev, lon_vertices, lat_vertices)

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
	zt = lev |> Array

    lat = lat |> Array
    lon = lon |> Array

    # same with lon_vertices
    lon_vertices = lon_vertices |> Array{Float64}
    lat_vertices = lat_vertices |> Array{Float64}

    # sort the vertices to mathc the default orientation
    vertexidx = vertexpermutation(lon_vertices, lat_vertices)
    lon_vertices = lon_vertices[vertexidx,:,:]
    lat_vertices = lat_vertices[vertexidx,:,:]

    C = CartesianIndices(size(lon))

    dirs = (:south, :east, :north, :west)
    edge_length_2D = Dict(d=>[verticalfacewidth(lon_vertices, lat_vertices, 𝑖.I[1], 𝑖.I[2], d) for 𝑖 in C] for d in dirs)
    distance_to_edge_2D = Dict(d=>[centroid2edgedistance(lon, lat, lon_vertices, lat_vertices, 𝑖.I[1], 𝑖.I[2], d) for 𝑖 in C] for d in dirs)

    arakawagrid = getgridtopology(lon_vertices, lat_vertices, zt)

	return (; area2D, v3D, thkcello, lon_vertices, lat_vertices, lon, lat, zt, edge_length_2D, distance_to_edge_2D, arakawagrid)
end

function makeindices(v3D)

    # LinearIndices and CartesianIndices required for building upwind operator
    nxyz = size(v3D)
    L = LinearIndices(nxyz)
    Lwet = L[.!isnan.(v3D)]
    N = length(Lwet)
    wet3D = falses(nxyz...)
    wet3D[Lwet] .= true
    Lwet3D = Array{Union{Int, Missing}, 3}(missing, nxyz...)
    Lwet3D[Lwet] .= 1:length(Lwet)
    C = CartesianIndices(nxyz)

    return (; wet3D, L, Lwet, N, Lwet3D, C)
end



# TODO generalized topology detection
# That is, how do I figure out how the boundaries are connected for other models than ACCESS?


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


