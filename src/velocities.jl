
"""
    velocity2fluxes(u, u_lon, u_lat, v, v_lon, v_lat, gridmetrics, ρ)

Return the fluxes (integrated over cell faces) given veloticies `u` and `v` and their (lon,lat).

The fluxes are calculated from the density ρ (number or 3D array) and from the cell face areas between the cells.
The mean density from the two cells that share the face is used.
The minimum thickness of the two cells that share the face is used.
"""
function velocity2fluxes(u, u_lon, u_lat, v, v_lon, v_lat, gridmetrics, ρ)

    # Unpack gridmetrics
    (; thkcello, edge_length_2D, v3D, lon_vertices, lat_vertices, zt) = gridmetrics

    # make indices
    indices = makeindices(v3D)

    # Interpolate to C-grid
    u, _, _, v, _, _ = interpolateontodefaultCgrid(u, u_lon, u_lat, v, v_lon, v_lat, gridmetrics)

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
    facefluxesfrommasstransport(; umo, vmo, gridmetrics, indices)

Return the fluxes integrated over each cell face (east, west, north, south, top, bottom)
given mass transport `umo` and `vmo` (fluxes across faces).

See also `facefluxes`.
"""
function facefluxesfrommasstransport(; umo, vmo, gridmetrics, indices)

    FillValue = umo.properties["_FillValue"]
    @assert isequal(FillValue, vmo.properties["_FillValue"])

    # Convert to in-memory Array to avoid slow getindexornan
    # Convert to Float64 for double-precision mass conservation
    umo = umo |> Array{Float64}
    vmo = vmo |> Array{Float64}

    return facefluxes(umo, vmo, gridmetrics, indices; FillValue)

end

"""
    facefluxesfromvelocities(; uo, uo_lon, uo_lat, vo, vo_lon, vo_lat, gridmetrics, indices, ρ)

Return the fluxes integrated over each cell face (east, west, north, south, top, bottom)
given either `uo`, `vo`, their lon/alt locations, and `gridmetrics` and `ρ`.

See also `facefluxes`.
"""
function facefluxesfromvelocities(; uo, uo_lon, uo_lat, vo, vo_lon, vo_lat, gridmetrics, indices, ρ)

    FillValue = uo.properties["_FillValue"]
    @assert isequal(FillValue, vo.properties["_FillValue"])

    # Convert to in-memory Array to avoid slow getindexornan
    # Convert to Float64 for double-precision mass conservation
    umo, vmo = velocity2fluxes(uo, uo_lon, uo_lat, vo, vo_lon, vo_lat, gridmetrics, ρ)

    return facefluxes(umo, vmo, gridmetrics, indices; FillValue)

end


function nofluxboundaries!(ϕᵢ, ϕⱼ, gridmetrics, indices)


    # unpack gridmetrics and indices
    (; C, wet3D) = indices
    (; gridtopology) = gridmetrics

    for i in eachindex(wet3D)

        E = i₊₁(C[i], gridtopology)
        N = j₊₁(C[i], gridtopology)

        # If i is land, just enforce zero east and north fluxes
        !wet3D[i] && (ϕᵢ[i] = ϕⱼ[i] = 0)

        # If eastern neighbor is nothing or land, enforce zero flux in the east direction
        (isnothing(E) || !wet3D[E]) && (ϕᵢ[i] = 0)

        # If northern neighbor is nothing or land, enforce zero flux in the north direction
        (isnothing(N) || !wet3D[N]) && (ϕⱼ[i] = 0)

    end

    return ϕᵢ, ϕⱼ

end

"""
    facefluxes(umo, vmo, gridmetrics, indices; FillValue)

Return the fluxes integrated over each cell face (east, west, north, south, top, bottom)
given the east and north fluxes, `umo` and `vmo`.

The west and south fluxes are computed by simple shift in coordinates.
The top and bottom fluxes are computed by mass conservation from the seafloor up.
"""
function facefluxes(umo, vmo, gridmetrics, indices; FillValue)

    umo, vmo = nofluxboundaries!(umo, vmo, gridmetrics, indices)

    # unpack gridmetrics and indices
    (; C) = indices
    (; gridtopology) = gridmetrics

    # Check that not all values are NaN or FillValue
    @assert !all(x -> isnan(x) || (x == FillValue), umo)
    @assert !all(x -> isnan(x) || (x == FillValue), vmo)

	@debug "Making ϕeast"
	ϕeast = replace(umo, NaN=>0.0, FillValue=>0.0) .|> Float64
	@debug "Making ϕwest"
    # ϕwest[i] is ϕeast[west of i]
    ϕwest = zeros(size(ϕeast))
    for i in eachindex(ϕwest)
        W = i₋₁(C[i], gridtopology)
        isnothing(W) && continue
        ϕwest[i] = ϕeast[W]
    end

	@debug "Making ϕnorth"
	# Check that south pole ϕnorth is zero (so that it causes no issues with circular shift)
	ϕnorth = replace(vmo, NaN=>0.0, FillValue=>0.0) .|> Float64
	@debug "Making ϕsouth"
    # ϕsouth[i] is ϕnorth[south of i]
    # Note from BP: This might break if there is a seam at the South Pole
	ϕsouth = zeros(size(ϕnorth))
    for i in eachindex(ϕsouth)
        S = j₋₁(C[i], gridtopology)
        isnothing(S) && continue
        ϕsouth[i] = ϕnorth[S]
    end

	@debug "Making ϕtop and ϕbottom"
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




