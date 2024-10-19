
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
    facefluxesfrommasstransport(; umo, vmo)

Return the fluxes integrated over each cell face (east, west, north, south, top, bottom)
given mass transport `umo` and `vmo` (fluxes across faces).

See also `facefluxes`.
"""
function facefluxesfrommasstransport(; umo, vmo)

    FillValue = umo.properties["_FillValue"]
    @assert isequal(FillValue, vmo.properties["_FillValue"])

    # Convert to in-memory Array to avoid slow getindexornan
    # Convert to Float64 for double-precision mass conservation
    umo = umo |> Array{Float64}
    vmo = vmo |> Array{Float64}

    return facefluxes(umo, vmo; FillValue)

end

"""
    facefluxesfromvelocities(; uo, uo_lon, uo_lat, vo, vo_lon, vo_lat, gridmetrics, ρ)

Return the fluxes integrated over each cell face (east, west, north, south, top, bottom)
given either `uo`, `vo`, their lon/alt locations, and `gridmetrics` and `ρ`.

See also `facefluxes`.
"""
function facefluxesfromvelocities(; uo, uo_lon, uo_lat, vo, vo_lon, vo_lat, gridmetrics, ρ)

    FillValue = uo.properties["_FillValue"]
    @assert isequal(FillValue, vo.properties["_FillValue"])

    # Convert to in-memory Array to avoid slow getindexornan
    # Convert to Float64 for double-precision mass conservation
    umo, vmo = velocity2fluxes(uo, uo_lon, uo_lat, vo, vo_lon, vo_lat, gridmetrics, ρ)

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




