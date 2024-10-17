
"""
    buildTadv(; ϕ, modelgrid, indices, ρ)

Build the advection operator Tadv.
"""
function buildTadv(; ϕ, modelgrid, indices, ρ)
    # default ρ = 1035 kg/m^3 is the value originally used by Chamberlain et al. (2019)

	@info "Building Tadv"
	𝑖s, 𝑗s, Tvals = upwind_advection_operator_sparse_entries(ϕ, modelgrid, indices, ρ)

    N = indices.N

	any(isnan.(Tvals)) && error("Tadv contains NaNs.")

	Tadv = sparse(𝑖s, 𝑗s, Tvals, N, N)

	return Tadv
end

"""
    buildTκH(; modelgrid, indices, ρ, κH)

Build the horizontal diffusivity operator TκH.
"""
function buildTκH(; modelgrid, indices, ρ, κH)

    N = indices.N

    # Wet mask for horizontal diffusivity
	ΩH = trues(N)

	@info "Building TκH"
	𝑖s, 𝑗s, Tvals = horizontal_diffusion_operator_sparse_entries(; modelgrid, indices, κH, ΩH)

	any(isnan.(Tvals)) && error("TκH contains NaNs.")

	TκH = sparse(𝑖s, 𝑗s, Tvals, N, N)

	return TκH
end


"""
    buildTκVML(; mlotst, modelgrid, indices, κVML)

Build the mixed layer diffusivity operator TκVML.
"""
function buildTκVML(; mlotst, modelgrid, indices, κVML)

    # Unpack model grid
    (; zt, ) = modelgrid

    # Unpack indices
    (; Lwet, N) = indices

	mlotst = mlotst |> Array # to prevent slow getindex for lazily loaded data?

	# Wet mask for mixed layer diffusivity
	Ω = replace(reshape(zt, 1, 1, length(zt)) .< mlotst, missing=>false)[Lwet]

	@info "Building TκVML "
	𝑖s, 𝑗s, Tvals = vertical_diffusion_operator_sparse_entries(; modelgrid, indices, κV = κVML, Ω)

	any(isnan.(Tvals)) && error("TκVML contains NaNs.")

	TκVML = sparse(𝑖s, 𝑗s, Tvals, N, N)

	return TκVML
end


"""
    buildTκVdeep(; mlotst, modelgrid, indices, κVdeep)

Build the deep diffusivity operator TκVdeep.
"""
function buildTκVdeep(; mlotst, modelgrid, indices, κVdeep)

    N = indices.N

	# Deep mask for vertical diffusivity
	Ω = trues(N) # TODO (maybe): make Ωdeep not overlap with ΩML at MLD?

	@info "Building TκVdeep"
	𝑖s, 𝑗s, Tvals = vertical_diffusion_operator_sparse_entries(; modelgrid, indices, κV = κVdeep, Ω)

	any(isnan.(Tvals)) && error("TκVdeep contains NaNs.")

	TκVdeep = sparse(𝑖s, 𝑗s, Tvals, N, N)

	return TκVdeep

end

"""
    transportmatrix(; ϕ, mlotst, modelgrid, indices, ρ, κH, κVML, κVdeep, Tadv, TκH, TκVML, TκVdeep)

Build the transport matrix, i.e., the flux-divergence operator T = Tadv + TκH + TκVML + TκVdeep,
and check divergence and mass conservation.
"""
function transportmatrix(; ϕ, mlotst, modelgrid, indices, ρ,
		κH = 500.0, # m^2/s,
		κVML = 0.1, # m^2/s,
		κVdeep = 1e-5, # m^2/s,
		Tadv = buildTadv(; ϕ, modelgrid, indices, ρ),
		TκH = buildTκH(; modelgrid, indices, ρ, κH),
		TκVML = buildTκVML(; mlotst, modelgrid, indices, κVML),
		TκVdeep = buildTκVdeep(; mlotst, modelgrid, indices, κVdeep),
	)

	@info "Building T"

	@time T = Tadv + TκH + TκVML + TκVdeep

	return (; T, Tadv, TκH, TκVML, TκVdeep)
end






# Some personal notes
# ϕ are water mass transports, in kg/s
# T[i,j] is the inverse timescale with which upwind j→i occurs such that
# The tendency of tracer χ at box i due to box j is given
# 	∂χ[i] = -T[i,j] * χ[j]    units:   1/s * mol/kg (or mol/m^3)
# It can be also computed as a upwind mass transfer if ϕ[j→i] ≥	0:
# 	∂χ[i] = ϕ[j→i] * χ[j] / m[i]     units    kg[j]/s * mol/kg[j] / kg[i] = mol/kg[i]
# and should also incur the opposite mass tendency at j:
# 	∂χ[j] = -ϕ[j→i] * χ[j] / m[j]     units    kg[j]/s * mol/kg[j] / kg[j] = mol/kg[j]
# Thus the matrix term should be constructed as
# 	T[i,j] = -ϕ[j→i] / m[i]             units = kg s⁻¹ / kg = s⁻¹
# and ϕ[j→i] / m[j] should be added to the diagonal T[j,j].

# Is this mass conserving?
# Only if
#   v[i] * T[i,i] + v[j] * T[j,i] ≈ 0
#   = vi ϕ[i→j] / m[i] + vj -ϕ[i→j] / m[j]
#   = ϕ[i→j] / ρ[i] - ϕ[i→j] / ρ[j]
# i.e., iff
#   ρ[i] = ρ[j]
# So I must use the mean density between facing cells.

# There are three ways to index:
# - Cartesian indices (i,j,k)
# - Linear index L𝑖
# - Wet linear index 𝑖 (that's the one we want to record for the matrix)
# So to fill T[𝑖,𝑗] -ϕ[𝑖→𝑗] / m[𝑖], I need be able to convert, in sequence:
# 𝑖 -> (i,j,k) -> neihghbour (i′,j′,k′) -> 𝑗
# The first 2 conversions are straightforard.
# For the last one, I make a 3D array filled with the wet linear indices

"""
    upwind_advection_operator_sparse_entries(ϕ, modelgrid, indices, ρ)

Return the sparse (i, j, v) for the upwind advection operator Tadv.
"""
function upwind_advection_operator_sparse_entries(ϕ, modelgrid, indices, ρ::Number)
	# If ρ is a scalar, broadcast it to the modelgrid size
	ρ = fill(ρ, size(modelgrid.v3D))
	return upwind_advection_operator_sparse_entries(ϕ, modelgrid, indices, ρ)
end
function upwind_advection_operator_sparse_entries(ϕ, modelgrid, indices, ρ)

    # Unpack model grid
    (; v3D, arakawagrid) = modelgrid
    # Unpack indices
    (; Lwet, Lwet3D, C) = indices

	any(isnan.(ρ[Lwet])) && error("ρ contains NaNs")

	𝑖s, 𝑗s, Tvals = Int[], Int[], Float64[]

    @time for 𝑖 in eachindex(Lwet)
		L𝑖 = Lwet[𝑖]
		C𝑖 = C[L𝑖]
		i, j, k = C𝑖.I
		v𝑖 = v3D[C𝑖]
		ρ𝑖 = ρ[C𝑖]
		# From West
		ϕwest = ϕ.west[C𝑖]
		if ϕwest > 0
			C𝑗 = i₋₁(C𝑖, arakawagrid)
			𝑗 = Lwet3D[C𝑗]
			v𝑗 = v3D[C𝑗]
			ρ𝑗 = ρ[C𝑗]
			pushTadvectionvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗, ϕwest, ρ𝑖, ρ𝑗, v𝑖, v𝑗)
		end
		# From East
		ϕeast = ϕ.east[C𝑖]
		if ϕeast < 0
			C𝑗 = i₊₁(C𝑖, arakawagrid)
			𝑗 = Lwet3D[C𝑗]
			v𝑗 = v3D[C𝑗]
			ρ𝑗 = ρ[C𝑗]
			pushTadvectionvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗, -ϕeast, ρ𝑖, ρ𝑗, v𝑖, v𝑗)
		end
		# From South
		ϕsouth = ϕ.south[C𝑖]
		if ϕsouth > 0
			C𝑗 = j₋₁(C𝑖, arakawagrid)
			𝑗 = Lwet3D[C𝑗]
			v𝑗 = v3D[C𝑗]
			ρ𝑗 = ρ[C𝑗]
			pushTadvectionvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗, ϕsouth, ρ𝑖, ρ𝑗, v𝑖, v𝑗)
		end
		# From North (Special case with north bipole)
		ϕnorth = ϕ.north[C𝑖]
		if ϕnorth < 0
			C𝑗 = j₊₁(C𝑖, arakawagrid)
			𝑗 = Lwet3D[C𝑗]
			v𝑗 = v3D[C𝑗]
			ρ𝑗 = ρ[C𝑗]
			pushTadvectionvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗, -ϕnorth, ρ𝑖, ρ𝑗, v𝑖, v𝑗)
		end
		# From Bottom
		ϕbottom = ϕ.bottom[C𝑖]
		if ϕbottom > 0
			C𝑗 = k₊₁(C𝑖, arakawagrid)
			𝑗 = Lwet3D[C𝑗]
			v𝑗 = v3D[C𝑗]
			ρ𝑗 = ρ[C𝑗]
			pushTadvectionvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗, ϕbottom, ρ𝑖, ρ𝑗, v𝑖, v𝑗)
		end
		# From Top
		ϕtop = ϕ.top[C𝑖]
		if ϕtop < 0 && k > 1 # Evaporation/precipitation -> no change to χ
			C𝑗 = k₋₁(C𝑖, arakawagrid)
			𝑗 = Lwet3D[C𝑗]
			v𝑗 = v3D[C𝑗]
			ρ𝑗 = ρ[C𝑗]
			pushTadvectionvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗, -ϕtop, ρ𝑖, ρ𝑗, v𝑖, v𝑗)
		end
	end
	return 𝑖s, 𝑗s, Tvals
end


"""
    pushTadvectionvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗, ϕ, ρ𝑖, ρ𝑗, v𝑖, v𝑗)

Pushes the sparse indices and values into (𝑖s, 𝑗s, Tvals) corresponding to the j→i advection.
"""
function pushTadvectionvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗, ϕ, ρ𝑖, ρ𝑗, v𝑖, v𝑗)
	ρ = (ρ𝑖 + ρ𝑗) / 2
	m𝑖 = ρ * v𝑖
	m𝑗 = ρ * v𝑗
	push!(𝑖s, 𝑖)
	push!(𝑗s, 𝑗)
	push!(Tvals, -ϕ / m𝑖)
	push!(𝑖s, 𝑗)
	push!(𝑗s, 𝑗)
	push!(Tvals, ϕ / m𝑗)
end



# Following Chamberlain et al. (2019)
#   ∂χ/∂t =   Δχ/d   *     κ       *  a   /   v
#         = gradient * diffusivity * area / volume
#
# units   = χ m⁻¹ ⋅ m² s⁻¹ ⋅ m² / m³ = χ s⁻¹
#
#   ∂χ/∂t[i] = κ a / d * (χ[j] - χ[i]) / v[i]
# so that
#   T[i,i] =  κ a / d / v[i]
#   T[i,j] = -κ a / d / v[i]

# Is this mass conserving?
# Only if
#   v[i] * T[i,i] + v[j] * T[j,i] ≈ 0
#   = vi κ aij / dij / vi - vj κ aji / dji / vj
# i.e., iff
#   aij / dij = aji / dji
#   and
#   v[i] * T[i,j] + v[j] * T[j,j] ≈ 0
#   = vi κ aij / dij / vi - vj κ aji / dji / vj.
# I.e., I must check that the distances and areas don't vary with orientation.
#



# Is this divergence free?
# Only if
#     T[i,j] + T[i,i] ≈ 0
# which is the case by construction




"""
    horizontal_diffusion_operator_sparse_entries(; modelgrid, indices, κH, ΩH)

Return the sparse (i, j, v) for the horizontal diffusion operator TκH.
"""
function horizontal_diffusion_operator_sparse_entries(; modelgrid, indices, κH, ΩH)

    # Unpack model grid
    (; v3D, edge_length_2D, lon, lat, thkcello, arakawagrid) = modelgrid
    # Unpack indices
    (; wet3D, Lwet, Lwet3D, C) = indices

	𝑖s, 𝑗s, Tvals = Int[], Int[], Float64[]

    ny = size(wet3D, 2) # Should not be needed once oppdir is dealt by topology functions

    @time for 𝑖 in eachindex(Lwet)
        ΩH[𝑖] || continue # only continue if inside ΩH
		L𝑖 = Lwet[𝑖]
		C𝑖 = C[L𝑖]
		i, j, k = C𝑖.I
		V = v3D[C𝑖]
		# From West
		C𝑗W = i₋₁(C𝑖, arakawagrid)
        if !isnothing(C𝑗W)
			𝑗W = Lwet3D[C𝑗W]
			if !ismissing(𝑗W) && ΩH[𝑗W]
				iW, jW, _ = C𝑗W.I
				# (𝑖 == 𝑗W) && @show(i, j, iW, jW)
				# I take the minimum area from both dirs (through which mixing goes through)
				aij = verticalfacearea(edge_length_2D, thkcello, i, j, k, :west)
				aji = verticalfacearea(edge_length_2D, thkcello, iW, jW, k, :east)
				a = min(aij, aji)
				# I take the mean distance from both dirs
				d = horizontalcentroiddistance(lon, lat, i, j, iW, jW)
				pushTmixingvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗W, κH, a, d, V)
			end
		end
        # From East
		C𝑗E = i₊₁(C𝑖, arakawagrid)
        if !isnothing(C𝑗E)
			𝑗E = Lwet3D[C𝑗E]
			if !ismissing(𝑗E) && ΩH[𝑗E]
				iE, jE, _ = C𝑗E.I
				# (𝑖 == 𝑗E) && @show(i, j, iE, jE)
				aij = verticalfacearea(edge_length_2D, thkcello, i, j, k, :east)
				aji = verticalfacearea(edge_length_2D, thkcello, iE, jE, k, :west)
				a = min(aij, aji)
				d = horizontalcentroiddistance(lon, lat, i, j, iE, jE)
				pushTmixingvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗E, κH, a, d, V)
			end
		end
        # From South
		C𝑗S = j₋₁(C𝑖, arakawagrid)
		if !isnothing(C𝑗S)
			𝑗S = Lwet3D[C𝑗S]
			if !ismissing(𝑗S) && ΩH[𝑗S]
				iS, jS, _ = C𝑗S.I
				# (𝑖 == 𝑗S) && @show(i, j, iS, jS)
				aij = verticalfacearea(edge_length_2D, thkcello, i, j, k, :south)
				aji = verticalfacearea(edge_length_2D, thkcello, iS, jS, k, :north)
				a = min(aij, aji)
				d = horizontalcentroiddistance(lon, lat, i, j, iS, jS)
				pushTmixingvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗S, κH, a, d, V)
			end
		end
        # From North
		C𝑗N = j₊₁(C𝑖, arakawagrid)
        if !isnothing(C𝑗N)
			𝑗N = Lwet3D[C𝑗N]
			if !ismissing(𝑗N) && ΩH[𝑗N]
				# (𝑖 == 𝑗N) && @show(i, j, iN, jN)
				iN, jN, _ = C𝑗N.I
				# Note that the opposite direction (oppdir) is still north at j == ny
				# TODO: implement this into a topology.jl function
				oppdir = (j == ny) ? :north : :south
				aij = verticalfacearea(edge_length_2D, thkcello, i, j, k, :north)
				aji = verticalfacearea(edge_length_2D, thkcello, iN, jN, k, oppdir)
				a = min(aij, aji)
				d = horizontalcentroiddistance(lon, lat, i, j, iN, jN)
				pushTmixingvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗N, κH, a, d, V)
			end
        end
	end

	return 𝑖s, 𝑗s, Tvals
end


"""
    pushTmixingvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗, κ, a, d, V)

Pushes the sparse indices and values into (𝑖s, 𝑗s, Tvals) corresponding to the j→i mixing.
"""
function pushTmixingvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗, κ, a, d, V)
	Tval = κ * a / (d * V)
    if d == 0 || V == 0 || isnan(Tval)
        @show 𝑖, 𝑗, κ, a, d, V
        error()
    end
    push!(𝑖s, 𝑖)
	push!(𝑗s, 𝑖)
	push!(Tvals, Tval)
	push!(𝑖s, 𝑖)
	push!(𝑗s, 𝑗)
	push!(Tvals, -Tval)
end





function vertical_diffusion_operator_sparse_entries(; modelgrid, indices, κV, Ω)

    # Unpack model grid
    (; v3D, area2D, zt, arakawagrid) = modelgrid
    # Unpack indices
    (; wet3D, Lwet, Lwet3D, C) = indices

    𝑖s, 𝑗s, Tvals = Int[], Int[], Float64[]
	nxyz = size(wet3D)
    _, _, nz = nxyz

    @time for 𝑖 in eachindex(Lwet)
        Ω[𝑖] || continue # only continue if inside Ω
		L𝑖 = Lwet[𝑖]
		C𝑖 = C[L𝑖]
		i, j, k = C𝑖.I
		V = v3D[C𝑖]
        a = area2D[i,j]
		# From Bottom
		C𝑗B = k₊₁(C𝑖, arakawagrid)
		if !isnothing(C𝑗B)
			𝑗B = Lwet3D[C𝑗B]
			if !ismissing(𝑗B) && Ω[𝑗B] # only continue if inside Ω
				_, _, k′ = C𝑗B.I
				d = abs(zt[k] - zt[k′])
				pushTmixingvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗B, κV, a, d, V)
			end
		end
		# From Top
		C𝑗T = k₋₁(C𝑖, arakawagrid)
		if !isnothing(C𝑗T)
			𝑗T = Lwet3D[C𝑗T]
			if !ismissing(𝑗T) && Ω[𝑗T] # only continue if inside Ω
				_, _, k′ = C𝑗T.I
				d = abs(zt[k] - zt[k′])
				pushTmixingvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗T, κV, a, d, V)
			end
		end
	end
	return 𝑖s, 𝑗s, Tvals
end