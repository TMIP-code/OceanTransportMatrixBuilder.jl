
"""
    buildTadv(; u, modelgrid, indices, ρ)

Build the advection operator Tadv.
"""
function buildTadv(; u, modelgrid, indices, ρ)
    # default ρ = 1035 kg/m^3 is the value originally used by Chamberlain et al. (2019)

	@info "Building Tadv"
	𝑖s, 𝑗s, Tvals = upwind_advection_operator_sparse_entries(; u, modelgrid, indices, ρ)

    N = indices.N

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
    (; Lwet) = indices

	# Wet mask for mixed layer diffusivity
	Ω = replace(reshape(zt, 1, 1, length(zt)) .< mlotst, missing=>false)[Lwet]

	@info "Building TκVML "
	𝑖s, 𝑗s, Tvals = vertical_diffusion_operator_sparse_entries(; modelgrid, indices, κV = κVML, Ω)

    N = indices.N

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

	TκVdeep = sparse(𝑖s, 𝑗s, Tvals, N, N)

	return TκVdeep

end

"""
    transportmatrix(; u, mlotst, modelgrid, indices, ρ, κH, κVML, κVdeep, Tadv, TκH, TκVML, TκVdeep)

Build the transport matrix, i.e., the flux-divergence operator T = Tadv + TκH + TκVML + TκVdeep,
and check divergence and mass conservation.
"""
function transportmatrix(; u, mlotst, modelgrid, indices, ρ,
		κH = 500.0, # m^2/s,
		κVML = 0.1, # m^2/s,
		κVdeep = 1e-5, # m^2/s,
		Tadv = buildTadv(; u, modelgrid, indices, ρ),
		TκH = buildTκH(; modelgrid, indices, ρ, κH),
		TκVML = buildTκVML(; mlotst, modelgrid, indices, κVML),
		TκVdeep = buildTκVdeep(; mlotst, modelgrid, indices, κVdeep),
	)

	@info "Building T"

	T = Tadv + TκH + TκVML + TκVdeep

	return (; T, Tadv, TκH, TκVML, TκVdeep)
end






# Some personal notes
# u and v are water mass transports, in kg/s
# We can convert them to m^3/s using ρ = 1035 kg/m^3
# but I don't think it matters much since since we normalize
# by mass (or v3D) to build T.
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

# There are three ways to index:
# - Cartesian indices (i,j,k)
# - Linear index Li
# - Wet linear index 𝑖 (that's the one we want to record for the matrix)
# So to fill T[𝑖,𝑗] -ϕ[𝑖→𝑗] / m[𝑖], I need be able to convert, in sequence:
# 𝑖 -> (i,j,k) -> neihghbour (i′,j′,k′) -> 𝑗
# The first 2 conversions are straightforard.
# For the last one, I make a 3D array filled with the wet linear indices

"""
    upwind_advection_operator_sparse_entries(; u, modelgrid, indices, ρ)

Return the sparse (i, j, v) for the upwind advection operator Tadv.
"""
function upwind_advection_operator_sparse_entries(; u, modelgrid, indices, ρ)

    # Unpack model grid
    (; v3D,) = modelgrid
    # Unpack indices
    (; wet3D, Lwet, Lwet3D, C) = indices


	𝑖s, 𝑗s, Tvals = Int[], Int[], Float64[]
	nxyz = size(wet3D)
    nx, ny, _ = nxyz

    @showprogress for 𝑖 in eachindex(Lwet)
		Li = Lwet[𝑖]
		i, j, k = C[Li].I
		m𝑖 = v3D[i,j,k] * ρ
		# From West
		ϕ = u.west[i,j,k]
		if ϕ > 0
			i′ = mod1(i - 1, nx)
			𝑗 = Lwet3D[i′,j,k]
			ismissing(𝑗) && @show(i, j, k, i′)
			m𝑗 = v3D[i′,j,k] * ρ
			pushTadvectionvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗, ϕ, m𝑖, m𝑗)
		end
		# From East
		ϕ = u.east[i,j,k]
		if ϕ < 0
			i′ = mod1(i + 1, nx)
			𝑗 = Lwet3D[i′,j,k]
			m𝑗 = v3D[i′,j,k] * ρ
			pushTadvectionvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗, -ϕ, m𝑖, m𝑗)
		end
		# From South
		ϕ = u.south[i,j,k]
		if ϕ > 0
			j′ = j - 1
			𝑗 = Lwet3D[i,j′,k]
			m𝑗 = v3D[i,j′,k] * ρ
			pushTadvectionvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗, ϕ, m𝑖, m𝑗)
		end
		# From North (Special case with north bipole)
		ϕ = u.north[i,j,k]
		if ϕ < 0
			if j == ny
				j′ = j
				i′ = nx - i + 1
			else
				j′ = j + 1
				i′ = i
			end
			𝑗 = Lwet3D[i′,j′,k]
			m𝑗 = v3D[i′,j′,k] * ρ
			pushTadvectionvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗, -ϕ, m𝑖, m𝑗)
		end
		# From Bottom
		ϕ = u.bottom[i,j,k]
		if ϕ > 0
			k′ = k + 1
			𝑗 = Lwet3D[i,j,k′]
			m𝑗 = v3D[i,j,k′] * ρ
			pushTadvectionvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗, ϕ, m𝑖, m𝑗)
		end
		# From Top
		ϕ = u.top[i,j,k]
		if ϕ < 0 && k > 1 # Evaporation/precipitation -> no change to χ
			k′ = k - 1
			𝑗 = Lwet3D[i,j,k′]
			m𝑗 = v3D[i,j,k′] * ρ
			pushTadvectionvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗, -ϕ, m𝑖, m𝑗)
		end
	end
	return 𝑖s, 𝑗s, Tvals
end


"""
    pushTadvectionvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗, ϕ, m𝑖, m𝑗)

Pushes the sparse indices and values into (𝑖s, 𝑗s, Tvals) corresponding to the j→i advection.
"""
function pushTadvectionvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗, ϕ, m𝑖, m𝑗)
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
    (; v3D, edge_length_2D, distance_to_edge_2D, DZT3d) = modelgrid
    # Unpack indices
    (; wet3D, Lwet, Lwet3D, C) = indices

	𝑖s, 𝑗s, Tvals = Int[], Int[], Float64[]
	nxyz = size(wet3D)
    nx, ny, _ = nxyz

    @showprogress for 𝑖 in eachindex(Lwet)
        ΩH[𝑖] || continue # only continue if inside ΩH
		Li = Lwet[𝑖]
		i, j, k = C[Li].I
		V = v3D[i,j,k]
		# From West
		iW, jW = mod1(i - 1, nx), j
		𝑗W = Lwet3D[iW, jW, k]
        # (𝑖 == 𝑗W) && @show(i, j, iW, jW)
        if !ismissing(𝑗W)
            # I take the minimum area from both dirs (through which mixing goes through)
            aij = verticalfacearea(edge_length_2D, DZT3d, i, j, k, :west)
            aji = verticalfacearea(edge_length_2D, DZT3d, iW, jW, k, :east)
            a = min(aij, aji)
            # I take the mean distance from both dirs
            dij = horizontalcentroiddistance(distance_to_edge_2D, i, j, iW, jW, :west)
            dji = horizontalcentroiddistance(distance_to_edge_2D, iW, jW, i, j, :east)
            d = (dij + dji) / 2
			pushTmixingvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗W, κH, a, d, V)
		end
        # From East
		iE, jE = mod1(i + 1, nx), j
		𝑗E = Lwet3D[iE, jE, k]
        # (𝑖 == 𝑗E) && @show(i, j, iE, jE)
        if !ismissing(𝑗E)
            aij = verticalfacearea(edge_length_2D, DZT3d, i, j, k, :east)
            aji = verticalfacearea(edge_length_2D, DZT3d, iE, jE, k, :west)
            a = min(aij, aji)
            dij = horizontalcentroiddistance(distance_to_edge_2D, i, j, iE, jE, :east)
            dji = horizontalcentroiddistance(distance_to_edge_2D, iE, jE, i, j, :west)
            d = (dij + dji) / 2
			pushTmixingvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗E, κH, a, d, V)
		end
        # From South
        if j > 1
            iS, jS = i, j - 1
            𝑗S = Lwet3D[iS, jS, k]
            # (𝑖 == 𝑗S) && @show(i, j, iS, jS)
            if !ismissing(𝑗S)
                aij = verticalfacearea(edge_length_2D, DZT3d, i, j, k, :south)
                aji = verticalfacearea(edge_length_2D, DZT3d, iS, jS, k, :north)
                a = min(aij, aji)
                dij = horizontalcentroiddistance(distance_to_edge_2D, i, j, iS, jS, :south)
                dji = horizontalcentroiddistance(distance_to_edge_2D, iS, jS, i, j, :north)
                d = (dij + dji) / 2
                pushTmixingvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗S, κH, a, d, V)
            end
        end
        # From North
        # Note that the opposite direction (oppdir) is still north at j == ny
        (iN, jN, oppdir) = (j == ny) ? (nx - i + 1, j, :north) : (i, j + 1, :south)
        𝑗N = Lwet3D[iN, jN, k]
        # (𝑖 == 𝑗N) && @show(i, j, iN, jN)
        if !ismissing(𝑗N)
            aij = verticalfacearea(edge_length_2D, DZT3d, i, j, k, :north)
            aji = verticalfacearea(edge_length_2D, DZT3d, iN, jN, k, oppdir)
            a = min(aij, aji)
            dij = horizontalcentroiddistance(distance_to_edge_2D, i, j, iN, jN, :north)
            dji = horizontalcentroiddistance(distance_to_edge_2D, iN, jN, i, j, oppdir)
            d = (dij + dji) / 2
            pushTmixingvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗N, κH, a, d, V)
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
    (; v3D, area2D, zt) = modelgrid
    # Unpack indices
    (; wet3D, Lwet, Lwet3D, C) = indices

    𝑖s, 𝑗s, Tvals = Int[], Int[], Float64[]
	nxyz = size(wet3D)
    _, _, nz = nxyz

    @showprogress for 𝑖 in eachindex(Lwet)
        Ω[𝑖] || continue # only continue if inside Ω
		Li = Lwet[𝑖]
		i, j, k = C[Li].I
		V = v3D[i,j,k]
        a = area2D[i,j]
		# From Bottom
        if k < nz
            k′ = k + 1
            𝑗B = Lwet3D[i,j,k′]
            if !ismissing(𝑗B) && Ω[𝑗B] # only continue if inside Ω
                d = abs(zt[k] - zt[k′])
                pushTmixingvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗B, κV, a, d, V)
            end
        end
		# From Top
        if k > 1
            k′ = k - 1
            𝑗T = Lwet3D[i,j,k′]
            if !ismissing(𝑗T) && Ω[𝑗T] # only continue if inside Ω
                d = abs(zt[k] - zt[k′])
                pushTmixingvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗T, κV, a, d, V)
            end
        end
	end
	return 𝑖s, 𝑗s, Tvals
end