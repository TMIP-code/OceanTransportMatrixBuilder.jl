# There are three ways to index:
# - Cartesian indices (i,j,k)
# - Linear index L𝑖
# - Wet linear index 𝑖 (that's the one we want to record for the matrix)
# So to fill T[𝑖,𝑗] -ϕ[𝑖→𝑗] / m[𝑖], I need be able to convert, in sequence:
# 𝑖 -> (i,j,k) -> neihghbour (i′,j′,k′) -> 𝑗
# The first 2 conversions are straightforard.
# For the last one, I make a 3D array filled with the wet linear indices

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

"""
    buildTadv(; ϕ, gridmetrics, indices, ρ, upwind = true)

Build the advection operator Tadv.
"""
function buildTadv(; ϕ, gridmetrics, indices, ρ, upwind = true)
    # default ρ = 1035 kg/m^3 is the value originally used by Chamberlain et al. (2019)

    @debug "Building Tadv"
    𝑖s, 𝑗s, Tvals = advection_operator_sparse_entries(ϕ, gridmetrics, indices, ρ; upwind)

    N = indices.N

    any(isnan.(Tvals)) && error("Tadv contains NaNs.")

    Tadv = sparse(𝑖s, 𝑗s, Tvals, N, N)

    return Tadv
end

"""
    buildTκH(; gridmetrics, indices, ρ, κH)

Build the horizontal diffusivity operator TκH.
"""
function buildTκH(; gridmetrics, indices, ρ, κH)

    N = indices.N

    # Wet mask for horizontal diffusivity
    ΩH = trues(N)

    @debug "Building TκH"
    𝑖s, 𝑗s, Tvals = horizontal_diffusion_operator_sparse_entries(; gridmetrics, indices, κH, ΩH)

    any(isnan.(Tvals)) && error("TκH contains NaNs.")

    TκH = sparse(𝑖s, 𝑗s, Tvals, N, N)

    return TκH
end


"""
    buildTκVML(; mlotst, gridmetrics, indices, κVML)

Build the mixed layer diffusivity operator TκVML.
"""
function buildTκVML(; mlotst, gridmetrics, indices, κVML)

    # Unpack model grid
    (; zt) = gridmetrics

    # Unpack indices
    (; Lwet, N) = indices

    mlotst = mlotst |> Array # to prevent slow getindexornan for lazily loaded data?

    # Wet mask for mixed layer diffusivity
    Ω = replace(reshape(zt, 1, 1, length(zt)) .< mlotst, missing => false)[Lwet]

    @debug "Building TκVML "
    𝑖s, 𝑗s, Tvals = vertical_diffusion_operator_sparse_entries(; gridmetrics, indices, κV = κVML, Ω)

    any(isnan.(Tvals)) && error("TκVML contains NaNs.")

    TκVML = sparse(𝑖s, 𝑗s, Tvals, N, N)

    return TκVML
end


"""
    buildTκVdeep(; mlotst, gridmetrics, indices, κVdeep)

Build the deep diffusivity operator TκVdeep.
"""
function buildTκVdeep(; mlotst, gridmetrics, indices, κVdeep)

    N = indices.N

    # κVdeep is a bit of a misnomer: It should be κVBG for "background".
    # And its mask is entire ocean, naturally.
    Ω = trues(N)

    @debug "Building TκVdeep"
    𝑖s, 𝑗s, Tvals = vertical_diffusion_operator_sparse_entries(; gridmetrics, indices, κV = κVdeep, Ω)

    any(isnan.(Tvals)) && error("TκVdeep contains NaNs.")

    TκVdeep = sparse(𝑖s, 𝑗s, Tvals, N, N)

    return TκVdeep

end

"""
    buildTGM(ρ, gridmetrics, indices; κGM = 600, maxslope = 0.01)

Build the Gent-McWilliams transport operator from density field.

# Arguments
- `ρ`: 3D density field (kg/m³) for GM calculation
- `gridmetrics`: Named tuple containing grid information
- `indices`: Named tuple containing index information
- `κGM`: Gent-McWilliams diffusivity coefficient (m²/s, default: 600)
- `maxslope`: Maximum density slope magnitude (default: 0.01)

# Returns
- `TGM`: Sparse matrix operator for GM bolus transport
"""
function buildTGM(ρ, gridmetrics, indices; κGM = 600, maxslope = 0.01)
    # Calculate GM bolus velocities from density field
    u_GM, v_GM = bolus_GM_velocity(ρ, gridmetrics, indices; κGM, maxslope)

    # For GM, we need to convert velocities to fluxes. Since GM velocities are
    # defined on the same grid as the density field, we can use the gridmetrics
    # to compute the fluxes directly.
    # Extract necessary grid information
    (; thkcello, edge_length_2D, v3D) = gridmetrics
    mean_ρ = mean(skipmissing(ρ))

    # Create flux fields with same size as velocity fields
    ϕᵢ_GM = zeros(size(u_GM))
    ϕⱼ_GM = zeros(size(v_GM))

    # Convert velocities to fluxes: flux = velocity * density * area * thickness
    # We need to handle the grid topology properly here
    for 𝑖 in eachindex(u_GM)
        i, j, k = indices.C[𝑖].I
        # East flux: u_GM * mean_ρ * edge_length * min_thickness
        ϕᵢ_GM[𝑖] = u_GM[𝑖] * mean_ρ * edge_length_2D[:east][i, j] * thkcello[i, j, k]
        # North flux: v_GM * mean_ρ * edge_length * min_thickness
        ϕⱼ_GM[𝑖] = v_GM[𝑖] * mean_ρ * edge_length_2D[:north][i, j] * thkcello[i, j, k]
    end

    # Build advection operator from GM fluxes
    # This creates a transport operator that represents the eddy-induced
    # bolus transport as an advection-like process
    return buildTadv(;
        ϕ = (east = ϕᵢ_GM, north = ϕⱼ_GM),
        gridmetrics, indices,
        ρ = mean_ρ,  # Use constant density for GM flux conversion
        upwind = true
    )
end

"""
    transportmatrix(; ϕ, mlotst, gridmetrics, indices, ρ, κH, κVML, κVdeep, Tadv, TκH, TκVML, TκVdeep, upwind)

Build the transport matrix, i.e., the flux-divergence operator T = Tadv + TκH + TκVML + TκVdeep,
and check divergence and mass conservation.
"""
function transportmatrix(;
        ϕ, mlotst, gridmetrics, indices, ρ,
        ρ_GM = nothing,  # Density field for Gent-McWilliams parameterization (optional)
        κH = 500.0, # m^2/s,
        κVML = 0.1, # m^2/s,
        κVdeep = 1.0e-5, # m^2/s,
        κGM = 600.0, # m^2/s - Gent-McWilliams diffusivity coefficient
        maxslope = 0.01, # Maximum density slope for GM slope limiting
        Tadv = nothing,
        TκH = nothing,
        TκVML = nothing,
        TκVdeep = nothing,
        TGM = nothing,  # Gent-McWilliams transport operator (optional)
        upwind = true,
    )

    # Build standard transport operators
    isnothing(Tadv) && (Tadv = buildTadv(; ϕ, gridmetrics, indices, ρ, upwind))
    isnothing(TκH) && (TκH = buildTκH(; gridmetrics, indices, ρ, κH))
    isnothing(TκVML) && (TκVML = buildTκVML(; mlotst, gridmetrics, indices, κVML))
    isnothing(TκVdeep) && (TκVdeep = buildTκVdeep(; mlotst, gridmetrics, indices, κVdeep))

    # Build Gent-McWilliams operator if density field is provided
    if !isnothing(ρ_GM)
        isnothing(TGM) && (TGM = buildTGM(ρ_GM, gridmetrics, indices; κGM, maxslope))
    else
        TGM = spzeros(indices.N, indices.N)  # Zero operator if no GM
    end

    @debug "Building transport matrix with GM contribution"

    T = Tadv + TκH + TκVML + TκVdeep + TGM

    return (; T, Tadv, TκH, TκVML, TκVdeep, TGM)
end


function preallocate_sparse_entries(sizehint)
    𝑖s = Int64[]
    𝑗s = Int64[]
    Tvals = Float64[]
    sizehint!(𝑖s, sizehint)
    sizehint!(𝑗s, sizehint)
    sizehint!(Tvals, sizehint)
    return 𝑖s, 𝑗s, Tvals
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

# For a centered scheme, where the concentration on the cell face is used,
# I still have the matrix operating as
# 	∂χ[i] = -T[i,j] * χ[j]    units:   1/s * mol/kg (or mol/m^3)
# But the centered mass transfer is (regardless of the sign of ϕ)
# 	∂χ[i] = ϕ[j→i] * (χ[j] + χ[i]) / 2m[i]     units    kg[j]/s * mol/kg[j] / kg[i] = mol/kg[i]
# and should also incur the opposite mass tendency at j:
# 	∂χ[j] = -ϕ[j→i] * (χ[j] + χ[i]) / 2m[j]     units    kg[j]/s * mol/kg[j] / kg[j] = mol/kg[j]
# Thus the matrix term should be constructed as
# 	T[i,j] = -ϕ[j→i] / 2m[i]             units = kg s⁻¹ / kg = s⁻¹
# and ϕ[j→i] / 2m[j] should be added to the diagonal T[j,j].
# So essentially just divide ϕ by 2 compared to upwind, but don't branch on the sign of ϕ.

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
    return nothing
end


# Is this mass conserving?
# Only if
#   v[i] * T[i,i] + v[j] * T[j,i] ≈ 0
#   = vi ϕ[i→j] / m[i] + vj -ϕ[i→j] / m[j]
#   = ϕ[i→j] / ρ[i] - ϕ[i→j] / ρ[j]
# i.e., iff
#   ρ[i] = ρ[j]
# So I must use the mean density between facing cells.

"""
    advection_operator_sparse_entries(ϕ, gridmetrics, indices, ρ; upwind = true)

Return the sparse (i, j, v) for the upwind advection operator Tadv.
"""
function advection_operator_sparse_entries(ϕ, gridmetrics, indices, ρ::Number; upwind = true)
    # If ρ is a scalar, broadcast it to the gridmetrics size
    ρ = fill(ρ, size(gridmetrics.v3D))
    return advection_operator_sparse_entries(ϕ, gridmetrics, indices, ρ; upwind)
end
function advection_operator_sparse_entries(ϕ, gridmetrics, indices, ρ; upwind = true)

    # Unpack model grid
    (; v3D, gridtopology) = gridmetrics
    # Unpack indices
    (; Lwet, Lwet3D, C, N) = indices

    any(isnan, ρ[Lwet]) && error("ρ contains NaNs")

    𝑖s, 𝑗s, Tvals = preallocate_sparse_entries(6N) # 6 directions for upwind advection

    for 𝑖 in eachindex(Lwet)
        L𝑖 = Lwet[𝑖]
        C𝑖 = C[L𝑖]
        i, j, k = C𝑖.I
        v𝑖 = v3D[C𝑖]
        ρ𝑖 = ρ[C𝑖]
        # From West
        ϕwest = upwind ? max(ϕ.west[C𝑖], 0) : ϕ.west[C𝑖] / 2
        if (ϕwest > 0) || (ϕwest < 0)
            C𝑗 = i₋₁(C𝑖, gridtopology)
            𝑗 = Lwet3D[C𝑗]
            v𝑗 = v3D[C𝑗]
            ρ𝑗 = ρ[C𝑗]
            pushTadvectionvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗, ϕwest, ρ𝑖, ρ𝑗, v𝑖, v𝑗)
        end
        # From East
        ϕeast = upwind ? min(ϕ.east[C𝑖], 0) : ϕ.east[C𝑖] / 2
        if (ϕeast > 0) || (ϕeast < 0)
            C𝑗 = i₊₁(C𝑖, gridtopology)
            𝑗 = Lwet3D[C𝑗]
            v𝑗 = v3D[C𝑗]
            ρ𝑗 = ρ[C𝑗]
            pushTadvectionvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗, -ϕeast, ρ𝑖, ρ𝑗, v𝑖, v𝑗)
        end
        # From South
        ϕsouth = upwind ? max(ϕ.south[C𝑖], 0) : ϕ.south[C𝑖] / 2
        if (ϕsouth > 0) || (ϕsouth < 0)
            C𝑗 = j₋₁(C𝑖, gridtopology)
            𝑗 = Lwet3D[C𝑗]
            v𝑗 = v3D[C𝑗]
            ρ𝑗 = ρ[C𝑗]
            pushTadvectionvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗, ϕsouth, ρ𝑖, ρ𝑗, v𝑖, v𝑗)
        end
        # From North (Special case with north bipole)
        ϕnorth = upwind ? min(ϕ.north[C𝑖], 0) : ϕ.north[C𝑖] / 2
        if (ϕnorth > 0) || (ϕnorth < 0)
            C𝑗 = j₊₁(C𝑖, gridtopology)
            𝑗 = Lwet3D[C𝑗]
            v𝑗 = v3D[C𝑗]
            ρ𝑗 = ρ[C𝑗]
            pushTadvectionvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗, -ϕnorth, ρ𝑖, ρ𝑗, v𝑖, v𝑗)
        end
        # From Bottom
        ϕbottom = upwind ? max(ϕ.bottom[C𝑖], 0) : ϕ.bottom[C𝑖] / 2
        if (ϕbottom > 0) || (ϕbottom < 0)
            C𝑗 = k₊₁(C𝑖, gridtopology)
            𝑗 = Lwet3D[C𝑗]
            v𝑗 = v3D[C𝑗]
            ρ𝑗 = ρ[C𝑗]
            pushTadvectionvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗, ϕbottom, ρ𝑖, ρ𝑗, v𝑖, v𝑗)
        end
        # From Top
        ϕtop = upwind ? min(ϕ.top[C𝑖], 0) : ϕ.top[C𝑖] / 2
        if (k > 1) && ((ϕtop > 0) || (ϕtop < 0)) # Evaporation/precipitation -> no change to χ
            C𝑗 = k₋₁(C𝑖, gridtopology)
            𝑗 = Lwet3D[C𝑗]
            v𝑗 = v3D[C𝑗]
            ρ𝑗 = ρ[C𝑗]
            pushTadvectionvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗, -ϕtop, ρ𝑖, ρ𝑗, v𝑖, v𝑗)
        end
    end
    return 𝑖s, 𝑗s, Tvals
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
    horizontal_diffusion_operator_sparse_entries(; gridmetrics, indices, κH, ΩH)

Return the sparse (i, j, v) for the horizontal diffusion operator TκH.
"""
function horizontal_diffusion_operator_sparse_entries(; gridmetrics, indices, κH, ΩH)

    # Unpack model grid
    (; v3D, edge_length_2D, lon, lat, thkcello, gridtopology, distance_to_neighbour_2D) = gridmetrics
    # Unpack indices
    (; wet3D, Lwet, Lwet3D, C, N) = indices

    𝑖s, 𝑗s, Tvals = preallocate_sparse_entries(8N) # 2 × 4 directions for horizontal diffusion

    ny = size(wet3D, 2) # Should not be needed once oppdir is dealt by topology functions

    for 𝑖 in eachindex(Lwet)
        ΩH[𝑖] || continue # only continue if inside ΩH
        L𝑖 = Lwet[𝑖]
        C𝑖 = C[L𝑖]
        C𝑖srf = horizontalindex(C𝑖)
        i, j, k = C𝑖.I
        V = v3D[C𝑖]
        # From West
        C𝑗W = i₋₁(C𝑖, gridtopology)
        if !isnothing(C𝑗W)
            𝑗W = Lwet3D[C𝑗W]
            if !ismissing(𝑗W) && ΩH[𝑗W]
                iW, jW, _ = C𝑗W.I
                # (𝑖 == 𝑗W) && @show(i, j, iW, jW)
                # I take the minimum area from both dirs (through which mixing goes through)
                aij = verticalfacearea(edge_length_2D, thkcello, i, j, k, :west)
                aji = verticalfacearea(edge_length_2D, thkcello, iW, jW, k, :east)
                a = min(aij, aji)
                d = distance_to_neighbour_2D[:west][C𝑖srf]
                pushTmixingvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗W, κH, a, d, V)
            end
        end
        # From East
        C𝑗E = i₊₁(C𝑖, gridtopology)
        if !isnothing(C𝑗E)
            𝑗E = Lwet3D[C𝑗E]
            if !ismissing(𝑗E) && ΩH[𝑗E]
                iE, jE, _ = C𝑗E.I
                # (𝑖 == 𝑗E) && @show(i, j, iE, jE)
                aij = verticalfacearea(edge_length_2D, thkcello, i, j, k, :east)
                aji = verticalfacearea(edge_length_2D, thkcello, iE, jE, k, :west)
                a = min(aij, aji)
                d = distance_to_neighbour_2D[:east][C𝑖srf]
                pushTmixingvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗E, κH, a, d, V)
            end
        end
        # From South
        C𝑗S = j₋₁(C𝑖, gridtopology)
        if !isnothing(C𝑗S)
            𝑗S = Lwet3D[C𝑗S]
            if !ismissing(𝑗S) && ΩH[𝑗S]
                iS, jS, _ = C𝑗S.I
                # (𝑖 == 𝑗S) && @show(i, j, iS, jS)
                aij = verticalfacearea(edge_length_2D, thkcello, i, j, k, :south)
                aji = verticalfacearea(edge_length_2D, thkcello, iS, jS, k, :north)
                a = min(aij, aji)
                d = distance_to_neighbour_2D[:south][C𝑖srf]
                pushTmixingvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗S, κH, a, d, V)
            end
        end
        # From North
        C𝑗N = j₊₁(C𝑖, gridtopology)
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
                d = distance_to_neighbour_2D[:north][C𝑖srf]
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
    push!(𝑖s, 𝑖)
    push!(𝑗s, 𝑖)
    push!(Tvals, Tval)
    push!(𝑖s, 𝑖)
    push!(𝑗s, 𝑗)
    push!(Tvals, -Tval)
    return nothing
end


function vertical_diffusion_operator_sparse_entries(; gridmetrics, indices, κV, Ω)

    # Unpack model grid
    (; v3D, area2D, zt, gridtopology) = gridmetrics
    # Unpack indices
    (; wet3D, Lwet, Lwet3D, C, N) = indices

    𝑖s, 𝑗s, Tvals = preallocate_sparse_entries(4 * N) # 2 × 2 directions for vertical diffusion

    nxyz = size(wet3D)
    _, _, nz = nxyz

    for 𝑖 in eachindex(Lwet)
        Ω[𝑖] || continue # only continue if inside Ω
        L𝑖 = Lwet[𝑖]
        C𝑖 = C[L𝑖]
        i, j, k = C𝑖.I
        V = v3D[C𝑖]
        a = area2D[i, j]
        # From Bottom
        C𝑗B = k₊₁(C𝑖, gridtopology)
        if !isnothing(C𝑗B)
            𝑗B = Lwet3D[C𝑗B]
            if !ismissing(𝑗B) && Ω[𝑗B] # only continue if inside Ω
                _, _, k′ = C𝑗B.I
                d = abs(zt[k] - zt[k′])
                pushTmixingvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗B, κV, a, d, V)
            end
        end
        # From Top
        C𝑗T = k₋₁(C𝑖, gridtopology)
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
