function localdensityslope(ρ, lon, lat, Z, I, gridtopology, dir)
    distances = verticalfacetriadgroupdistances(lon, lat, Z, I, gridtopology, dir)
    localρ = verticalfacetriadgroupvalues(ρ, I, gridtopology, dir)
    return localtriadderivative(localρ, distances)
end
function globaldensityslope(ρ, gridmetrics, indices, dir)
    (; lon, lat, Z, gridtopology) = gridmetrics
    (; Lwet, C) = indices
    S = fill(NaN, size(ρ))
    for 𝑖 in eachindex(Lwet)
        L𝑖 = Lwet[𝑖]
        C𝑖 = C[L𝑖]
        S[C𝑖] = localdensityslope(ρ, lon, lat, Z, C𝑖, gridtopology, dir)
    end
    return S
end
function localpotentialdensityslope(gsw_rho, so, ct, lon, lat, Z, I, gridtopology, dir)
    zref = Z[I]
    distances = verticalfacetriadgroupdistances(lon, lat, Z, I, gridtopology, dir)
    localso = verticalfacetriadgroupvalues(so, I, gridtopology, dir)
    localct = verticalfacetriadgroupvalues(ct, I, gridtopology, dir)
    localρθ = gsw_rho.(localso, localct, zref)
    return localtriadderivative(localρθ, distances)
end
function globalpotentialdensityslope(gsw_rho, so, ct, gridmetrics, indices, dir)
    (; lon, lat, Z, gridtopology) = gridmetrics
    (; Lwet, C) = indices
    S = fill(NaN, size(so))
    for 𝑖 in eachindex(Lwet)
        L𝑖 = Lwet[𝑖]
        C𝑖 = C[L𝑖]
        S[C𝑖] = localpotentialdensityslope(gsw_rho, so, ct, lon, lat, Z, C𝑖, gridtopology, dir)
    end
    return S
end


"""
    bolus_GM_velocity(ρ, gridmetrics, indices; κGM = 600, maxslope = 0.01)

Compute Gent-McWilliams bolus velocities from density field ρ.

The function calculates the eddy-induced bolus velocities that represent the effect
of mesoscale eddies on large-scale ocean circulation. This implementation uses density
slopes computed via triad derivatives and applies slope limiting and tapering for
numerical stability.

# Arguments
- `ρ`: 3D density field (kg/m³)
- `gridmetrics`: Named tuple containing grid information (lon, lat, Z, gridtopology, etc.)
- `indices`: Named tuple containing index information (Lwet, C, etc.)
- `κGM`: Gent-McWilliams diffusivity coefficient (m²/s, default: 600)
- `maxslope`: Maximum density slope magnitude (default: 0.01)

# Returns
- `u_GM, v_GM`: Bolus velocities in x and y directions (m/s)

# Mathematical Formulation
The GM bolus velocity is computed as:
1. Calculate density slopes: Sᵢ = ∂ρ/∂x, Sⱼ = ∂ρ/∂y
2. Apply slope limiting: S = clamp(S, -maxslope, maxslope)
3. Apply taper function: taper = 0.5*(1 + tanh((Sc - |S|)/Sd))
4. Compute bolus velocities: u_GM = -κGM * ∂/∂z (taper * Sᵢ)

# References
- Gent, P. R., and J. C. McWilliams, 1990: Isopycnal mixing in ocean circulation models.
  J. Phys. Oceanogr., 20, 150-155.
"""
function bolus_GM_velocity(ρ, gridmetrics, indices; κGM = 600, maxslope = 0.01)
    # Calculate density slopes in i and j directions using the existing slope functions
    # These represent the horizontal density gradients: Sᵢ ≈ ∂ρ/∂x, Sⱼ ≈ ∂ρ/∂y
    Sᵢ = globaldensityslope(ρ, gridmetrics, indices, Icoord())
    Sⱼ = globaldensityslope(ρ, gridmetrics, indices, Jcoord())

    # Apply slope limiting to prevent numerical instability
    # This caps the magnitude of density slopes to physically reasonable values
    Sᵢ = clamp.(Sᵢ, -maxslope, maxslope)
    Sⱼ = clamp.(Sⱼ, -maxslope, maxslope)

    # Apply taper function to smoothly reduce mixing in regions of strong slopes
    # This prevents excessive mixing where isopycnals are steeply sloped
    Sc = 0.004  # Critical slope - where taper starts to reduce mixing
    Sd = 0.001  # Taper width - controls how quickly mixing is reduced
    slope_magnitude = sqrt.(Sᵢ.^2 + Sⱼ.^2)
    taper = @. 0.5 * (1 + tanh((Sc - slope_magnitude) / Sd))
    Sᵢ = taper .* Sᵢ
    Sⱼ = taper .* Sⱼ

    # Calculate bolus velocities from density slopes
    # u_GM = -κGM * ∂/∂z (Sᵢ), v_GM = -κGM * ∂/∂z (Sⱼ)
    # The negative sign follows the conventional formulation where
    # eddy-induced transport opposes the mean density gradients
    u_GM = -globalverticaldyadderivative(κGM .* Sᵢ, gridmetrics, indices)
    v_GM = -globalverticaldyadderivative(κGM .* Sⱼ, gridmetrics, indices)

    return u_GM, v_GM
end
