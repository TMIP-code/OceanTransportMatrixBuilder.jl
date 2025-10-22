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
    bolus_GM_velocity(σ, gridmetrics; κGM = 600, maxslope = 0.01)

Returns the bolus velocity field due to the Gent-McWilliams parameterization,
computed from the neutral density field `σ` (or potential density, ρθ in kg/m³).

Note: This is experimental at this stage.
"""
function bolus_GM_velocity(ρ, gridmetrics, indices; κGM = 600, maxslope = 0.01)
    # TODO: implement with neutral density (or potential density, ρθ in kg/m³)
    # function bolus_GM_velocity(so, thetao, gridmetrics; κGM = 600, maxslope = 0.01)

    # Sᵢ = globalpotentialdensityslope(gsw_rho, so, ct, gridmetrics, indices, Icoord())
    # Sⱼ = globalpotentialdensityslope(gsw_rho, so, ct, gridmetrics, indices, Jcoord())
    Sᵢ = globalverticalfacetriadderivative(ρ, gridmetrics, indices, Icoord())
    Sⱼ = globalverticalfacetriadderivative(ρ, gridmetrics, indices, Jcoord())

    # cap the slope
    Sᵢ = clamp.(Sᵢ, -maxslope, maxslope)
    Sⱼ = clamp.(Sⱼ, -maxslope, maxslope)

    Sc = 0.004
    Sd = 0.001
    # that should not work since Sᵢ and Sⱼ are not colocated
    taper = @. 0.5 * (1 + tanh((Sc - sqrt(Sᵢ^2 + Sⱼ^2)) / Sd))
    Sᵢ = taper .* Sᵢ
    Sⱼ = taper .* Sⱼ

    # TODO Include Z into gridmetrics as Z3D so it's clear what it is
    # TODO Apply vertical dyad derivative to the density Slopes
    # TODO Add some tests on the derivatives. Maybe set some fields like
    # χ(x,y,z) = x^2 + y^2 + z^2 and check the derivatives at some points. Not sure. Think more about it.

    # TODO clean up this code and maybe split it into indexing (including dyad and triad indexing)
    # and derivatives (including forward, backward, centered, and dyad and triad derivatives)

    # Take the vertical derivative of the density slope in x
    u = globalverticaldyadderivative(κGM .* Sᵢ, gridmetrics, indices)
    v = globalverticaldyadderivative(κGM .* Sⱼ, gridmetrics, indices)

    return u, v
end
