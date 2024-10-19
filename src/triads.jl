




###### TRIADS ######

"""
TriadGroup

I'm sort of following the POP manual cardinal orientation here for naming convention
because it looks simpler (except I use uppercase).

```
NW───N───NE
│    │    │
W────C────E
│    │    │
SW───S───SE
```

Almost quoted from said manual (I modified center "C" and use uppercase):

> The central tracer point (T-point) is denoted by "C"
    and the surrounding T-points by "N", "S", etc. (this is
    for notation only, points "N" and "S" are displaced in the vertical direction.)
"""
struct CenteredTriadGroupValues <: FieldVector{5, Float64}
    C
    E
    W
    N
    S
end
struct CenteredTriadGroupDistances <: FieldVector{4, Float64}
    CW
    CE
    CN
    CS
end
# TODO fix centeredtriadgroup functions
function centeredtriadgroup(χ, I, lon, lat, Z, gridtopology, dir)
    E = idx₊₁(dir)(I, gridtopology)
    W = idx₋₁(dir)(I, gridtopology)
    N = k₋₁(I, gridtopology)
    S = k₊₁(I, gridtopology)
    vals = CenteredTriadGroupValues(
        getindexornan(χ, I),
        getindexornan(χ, E),
        getindexornan(χ, W),
        getindexornan(χ, N),
        getindexornan(χ, S),
    )
    distances = CenteredTriadGroupDistances(
        horizontaldistance(lon, lat, I, E),
        horizontaldistance(lon, lat, I, W),
        verticaldistance(Z, I, N),
        verticaldistance(Z, I, S),
    )
    return vals, distances
end
function localtriadderivative(vals::CenteredTriadGroupValues, distances::CenteredTriadGroupDistances)
    Δvals = (
        (vals.E - vals.C) / distances.CE,
        (vals.C - vals.W) / distances.CW,
        (vals.N - vals.C) / distances.CN,
        (vals.C - vals.S) / distances.CS,
    )
    weights = (!isnan(Δval) for Δval in Δvals)
    return sum(w * Δval for (w, Δval) in zip(weights, Δvals)) / sum(weights)
end


struct VerticalFaceTriadGroupValues <: FieldVector{6, Float64}
    C
    N
    S
    E
    NE
    SE
end
struct VerticalFaceTriadGroupDistances <: FieldVector{5, Float64}
    CN
    CS
    CE
    ENE
    ESE
end
function verticalfacetriadgroupindices(I, gridtopology, dir)
    N = k₋₁(I, gridtopology)
    S = k₊₁(I, gridtopology)
    E = idx₊₁(dir)(I, gridtopology)
    NE = k₋₁(E, gridtopology)
    SE = k₊₁(E, gridtopology)
    return N, S, E, NE, SE
end
function verticalfacetriadgroupvalues(χ, I, gridtopology, dir)
    N, S, E, NE, SE = verticalfacetriadgroupindices(I, gridtopology, dir)
    return VerticalFaceTriadGroupValues(
        getindexornan(χ, I),
        getindexornan(χ, N),
        getindexornan(χ, S),
        getindexornan(χ, E),
        getindexornan(χ, NE),
        getindexornan(χ, SE),
    )
end
function verticalfacetriadgroupdistances(lon, lat, Z, I, gridtopology, dir)
    N, S, E, NE, SE = verticalfacetriadgroupindices(I, gridtopology, dir)
    return VerticalFaceTriadGroupDistances(
        verticaldistance(Z, I, N),
        verticaldistance(Z, I, S),
        horizontaldistance(lon, lat, I, E),
        verticaldistance(Z, E, NE),
        verticaldistance(Z, E, SE),
    )
end

function localtriadderivative(vals::VerticalFaceTriadGroupValues, distances::VerticalFaceTriadGroupDistances)
    Δvals = (
        (vals.N - vals.C) / distances.CN,
        (vals.C - vals.S) / distances.CS,
        (vals.E - vals.C) / distances.CE,
        (vals.NE - vals.E) / distances.ENE,
        (vals.E - vals.SE) / distances.ESE,
    )
    weights = (!isnan(Δval) for Δval in Δvals)
    return sum(w * Δval for (w, Δval) in zip(weights, Δvals)) / sum(weights)
end
function globalverticalfacetriadderivative(χ, gridmetrics, indices, dir)
    (; lon, lat, Z3D, gridtopology) = gridmetrics
    ∂χ = fill(NaN, size(χ))
    (; Lwet, C) = indices
    for 𝑖 in eachindex(Lwet)
        L𝑖 = Lwet[𝑖]
        C𝑖 = C[L𝑖]
        vals = verticalfacetriadgroupvalues(χ, C𝑖, gridtopology, dir)
        distances = verticalfacetriadgroupdistances(lon, lat, Z3D, C𝑖, gridtopology, dir)
        ∂χ[C𝑖] = localtriadderivative(vals, distances)
    end
    return ∂χ
end
