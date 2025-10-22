#########
# DYADS #
#########

struct CenteredDyadGroup end
struct VerticalFaceDyadGroup end
# struct HorizontalFaceDyadGroup end

# Just like the TRIADS below, but for a single variable
# The idea is to take the centered difference if possible (no NaNs),
# or the forward or the backward derivative otherwise.
# My understanding is that this is just what the triads do
# by taking the average of all possible combinations.


"""
DyadGroup

Same naming convention as for the triads.

```
     N
     │
W────C────E
     │
     S
```
"""
struct VerticalDyadGroupValues <: FieldVector{3, Float64}
    C
    N
    S
end
struct VerticalDyadGroupDistances <: FieldVector{2, Float64}
    CN
    CS
end
function verticaldyadgroupindices(I, gridtopology)
    N = k₋₁(I, gridtopology)
    S = k₊₁(I, gridtopology)
    return N, S
end
function verticaldyadgroupvalues(χ, I, gridtopology)
    N, S = verticaldyadgroupindices(I, gridtopology)
    return VerticalDyadGroupValues(
        getindexornan(χ, I),
        getindexornan(χ, N),
        getindexornan(χ, S),
    )
end
function verticaldyadgroupdistances(Z, I, gridtopology)
    N, S = verticaldyadgroupindices(I, gridtopology)
    return VerticalDyadGroupDistances(
        verticaldistance(Z, I, N),
        verticaldistance(Z, I, S),
    )
end
function localdyadderivative(vals::VerticalDyadGroupValues, distances::VerticalDyadGroupDistances)
    Δvals = (
        (vals.N - vals.C) / distances.CN,
        (vals.C - vals.S) / distances.CS,
    )
    weights = (!isnan(Δval) for Δval in Δvals)
    return sum(w * Δval for (w, Δval) in zip(weights, Δvals)) / sum(weights)
end
function globalverticaldyadderivative(χ, gridmetrics, indices)
    (; Z3D, gridtopology) = gridmetrics
    ∂χ = fill(NaN, size(χ))
    (; Lwet, C) = indices
    for 𝑖 in eachindex(Lwet)
        L𝑖 = Lwet[𝑖]
        C𝑖 = C[L𝑖]
        vals = verticaldyadgroupvalues(χ, C𝑖, gridtopology)
        distances = verticaldyadgroupdistances(Z3D, C𝑖, gridtopology)
        ∂χ[C𝑖] = localdyadderivative(vals, distances)
    end
    return ∂χ
end
