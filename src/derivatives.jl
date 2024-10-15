



"""
    function δᵢ(g, )

Returns the discrete difference
"""


# forward derivative in the i direction for centered data (A-grid)
function ∂ᵢ₊(χ, modelgrid)
    χ = χ |> Array
    (; lon, lat, gridtype) = modelgrid
    ∂ᵢ₊χ = zeros(size(χ))
    for I in CartesianIndices(lon)
        P = (lon[I], lat[I])
        Iᵢ₊₁ = i₊₁(I, gridtype)
        Pᵢ₊₁ = (lon[Iᵢ₊₁], lat[Iᵢ₊₁])
        d = haversine(P, Pᵢ₊₁)
        ∂ᵢ₊χ[I, :] = (χ[Iᵢ₊₁, :] - χ[I, :]) / d
    end
    return ∂ᵢ₊χ
end
# backward derivative in the i direction for centered data (A-grid)
function ∂ᵢ₋(χ, modelgrid)
    χ = χ |> Array
    (; lon, lat, gridtype) = modelgrid
    ∂ᵢ₋χ = zeros(size(χ))
    for I in CartesianIndices(lon)
        P = (lon[I], lat[I])
        Iᵢ₋₁ = i₋₁(I, gridtype)
        Pᵢ₋₁ = (lon[Iᵢ₋₁], lat[Iᵢ₋₁])
        d = haversine(P, Pᵢ₋₁)
        ∂ᵢ₋χ[I, :] = (χ[I, :] - χ[Iᵢ₋₁, :]) / d
    end
    return ∂ᵢ₋χ
end
# forward derivative in the j direction
function ∂ⱼ₊(χ, modelgrid)
    χ = χ |> Array
    (; lon, lat, gridtype) = modelgrid
    ∂ⱼ₊χ = zeros(size(χ))
    for I in CartesianIndices(lon)
        P = (lon[I], lat[I])
        Iᵢ₊₁ = j₊₁(I, gridtype)
        Pᵢ₊₁ = (lon[Iᵢ₊₁], lat[Iᵢ₊₁])
        d = haversine(P, Pᵢ₊₁)
        ∂ⱼ₊χ[I, :] = (χ[Iᵢ₊₁, :] - χ[I, :]) / d
    end
    return ∂ⱼ₊χ
end
# backward derivative in the j direction
function ∂ⱼ₋(χ, modelgrid)
    χ = χ |> Array
    (; lon, lat, gridtype) = modelgrid
    ∂ⱼ₋χ = zeros(size(χ))
    for I in CartesianIndices(lon)
        P = (lon[I], lat[I])
        Iᵢ₋₁ = j₋₁(I, gridtype)
        Pᵢ₋₁ = (lon[Iᵢ₋₁], lat[Iᵢ₋₁])
        d = haversine(P, Pᵢ₋₁)
        ∂ⱼ₋χ[I, :] = (χ[I, :] - χ[Iᵢ₋₁, :]) / d
    end
    return ∂ⱼ₋χ
end
# derivative in the k direction
function ∂ₖ₊(χ, modelgrid)
    χ = χ |> Array
    (; zt, DZT3d, gridtype) = modelgrid
    ∂ₖ₊χ = zeros(size(χ))
    for I in CartesianIndices(χ)
        Iᵢ₊₁ = k₊₁(I, gridtype)
        isnothing(Iᵢ₊₁) && continue
        h = (DZT3d[I] + DZT3d[Iᵢ₊₁]) / 2
        ∂ₖ₊χ[I] = (χ[Iᵢ₊₁] - χ[I]) / h
    end
    return ∂ₖ₊χ
end
# backward derivative in the k direction
function ∂ₖ₋(χ, modelgrid)
    χ = χ |> Array
    (; zt, DZT3d, gridtype) = modelgrid
    ∂ₖ₋χ = zeros(size(χ))
    for I in CartesianIndices(χ)
        Iᵢ₋₁ = k₋₁(I, gridtype)
        isnothing(Iᵢ₋₁) && continue
        h = (DZT3d[I] + DZT3d[Iᵢ₋₁]) / 2
        ∂ₖ₋χ[I] = (χ[I] - χ[Iᵢ₋₁]) / h
    end
    return ∂ₖ₋χ
end

# interpolation in the i direction
function itpₖ₊(χ, modelgrid)
    (; gridtype) = modelgrid
    χ = χ |> Array
    itpₖ₊χ = zeros(size(χ))
    for I in CartesianIndices(χ)
        Iᵢ₊₁ = k₊₁(I, gridtype)
        isnothing(Iᵢ₊₁) && continue
        itpₖ₊χ[I] = (χ[Iᵢ₊₁] + χ[I]) / 2
    end
    return itpₖ₊χ
end
function itpₖ₋(χ, modelgrid)
    (; gridtype) = modelgrid
    χ = χ |> Array
    itpₖ₋χ = zeros(size(χ))
    for I in CartesianIndices(χ)
        Iᵢ₋₁ = k₋₁(I, gridtype)
        isnothing(Iᵢ₋₁) && continue
        itpₖ₋χ[I] = (χ[Iᵢ₋₁] + χ[I]) / 2
    end
    return itpₖ₋χ
end
function itpᵢ₊(χ, modelgrid)
    (; gridtype) = modelgrid
    χ = χ |> Array
    itpᵢ₊χ = zeros(size(χ))
    for I in CartesianIndices(χ)
        Iᵢ₊₁ = i₊₁(I, gridtype)
        itpᵢ₊χ[I] = (χ[Iᵢ₊₁] + χ[I]) / 2
    end
    return itpᵢ₊χ
end
function itpⱼ₊(χ, modelgrid)
    (; gridtype) = modelgrid
    χ = χ |> Array
    itpⱼ₊χ = zeros(size(χ))
    for I in CartesianIndices(χ)
        Iᵢ₊₁ = j₊₁(I, gridtype)
        itpⱼ₊χ[I] = (χ[Iᵢ₊₁] + χ[I]) / 2
    end
    return itpⱼ₊χ
end


function bolus_GM_velocity(σ, modelgrid; κGM = 600, maxslope = 0.01)
    # σ is neutral density (or potential density, ρθ in kg/m³)
    σ = replace(σ, missing => 0, NaN => 0)
    # σ is cell centered (A-grid) at [i, j, k]
    ∂ᵢσ = ∂ᵢ₊(σ, modelgrid) # at [i+½, j, k]
    ∂ⱼσ = ∂ⱼ₊(σ, modelgrid) # at [i, j+½, k]
    ∂ₖσ = ∂ₖ₊(σ, modelgrid) # at [i, j, k+½]
    # Interpolate to matching B-grid corners
    ∂ᵢσ = itpₖ₊(∂ᵢσ, modelgrid) # at [i+½, j, k+½]
    ∂ⱼσ = itpₖ₊(∂ⱼσ, modelgrid) # at [i, j+½, k+½]
    ∂ₖσᵢ = itpᵢ₊(∂ₖσ, modelgrid) # at [i+½, j, k+½]
    ∂ₖσⱼ = itpⱼ₊(∂ₖσ, modelgrid) # at [i, j+½, k+½]
    # Compute the slope of the density field
    Sᵢ = ∂ᵢσ ./ ∂ₖσᵢ # at [i+½, j, k+½]
    Sⱼ = ∂ⱼσ ./ ∂ₖσⱼ # at [i, j+½, k+½]vghtk88VK f
    # cap the slope
    Sᵢ = clamp.(Sᵢ, -maxslope, maxslope)
    Sⱼ = clamp.(Sⱼ, -maxslope, maxslope)
    # Take the vertical derivative of the density slope in x
    u = +∂ₖ₋(κGM * Sᵢ, modelgrid) # at [i+½, j, k] # +sign because z is positive downwards
    v = +∂ₖ₋(κGM * Sⱼ, modelgrid) # at [i, j+½, k] # +sign because z is positive downwards
    u = replace(u, NaN => 0)
    v = replace(v, NaN => 0)
    return u, v
end


"""
    dᵢ₊(lon, lat)

Returns the distance between neighbors in the i direction.
"""
function dᵢ₊(lon, lat)
    d = zeros(size(lon))
    for I in CartesianIndices(lon)
        P = (lon[I], lat[I])
        Iᵢ₊₁ = i₊₁(I, modelgrid)
        Pᵢ₊₁ = (lon[Iᵢ₊₁], lat[Iᵢ₊₁])
        d[I] = haversine(P, Pᵢ₊₁)
    end
    return d
end
function dⱼ₊(lon, lat)
    d = zeros(size(lon))
    for I in CartesianIndices(lon)
        P = (lon[I], lat[I])
        Iᵢ₊₁ = j₊₁(I, modelgrid)
        Pᵢ₊₁ = (lon[Iᵢ₊₁], lat[Iᵢ₊₁])
        d[I] = haversine(P, Pᵢ₊₁)
    end
    return d
end



# function ∂ᵢ(χ, modelgrid, shift=0,  d=:f)
#     (;lon, lat, gridtype) = modelgrid
#     m = 1
#     ∂ᵢχ = zeros(size(χ))
#     # Draw a line along coordinate i and compute the distances
#     d = dᵢ₊(lon, lat)
#     for I in CartesianIndices(lon)
#         movingI = SVector(ishift(I, gridtype, xxx))
#         i, j = I.I
#         x = SVector(0, d[i])
#         w = stencil(x, 0, χ)
#     end
#     for (j, k) in axes(χ, 2), axes(χ, 3)
#         x = [χ[i, j, k] for i in axes(χ, 1)]
#         ∂ᵢχ[:, j, k] = stencil(x, shift, m)
#     end
#     return ∂ᵢχ
# end


"""
    stencil(x, x₀, m)

returns the stencil weights wₖ such that

f⁽ᵐ⁾(x₀) ≈ ∑ₖ₌₁ⁿ wₖ f(xₖ)

from https://discourse.julialang.org/t/generating-finite-difference-stencils/85876/5
"""
function stencil(x::AbstractVector{<:Real}, x₀::Real, m::Integer)
    ℓ = 0:length(x)-1
    m in ℓ || throw(ArgumentError("order $m ∉ $ℓ"))
    A = @. (x' - x₀)^ℓ / factorial(ℓ)
    return A \ (ℓ .== m) # vector of weights w
end