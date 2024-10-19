
####################################
# Forward and backward derivatives #
####################################

# Derivative functions (forward and backward)
∂(::Icoord) = (∂ᵢ₊, ∂ᵢ₋)
∂(::Jcoord) = (∂ⱼ₊, ∂ⱼ₋)
∂(::Kcoord) = (∂ₖ₊, ∂ₖ₋)

# local horizontal derivatives ("local" means "at index I")
function ∂ₕ₊(χ, lon, lat, gridtopology, I, dir)
    J = idx₊₁(dir)(I, gridtopology)
    d = horizontaldistance(lon, lat, I, J)
    return (χ[J] - χ[I]) / d
end
function ∂ₕ₋(χ, lon, lat, gridtopology, I, dir)
    J = idx₋₁(dir)(I, gridtopology)
    d = horizontaldistance(lon, lat, I, J)
    return (χ[I] - χ[J]) / d
end
∂ᵢ₊(args...) = ∂ₕ₊(args..., Icoord())
∂ᵢ₋(args...) = ∂ₕ₋(args..., Icoord())
∂ⱼ₊(args...) = ∂ₕ₊(args..., Jcoord())
∂ⱼ₋(args...) = ∂ₕ₋(args..., Jcoord())

# local vertical derivatives
function ∂ₖ₊(χ, Z, gridtopology, I)
    J = k₊₁(I, gridtopology)
    d = verticaldistance(Z, I, J)
    return (χ[J] - χ[I]) / d
end
function ∂ₖ₋(χ, Z, gridtopology, I)
    J = k₋₁(I, gridtopology)
    d = verticaldistance(Z, I, J)
    return (χ[I] - χ[J]) / d
end

# global forward derivative ("global" means "for all indices")
function ∂ₕ₊(χ, gridmetrics, dir)
    # χ = χ |> Array # TODO: remove this line if not required
    (; lon, lat, gridtopology) = gridmetrics
    ∂ₕ₊χ = zeros(size(χ))
    for I in CartesianIndices(∂ₕ₊χ)
        ∂ₕ₊χ[I] = ∂ₕ₊(χ, lon, lat, gridtopology, I, dir)
    end
    return ∂ₕ₊χ
end
# global backward derivative
function ∂ₕ₋(χ, gridmetrics, dir)
    # χ = χ |> Array # TODO: remove this line if not required
    (; lon, lat, gridtopology) = gridmetrics
    ∂ₕ₋χ = zeros(size(χ))
    for I in CartesianIndices(∂ₕ₋χ)
        ∂ₕ₋χ[I] = ∂ₕ₋(χ, lon, lat, gridtopology, I, dir)
    end
    return ∂ₕ₋χ
end
# global forard derivative in the k direction
function ∂ₖ₊(χ, gridmetrics)
    # χ = χ |> Array # TODO: remove this line if not required
    (; Z3D, thkcello, gridtopology) = gridmetrics
    ∂ₖ₊χ = zeros(size(χ))
    for I in CartesianIndices(∂ₖ₊χ)
        ∂ₖ₊χ[I] = ∂ₖ₊(χ, Z3D, gridtopology, I)
    end
    return ∂ₖ₊χ
end
# backward derivative in the k direction
function ∂ₖ₋(χ, gridmetrics)
    # χ = χ |> Array # TODO: remove this line if not required
    (; Z3D, thkcello, gridtopology) = gridmetrics
    ∂ₖ₋χ = zeros(size(χ))
    for I in CartesianIndices(χ)
        ∂ₖ₋χ[I] = ∂ₖ₋(χ, Z3D, gridtopology, I)
    end
    return ∂ₖ₋χ
end






# TODO: Check that I am properly dealing with unevenly spaced grids
# The commented code below could help with that, as it provides the weights
# for a finite difference stencil for any arrangement of grid points
# """
#     stencil(x, x₀, m)

# returns the stencil weights wₖ such that

# f⁽ᵐ⁾(x₀) ≈ ∑ₖ₌₁ⁿ wₖ f(xₖ)

# from https://discourse.julialang.org/t/generating-finite-difference-stencils/85876/5
# """
# function stencil(x::AbstractVector{<:Real}, x₀::Real, m::Integer)
#     ℓ = 0:length(x)-1
#     m in ℓ || throw(ArgumentError("order $m ∉ $ℓ"))
#     A = @. (x' - x₀)^ℓ / factorial(ℓ)
#     return A \ (ℓ .== m) # vector of weights w
# end


