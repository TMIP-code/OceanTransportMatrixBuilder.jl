"""
    LUMP, SPRAY = lump_and_spray(wet3D, volume; di=2, dj=2, dk=1)

returns the LUMP and SPRAY matrices that go with the
coarsened grid (also returned).

To get the coarsened vector from fine vector x, use

    LUMP * x

To get the coarsened operator from fine operator T, use

    LUMP * T * SPRAY

The di, dj, and dk options (default di = dj = 2 and dk = 1)
set the size of the coarsening.
Let me explain with an example:

    lump_and_spray(grd, di=3, dj=4)

will lump:
- every 3 cells in the x direction (lon),
- every 4 cells in the y direction (lat),
- and every 1 cell (default) in the z direction (depth).

Expected grid arrangement is OCIM2 like, i.e.,
lat × lon × depth.

You can also provide a vector of indices instead of a scalar
in order to customize the coarsening. For example,

    lump_and_spray(grd, di=[1 1 1 1 1 1 1 2 2 2 ... n])

will lump all the boxes in the x dimension that are marked
with a 1, then all those marked with a 2, and so on.
Confusing? Let me give a practical example. Say you want to
coarsen the ACCESS grid in the vertical to match the OCIM grid.
The you can do

    [~, zidx] = min(abs(ACCESSgrd.zt - OCIMgrd.zt'), [], 1)

to get the OCIM2 z-indices that are closest to the
ACCESS z-indices. And then you can pass it to this function via

    lump_and_spray(grd, dk=zidx)
"""
lump_and_spray(wet3D, volume; di=2, dj=2, dk=1) = _lump_and_spray(wet3D, volume, di, dj, dk)
function _lump_and_spray(wet3D, volume, di::Int, dj::Int, dk::Int)
    # grd wet array and vector and sizes
    nx, ny, nz = size(wet3D)
    # Convert scalar lumping options into lumping indices vector
    # So that syntax like `di=2` works
    vi = repeat(1:nx, inner=di)[1:nx]
    vj = repeat(1:ny, inner=dj)[1:ny]
    vk = repeat(1:nz, inner=dk)[1:nz]

    return _lump_and_spray(wet3D, volume, vi, vj, vk)
end
function _lump_and_spray(wet3D, volume, vi, vj, vk)
    wet = wet3D[:]
    nx, ny, nz = size(wet3D)
    # Create lumping matrices in each dimension
    LUMPx = sparse(vi, 1:nx, true)
    LUMPy = sparse(vj, 1:ny, true)
    LUMPz = sparse(vk, 1:nz, true)
    # kron each dimension to build whole LUMP matrix
    LUMP = kron(LUMPz, kron(LUMPy, LUMPx))

    # Find wet points in coarsened grid
    wet_c = LUMP * wet .> 0
    nx_c = vi[end]
    ny_c = vj[end]
    nz_c = vk[end]
    wet3D_c = fill(false, nx_c, ny_c, nz_c)
    wet3D_c[wet_c] .= true

    # Extract only indices of wet grd points
    LUMP = LUMP[wet_c, wet]

    # Make the LUMP operator volume-conserving
    # by volume integrating on the right and dividing by the coarse
    # volume on the left
    volume_c = LUMP * volume
    LUMP = sparse(Diagonal(1 ./ volume_c)) * LUMP * sparse(Diagonal(volume))

    # The SPRAY operator just copies the values back
    # so it is sinply 1's with the transposed sparsity structure
    SPRAY = copy(LUMP')
    SPRAY.nzval .= 1

    return LUMP, SPRAY, wet3D_c, volume_c
end



"""
	edges = edges_from_midpoints(midpoints, lims)

returns the edges given the midpoints and limits
using backslash. That is, the least squares solution of

	midpoints = (edges[1:end-1] + edges[2:end]) / 2
	edges[1] = lims[1]
	edges[end] = lims[2]
    diff()
"""
function edges_from_midpoints(midpoints, lims)
    # TODO this creates noisy diffs so maybe worth
    # including a smoothing parameter for the diff?
    # Also not sure I actually really need these at all
    N = length(midpoints)
    I = sparse(Diagonal(ones(N)))
    I_left = [I zeros(N)]
    I_right = [zeros(N) I]
    A = (I_left + I_right) / 2
    M = [A; 1 zeros(1, N); zeros(1, N) 1]
	return M \ [midpoints; collect(lims)]
end




# Some personal notes
# u and v are water mass transports, in kg/s
# We can convert them to m^3/s using ρ = 1035 kg/m^3
# but I don't think it matters much since since we normalize
# by mass (or volume3D) to build T.
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

function pushTadvectionvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗, ϕ, m𝑖, m𝑗)
	push!(𝑖s, 𝑖)
	push!(𝑗s, 𝑗)
	push!(Tvals, -ϕ / m𝑖)
	push!(𝑖s, 𝑗)
	push!(𝑗s, 𝑗)
	push!(Tvals, ϕ / m𝑗)
end

function build_upwind_advection_operator(uwest::Array{T,3}, ueast, vsouth, vnorth, wbottom, wtop, Lwet, Lwet3D, C, volume3D, ρ) where T
	𝑖s, 𝑗s, Tvals = Int[], Int[], nonmissingtype(T)[]
	nxyz = size(uwest)
    nx, ny, nz = nxyz
    @showprogress for 𝑖 in eachindex(Lwet)
		Li = Lwet[𝑖]
		i, j, k = C[Li].I
		m𝑖 = volume3D[i,j,k] * ρ
		# From West
		ϕ = uwest[i,j,k]
		if ϕ > 0
			i′ = mod1(i - 1, nx)
			𝑗 = Lwet3D[i′,j,k]
			m𝑗 = volume3D[i′,j,k] * ρ
			pushTadvectionvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗, ϕ, m𝑖, m𝑗)
		end
		# From East
		ϕ = ueast[i,j,k]
		if ϕ < 0
			i′ = mod1(i + 1, nx)
			𝑗 = Lwet3D[i′,j,k]
			m𝑗 = volume3D[i′,j,k] * ρ
			pushTadvectionvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗, -ϕ, m𝑖, m𝑗)
		end
		# From South
		ϕ = vsouth[i,j,k]
		if ϕ > 0
			j′ = j - 1
			𝑗 = Lwet3D[i,j′,k]
			m𝑗 = volume3D[i,j′,k] * ρ
			pushTadvectionvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗, ϕ, m𝑖, m𝑗)
		end
		# From North (Special case with north bipole)
		ϕ = vnorth[i,j,k]
		if ϕ < 0
			if j == ny
				j′ = j
				i′ = nx - i + 1
			else
				j′ = j + 1
				i′ = i
			end
			𝑗 = Lwet3D[i′,j′,k]
			m𝑗 = volume3D[i′,j′,k] * ρ
			pushTadvectionvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗, -ϕ, m𝑖, m𝑗)
		end
		# From Bottom
		ϕ = wbottom[i,j,k]
		if ϕ > 0
			k′ = k + 1
			𝑗 = Lwet3D[i,j,k′]
			m𝑗 = volume3D[i,j,k′] * ρ
			pushTadvectionvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗, ϕ, m𝑖, m𝑗)
		end
		# From Top
		ϕ = wtop[i,j,k]
		if ϕ < 0 && k > 1 # Evaporation/precipitation -> no change to χ
			k′ = k - 1
			𝑗 = Lwet3D[i,j,k′]
			m𝑗 = volume3D[i,j,k′] * ρ
			pushTadvectionvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗, -ϕ, m𝑖, m𝑗)
		end
	end
	return 𝑖s, 𝑗s, Tvals
end



# View form top for vlon and vlat vertices
#    4 ────┐ 3
#          │
#    1 ────┘ 2
function vertexindices(dir)
    dir == :south ? (1, 2) :
    dir == :east ? (2, 3) :
    dir == :north ? (3, 4) :
    dir == :west ? (1, 4) :
    error()
end
vertexpoint(vlon, vlat, i, j, vertexidx) = (vlon[vertexidx,i,j], vlat[vertexidx,i,j])
function verticalfacewidth(vlon, vlat, i, j, dir)
    a, b = vertexindices(dir)
    A = vertexpoint(vlon, vlat, i, j, a)
    B = vertexpoint(vlon, vlat, i, j, b)
    haversine(A, B)
end
verticalfacewidth(edge_length_2D, i, j, dir) = edge_length_2D[dir][i, j]

function verticalfacearea(vlon, vlat, lev_bnds_or_thkcello, i, j, k, dir)
    height = cellthickness(lev_bnds_or_thkcello, i, j, k)
    width = verticalfacewidth(vlon, vlat, i, j, dir)
    height * width
end
function verticalfacearea(edge_length_2D, lev_bnds_or_thkcello, i, j, k, dir)
    height = cellthickness(lev_bnds_or_thkcello, i, j, k)
    width = verticalfacewidth(edge_length_2D, i, j, dir)
    height * width
end

cellthickness(lev_bnds::Matrix, i, j, k) = abs(lev_bnds[2,k] - lev_bnds[1,k])
cellthickness(thkcello, i, j, k) = thkcello[i, j, k]


function centroid2edgedistance(lon, lat, vlon, vlat, i, j, dir)
    a, b = vertexindices(dir)
    C = (lon[i, j], lat[i, j])
    A = vertexpoint(vlon, vlat, i, j, a)
    B = vertexpoint(vlon, vlat, i, j, b)
    M = (A .+ B) ./ 2
    haversine(C, M)
end


function horizontalcentroiddistance(distance_to_edge_2D, iA, jA, iB, jB, dir)
    if dir == :south
        distance_to_edge_2D[:south][iA, jA] + distance_to_edge_2D[:north][iB, jB]
    elseif dir == :north
        if jA == jB # if on North wall, A and B both connect via north
            distance_to_edge_2D[:north][iA, jA] + distance_to_edge_2D[:north][iB, jB]
        else
            distance_to_edge_2D[:north][iA, jA] + distance_to_edge_2D[:south][iB, jB]
        end
    elseif dir == :west
        distance_to_edge_2D[:west][iA, jA] + distance_to_edge_2D[:east][iB, jB]
    elseif dir == :east
        distance_to_edge_2D[:east][iA, jA] + distance_to_edge_2D[:west][iB, jB]
    end
end



"""
as2D(x, wet3D)

`x` must be a vector with `length(x) == sum(wet3D[:,:,1])`
"""
function as2D(x, wet3D)
    x2D = fill(NaN, size(wet3D)[1:2])
    @views x2D[wet3D[:,:,1]] .= x
    x2D
end

"""
as3D(x, wet3D)

`x` must be a vector with `length(x) == sum(wet3D)`
"""
function as3D(x, wet3D)
    x3D = fill(NaN, size(wet3D))
    @views x3D[wet3D] .= x
    x3D
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


function areafun3D(vlon, vlat, vlev_bnds_or_thkcello, C, Lwet, dir)
    area1D = [verticalfacearea(vlon, vlat, vlev_bnds_or_thkcello, I.I[1], I.I[2], I.I[3], dir) for I in C[Lwet]]
    area3D = fill(NaN, size(C))
    area3D[Lwet] .= area1D
    return area3D
end


function build_horizontal_diffusivity_operator(κ::T, Lwet, Lwet3D, C, volume3D, Ω,
    lon, lat, vertices_longitude, vertices_latitude, lev_bnds_or_thkcello) where T

    vlon = vertices_longitude
    vlat = vertices_latitude
    dirs = (:south, :east, :north, :west)
    edge_length_2D = Dict(dir=>[verticalfacewidth(vlon, vlat, 𝑖.I[1], 𝑖.I[2], dir) for 𝑖 in C[:,:,1]] for dir in dirs)
    distance_to_edge_2D = Dict(dir=>[centroid2edgedistance(lon, lat, vlon, vlat, 𝑖.I[1], 𝑖.I[2], dir) for 𝑖 in C[:,:,1]] for dir in dirs)

    return build_horizontal_diffusivity_operator(κ::T, Lwet, Lwet3D, C, volume3D, Ω,
        edge_length_2D, distance_to_edge_2D, lev_bnds_or_thkcello)
end

function build_horizontal_diffusivity_operator(κ::T, Lwet, Lwet3D, C, volume3D, Ω,
    edge_length_2D, distance_to_edge_2D, lev_bnds_or_thkcello) where T
	𝑖s, 𝑗s, Tvals = Int[], Int[], T[]
	nxyz = size(Lwet3D)
    nx, ny, nz = nxyz

    @showprogress for 𝑖 in eachindex(Lwet)
        Ω[𝑖] || continue # only continue if inside Ω
		Li = Lwet[𝑖]
		i, j, k = C[Li].I
		V = volume3D[i,j,k]
		# From West
		iW, jW = mod1(i - 1, nx), j
		𝑗W = Lwet3D[iW, jW, k]
        # (𝑖 == 𝑗W) && @show(i, j, iW, jW)
        if !ismissing(𝑗W)
            # I take the minimum area from both dirs (through which mixing goes through)
            aij = verticalfacearea(edge_length_2D, lev_bnds_or_thkcello, i, j, k, :west)
            aji = verticalfacearea(edge_length_2D, lev_bnds_or_thkcello, iW, jW, k, :east)
            a = min(aij, aji)
            # I take the mean distance from both dirs
            dij = horizontalcentroiddistance(distance_to_edge_2D, i, j, iW, jW, :west)
            dji = horizontalcentroiddistance(distance_to_edge_2D, iW, jW, i, j, :east)
            d = (dij + dji) / 2
			pushTmixingvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗W, κ, a, d, V)
		end
        # From East
		iE, jE = mod1(i + 1, nx), j
		𝑗E = Lwet3D[iE, jE, k]
        # (𝑖 == 𝑗E) && @show(i, j, iE, jE)
        if !ismissing(𝑗E)
            aij = verticalfacearea(edge_length_2D, lev_bnds_or_thkcello, i, j, k, :east)
            aji = verticalfacearea(edge_length_2D, lev_bnds_or_thkcello, iE, jE, k, :west)
            a = min(aij, aji)
            dij = horizontalcentroiddistance(distance_to_edge_2D, i, j, iE, jE, :east)
            dji = horizontalcentroiddistance(distance_to_edge_2D, iE, jE, i, j, :west)
            d = (dij + dji) / 2
			pushTmixingvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗E, κ, a, d, V)
		end
        # From South
        if j > 1
            iS, jS = i, j - 1
            𝑗S = Lwet3D[iS, jS, k]
            # (𝑖 == 𝑗S) && @show(i, j, iS, jS)
            if !ismissing(𝑗S)
                aij = verticalfacearea(edge_length_2D, lev_bnds_or_thkcello, i, j, k, :south)
                aji = verticalfacearea(edge_length_2D, lev_bnds_or_thkcello, iS, jS, k, :north)
                a = min(aij, aji)
                dij = horizontalcentroiddistance(distance_to_edge_2D, i, j, iS, jS, :south)
                dji = horizontalcentroiddistance(distance_to_edge_2D, iS, jS, i, j, :north)
                d = (dij + dji) / 2
                pushTmixingvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗S, κ, a, d, V)
            end
        end
        # From North
        # Note that the opposite direction (oppdir) is still north at j == ny
        (iN, jN, oppdir) = (j == ny) ? (nx - i + 1, j, :north) : (i, j + 1, :south)
        𝑗N = Lwet3D[iN, jN, k]
        # (𝑖 == 𝑗N) && @show(i, j, iN, jN)
        if !ismissing(𝑗N)
            aij = verticalfacearea(edge_length_2D, lev_bnds_or_thkcello, i, j, k, :north)
            aji = verticalfacearea(edge_length_2D, lev_bnds_or_thkcello, iN, jN, k, oppdir)
            a = min(aij, aji)
            dij = horizontalcentroiddistance(distance_to_edge_2D, i, j, iN, jN, :north)
            dji = horizontalcentroiddistance(distance_to_edge_2D, iN, jN, i, j, oppdir)
            d = (dij + dji) / 2
            pushTmixingvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗N, κ, a, d, V)
        end
	end

	return 𝑖s, 𝑗s, Tvals
end




# Note that in the case of horizontal neighbors, a = v / Δz
function build_vertical_diffusivity_operator(κ::T, Lwet, Lwet3D, C, volume3D, area2D, Ω, lev) where T
	𝑖s, 𝑗s, Tvals = Int[], Int[], T[]
	nxyz = size(Lwet3D)
    _, _, nz = nxyz

    @showprogress for 𝑖 in eachindex(Lwet)
        Ω[𝑖] || continue # only continue if inside Ω
		Li = Lwet[𝑖]
		i, j, k = C[Li].I
		V = volume3D[i,j,k]
        a = area2D[i,j]
		# From Bottom
        if k < nz
            k′ = k + 1
            𝑗B = Lwet3D[i,j,k′]
            if !ismissing(𝑗B) && Ω[𝑗B] # only continue if inside Ω
                d = abs(lev[k] - lev[k′])
                pushTmixingvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗B, κ, a, d, V)
            end
        end
		# From Top
        if k > 1
            k′ = k - 1
            𝑗T = Lwet3D[i,j,k′]
            if !ismissing(𝑗T) && Ω[𝑗T] # only continue if inside Ω
                d = abs(lev[k] - lev[k′])
                pushTmixingvalues!(𝑖s, 𝑗s, Tvals, 𝑖, 𝑗T, κ, a, d, V)
            end
        end
	end
	return 𝑖s, 𝑗s, Tvals
end



"""
tests if diagonal elements are > 0 and off-diagonal are < 0.
"""
function isdivergence(T)
    diagT = sparse(Diagonal(T))

    @show posdiag = all(diagT.nzval .> 0)
    @show negoffdiag = all((T - diagT).nzval .< 0)

    posdiag & negoffdiag
end






# Figuring out the "topology" of the grid
# Right now I think I can only use standard grids or tripolar grids.
function gridtopology(volumedata)
    !haskey(volumedata, "vertices_longitude")
    if !haskey(volumedata, "vertices_longitude") ||
        !haskey(volumedata, "vertices_latitude")
        return nothing
    end

    vlonnorth = volumedata["vertices_longitude"][:,:,end]
    vlatnorth = volumedata["vertices_latitude"][:,:,end]

    # The longitude of the two north poles
    lonpoles = unique(vlonnorth[end,:])
    length(lonpoles) ≠ 2 && return nothing

    # The i index (lon) of the box[i,end] that connects with box[1,end]
    NEpoints = collect(zip(vlonnorth[3,:,end], vlatnorth[3,:,end]))
    NWpoints = collect(zip(vlonnorth[4,:,end], vlatnorth[4,:,end]))
    ilon1 = only(findall(([NEpoints[1]] .== NWpoints) .& ([NWpoints[1]] .== NEpoints)))
    nx = length(NEpoints)

    return i -> mod1(ilon1 - i + 1, nx)

end









