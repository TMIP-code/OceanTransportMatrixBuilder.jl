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





function areafun3D(vlon, vlat, vlev_bnds_or_thkcello, C, Lwet, dir)
    area1D = [verticalfacearea(vlon, vlat, vlev_bnds_or_thkcello, I.I[1], I.I[2], I.I[3], dir) for I in C[Lwet]]
    area3D = fill(NaN, size(C))
    area3D[Lwet] .= area1D
    return area3D
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









