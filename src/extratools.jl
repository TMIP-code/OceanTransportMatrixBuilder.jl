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

    lump_and_spray(wet3D, volume; di=3, dj=4)

will lump:
- every 3 cells in the x direction (lon),
- every 4 cells in the y direction (lat),
- and every 1 cell (default) in the z direction (depth).

Expected grid arrangement is OCIM2 like, i.e.,
lat × lon × depth.

You can also provide a vector of indices instead of a scalar
in order to customize the coarsening. For example,

    lump_and_spray(wet3D, volume; di=[1 1 1 1 1 1 1 2 2 2 ... n])

will lump all the boxes in the x dimension that are marked
with a 1, then all those marked with a 2, and so on.
Confusing? Let me give a practical example. Say you want to
coarsen the ACCESS grid in the vertical to match the OCIM grid.
The you can do

    [~, zidx] = min(abs(ACCESSgrd.zt - OCIMgrd.zt'), [], 1)

to get the OCIM2 z-indices that are closest to the
ACCESS z-indices. And then you can pass it to this function via

    lump_and_spray(wet3D, volume; dk=zidx)

TODO: Fix these docs that have not been fully translated from MATLAB version.
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





















