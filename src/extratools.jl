"""
    LUMP, SPRAY, v_c = lump_and_spray(wet3D, volume; f=Returns(true),di=2, dj=2, dk=1)

returns the LUMP and SPRAY matrices, and the coarsened volume vector.

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

The optional function argument `f` can be used to specify a region
where the lumping should occur.
`f` should be a function of the indices, e.g.,

    f(i,j,k) = lat[i,j] > -40 # only lump north of 40°S

Outside of this region, no lumping.
Default is `f=Returns(true)`, i.e. lump everywhere.
"""
function lump_and_spray(wet3D, volume; f=Returns(true), di=2, dj=2, dk=1)

    # extend the grid to avoid lumping cells outside of bounds
    nxyz = size(wet3D)
    LUMPidx = zeros(Int, nxyz .+ (di - 1, dj - 1, dk - 1))

    # loop once over all grid cells and assign lumped indices
    c = 1 # index / counter
    neighbours = CartesianIndices((0:di-1, 0:dj-1, 0:dk-1))
    C = CartesianIndices(nxyz)
    for 𝑖 in eachindex(C)
        C𝑖 = C[𝑖]
        i, j, k = C𝑖.I
        # assign index if not indexed yet
        if LUMPidx[C𝑖] == 0
            if f(i,j,k) # to all lumped cells if in region
                LUMPidx[C𝑖 .+ neighbours] .= c
            else # to only this cell if outside region
                LUMPidx[C𝑖] = c
            end
            c += 1
        end
    end
    LUMP = sparse(LUMPidx[:], 1:length(LUMPidx), 1)

    # Find wet points in coarsened grid
    wet = wet3D[:]
    wet_c = LUMP * wet .> 0

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

    return LUMP, SPRAY, volume_c
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





















