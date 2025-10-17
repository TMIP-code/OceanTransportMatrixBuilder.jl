"""
    LUMP, SPRAY, v_c = lump_and_spray(wet3D, vol, T, mask=trues(size(wet3D)); di=2, dj=2, dk=1)

returns the LUMP and SPRAY matrices, and the coarsened volume vector.

The `vol` vector of volumes is used to make the LUMP operator volume-conserving.
The `T` matrix is used to determine connectivity of the grid cells, so as to avoid lumping disconnected cells.
The `mask` boolean array (same size as `wet3D`) can be used to restrict the lumping to a certain region.

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
function lump_and_spray(wet3D, vol, T, mask=trues(size(wet3D)); di=2, dj=2, dk=1)

    # extend the grid to avoid lumping cells outside of bounds
    nxyz = size(wet3D)
    C = CartesianIndices(nxyz)
    LUMPidx = zeros(Int, nxyz .+ (di - 1, dj - 1, dk - 1))

    wet = wet3D[:]
    i, j = findnz(T)
    L = LinearIndices(nxyz .+ (di - 1, dj - 1, dk - 1))
    wet3Dext = falses(nxyz .+ (di - 1, dj - 1, dk - 1))
    wet3Dext[C[wet3D]] .= true
    Lwet = L[wet3Dext] # careful here, as L is larger than wet3D

    connectivitymatrix = sparse(Lwet[i], Lwet[j], true, length(LUMPidx), length(LUMPidx))

    # loop once over all grid cells and assign lumped indices
    c = 2 # index / counter (we start at 2 because 1 is reserved for dry cells)
    neighbours = CartesianIndices((0:di-1, 0:dj-1, 0:dk-1))
    for 𝑖 in eachindex(C)
        C𝑖 = C[𝑖]
        # skip if index already assigned and in mask
        #(if assigned but not in mask, must reassign)
        LUMPidx[C𝑖] > 0 && mask[C𝑖] && continue
        if mask[C𝑖] # if in region, assign all neighbours also in region
            # list of neighbours
            L𝑖 = L[C𝑖 .+ neighbours][:]
            # First, assign dry index to dry cells
            localwet = wet3Dext[L𝑖]
            dryidx = L𝑖[.!localwet]
            LUMPidx[dryidx] .= 1
            # Then build a simple connectivity graph of the remaining cells to be lumped
            wetidx = L𝑖[localwet]
            localconnectivitymatrix = view(connectivitymatrix, wetidx, wetidx)
            g = SimpleGraph(localconnectivitymatrix)
            # and split it in connected components
            for comp in connected_components(g)
                LUMPidx[wetidx[comp]] .= c
                c += 1
            end
        else # if outside region, assign new index
            LUMPidx[C𝑖] = c
            c += 1
        end
    end

    # Remove ghost cells outside of original grid
    LUMP = sparse(LUMPidx[C][:], 1:length(C), 1)

    # Find wet points in coarsened grid
    wet_c = LUMP * wet .> 0

    # Extract only indices of wet grd points
    LUMP = LUMP[wet_c, wet]

    # Make the LUMP operator volume-conserving
    # by volume integrating on the right and dividing by the coarse
    # volume on the left
    vol_c = LUMP * vol
    LUMP = sparse(Diagonal(1 ./ vol_c)) * LUMP * sparse(Diagonal(vol))

    # The SPRAY operator just copies the values back
    # so it is sinply 1's with the transposed sparsity structure
    SPRAY = copy(LUMP')
    SPRAY.nzval .= 1

    # Print some info
    nwet = size(LUMP, 2)
    nwet_c = size(LUMP, 1)
    @info """LUMP and SPRAY:
    Matrix size reduction: $(round(100(1 - nwet_c/nwet)))% ($nwet -> $nwet_c)
    """

    return LUMP, SPRAY, vol_c
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





















