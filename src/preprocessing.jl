

function makeualldirections(; umo, vmo)

	iseastborder = isnan.(umo[[2:end;1],:,:]) .& .!isnan.(umo)
	isnorthborder = isnan.([vmo[:,2:end,:] vmo[end:-1:1,end:end,:]]) .& .!isnan.(vmo)
	umo[iseastborder] .= 0
	vmo[isnorthborder] .= 0

	@info "Making ueast and uwest"
	ueast = replace(umo, NaN=>0.0) .|> Float64
	uwest = ueast[[end;1:end-1],:,:] .|> Float64

	@info "Making vnorth and vsouth"
	# Check that south pole vnorth is zero (so that it causes no issues with circular shift)
	vnorth = replace(vmo, NaN=>0.0) .|> Float64
	vsouth = vnorth[:,[end;1:end-1],:] .|> Float64

	@info "Making wtop and wbottom"
	# Then build wtop and wbottom from bottom up
	# Mass conservation implies that
	# 	uwest + vsouth + wbottom - ueast - vnorth - wtop = 0
	# except at the top.
	# We could build w from the bottom up, by looking at each
	# water column's sea floor, but it's simpler to go from the top down,
	# and then remove wbottom
	wbottom = similar(ueast)
	wtop = similar(ueast)
	for k in reverse(eachindex(axes(ueast, 3)))
		if k == lastindex(axes(ueast, 3))
			@views @. wbottom[:,:,k] = 0 # seafloor wbottom is zero
		else
			@views @. wbottom[:,:,k] = wtop[:,:,k+1] # otherwise it's wtop from below
		end
		@views @. wtop[:,:,k] = wbottom[:,:,k] + uwest[:,:,k] + vsouth[:,:,k] - ueast[:,:,k] - vnorth[:,:,k]
	end

    ualldirs = (
        east = ueast,
        west = uwest,
        north = vnorth,
        south = vsouth,
        top = wtop,
        bottom = wbottom
    )

	return ualldirs
end


function makemodelgrid(; areacello_ds, volcello_ds, mlotst_ds)

    areacello = areacello_ds["areacello"] |> Array{Float64}
    volcello = volcello_ds["volcello"] |> Array{Union{Missing, Float64}}

	# area and volume
	v3D = replace(volcello, missing => NaN, 0 => NaN)
    area2D = replace(areacello, missing => NaN, 0 => NaN)

	# depth and cell height
	DZT3d = v3D ./ area2D
	zt = volcello_ds["lev"] |> Array
    lat = volcello_ds["lat"] |> Array
    lon = volcello_ds["lon"] |> Array
    lon_vertices = mlotst_ds["lon_verticies"] |> Array # <- typo from xmip (https://github.com/jbusecke/xMIP/issues/369)
    lat_vertices = mlotst_ds["lat_verticies"] |> Array # <- typo from xmip (https://github.com/jbusecke/xMIP/issues/369)

    C = CartesianIndices(size(lon))

    dirs = (:south, :east, :north, :west)
    edge_length_2D = Dict(d=>[verticalfacewidth(lon_vertices, lat_vertices, 𝑖.I[1], 𝑖.I[2], d) for 𝑖 in C] for d in dirs)
    distance_to_edge_2D = Dict(d=>[centroid2edgedistance(lon, lat, lon_vertices, lat_vertices, 𝑖.I[1], 𝑖.I[2], d) for 𝑖 in C] for d in dirs)

	return (; area2D, v3D, DZT3d, lon_vertices, lat_vertices, lon, lat, zt, edge_length_2D, distance_to_edge_2D)
end

function makeindices(v3D)

    # LinearIndices and CartesianIndices required for building upwind operator
    nxyz = size(v3D)
    L = LinearIndices(nxyz)
    Lwet = L[.!isnan.(v3D)]
    N = length(Lwet)
    wet3D = falses(nxyz...)
    wet3D[Lwet] .= true
    Lwet3D = Array{Union{Int, Missing}, 3}(missing, nxyz...)
    Lwet3D[Lwet] .= 1:length(Lwet)
    C = CartesianIndices(nxyz)

    return (; wet3D, L, Lwet, N, Lwet3D, C)
end



# TODO generalized topology detection
# That is, how do I figure out how the boundaries are connected for other models than ACCESS?




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
