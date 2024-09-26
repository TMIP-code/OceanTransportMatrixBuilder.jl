


@testmodule BuiltMatrix begin

    using Test
    using OceanTransportMatrixBuilder
    using NetCDF
    using YAXArrays
    using SparseArrays
    using LinearAlgebra

    # My local directory for input files
    inputdir = "/Users/benoitpasquier/Data/TMIP/data/ACCESS-ESM1-5/historical/r1i1p1f1/Jan1990-Dec1999"

    # Load umo, vmo, mlotst, volcello, and areacello
    umo_ds = open_dataset(joinpath(inputdir, "umo.nc"))
    vmo_ds = open_dataset(joinpath(inputdir, "vmo.nc"))
    mlotst_ds = open_dataset(joinpath(inputdir, "mlotst.nc"))
    volcello_ds = open_dataset(joinpath(inputdir, "volcello.nc"))
    areacello_ds = open_dataset(joinpath(inputdir, "areacello.nc"))

    mlotst = mlotst_ds["mlotst"] |> Array{Float64}

    # Make ualldirs
    u = makeualldirections(; umo_ds, vmo_ds)

    # Make makemodelgrid
    modelgrid = makemodelgrid(; areacello_ds, volcello_ds, mlotst_ds)

    # Make indices
    indices = makeindices(modelgrid.v3D)

    # Make transport matrix
    (; T, Tadv, TκH, TκVML, TκVdeep) = transportmatrix(; u, mlotst, modelgrid, indices,
        ρ = 1025.0,
        κH = 500.0, # m^2/s
        κVML = 0.1, # m^2/s
        κVdeep = 1e-5, # m^2/s
    )

    Tsyms = (:T, :Tadv, :TκH, :TκVML, :TκVdeep)
	for Ttest in (T, Tadv, TκH, TκVML, TκVdeep)
        @test Ttest isa SparseMatrixCSC{Float64, Int}
    end

end


@testmodule BuiltACCESSModelGrids begin

    using Test
    using OceanTransportMatrixBuilder
    using NetCDF
    using YAXArrays
    using SparseArrays
    using LinearAlgebra

    modelgrids = Dict()
    indicess = Dict()
    mlotsts = Dict()

    models = ("ACCESS-ESM1-5", "ACCESS-CM2", "ACCESS1-3")
    members = ("r1i1p1f1", "r1i1p1f1", "r1i1p1")

    for (model, member) in zip(models, members)

        # My local directory for input files
        inputdir = "/Users/benoitpasquier/Data/TMIP/data/$model/historical/$member/Jan1990-Dec1999"

        # Load umo, vmo, mlotst, volcello, and areacello
        mlotst_ds = open_dataset(joinpath(inputdir, "mlotst.nc"))
        volcello_ds = open_dataset(joinpath(inputdir, "volcello.nc"))
        areacello_ds = open_dataset(joinpath(inputdir, "areacello.nc"))

        # Make makemodelgrid
        modelgrids[model] = makemodelgrid(; areacello_ds, volcello_ds, mlotst_ds)

        # Make indices
        indicess[model] = makeindices(modelgrids[model].v3D)

        mlotsts[model] = mlotst_ds["mlotst"] |> Array{Float64}

    end

end

@testitem "Timescales (divergence and mass conservation)" setup=[BuiltMatrix] tags=[:skipci] begin

    using Unitful
    using Unitful: s, Myr
    using LinearAlgebra

    # unpack transport matrices
    (; T, Tadv, TκH, TκVML, TκVdeep, Tsyms) = BuiltMatrix
    # unpack model grid
    (; v3D,) = BuiltMatrix.modelgrid
    # unpack indices
    (; wet3D, N) = BuiltMatrix.indices

	e1 = ones(N)
	v = v3D[wet3D]
	for (Ttest, Tsym) in zip((T, Tadv, TκH, TκVML, TκVdeep), Tsyms)
		@info "  $Tsym:"

        τdiv = ustrip(Myr, norm(e1) / norm(Ttest * e1) * s)
        Tsym ∉ (:T, :Tadv) && @test τdiv > 1e6
		@info "    div: $(round(τdiv, sigdigits=2)) Myr"

		τvol = ustrip(Myr, norm(v) / norm(Ttest' * v) * s)
        @test τvol > 1e6
		@info "    vol: $(round(τvol, sigdigits=2)) Myr"
	end
end

@testitem "Test if flux divergence (not convergence)" setup=[BuiltMatrix] tags=[:skipci] begin

    using SparseArrays
    using LinearAlgebra

    # unpack transport matrices
    (; T) = BuiltMatrix

    # tests if diagonal elements are > 0 and off-diagonal are < 0.
    diagT = sparse(Diagonal(T))
    @test all(diagT.nzval .> 0)
    @test all((T - diagT).nzval .< 0)

end

@testitem "Ideal age (coarsened)" setup=[BuiltMatrix] tags=[:skipci] begin

    using SparseArrays
    using LinearAlgebra
    using Unitful
    using Unitful: s, yr


    # unpack model grid
    (; v3D,) = BuiltMatrix.modelgrid
    # unpack indices
    (; wet3D, N) = BuiltMatrix.indices

    v = v3D[wet3D]

    @info "coarsening grid"
    LUMP, SPRAY, wet3D_c, v_c = OceanTransportMatrixBuilder.lump_and_spray(wet3D, v; di=2, dj=2, dk=1)

    # unpack transport matrices
    (; T) = BuiltMatrix

    # surface mask
    issrf3D = copy(wet3D)
    issrf3D[:,:,2:end] .= false
    issrf = issrf3D[wet3D]
    # Ideal mean age Γ is governed by
    # 	∂Γ/∂t + T Γ = 1 - M Γ
    # where M is matrix mask of surface with short timescale (1s)
    sΓ = ones(size(v))
    T_c = LUMP * T * SPRAY
    issrf_c = LUMP * issrf .> 0
    M_c = sparse(Diagonal(issrf_c))
    sΓ_c = LUMP * sΓ
    @info "Solving ideal mean age"
    Γ_c = (T_c + M_c) \ sΓ_c
    Γ = SPRAY * Γ_c
    Γyr = ustrip.(yr, Γ .* s)
    Γ3D = OceanTransportMatrixBuilder.as3D(Γyr, wet3D)

    @test 0 < (v' * Γyr) / sum(v) < 2000

end


@testitem "grid checks" setup=[BuiltACCESSModelGrids] tags=[:skipci] begin

    using GLMakie
    using Unitful
    using Unitful: m, km

    function customdecorations!(ax; i = 1, j = 1, imax = 1)
        hidexdecorations!(ax,
            label = i < imax, ticklabels = i < imax, ticks = i < imax, grid = false)
        hideydecorations!(ax,
            label = j > 1, ticklabels = j > 1, ticks = j > 1, grid = false)
    end
    m3 = rich("m", superscript("3"))
    km3 = rich("km", superscript("3"))
    m2 = rich("m", superscript("2"))
    km2 = rich("km", superscript("2"))

    model0 = "ACCESS-ESM1-5"
    modelgrid0 = BuiltACCESSModelGrids.modelgrids[model0]
    indices0 = BuiltACCESSModelGrids.indicess[model0]

    for model in keys(BuiltACCESSModelGrids.modelgrids)

        @info "Checking grid for $model"
        modelgrid = BuiltACCESSModelGrids.modelgrids[model]
        indices = BuiltACCESSModelGrids.indicess[model]

        (; wet3D, L, Lwet, N, Lwet3D, C) = indices

        (; lon_vertices, lat_vertices, edge_length_2D, distance_to_edge_2D, lon, lat, zt, area2D, v3D, DZT3d) = modelgrid

        @test isequal(wet3D, indices0.wet3D)
        @test isequal(L, indices0.L)
        @test isequal(Lwet, indices0.Lwet)
        @test isequal(N, indices0.N)
        @test isequal(Lwet3D, indices0.Lwet3D)
        @test isequal(C, indices0.C)

        @test isequal(zt, modelgrid0.zt)
        @test isequal(lon, modelgrid0.lon)
        @test isequal(lat, modelgrid0.lat)
        @test isequal(area2D, modelgrid0.area2D)

        @test isapprox(v3D, modelgrid0.v3D, nans = true, rtol = 0.15)
        @test isapprox(DZT3d, modelgrid0.DZT3d, nans = true, rtol = 0.15)

        @test isapprox(mod.(lon_vertices, 360), mod.(modelgrid0.lon_vertices, 360), nans = true)
        @test isapprox(lat_vertices, modelgrid0.lat_vertices, nans = true)
        for dir in (:east, :west, :north, :south)
            @test isapprox(edge_length_2D[dir], modelgrid0.edge_length_2D[dir], nans = true)
            @test isapprox(distance_to_edge_2D[dir], modelgrid0.distance_to_edge_2D[dir], nans = true)
        end

        colorrange = (0, 120)
        fig = Figure(size=(1200, 500));
        for (j, dir) in enumerate((:east, :west, :north, :south))
            i = 1
            ax = Axis(fig[i,j], xlabel = "i", ylabel = "j")
            heatmap!(ax, ustrip.(km, edge_length_2D[dir] * m); colorrange)
            customdecorations!(ax; i, j, imax = 2)

            i = 2
            ax = Axis(fig[i,j], xlabel = "i", ylabel = "j")
            heatmap!(ax, ustrip.(km, distance_to_edge_2D[dir] * m); colorrange)
            customdecorations!(ax; i, j, imax = 2)

            Label(fig[0,j], text=string(dir), tellwidth=false)
        end
        Label(fig[1,0], text="edge_length_2D", tellheight=false, rotation=π/2)
        Label(fig[2,0], text="distance_to_edge_2D", tellheight=false, rotation=π/2)
        Label(fig[-1,1:4], text="Distances check $model model", tellwidth=false, fontsize=20)
        cb = Colorbar(fig[1:2,5]; limits=colorrange, label="km")
        cb.height = Relative(0.8)
        fig
        outputfile = "plots/distances_check_$model.png"
        @info "Saving distances check as image file:\n  $(outputfile)"
        save(outputfile, fig)



        fig = Figure()
        ax = Axis(fig[1,1], xlabel = "longitude (°)", ylabel = "latitude (°)")
        aligns = [(:left, :bottom), (:right, :bottom), (:right, :top), (:left, :top)]
        is = 220 .+ [1, 2, 1]
        js = 270 .+ [1, 1, 2]
        itxts = ["i", "i+1", "i"]
        jtxts = ["j", "j", "j+1"]
        colors = Makie.wong_colors(3)
        for (i, j, itxt, jtxt, color, align) in zip(is, js, itxts, jtxts, colors, aligns)
            scatterlines!(ax, lon_vertices[:,i,j], lat_vertices[:,i,j]; marker=:circle, color)
            text!(ax, collect(zip(lon_vertices[:,i,j], lat_vertices[:,i,j])); align=align, text=string.(1:4), color)
            text!(ax, lon[i, j], lat[i, j]; text="($itxt, $jtxt)", color)
        end
        Label(fig[0,1], text="Vertices check $model model", tellwidth = false, fontsize = 20)
        fig
        outputfile = "plots/vertices_check_local_$model.png"
        @info "Saving vertices check as image file:\n  $(outputfile)"
        save(outputfile, fig)

        # Check that the correct vertex order was applied to all cells
        @test lon_vertices[2,1:end-1,:] == lon_vertices[1,2:end,:]
        @test lon_vertices[3,1:end-1,:] == lon_vertices[4,2:end,:]
        @test lat_vertices[2,1:end-1,:] == lat_vertices[1,2:end,:]
        @test lat_vertices[3,1:end-1,:] == lat_vertices[4,2:end,:]
        @test lon_vertices[3,:,1:end-1] == lon_vertices[2,:,2:end]
        @test lon_vertices[4,:,1:end-1] == lon_vertices[1,:,2:end]
        @test lat_vertices[3,:,1:end-1] == lat_vertices[2,:,2:end]
        @test lat_vertices[4,:,1:end-1] == lat_vertices[1,:,2:end]



        # plot lon lat grid
        fig = Figure(size=(2000, 1000))

        loncut = (lon[1,1] + lon[end,1]) / 2
        lon = mod.(lon .- loncut, 360) .+ loncut

        ax = Axis(fig[1,1], xlabel = "longitude (°)", ylabel = "latitude (°)")
        lonnorthsouth = reduce(vcat, [[col; NaN] for col in eachcol(lon)])
        latnorthsouth = reduce(vcat, [[col; NaN] for col in eachcol(lat)])
        lines!(ax, lonnorthsouth, latnorthsouth, linewidth = 0.5)
        loneastwest = reduce(vcat, [[row; NaN] for row in eachrow(lon)])
        lateastwest = reduce(vcat, [[row; NaN] for row in eachrow(lat)])
        lines!(ax, loneastwest, lateastwest, linewidth = 0.5)

        Label(fig[0,1], text="Longitude/latitude grid check $model model", tellwidth = false, fontsize = 20)
        fig
        outputfile = "plots/lonlat_check_local_$model.png"
        @info "Saving lon/lat check as image file:\n  $(outputfile)"
        save(outputfile, fig)



        # plot area
        fig = Figure(size=(1000, 500))
        ax = Axis(fig[1,1], xlabel = "i", ylabel = "j")
        cf = contourf!(ax, ustrip.(km^2, area2D * m^2))
        cb = Colorbar(fig[1,2], cf; label = km2)
        Label(fig[0,1], text="Area grid check $model model", tellwidth = false, fontsize = 20)
        fig
        outputfile = "plots/area_check_local_$model.png"
        @info "Saving area check as image file:\n  $(outputfile)"
        save(outputfile, fig)




        # plot v3D every 1000m
        fig = Figure(size=(1000, 750))
        zs = [0 1000; 2000 3000; 4000 5000]
        colorrange = (0, 3000) # km^3
        hm = nothing
        for idx in eachindex(IndexCartesian(), zs)
            i, j = Tuple(idx)
            ax = Axis(fig[i,j], xlabel = "i", ylabel = "j")
            k = findfirst(zt .> zs[i,j])
            hm = heatmap!(ax, ustrip.(km^3, v3D[:,:,k] * m^3); colorrange)
            customdecorations!(ax; i, j, imax = size(zs, 1))
            text!(ax, 0.5, 1, text = "$(zs[i,j])m", align = (:center, :top), offset = (0, -2), space = :relative)
        end
        cb = Colorbar(fig[1:3,end+1], hm; label = km3)
        cb.height = Relative(0.666)
        Label(fig[0,1:2], text="Volume grid check $model model", tellwidth = false, fontsize = 20)
        fig
        outputfile = "plots/volume_check_local_$model.png"
        @info "Saving volume check as image file:\n  $(outputfile)"
        save(outputfile, fig)



        # plot DZT3d every 1000m
        fig = Figure(size=(1000, 750))
        zs = [0 1000; 2000 3000; 4000 5000]
        colorrange = (0, 500) # m
        hm = nothing
        for idx in eachindex(IndexCartesian(), zs)
            i, j = Tuple(idx)
            ax = Axis(fig[i,j], xlabel = "i", ylabel = "j")
            k = findfirst(zt .> zs[i,j])
            hm = heatmap!(ax, ustrip.(m, DZT3d[:,:,k] * m); colorrange)
            customdecorations!(ax; i, j, imax = size(zs, 1))
            text!(ax, 0.5, 1, text = "$(zs[i,j])m", align = (:center, :top), offset = (0, -2), space = :relative)
        end
        cb = Colorbar(fig[1:3,end+1], hm; label = "m")
        cb.height = Relative(0.666)
        Label(fig[0,1:2], text="Cell thickness grid check $model model", tellwidth = false, fontsize = 20)
        fig
        outputfile = "plots/cell_thickness_check_local_$model.png"
        @info "Saving cell thickness check as image file:\n  $(outputfile)"
        save(outputfile, fig)




        # Compare volume
        fig = Figure()
        ax = Axis(fig[1,1], xlabel = rich("$model0 cell volumes (", km3, ")"), ylabel = rich("$model cell volumes (", km3, ")"))
        X = ustrip.(km^3, modelgrid0.v3D[indices0.wet3D] * m^3)
        Y = ustrip.(km^3, v3D[wet3D] * m^3)
        scatter!(ax, X, Y, markersize = 0.5)
        ax = Axis(fig[1,2], xlabel = rich("$model0 cell volumes (", km3, ")"), ylabel = "relative error (%)")
        scatter!(ax, X, 100 * abs.(Y .- X) ./ X, markersize = 0.5)
        outputfile = "plots/volume_comparison_$(model0)_vs_$(model).png"
        @info "Saving volume comparison as image file:\n  $(outputfile)"
        save(outputfile, fig)




        # plot mixed layer depths
        fig = Figure(size=(1000, 750))
        colorrange = (50, 5000) # m
        colorscale = Makie.pseudolog10
        hm = nothing
        ax = Axis(fig[1,1], xlabel = "i", ylabel = "j")
        hm = heatmap!(ax, BuiltACCESSModelGrids.mlotsts[model]; colorrange, colorscale)
        cb = Colorbar(fig[1,end+1], hm; label = "MLD (m)", ticks=[50, 100, 200, 500, 1000, 2000, 5000])
        cb.height = Relative(0.666)
        Label(fig[0,1], text="$model mixed-layer depth", tellwidth = false, fontsize = 20)
        fig
        outputfile = "plots/MLD_check_local_$model.png"
        @info "Saving MLD check as image file:\n  $(outputfile)"
        save(outputfile, fig)

    end



end