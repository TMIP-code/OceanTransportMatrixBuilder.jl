


@testitem "grid checks" tags = [:skipci] begin

    using Test
    using OceanTransportMatrixBuilder
    using NetCDF
    using YAXArrays
    using SparseArrays
    using LinearAlgebra
    using GLMakie
    using Unitful
    using Unitful: m, km

    # TODO fix this so that I don't have to track all members
    models = ("ACCESS-ESM1-5", "ACCESS-CM2", "ACCESS1-3", "ACCESS1-0")
    members = ("r1i1p1f1", "r1i1p1f1", "r1i1p1", "r1i1p1")

    # models_untested = ("CMCC-CM2-HR4",)
    # members_untested = ("r1i1p1f1",)

    @show version = "v$(pkgversion(OceanTransportMatrixBuilder))"
    outputdir = joinpath("plots", version)
    mkpath(outputdir)


    function customdecorations!(ax; i=1, j=1, imax=1)
        hidexdecorations!(ax,
            label=i < imax, ticklabels=i < imax, ticks=i < imax, grid=false)
        hideydecorations!(ax,
            label=j > 1, ticklabels=j > 1, ticks=j > 1, grid=false)
    end
    m3 = rich("m", superscript("3"))
    km3 = rich("km", superscript("3"))
    m2 = rich("m", superscript("2"))
    km2 = rich("km", superscript("2"))

    for (model, member) in zip(models, members)

        # My local directory for input files
        inputdir = "/Users/benoitpasquier/Data/TMIP/data/$model/historical/$member/Jan1990-Dec1999"

        # Load datasets
        volcello_ds = open_dataset(joinpath(inputdir, "volcello.nc"))
        areacello_ds = open_dataset(joinpath(inputdir, "areacello.nc"))
        mlotst_ds = open_dataset(joinpath(inputdir, "mlotst.nc"))

        # Load variables in memory
        mlotst = readcubedata(mlotst_ds.mlotst)
        areacello = readcubedata(areacello_ds.areacello)
        volcello = readcubedata(volcello_ds.volcello)
        lon = readcubedata(volcello_ds.lon)
        lat = readcubedata(volcello_ds.lat)
        lev = volcello_ds.lev
        # Identify the vertices keys (vary across CMIPs / models)
        volcello_keys = propertynames(volcello_ds)
        lon_vertices_key = volcello_keys[findfirst(x -> occursin("lon", x) & occursin("vert", x), string.(volcello_keys))]
        lat_vertices_key = volcello_keys[findfirst(x -> occursin("lat", x) & occursin("vert", x), string.(volcello_keys))]
        lon_vertices = readcubedata(getproperty(volcello_ds, lon_vertices_key))
        lat_vertices = readcubedata(getproperty(volcello_ds, lat_vertices_key))


        # Make makegridmetrics
        gridmetrics = makegridmetrics(; areacello, volcello, lon, lat, lev, lon_vertices, lat_vertices)
        (; lon_vertices, lat_vertices, edge_length_2D, distance_to_edge_2D, lon, lat, zt, area2D, v3D, thkcello) = gridmetrics

        # Make indices
        indices = makeindices(v3D)
        (; wet3D, L, Lwet, N, Lwet3D, C) = indices

        @info "Checking grid for $model"


        colorrange = (0, 120)
        fig = Figure(size=(1200, 500))
        for (j, dir) in enumerate((:east, :west, :north, :south))
            i = 1
            ax = Axis(fig[i, j], xlabel="i", ylabel="j")
            heatmap!(ax, ustrip.(km, edge_length_2D[dir] * m); colorrange)
            customdecorations!(ax; i, j, imax=2)

            i = 2
            ax = Axis(fig[i, j], xlabel="i", ylabel="j")
            heatmap!(ax, ustrip.(km, distance_to_edge_2D[dir] * m); colorrange)
            customdecorations!(ax; i, j, imax=2)

            Label(fig[0, j], text=string(dir), tellwidth=false)
        end
        Label(fig[1, 0], text="edge_length_2D", tellheight=false, rotation=π / 2)
        Label(fig[2, 0], text="distance_to_edge_2D", tellheight=false, rotation=π / 2)
        Label(fig[-1, 1:4], text="Distances check $model model", tellwidth=false, fontsize=20)
        cb = Colorbar(fig[1:2, 5]; limits=colorrange, label="km")
        cb.height = Relative(0.8)
        fig
        outputfile = joinpath(outputdir, "distances_check_$model.png")
        @info "Saving distances check as image file:\n  $(joinpath("test", outputfile))"
        save(outputfile, fig)



        fig = Figure()
        ax = Axis(fig[1, 1], xlabel="longitude (°)", ylabel="latitude (°)")
        aligns = [(:left, :bottom), (:right, :bottom), (:right, :top), (:left, :top)]
        is = 220 .+ [1, 2, 1]
        js = 270 .+ [1, 1, 2]
        itxts = ["i", "i+1", "i"]
        jtxts = ["j", "j", "j+1"]
        colors = Makie.wong_colors(3)
        for (i, j, itxt, jtxt, color, align) in zip(is, js, itxts, jtxts, colors, aligns)
            scatterlines!(ax, lon_vertices[:, i, j], lat_vertices[:, i, j]; marker=:circle, color)
            text!(ax, collect(zip(lon_vertices[:, i, j], lat_vertices[:, i, j])); align=align, text=string.(1:4), color)
            text!(ax, lon[i, j], lat[i, j]; text="($itxt, $jtxt)", color)
        end
        Label(fig[0, 1], text="Vertices check $model model", tellwidth=false, fontsize=20)
        fig
        outputfile = joinpath(outputdir, "vertices_check_local_$model.png")
        @info "Saving vertices check as image file:\n  $(joinpath("test", outputfile))"
        save(outputfile, fig)

        # Check that the correct vertex order was applied to all cells
        @test lon_vertices[2, 1:end-1, :] == lon_vertices[1, 2:end, :]
        @test lon_vertices[3, 1:end-1, :] == lon_vertices[4, 2:end, :]
        @test lat_vertices[2, 1:end-1, :] == lat_vertices[1, 2:end, :]
        @test lat_vertices[3, 1:end-1, :] == lat_vertices[4, 2:end, :]
        @test lon_vertices[3, :, 1:end-1] == lon_vertices[2, :, 2:end]
        @test lon_vertices[4, :, 1:end-1] == lon_vertices[1, :, 2:end]
        @test lat_vertices[3, :, 1:end-1] == lat_vertices[2, :, 2:end]
        @test lat_vertices[4, :, 1:end-1] == lat_vertices[1, :, 2:end]



        # plot lon lat grid edges
        fig = Figure(size=(2000, 1000))
        loncut = lon_vertices[1, 1, 1]
        # TODO: this is currently not working for plotting cells that wrap around the dateline
        # It's OK for longitudes to be in the same window for checking connections,
        # but it would be good to implement a simple function to fix that for plotting.
        lon_edges = [
            lon_vertices[1, :, :] lon_vertices[4, :, end]
            lon_vertices[2, end, :]' lon_vertices[3, end, end]
        ]
        lon_edges[1:end-1, :] .= mod.(lon_edges[1:end-1, :] .- loncut, 360) .+ loncut
        lon_edges[end, :] .= mod1.(lon_edges[end, :] .- loncut, 360) .+ loncut
        lat_edges = [
            lat_vertices[1, :, :] lat_vertices[4, :, end]
            lat_vertices[2, end, :]' lat_vertices[3, end, end]
        ]
        ax = Axis(fig[1, 1], xlabel="longitude (°)", ylabel="latitude (°)")
        lonnorthsouth = reduce(vcat, [[col; NaN] for col in eachcol(lon_edges)])
        latnorthsouth = reduce(vcat, [[col; NaN] for col in eachcol(lat_edges)])
        loneastwest = reduce(vcat, [[row; NaN] for row in eachrow(lon_edges)])
        lateastwest = reduce(vcat, [[row; NaN] for row in eachrow(lat_edges)])
        lines!(ax, [lonnorthsouth; NaN; loneastwest], [latnorthsouth; NaN; lateastwest], linewidth=0.5)
        # Add boarders
        lines!(ax, lon_edges[1, :], lat_edges[1, :], linewidth=2, color=:black)
        lines!(ax, lon_edges[end, :], lat_edges[end, :], linewidth=2, color=:red)
        lines!(ax, lon_edges[:, 1], lat_edges[:, 1], linewidth=2, color=:black)
        lines!(ax, lon_edges[:, end], lat_edges[:, end], linewidth=2, color=:green)
        Label(fig[0, 1], text="Grid cell edges check $model model", tellwidth=false, fontsize=20)
        ylims!(ax, (50, 91))
        fig
        outputfile = joinpath(outputdir, "grid_cell_edges_check_local_$model.png")
        @info "Saving lon/lat check as image file:\n  $(outputfile)"
        save(outputfile, fig)



        # plot lon lat grid centers
        fig = Figure(size=(2000, 1000))
        loncut = (lon[1, 1] + lon[end, 1]) / 2
        lon2 = mod.(lon .- loncut, 360) .+ loncut
        ax = Axis(fig[1, 1], xlabel="longitude (°)", ylabel="latitude (°)")
        lonnorthsouth = reduce(vcat, [[col; NaN] for col in eachcol(lon2)[1:20:end]])
        latnorthsouth = reduce(vcat, [[col; NaN] for col in eachcol(lat)[1:20:end]])
        lines!(ax, lonnorthsouth, latnorthsouth, linewidth=0.5)
        loneastwest = reduce(vcat, [[row; NaN] for row in eachrow(lon2)[1:20:end]])
        lateastwest = reduce(vcat, [[row; NaN] for row in eachrow(lat)[1:20:end]])
        lines!(ax, loneastwest, lateastwest, linewidth=0.5)
        Label(fig[0, 1], text="Grid cell centers check $model model", tellwidth=false, fontsize=20)
        fig
        outputfile = joinpath(outputdir, "grid_cell_centers_check_local_$model.png")
        @info "Saving lon/lat check as image file:\n  $(joinpath("test", outputfile))"
        save(outputfile, fig)



        # plot area
        fig = Figure(size=(1000, 500))
        ax = Axis(fig[1, 1], xlabel="i", ylabel="j")
        cf = contourf!(ax, ustrip.(km^2, area2D * m^2))
        cb = Colorbar(fig[1, 2], cf; label=km2)
        Label(fig[0, 1], text="Area grid check $model model", tellwidth=false, fontsize=20)
        fig
        outputfile = joinpath(outputdir, "area_check_local_$model.png")
        @info "Saving area check as image file:\n  $(joinpath("test", outputfile))"
        save(outputfile, fig)




        # plot v3D every 1000m
        fig = Figure(size=(1000, 750))
        zs = [0 1000; 2000 3000; 4000 5000]
        colorrange = (0, 3000) # km^3
        hm = nothing
        for idx in eachindex(IndexCartesian(), zs)
            i, j = Tuple(idx)
            ax = Axis(fig[i, j], xlabel="i", ylabel="j")
            k = findfirst(zt .> zs[i, j])
            hm = heatmap!(ax, ustrip.(km^3, v3D[:, :, k] * m^3); colorrange)
            customdecorations!(ax; i, j, imax=size(zs, 1))
            text!(ax, 0.5, 1, text="$(zs[i,j])m", align=(:center, :top), offset=(0, -2), space=:relative)
        end
        cb = Colorbar(fig[1:3, end+1], hm; label=km3)
        cb.height = Relative(0.666)
        Label(fig[0, 1:2], text="Volume grid check $model model", tellwidth=false, fontsize=20)
        fig
        outputfile = joinpath(outputdir, "volume_check_local_$model.png")
        @info "Saving volume check as image file:\n  $(joinpath("test", outputfile))"
        save(outputfile, fig)



        # plot thkcello every 1000m
        fig = Figure(size=(1000, 750))
        zs = [0 1000; 2000 3000; 4000 5000]
        colorrange = (0, 500) # m
        hm = nothing
        for idx in eachindex(IndexCartesian(), zs)
            i, j = Tuple(idx)
            ax = Axis(fig[i, j], xlabel="i", ylabel="j")
            k = findfirst(zt .> zs[i, j])
            hm = heatmap!(ax, ustrip.(m, thkcello[:, :, k] * m); colorrange)
            customdecorations!(ax; i, j, imax=size(zs, 1))
            text!(ax, 0.5, 1, text="$(zs[i,j])m", align=(:center, :top), offset=(0, -2), space=:relative)
        end
        cb = Colorbar(fig[1:3, end+1], hm; label="m")
        cb.height = Relative(0.666)
        Label(fig[0, 1:2], text="Cell thickness grid check $model model", tellwidth=false, fontsize=20)
        fig
        outputfile = joinpath(outputdir, "cell_thickness_check_local_$model.png")
        @info "Saving cell thickness check as image file:\n  $(joinpath("test", outputfile))"
        save(outputfile, fig)



        # plot mixed layer depths
        fig = Figure(size=(1000, 750))
        colorrange = (50, 5000) # m
        colorscale = Makie.pseudolog10
        hm = nothing
        ax = Axis(fig[1, 1], xlabel="i", ylabel="j")
        hm = heatmap!(ax, mlotst; colorrange, colorscale)
        cb = Colorbar(fig[1, end+1], hm; label="MLD (m)", ticks=[50, 100, 200, 500, 1000, 2000, 5000])
        cb.height = Relative(0.666)
        Label(fig[0, 1], text="$model mixed-layer depth", tellwidth=false, fontsize=20)
        fig
        outputfile = joinpath(outputdir, "MLD_check_local_$model.png")
        @info "Saving MLD check as image file:\n  $(joinpath("test", outputfile))"
        save(outputfile, fig)

    end



end