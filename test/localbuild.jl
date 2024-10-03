


@testmodule LocalBuiltMatrix begin

    using Test
    using OceanTransportMatrixBuilder
    using NetCDF
    using YAXArrays

    # stdlib
    using SparseArrays
    using LinearAlgebra

    # My local directory for input files
    inputdir = "/Users/benoitpasquier/Data/TMIP/data/ACCESS-ESM1-5/historical/r1i1p1f1/Jan1990-Dec1999"

    # Load datasets
    umo_ds = open_dataset(joinpath(inputdir, "umo.nc"))
    vmo_ds = open_dataset(joinpath(inputdir, "vmo.nc"))
    uo_ds = open_dataset(joinpath(inputdir, "uo.nc"))
    vo_ds = open_dataset(joinpath(inputdir, "vo.nc"))
    mlotst_ds = open_dataset(joinpath(inputdir, "mlotst.nc"))
    volcello_ds = open_dataset(joinpath(inputdir, "volcello.nc"))
    areacello_ds = open_dataset(joinpath(inputdir, "areacello.nc"))

    # Load variables
    umo = umo_ds.umo
    vmo = vmo_ds.vmo
    uo = uo_ds.uo
    vo = vo_ds.vo
    mlotst = mlotst_ds.mlotst
    areacello = areacello_ds.areacello
    volcello = volcello_ds.volcello
    lon = volcello_ds.lon
    lat = volcello_ds.lat
    lev = volcello_ds.lev
    lon_vertices = volcello_ds.lon_verticies # xmip issue: https://github.com/jbusecke/xMIP/issues/369
    lat_vertices = volcello_ds.lat_verticies # xmip issue: https://github.com/jbusecke/xMIP/issues/369

    # Plot location of cell center for volcello, umo, vmo, uo, vo
    using GLMakie
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "lon", ylabel = "lat")
    lines!(ax, lon_vertices[:,1,1] |> Vector, lat_vertices[:,1,1] |> Vector)
    text!(ax, lon[1], lat[1]; text="area (i,j)", align = (:center, :bottom))
    text!(ax, umo_ds.lon[1], umo_ds.lat[1]; text="umo (i,j)", align = (:left, :center))
    text!(ax, uo_ds.lon[1], uo_ds.lat[1]; text="uo (i,j)", align = (:left, :center))
    text!(ax, vmo_ds.lon[1], vmo_ds.lat[1]; text="vmo (i,j)", align = (:center, :bottom))
    text!(ax, vo_ds.lon[1], vo_ds.lat[1]; text="vo (i,j)", align = (:center, :bottom))
    fig

    # Some parameter values
    ρ = 1035.0    # kg/m^3
    κH = 500.0    # m^2/s
    κVML = 0.1    # m^2/s
    κVdeep = 1e-5 # m^2/s

    # Make makemodelgrid
    modelgrid = makemodelgrid(; areacello, volcello, lon, lat, lev, lon_vertices, lat_vertices)

    # Make ualldirs
    ϕ = facefluxesfrommasstransport(; umo, vmo)
    ϕ_bis = facefluxesfromvelocities(; uo, vo, modelgrid, ρ)

    # Make indices
    indices = makeindices(modelgrid.v3D)

    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "i", ylabel = "j")
    hm = heatmap!(ax, indices.wet3D[:,:,1])
    fig

    # Make transport matrix
    (; T, Tadv, TκH, TκVML, TκVdeep) = transportmatrix(; ϕ, mlotst, modelgrid, indices, ρ, κH, κVML, κVdeep)

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

        # Load datasets
        volcello_ds = open_dataset(joinpath(inputdir, "volcello.nc"))
        areacello_ds = open_dataset(joinpath(inputdir, "areacello.nc"))
        mlotst_ds = open_dataset(joinpath(inputdir, "mlotst.nc"))
        mlotsts[model] = mlotst_ds.mlotst

        areacello = areacello_ds.areacello
        volcello = volcello_ds.volcello
        lon = volcello_ds.lon
        lat = volcello_ds.lat
        lev = volcello_ds.lev
        if model ∈ ("ACCESS-ESM1-5", "ACCESS-CM2")
            lon_vertices = volcello_ds.lon_verticies # xmip issue: https://github.com/jbusecke/xMIP/issues/369
            lat_vertices = volcello_ds.lat_verticies # xmip issue: https://github.com/jbusecke/xMIP/issues/369
        elseif model ∈ ("ACCESS1-3", "ACCESS1-0")
            lon_vertices = volcello_ds.lon_vertices # no xmip so default key
            lat_vertices = volcello_ds.lat_vertices # no xmip so default key
        else
            error("Need to hardcode the lon_vertices and lat_vertices keys for $model in the tests")
        end

        # Make makemodelgrid
        modelgrids[model] = makemodelgrid(; areacello, volcello, lon, lat, lev, lon_vertices, lat_vertices)

        # Make indices
        indicess[model] = makeindices(modelgrids[model].v3D)

    end

end

@testitem "Timescales (divergence and mass conservation)" setup=[LocalBuiltMatrix] tags=[:skipci] begin

    using Unitful
    using Unitful: s, Myr
    using LinearAlgebra

    # unpack transport matrices
    (; T, Tadv, TκH, TκVML, TκVdeep, Tsyms) = LocalBuiltMatrix
    # unpack model grid
    (; v3D,) = LocalBuiltMatrix.modelgrid
    # unpack indices
    (; wet3D, N) = LocalBuiltMatrix.indices

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

@testitem "Test if flux divergence (not convergence)" setup=[LocalBuiltMatrix] tags=[:skipci] begin

    using SparseArrays
    using LinearAlgebra

    # unpack transport matrices
    (; T) = LocalBuiltMatrix

    # tests if diagonal elements are > 0 and off-diagonal are < 0.
    diagT = sparse(Diagonal(T))
    @test all(diagT.nzval .> 0)
    @test all((T - diagT).nzval .< 0)

end

@testitem "Ideal age (coarsened)" setup=[LocalBuiltMatrix] tags=[:skipci] begin

    using SparseArrays
    using LinearAlgebra
    using Unitful
    using Unitful: s, yr


    # unpack model grid
    (; v3D,) = LocalBuiltMatrix.modelgrid
    # unpack indices
    (; wet3D, N) = LocalBuiltMatrix.indices

    v = v3D[wet3D]

    @info "coarsening grid"
    LUMP, SPRAY, wet3D_c, v_c = OceanTransportMatrixBuilder.lump_and_spray(wet3D, v; di=2, dj=2, dk=1)

    # unpack transport matrices
    (; T) = LocalBuiltMatrix

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

@testitem "mass transport vs velocity checks" setup=[LocalBuiltMatrix] tags=[:skipci] begin

    using GLMakie
    using Makie.StructArrays

    (; ϕ, ϕ_bis) = LocalBuiltMatrix
    (; wet3D) = LocalBuiltMatrix.indices

    kgs⁻¹ = rich("kg s", superscript("−1"))
    fig = Figure(size = (1000, 600))
    uvw = [:u, :v, :w]
    dirs = [:east :north :top; :west :south :bottom]
    colormap = to_colormap(:BuPu_9)
    colormap[1] = RGBAf(1, 1, 1, 1)
    for idx in eachindex(IndexCartesian(), dirs)
        i, j = Tuple(idx)
        dir = dirs[idx]
        xlabel = rich("mass transport ϕ (", kgs⁻¹, ")")
        ylabel = rich("velocity ϕ (", kgs⁻¹, ")")
        ax = Axis(fig[i,j]; xlabel, ylabel, aspect = 1)
        imax = size(dirs, 1)
        hidexdecorations!(ax, label = i < imax, ticklabels = i < imax, ticks = i < imax, grid = false)
        hideydecorations!(ax, label = j > 1, ticklabels = false, ticks = false, grid = false)
        X = getproperty(ϕ, dir)[wet3D]
        Y = getproperty(ϕ_bis, dir)[wet3D]
        XYmax = extrema([X; Y])
        scatter!(ax, X, Y, markersize = 0.5)
        # points = Point2f.(X, Y)
        # plt = datashader!(ax, points; local_operation = x -> log10(x + 1), colormap)
        # translate!(plt, 0, 0, -100) # move the plot behind the grid
        xlims!(ax, XYmax)
        ylims!(ax, XYmax)
        text!(ax, 0, 1, text = string(dir), align = (:left, :top), offset = (5, -5), space = :relative)
        # (i == imax) && linkaxes!(contents(fig[1:imax,j])...)
    end
    # Labels
    for j in axes(dirs, 2)
        text = "$(uvw[j])o vs $(uvw[j])mo"
        Label(fig[0,j]; text, tellwidth = false)
    end
    Label(fig[-1,:], text = "Cell area fluxes from mass transport vs velocities", tellwidth = false, fontsize = 20)
    fig
    outputfile = "plots/cell_face_fluxes_check_local_ACCESS-ESM1-5.png"
    @info "Saving cell face fluxes check as image file:\n $(joinpath("test", outputfile))"
    save(outputfile, fig)

    for dir in dirs
        @test_broken isapprox(tovec(ϕ, dir), tovec(ϕ_bis, dir), rtol = 0.1)
    end

    # Difference between the fluxes from velocities and mass transport
    (; uo_ds, vo_ds, umo_ds, vmo_ds, modelgrid, ρ) = LocalBuiltMatrix
    uo = uo_ds.uo |> Array{Float64}
    vo = vo_ds.vo |> Array{Float64}
    umo = umo_ds.umo |> Array{Float64}
    vmo = vmo_ds.vmo |> Array{Float64}
    umo_bis, vmo_bis = OceanTransportMatrixBuilder.velocity2fluxes(; uo, vo, modelgrid, ρ)
    colorrange = 1e9 .* (-1, 1)
    colormap = cgrad(:RdBu, rev=true)
    Δcolorrange = (-100, 100)
    Δcolormap = cgrad(:PRGn, rev=true)
    # for k in axes(umo, 3)
    for k in [1]
        local fig = Figure(size=(1200, 500))

        local ax = Axis(fig[1,1], xlabel = "i", ylabel = "j")
        heatmap!(ax, umo_bis[:,:,k]; colormap, colorrange)
        text!(ax, 0, 1; text = "umo_bis", align = (:left, :top), offset = (5, -5), space = :relative)
        hidedecorations!(ax)

        ax = Axis(fig[2,1], xlabel = "i", ylabel = "j")
        heatmap!(ax, vmo_bis[:,:,k]; colormap, colorrange)
        text!(ax, 0, 1; text = "vmo_bis", align = (:left, :top), offset = (5, -5), space = :relative)
        hidedecorations!(ax)

        ax = Axis(fig[1,2], xlabel = "i", ylabel = "j")
        heatmap!(ax, umo[:,:,k]; colormap, colorrange)
        text!(ax, 0, 1; text = "umo", align = (:left, :top), offset = (5, -5), space = :relative)
        hidedecorations!(ax)

        ax = Axis(fig[2,2], xlabel = "i", ylabel = "j")
        hm = heatmap!(ax, vmo[:,:,k]; colormap, colorrange)
        text!(ax, 0, 1; text = "vmo", align = (:left, :top), offset = (5, -5), space = :relative)
        hidedecorations!(ax)

        ax = Axis(fig[1,3], xlabel = "i", ylabel = "j")
        heatmap!(ax, 100((umo_bis - umo) ./ umo)[:,:,k]; colormap = Δcolormap, colorrange = Δcolorrange)
        text!(ax, 0, 1; text = "(umo_bis − umo) / umo", align = (:left, :top), offset = (5, -5), space = :relative)
        hidedecorations!(ax)

        ax = Axis(fig[2,3], xlabel = "i", ylabel = "j")
        Δhm = heatmap!(ax, 100((vmo_bis - vmo) ./ vmo)[:,:,k]; colormap = Δcolormap, colorrange = Δcolorrange)
        text!(ax, 0, 1; text = "(vmo_bis − vmo) / vmo", align = (:left, :top), offset = (5, -5), space = :relative)
        hidedecorations!(ax)

        cb = Colorbar(fig[3,1:2], hm; label = "kg s⁻¹", vertical = false, flipaxis = false)
        cb.width = Relative(0.666)

        Δcb = Colorbar(fig[3,3], Δhm; label = "%", vertical = false, flipaxis = false)
        Δcb.width = Relative(1)

        Label(fig[0,:], text = "ϕ_bis - ϕ (k=$k)", tellwidth = false, fontsize = 20)
        outputdir = "plots/fluxes_from_velocity"
        mkpath(outputdir)
        local outputfile = joinpath(outputdir, "k=$(k)_check_local_ACCESS-ESM1-5.png")
        @info "Saving fluxes comparison k=$k check as image file:\n $(joinpath("test", outputfile))"
        save(outputfile, fig)
    end

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

            Label(fig[0,j], text = string(dir), tellwidth=false)
        end
        Label(fig[1,0], text = "edge_length_2D", tellheight=false, rotation=π/2)
        Label(fig[2,0], text = "distance_to_edge_2D", tellheight=false, rotation=π/2)
        Label(fig[-1,1:4], text = "Distances check $model model", tellwidth=false, fontsize=20)
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
            text!(ax, collect(zip(lon_vertices[:,i,j], lat_vertices[:,i,j])); align = align, text = string.(1:4), color)
            text!(ax, lon[i, j], lat[i, j]; text = "($itxt, $jtxt)", color)
        end
        Label(fig[0,1], text = "Vertices check $model model", tellwidth = false, fontsize = 20)
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

        Label(fig[0,1], text = "Longitude/latitude grid check $model model", tellwidth = false, fontsize = 20)
        fig
        outputfile = "plots/lonlat_check_local_$model.png"
        @info "Saving lon/lat check as image file:\n  $(outputfile)"
        save(outputfile, fig)



        # plot area
        fig = Figure(size=(1000, 500))
        ax = Axis(fig[1,1], xlabel = "i", ylabel = "j")
        cf = contourf!(ax, ustrip.(km^2, area2D * m^2))
        cb = Colorbar(fig[1,2], cf; label = km2)
        Label(fig[0,1], text = "Area grid check $model model", tellwidth = false, fontsize = 20)
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
        Label(fig[0,1:2], text = "Volume grid check $model model", tellwidth = false, fontsize = 20)
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
        Label(fig[0,1:2], text = "Cell thickness grid check $model model", tellwidth = false, fontsize = 20)
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
        Label(fig[0,1], text = "$model mixed-layer depth", tellwidth = false, fontsize = 20)
        fig
        outputfile = "plots/MLD_check_local_$model.png"
        @info "Saving MLD check as image file:\n  $(outputfile)"
        save(outputfile, fig)

    end



end