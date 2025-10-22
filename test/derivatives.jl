@testitem "derivatives" setup = [LocalBuiltMatrix] tags = [:skipci] begin

    using NaNStatistics
    using CairoMakie

    (;
        gridmetrics, indices, ρ, ρθ, v3D,
        lat, lon, zt, uo, vo, uo_lon, uo_lat, vo_lon, vo_lat,
        lon_vertices, lat_vertices, indices,
        umo, vmo, umo_lon, umo_lat, vmo_lon, vmo_lat,
        model, member, outputdir,
    ) = LocalBuiltMatrix
    (; Z3D) = gridmetrics

    # fake z data that looks like z^2
    Zsquared = Z3D .^ 2

    ∂Zsquared_∂z = OceanTransportMatrixBuilder.globalverticaldyadderivative(Zsquared, gridmetrics, indices)

    # plot for sanity check
    begin # plot zonal average density
        Zsquared2D = dropdims(nansum(Zsquared .* v3D, dims = 1) ./ nansum(v3D, dims = 1), dims = 1)
        ∂Zsquared_∂z2D = dropdims(nansum(∂Zsquared_∂z .* v3D, dims = 1) ./ nansum(v3D, dims = 1), dims = 1)
        fig = Figure()
        ax = Axis(fig[1, 1], xlabel = "latitude (°)", ylabel = "depth (m)")
        # levels = 25:0.1:30
        colormap = :viridis
        levels = (0:5:35) .* 1.0e6
        # co = contourf!(ax, dropdims(maximum(lat |> Array, dims=1), dims=1), zt |> Array, Γ2D; levels, colormap)
        co = contourf!(ax, dropdims(maximum(lat |> Array, dims = 1), dims = 1), zt |> Array, Zsquared2D; levels, colormap)
        Zsquaredunit = rich("m", superscript("2"))
        cb = Colorbar(fig[1, 2], co; label = rich("Zsquared (", Zsquaredunit, ")"), tellheight = false)
        cb.height = Relative(4 / 5)
        ylims!(ax, (6000, 0))

        ax = Axis(fig[2, 1], xlabel = "latitude (°)", ylabel = "depth (m)")
        # levels = 25:0.1:30
        colormap = :magma
        levels = (-6:1:0) .* 2.0e3
        # co = contourf!(ax, dropdims(maximum(lat |> Array, dims=1), dims=1), zt |> Array, Γ2D; levels, colormap)
        co = contourf!(ax, dropdims(maximum(lat |> Array, dims = 1), dims = 1), zt |> Array, ∂Zsquared_∂z2D; levels, colormap)
        ∂Zsquared_∂zunit = rich("m")
        cb = Colorbar(fig[2, 2], co; label = rich("∂Zsquared/∂z (", ∂Zsquared_∂zunit, ")"), tellheight = false)
        cb.height = Relative(4 / 5)
        ylims!(ax, (6000, 0))
        fig
    end
    outputfile = joinpath(outputdir, "derivatives_z_dyads_$model.png")
    @info "Saving dyads z derivative zonal average as image file:\n $(joinpath("test", outputfile))"
    save(outputfile, fig)


    # test slopes on density
    ∂z∂yρθ = OceanTransportMatrixBuilder.globalverticalfacetriadderivative(ρθ, gridmetrics, indices, OceanTransportMatrixBuilder.Jcoord())

    # plot for sanity check
    begin # plot zonal average density
        # ρθ2D = dropdims(nansum(ρθ .* v3D, dims = 1) ./ nansum(v3D, dims = 1), dims = 1)
        # ∂z∂yρθ2D = dropdims(nansum(∂z∂yρθ .* v3D, dims = 1) ./ nansum(v3D, dims = 1), dims = 1)
        ilon = 130 # in the Pacific
        ρθ2D = ρθ[ilon, :, :]
        ∂z∂yρθ2D = ∂z∂yρθ[ilon, :, :]
        fig = Figure()
        ax = Axis(fig[1, 1], xlabel = "latitude (°)", ylabel = "depth (m)")
        levels = -100:10:100
        colormap = :balance
        # Filled contour for the slope
        cof = contourf!(ax, dropdims(maximum(lat |> Array, dims = 1), dims = 1), zt |> Array, clamp.(100∂z∂yρθ2D, -99.9, 99.9); levels, colormap)
        # Add density contours
        levels = 1000 .+ (0:0.1:100)
        co = contour!(ax, dropdims(maximum(lat |> Array, dims = 1), dims = 1), zt |> Array, ρθ2D; levels, color = :black, linewidth = 0.5)
        translate!(co, 0, 0, 100)
        ylims!(ax, (6000, 0))
        ∂z∂yρθunit = rich("%")
        cb = Colorbar(fig[1, 2], cof; label = rich("(∂ρθ/∂y) / (∂ρθ/∂z) (", ∂z∂yρθunit, ")"), tellheight = false)
        cb.height = Relative(4 / 5)
        fig
    end
    outputfile = joinpath(outputdir, "derivatives_zy_triads_$model.png")
    @info "Saving triads zy derivative meridional transect as image file:\n $(joinpath("test", outputfile))"
    save(outputfile, fig)

end
