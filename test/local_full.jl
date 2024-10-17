


@testmodule LocalBuiltMatrix begin

    using Test
    using OceanTransportMatrixBuilder
    using NetCDF
    using YAXArrays
    using GibbsSeaWater
    using DimensionalData
    using NaNStatistics

    # stdlib
    using SparseArrays
    using LinearAlgebra

    # My local directory for input files
    model = "ACCESS-ESM1-5"
    member = "r1i1p1f1"
    inputdir = "/Users/benoitpasquier/Data/TMIP/data/$model/historical/$member/Jan1990-Dec1999"

    # and for output files
    @show version = "v$(pkgversion(OceanTransportMatrixBuilder))"
    outputdir = joinpath("plots", version)
    mkpath(outputdir)

    # Load datasets
    umo_ds = open_dataset(joinpath(inputdir, "umo.nc"))
    vmo_ds = open_dataset(joinpath(inputdir, "vmo.nc"))
    uo_ds = open_dataset(joinpath(inputdir, "uo.nc"))
    vo_ds = open_dataset(joinpath(inputdir, "vo.nc"))
    mlotst_ds = open_dataset(joinpath(inputdir, "mlotst.nc"))
    volcello_ds = open_dataset(joinpath(inputdir, "volcello.nc"))
    areacello_ds = open_dataset(joinpath(inputdir, "areacello.nc"))


    # Load variables in memory
    umo = readcubedata(umo_ds.umo)
    vmo = readcubedata(vmo_ds.vmo)
    umo_lon = readcubedata(umo_ds.lon)
    umo_lat = readcubedata(umo_ds.lat)
    vmo_lon = readcubedata(vmo_ds.lon)
    vmo_lat = readcubedata(vmo_ds.lat)
    mlotst = readcubedata(mlotst_ds.mlotst)
    areacello = readcubedata(areacello_ds.areacello)
    volcello = readcubedata(volcello_ds.volcello)
    lon = readcubedata(volcello_ds.lon)
    lat = readcubedata(volcello_ds.lat)
    lev = volcello_ds.lev
    lon_vertices = readcubedata(volcello_ds.lon_verticies) # xmip issue: https://github.com/jbusecke/xMIP/issues/369
    lat_vertices = readcubedata(volcello_ds.lat_verticies) # xmip issue: https://github.com/jbusecke/xMIP/issues/369

    # Make makemodelgrid
    modelgrid = makemodelgrid(; areacello, volcello, lon, lat, lev, lon_vertices, lat_vertices)
    (; lon_vertices, lat_vertices, v3D, zt, thkcello) = modelgrid

    uo = readcubedata(uo_ds.uo)
    vo = readcubedata(vo_ds.vo)
    uo_lon = readcubedata(uo_ds.lon)
    uo_lat = readcubedata(uo_ds.lat)
    vo_lon = readcubedata(vo_ds.lon)
    vo_lat = readcubedata(vo_ds.lat)

    # Load thetato and so to compute density
    thetao_ds = open_dataset(joinpath(inputdir, "thetao.nc"))
    so_ds = open_dataset(joinpath(inputdir, "so.nc"))
    # Load variables in memory
    thetao = readcubedata(thetao_ds.thetao)
    @test 0 < nanmean(thetao) < 20
    so = readcubedata(so_ds.so)
    @show 30 < nanmean(so) < 40
    # Convert thetao and so to density
    ct = gsw_ct_from_pt.(so, thetao)
    @show nanmean(ct)
    ρθ = gsw_rho.(so, ct, 0)
    @show nanmean(ρθ)
    # from MATLAB GSW toolbox:
    # gsw_rho.(so, ct, p)
    # so = Absolute Salinity (g/kg)
    # ct = Conservative Temperature (ITS-90) (°C)
    # p = sea pressure (dbar) (here using 0 pressure to get potential density
    # TODO: CHECK IF THIS IS CORRECT

    # Some parameter values
    ρ = 1035.0    # kg/m^3
    # Alternatively, rebuild density from thetao, so, and depth as approximate pressure
    ZBOT3D = cumsum(thkcello, dims = 3)
    Z3D = ZBOT3D - 0.5 * thkcello
    ρ = gsw_rho.(so, ct, Z3D)

    # Diffusivites
    κH = 500.0    # m^2/s
    κVML = 0.1    # m^2/s
    κVdeep = 1e-5 # m^2/s

    # Make fuxes from all directions
    ϕ = facefluxesfrommasstransport(; umo, vmo)

    # Make fuxes from all directions from velocities
    ϕ_bis = facefluxesfromvelocities(; uo, uo_lon, uo_lat, vo, vo_lon, vo_lat, modelgrid, ρ)

    # Make indices
    indices = makeindices(modelgrid.v3D)

    @test all(.!isnan.(ρ[indices.wet3D])) == true

    # Make transport matrix
    (; T, Tadv, TκH, TκVML, TκVdeep) = transportmatrix(; ϕ, mlotst, modelgrid, indices, ρ, κH, κVML, κVdeep)

    Tsyms = (:T, :Tadv, :TκH, :TκVML, :TκVdeep)
	for Ttest in (T, Tadv, TκH, TκVML, TκVdeep)
        @test Ttest isa SparseMatrixCSC{Float64, Int}
    end


end

@testitem "velocities and mass transports" setup=[LocalBuiltMatrix] tags=[:skipci] begin

    using NaNStatistics
    using GLMakie

    (; modelgrid, indices, ρθ, v3D, lat, lon, zt, uo, vo, uo_lon, uo_lat, vo_lon, vo_lat,
    lon_vertices, lat_vertices, indices,
    umo, vmo, umo_lon, umo_lat, vmo_lon, vmo_lat, model, member, outputdir) = LocalBuiltMatrix



    # plot for sanity check
    begin # plot zonal average
        ρθ2D = dropdims(nansum(ρθ .* v3D, dims = 1) ./ nansum(v3D, dims = 1), dims = 1)
        fig = Figure()
        ax = Axis(fig[1,1], xlabel = "latitude (°)", ylabel = "depth (m)")
        # levels = 25:0.1:30
        colormap = :viridis
        # co = contourf!(ax, dropdims(maximum(lat |> Array, dims=1), dims=1), zt |> Array, Γ2D; levels, colormap)
        co = contourf!(ax, dropdims(maximum(lat |> Array, dims=1), dims=1), zt |> Array, ρθ2D; colormap)
        cb = Colorbar(fig[1, 2], co; label = "Potential density (?)", tellheight = false)
        cb.height = Relative(2/3)
        ylims!(ax, (6000, 0))
        Label(fig[0,1], text = "$model $member Potential density", tellwidth = false)
        fig
    end
    outputfile = joinpath(outputdir, "potential_density_$model.png")
    @info "Saving ideal age as image file:\n $(joinpath("test", outputfile))"
    save(outputfile, fig)

    κGM = 600 # m^2/s
    maxslope = 0.01
    uGM, vGM = OceanTransportMatrixBuilder.bolus_GM_velocity(ρθ, modelgrid; κGM, maxslope)

    # Plot location of cell center for volcello, umo, vmo, uo, vo
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "lon", ylabel = "lat")
    lines!(ax, lon_vertices[:,1,1] |> Vector, lat_vertices[:,1,1] |> Vector)
    text!(ax, lon[1], lat[1]; text="area (i,j)", align = (:center, :bottom))
    text!(ax, umo_lon[1], umo_lat[1]; text="umo (i,j)", align = (:left, :center))
    text!(ax, uo_lon[1], uo_lat[1]; text="uo (i,j)", align = (:left, :center))
    text!(ax, umo_lon[1], umo_lat[1]; text="vmo (i,j)", align = (:center, :bottom))
    text!(ax, vo_lon[1], vo_lat[1]; text="vo (i,j)", align = (:center, :bottom))
    fig

    uo2, uo2_lon, uo2_lat, vo2, vo2_lon, vo2_lat = OceanTransportMatrixBuilder.interpolateontodefaultCgrid(uo, uo_lon, uo_lat, vo, vo_lon, vo_lat, modelgrid)
    fig = Figure(size=(1500, 800))
    ax = Axis(fig[1,1], xlabel = "lon", ylabel = "lat", xgridvisible = false, ygridvisible = false)
    i = 85 .+ (-1:5); j = 2 .+ (-1:4); k = 1
    [poly!(ax, lon_vertices[:,i,j] |> vec, lat_vertices[:,i,j] |> vec, color = (:blue, 0.1)) for i in i for j in j if indices.wet3D[i,j,1]];
    [lines!(ax, lon_vertices[:,i,j] |> vec, lat_vertices[:,i,j] |> vec, color = (:black, 0.1)) for i in i for j in j];
    # lines!(ax, lon_vertices[1,i,j] |> Array, lat_vertices[1,i,j] |> Array, color = :black)
    arrows!(ax, uo_lon[i,j,k] |> vec, uo_lat[i,j,k] |> vec, uo[i,j,k] |> vec, 0uo[i,j,k] |> vec, arrowsize = 10, lengthscale = 10, arrowcolor = :darkblue, linecolor = :darkblue)
    arrows!(ax, vo_lon[i,j,k] |> vec, vo_lat[i,j,k] |> vec, 0vo[i,j,k] |> vec, vo[i,j,k] |> vec, arrowsize = 10, lengthscale = 10, arrowcolor = :darkred, linecolor = :darkred)
    arrows!(ax, uo2_lon[i,j,k] |> vec, uo2_lat[i,j,k] |> vec, uo2[i,j,k] |> vec, 0uo2[i,j,k] |> vec, arrowsize = 10, lengthscale = 10, arrowcolor = :blue, linecolor = :blue)
    arrows!(ax, vo2_lon[i,j,k] |> vec, vo2_lat[i,j,k] |> vec, 0vo2[i,j,k] |> vec, vo2[i,j,k] |> vec, arrowsize = 10, lengthscale = 10, arrowcolor = :red, linecolor = :red)
    arrows!(ax, uo2_lon[i,j,k] |> vec, uo2_lat[i,j,k] .+ 0.02 |> vec, umo[i,j,k] |> vec, 0umo[i,j,k] |> vec, arrowsize = 10, lengthscale = 3e-8, arrowcolor = (:blue, 0.3), linecolor = (:blue, 0.3))
    arrows!(ax, vo2_lon[i,j,k] .+ 0.02 |> vec, vo2_lat[i,j,k] |> vec, 0vmo[i,j,k] |> vec, vmo[i,j,k] |> vec, arrowsize = 10, lengthscale = 3e-8, arrowcolor = (:red, 0.3), linecolor = (:red, 0.3))
    fig
    ax = Axis(fig[1,2], xlabel = "lon", ylabel = "lat", xgridvisible = false, ygridvisible = false)
    i = 160 .+ (-3:3); j = 200 .+ (-2:2); k = 1
    [poly!(ax, lon_vertices[:,i,j] |> vec, lat_vertices[:,i,j] |> vec, color = (:blue, 0.1)) for i in i for j in j if indices.wet3D[i,j,1]];
    [lines!(ax, lon_vertices[:,i,j] |> vec, lat_vertices[:,i,j] |> vec, color = (:black, 0.1)) for i in i for j in j];
    # lines!(ax, lon_vertices[1,i,j] |> Array, lat_vertices[1,i,j] |> Array, color = :black)
    arrows!(ax, uo_lon[i,j,k] |> vec, uo_lat[i,j,k] |> vec, uo[i,j,k] |> vec, 0uo[i,j,k] |> vec, arrowsize = 10, lengthscale = 10, arrowcolor = :darkblue, linecolor = :darkblue)
    arrows!(ax, vo_lon[i,j,k] |> vec, vo_lat[i,j,k] |> vec, 0vo[i,j,k] |> vec, vo[i,j,k] |> vec, arrowsize = 10, lengthscale = 10, arrowcolor = :darkred, linecolor = :darkred)
    arrows!(ax, uo2_lon[i,j,k] |> vec, uo2_lat[i,j,k] |> vec, uo2[i,j,k] |> vec, 0uo2[i,j,k] |> vec, arrowsize = 10, lengthscale = 10, arrowcolor = :blue, linecolor = :blue)
    arrows!(ax, vo2_lon[i,j,k] |> vec, vo2_lat[i,j,k] |> vec, 0vo2[i,j,k] |> vec, vo2[i,j,k] |> vec, arrowsize = 10, lengthscale = 10, arrowcolor = :red, linecolor = :red)
    arrows!(ax, uo2_lon[i,j,k] |> vec, uo2_lat[i,j,k] .+ 0.02 |> vec, umo[i,j,k] |> vec, 0umo[i,j,k] |> vec, arrowsize = 10, lengthscale = 1e-8, arrowcolor = (:blue, 0.3), linecolor = (:blue, 0.3))
    arrows!(ax, vo2_lon[i,j,k] .+ 0.02 |> vec, vo2_lat[i,j,k] |> vec, 0vmo[i,j,k] |> vec, vmo[i,j,k] |> vec, arrowsize = 10, lengthscale = 1e-8, arrowcolor = (:red, 0.3), linecolor = (:red, 0.3))
    fig

    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "i", ylabel = "j")
    hm = heatmap!(ax, indices.wet3D[:,:,1])
    translate!(hm, 0, 0, -100) # move the plot behind the grid
    fig
end

@testitem "Divergence and mass conservation" setup=[LocalBuiltMatrix] tags=[:skipci] begin

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

@testitem "Test flux divergence" setup=[LocalBuiltMatrix] tags=[:skipci] begin

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
    using NaNStatistics
    using GLMakie

    (; modelgrid, indices, T, model, member, outputdir) = LocalBuiltMatrix

    # unpack model grid
    (; v3D, lat, zt) = modelgrid
    # unpack indices
    (; wet3D, N) = indices

    v = v3D[wet3D]

    @info "coarsening grid"
    LUMP, SPRAY, wet3D_c, v_c = OceanTransportMatrixBuilder.lump_and_spray(wet3D, v; di=2, dj=2, dk=1)

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
    begin # plot zonal average
        Γ2D = dropdims(nansum(Γ3D .* v3D, dims = 1) ./ nansum(v3D, dims = 1), dims = 1)
        fig = Figure()
        ax = Axis(fig[1,1], xlabel = "latitude (°)", ylabel = "depth (m)")
        levels = 0:100:2000
        colormap = :viridis
        co = contourf!(ax, dropdims(maximum(lat |> Array, dims=1), dims=1), zt |> Array, Γ2D; levels, colormap)
        cb = Colorbar(fig[1, 2], co; label = "Ideal mean age (yr)", tellheight = false)
        cb.height = Relative(2/3)
        ylims!(ax, (6000, 0))
        Label(fig[0,1], text = "$model $member ideal mean age", tellwidth = false)
        fig
    end
    outputfile = joinpath(outputdir, "ideal_age_$model.png")
    @info "Saving ideal age as image file:\n $(joinpath("test", outputfile))"
    save(outputfile, fig)
    @test 0 < (v' * Γyr) / sum(v) < 2000

end

@testitem "mass transport vs velocity checks" setup=[LocalBuiltMatrix] tags=[:skipci] begin

    using GLMakie
    using Makie.StructArrays

    outputdir = LocalBuiltMatrix.outputdir
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
    outputfile = joinpath(outputdir, "cell_face_fluxes_check_local_ACCESS-ESM1-5.png")
    @info "Saving cell face fluxes check as image file:\n $(joinpath("test", outputfile))"
    save(outputfile, fig)

    for dir in dirs
        @test_broken isapprox(tovec(ϕ, dir), tovec(ϕ_bis, dir), rtol = 0.1)
    end

    # Difference between the fluxes from velocities and mass transport
    (; uo, vo, umo, vmo, uo_lon, uo_lat, vo_lon, vo_lat, modelgrid, ρ) = LocalBuiltMatrix
    umo_bis, vmo_bis = velocity2fluxes(uo, uo_lon, uo_lat, vo, vo_lon, vo_lat, modelgrid, ρ)
    colorrange = 1e9 .* (-1, 1)
    colormap = cgrad(:RdBu, rev=true)
    Δcolorrange = (-5, 5)
    Δcolormap = cgrad(:PRGn, rev=true)
    for k in axes(umo, 3)
    # for k in [1]
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
        outputdir2 = joinpath(outputdir, "fluxes_from_velocity")
        mkpath(outputdir2)
        local outputfile = joinpath(outputdir2, "k=$(k)_check_local_ACCESS-ESM1-5.png")
        @info "Saving fluxes comparison k=$k check as image file:\n $(joinpath("test", outputfile))"
        save(outputfile, fig)
    end

end

