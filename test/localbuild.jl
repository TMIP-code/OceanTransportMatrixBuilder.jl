


@testmodule LocalBuild begin

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

@testitem "Timescales (divergence and mass conservation)" setup=[LocalBuild] tags=[:skipci] begin

    using Unitful
    using Unitful: s, Myr
    using LinearAlgebra

    # unpack transport matrices
    (; T, Tadv, TκH, TκVML, TκVdeep, Tsyms) = LocalBuild
    # unpack model grid
    (; v3D,) = LocalBuild.modelgrid
    # unpack indices
    (; wet3D, N) = LocalBuild.indices

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

@testitem "Test if flux divergence (not convergence)" setup=[LocalBuild] tags=[:skipci] begin

    using SparseArrays
    using LinearAlgebra

    # unpack transport matrices
    (; T) = LocalBuild

    # tests if diagonal elements are > 0 and off-diagonal are < 0.
    diagT = sparse(Diagonal(T))
    @test all(diagT.nzval .> 0)
    @test all((T - diagT).nzval .< 0)

end

@testitem "Ideal age (coarsened)" setup=[LocalBuild] tags=[:skipci] begin

    using SparseArrays
    using LinearAlgebra
    using Unitful
    using Unitful: s, yr


    # unpack model grid
    (; v3D,) = LocalBuild.modelgrid
    # unpack indices
    (; wet3D, N) = LocalBuild.indices

    v = v3D[wet3D]

    @info "coarsening grid"
    LUMP, SPRAY, wet3D_c, v_c = OceanTransportMatrixBuilder.lump_and_spray(wet3D, v; di=2, dj=2, dk=1)

    # unpack transport matrices
    (; T) = LocalBuild

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


@testitem "grid checks" setup=[LocalBuild] tags=[:skipci] begin

    using GLMakie

    (; lon_vertices, lat_vertices, edge_length_2D, distance_to_edge_2D, lon, lat) = LocalBuild.modelgrid

    colorrange = (0, 1e5)
    fig = Figure(size=(1200, 500));
    for (i, dir) in enumerate((:east, :west, :north, :south))
        ax = Axis(fig[1,i])
        heatmap!(ax, edge_length_2D[dir]; colorrange)
        ax = Axis(fig[2,i])
        heatmap!(ax, distance_to_edge_2D[dir]; colorrange)
        Label(fig[0,i], text=string(dir), tellwidth=false)
    end
    Label(fig[1,0], text="edge_length_2D", tellheight=false, rotation=π/2)
    Label(fig[2,0], text="distance_to_edge_2D", tellheight=false, rotation=π/2)
    fig
    outputfile = "plots/distances_check.png"
    @info "Saving ideal mean age as image file:\n  $(outputfile)"
    save(outputfile, fig)



    fig = Figure()
    ax = Axis(fig[1,1])
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
    fig
    outputfile = "plots/vertices_check.png"
    @info "Saving ideal mean age as image file:\n  $(outputfile)"
    save(outputfile, fig)

end