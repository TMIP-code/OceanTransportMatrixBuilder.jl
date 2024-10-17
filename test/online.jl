

@testitem "ACCESS-ESM1-5" begin

    # using Test
    # using OceanTransportMatrixBuilder
    using NetCDF
    using YAXArrays
    using Zarr
    using DataFrames
    using CSV
    using Unitful
    using Unitful: s, Myr
    using Downloads: download

    # stdlib
    using SparseArrays
    using LinearAlgebra

    @info "Loading catalog"
    catalog = CSV.read(download("https://storage.googleapis.com/cmip6/cmip6-zarr-consolidated-stores.csv"), DataFrame)

    @info "Finding ztores for uo, vo, umo, vmo, mlotst, areacello, and volcello"
    common_filter(row) = row.source_id == "ACCESS-ESM1-5" && row.member_id == "r1i1p1f1" && row.experiment_id == "historical"
    @show areacello_store = filter(row -> row.variable_id == "areacello" && row.table_id == "Ofx" && common_filter(row), catalog).zstore[1]
    @show volcello_store = filter(row -> row.variable_id == "volcello" && row.table_id == "Ofx" && common_filter(row), catalog).zstore[1]
    @show mlotst_store = filter(row -> row.variable_id == "mlotst" && row.table_id == "Omon" && common_filter(row), catalog).zstore[1]
    @show umo_store = filter(row -> row.variable_id == "umo" && row.table_id == "Omon" && common_filter(row), catalog).zstore[1]
    @show vmo_store = filter(row -> row.variable_id == "vmo" && row.table_id == "Omon" && common_filter(row), catalog).zstore[1]
    @show uo_store = filter(row -> row.variable_id == "uo" && row.table_id == "Omon" && common_filter(row), catalog).zstore[1]
    @show vo_store = filter(row -> row.variable_id == "vo" && row.table_id == "Omon" && common_filter(row), catalog).zstore[1]

    @info "Lazily load data using Zarr"
    volcello_zarr = zopen(volcello_store, consolidated=true)
    areacello_zarr = zopen(areacello_store, consolidated=true)
    mlotst_zarr = zopen(mlotst_store, consolidated=true)
    umo_zarr = zopen(umo_store, consolidated=true)
    vmo_zarr = zopen(vmo_store, consolidated=true)
    uo_zarr = zopen(uo_store, consolidated=true)
    vo_zarr = zopen(vo_store, consolidated=true)

    @info "Select first time step to limit download size"
    volcello_ds = open_dataset(volcello_zarr)
    areacello_ds = open_dataset(areacello_zarr)
    mlotst_ds = open_dataset(mlotst_zarr)[time = 1]
    umo_ds = open_dataset(umo_zarr)[time = 1]
    vmo_ds = open_dataset(vmo_zarr)[time = 1]
    uo_ds = open_dataset(uo_zarr)[time = 1]
    vo_ds = open_dataset(vo_zarr)[time = 1]

    @info "Load variables in memory"
    umo = readcubedata(umo_ds.umo)
    vmo = readcubedata(vmo_ds.vmo)
    uo = readcubedata(uo_ds.uo)
    vo = readcubedata(vo_ds.vo)
    uo_lon = readcubedata(uo_ds.longitude)
    uo_lat = readcubedata(uo_ds.latitude)
    vo_lon = readcubedata(vo_ds.longitude)
    vo_lat = readcubedata(vo_ds.latitude)
    mlotst = readcubedata(mlotst_ds.mlotst)
    areacello = readcubedata(areacello_ds.areacello)
    volcello = readcubedata(volcello_ds.volcello)
    lon = readcubedata(volcello_ds.longitude)
    lat = readcubedata(volcello_ds.latitude)
    lev = volcello_ds.lev
    lon_vertices = readcubedata(volcello_ds.vertices_longitude) # no xmip so must use default dataset propery names
    lat_vertices = readcubedata(volcello_ds.vertices_latitude) # no xmip so must use default dataset propery names

    # Some parameter values
    ρ = 1035.0    # kg/m^3
    κH = 500.0    # m^2/s
    κVML = 0.1    # m^2/s
    κVdeep = 1e-5 # m^2/s

    # Make makemodelgrid
    modelgrid = makemodelgrid(; areacello, volcello, lon, lat, lev, lon_vertices, lat_vertices)

    # Make fuxes from all directions
    ϕ = facefluxesfrommasstransport(; umo, vmo)

    # Make fuxes from all directions
    ϕ_bis = facefluxesfromvelocities(; uo, uo_lon, uo_lat, vo, vo_lon, vo_lat, modelgrid, ρ)

    for dir in (:east, :west, :north, :south, :top, :bottom)
        @test_broken isapprox(getpropery(ϕ, dir), getpropery(ϕ_bis, dir), rtol = 0.1)
    end

    # Make indices
    indices = makeindices(modelgrid.v3D)

    # Make transport matrix
    (; T, Tadv, TκH, TκVML, TκVdeep) = transportmatrix(; ϕ, mlotst, modelgrid, indices, ρ, κH, κVML, κVdeep)

    Tsyms = (:T, :Tadv, :TκH, :TκVML, :TκVdeep)
	for Ttest in (T, Tadv, TκH, TκVML, TκVdeep)
        @test Ttest isa SparseMatrixCSC{Float64, Int}
    end

	# Check divergence and mass conservation

    # unpack model grid
    (; v3D,) = modelgrid
    # unpack indices
    (; wet3D, N) = indices
	e1 = ones(N)
	v = v3D[wet3D]
	Tsyms = (:T, :Tadv, :TκH, :TκVML, :TκVdeep)
	@info "Timescales (divergence and mass conservation)"
	for (Ttest, Tsym) in zip((T, Tadv, TκH, TκVML, TκVdeep), Tsyms)
		@info "  $Tsym:"

        τdiv = ustrip(Myr, norm(e1) / norm(Ttest * e1) * s)
        Tsym ∉ (:T, :Tadv) && @test τdiv > 1e6
		@info "    div: $(round(τdiv, sigdigits=2)) Myr"

		τvol = ustrip(Myr, norm(v) / norm(Ttest' * v) * s)
        @test τvol > 1e6
		@info "    vol: $(round(τvol, sigdigits=2)) Myr"
	end

    # tests if diagonal elements are > 0 and off-diagonal are < 0.
    diagT = sparse(Diagonal(T))

    @test all(diagT.nzval .> 0)
    @test all((T - diagT).nzval .< 0)

end