

@testitem "building matrices" begin

    using NetCDF
    using YAXArrays
    using Zarr
    using DataFrames
    using CSV
    using Unitful
    using Unitful: cm, m, km, kg, s, d, yr, Myr, hr # Some handy unit shortcuts

    # # stdlib
    using SparseArrays
    using LinearAlgebra

    @info "Loading catalog"
    catalog = CSV.read(download("https://storage.googleapis.com/cmip6/cmip6-zarr-consolidated-stores.csv"), DataFrame)
    store = filter(catalog) do row
        row.source_id == "ACCESS-ESM1-5" && row.member_id=="r1i1p1f1" && row.variable_id∈("areacello","volcello") && row.experiment_id=="historical"
    end
    @info "Finding ztores for umo, vmo, mlotst, areacello, and volcello"
    common_filter(row) = row.source_id == "ACCESS-ESM1-5" && row.member_id == "r1i1p1f1" && row.experiment_id == "historical"
    @show areacello_store = filter(row -> row.variable_id == "areacello" && row.table_id == "Ofx" && common_filter(row), catalog).zstore[1]
    @show volcello_store = filter(row -> row.variable_id == "volcello" && row.table_id == "Ofx" && common_filter(row), catalog).zstore[1]
    @show mlotst_store = filter(row -> row.variable_id == "mlotst" && row.table_id == "Omon" && common_filter(row), catalog).zstore[1]
    @show umo_store = filter(row -> row.variable_id == "umo" && row.table_id == "Omon" && common_filter(row), catalog).zstore[1]
    @show vmo_store = filter(row -> row.variable_id == "vmo" && row.table_id == "Omon" && common_filter(row), catalog).zstore[1]

    @info "Lazily load data using Zarr"
    umo_zarr = zopen(umo_store, consolidated=true)
    vmo_zarr = zopen(vmo_store, consolidated=true)
    mlotst_zarr = zopen(mlotst_store, consolidated=true)
    volcello_zarr = zopen(volcello_store, consolidated=true)
    areacello_zarr = zopen(areacello_store, consolidated=true)

    @info "Loading umo, vmo, and mlotst at first time step (to limit download size)"
    umo_ds = open_dataset(umo_zarr)[time = 1]
    display(umo_ds)
    vmo_ds = open_dataset(vmo_zarr)[time = 1]
    display(vmo_ds)
    mlotst_ds = open_dataset(mlotst_zarr)[time = 1]
    display(mlotst_ds)
    @info "Loading fixed volcello and areacello"
    volcello_ds = open_dataset(volcello_zarr)
    display(volcello_ds)
    areacello_ds = open_dataset(areacello_zarr)
    display(areacello_ds)

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
        @test Ttest isa SparseMatrixCSC{Float64, Int}
		@info "  $Tsym:"
		τdiv = ustrip(Myr, norm(e1) / norm(Ttest * e1) * s)
		@info "    div: $(round(τdiv, sigdigits=2)) Myr"
		τvol = ustrip(Myr, norm(v) / norm(Ttest' * v) * s)
		@info "    vol: $(round(τvol, sigdigits=2)) Myr"
	end

    # tests if diagonal elements are > 0 and off-diagonal are < 0.
    diagT = sparse(Diagonal(T))

    @test all(diagT.nzval .> 0)
    @test all((T - diagT).nzval .< 0)



end