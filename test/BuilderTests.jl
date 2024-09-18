

@testitem "building matrices" default_imports=false begin
    using Pkg
    Pkg.activate(@__DIR__)
    Pkg.instantiate()
    Pkg.resolve()
    using OceanTransportMatrixBuilder
    using Test
    using NetCDF
    using YAXArrays
    using Unitful
    using Unitful: cm, m, km, kg, s, d, yr, Myr, hr # Some handy unit shortcuts
    # using JLD2
    # using MAT

    # # should go into OceanTransportMatrixBuilder
    # using NaNStatistics
    # using ProgressMeter
    # using Format
    # using UnicodePlots
    # using CairoMakie
    # using GeoMakie

    # # stdlib
    # using Statistics
    # using StatsBase
    # using LinearAlgebra
    # using Dates
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

    umo = umo_ds["umo"] |> Array{Float64}
    vmo = vmo_ds["vmo"] |> Array{Float64}
    mlotst = mlotst_ds["mlotst"] |> Array{Float64}
    volcello = volcello_ds["volcello"] |> Array{Union{Missing, Float64}}
    areacello = areacello_ds["areacello"] |> Array{Float64}

    # Make ualldirs
    u = makeualldirections(; umo, vmo)

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

end