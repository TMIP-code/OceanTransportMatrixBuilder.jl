using Pkg
Pkg.activate(".")
Pkg.instantiate()
using TestEnv
TestEnv.activate();

using Revise

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
inputdir = "$(ENV["HOME"])/Data/TMIP/data/$model/historical/$member/Jan1990-Dec1999"

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

# Make makegridmetrics
gridmetrics = makegridmetrics(; areacello, volcello, lon, lat, lev, lon_vertices, lat_vertices)
(; lon_vertices, lat_vertices, Z3D, v3D, zt, thkcello) = gridmetrics

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
ct = map(gsw_ct_from_pt, so, thetao)
ρ = map(gsw_rho, so, ct, Z3D)
# TODO: Check if this is correct usage of gsw functions!
# Alternatively use a fixed density:
# ρ = 1035.0    # kg/m^3

# Below is commented out but should eventually be teseted for neutral/potential density
# @show nanmean(ct)
ρθ = gsw_rho.(so, ct, 1000)
# @show nanmean(ρθ)
# from MATLAB GSW toolbox:
# gsw_rho.(so, ct, p)
# so = Absolute Salinity (g/kg)
# ct = Conservative Temperature (ITS-90) (°C)
# p = sea pressure (dbar) (here using 0 pressure to get potential density

# Diffusivites
κH = 500.0    # m^2/s
κVML = 0.1    # m^2/s
κVdeep = 1e-5 # m^2/s

# Make indices
indices = makeindices(gridmetrics.v3D)

@test all(.!isnan.(ρ[indices.wet3D])) == true

# Make fuxes from all directions
ϕ = facefluxesfrommasstransport(; umo, vmo, gridmetrics, indices)

# Make fuxes from all directions from velocities
ϕ_bis = facefluxesfromvelocities(; uo, uo_lon, uo_lat, vo, vo_lon, vo_lat, gridmetrics, indices, ρ)

# Make transport matrix
(; T, Tadv, TκH, TκVML, TκVdeep) = transportmatrix(; ϕ, mlotst, gridmetrics, indices, ρ, κH, κVML, κVdeep)

Tsyms = (:T, :Tadv, :TκH, :TκVML, :TκVdeep)
for Ttest in (T, Tadv, TκH, TκVML, TκVdeep)
    @test Ttest isa SparseMatrixCSC{Float64, Int}
end

# Make transport matrix without upwind this time
(; T, Tadv, TκH, TκVML, TκVdeep) = transportmatrix(; ϕ, mlotst, gridmetrics, indices, ρ, κH, κVML, κVdeep, upwind = false)


# @profview transportmatrix(; ϕ, mlotst, gridmetrics, indices, ρ, κH, κVML, κVdeep)
# @profview_allocs transportmatrix(; ϕ, mlotst, gridmetrics, indices, ρ, κH, κVML, κVdeep) sample_rate=0.9
