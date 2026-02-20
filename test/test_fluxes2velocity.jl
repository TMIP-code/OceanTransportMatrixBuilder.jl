# Standalone round-trip test for fluxes2velocity
# Run with: julia --project test/test_fluxes2velocity.jl

using Test
using OceanTransportMatrixBuilder
using NetCDF
using YAXArrays
using GibbsSeaWater
using DimensionalData
using NaNStatistics

model = "ACCESS-ESM1-5"
member = "r1i1p1f1"
inputdir = "$(ENV["HOME"])/Data/TMIP/data/$model/historical/$member/Jan1990-Dec1999"

# Load datasets
uo_ds = open_dataset(joinpath(inputdir, "uo.nc"))
vo_ds = open_dataset(joinpath(inputdir, "vo.nc"))
volcello_ds = open_dataset(joinpath(inputdir, "volcello.nc"))
areacello_ds = open_dataset(joinpath(inputdir, "areacello.nc"))
thetao_ds = open_dataset(joinpath(inputdir, "thetao.nc"))
so_ds = open_dataset(joinpath(inputdir, "so.nc"))

areacello = readcubedata(areacello_ds.areacello)
volcello = readcubedata(volcello_ds.volcello)
lon = readcubedata(volcello_ds.lon)
lat = readcubedata(volcello_ds.lat)
lev = volcello_ds.lev
lon_vertices = readcubedata(volcello_ds.lon_verticies)
lat_vertices = readcubedata(volcello_ds.lat_verticies)
uo = readcubedata(uo_ds.uo)
vo = readcubedata(vo_ds.vo)
uo_lon = readcubedata(uo_ds.lon)
uo_lat = readcubedata(uo_ds.lat)
vo_lon = readcubedata(vo_ds.lon)
vo_lat = readcubedata(vo_ds.lat)
thetao = readcubedata(thetao_ds.thetao) |> Array
so = readcubedata(so_ds.so) |> Array

gridmetrics = makegridmetrics(; areacello, volcello, lon, lat, lev, lon_vertices, lat_vertices)
(; Z3D, v3D) = gridmetrics
ct = gsw_ct_from_pt.(so, thetao)
ρ = gsw_rho.(so, ct, Z3D)
indices = makeindices(v3D)
(; wet3D) = indices

# Round-trip: velocity → flux → velocity
u_cgrid, _, _, v_cgrid, _, _ = OceanTransportMatrixBuilder.interpolateontodefaultCgrid(uo, uo_lon, uo_lat, vo, vo_lon, vo_lat, gridmetrics)
umo_bis, vmo_bis = velocity2fluxes(uo, uo_lon, uo_lat, vo, vo_lon, vo_lat, gridmetrics, ρ)
u_rec, v_rec = fluxes2velocity(umo_bis, vmo_bis, gridmetrics, ρ)

@test u_rec[wet3D] ≈ u_cgrid[wet3D]
@test v_rec[wet3D] ≈ v_cgrid[wet3D]

println("All tests passed!")
