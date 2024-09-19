# OceanTransportMatrixBuilder

*A Julia package to build ocean transport matrices from CMIP model output.*

[![Build Status](https://github.com/TMIP-code/OceanTransportMatrixBuilder.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/TMIP-code/OceanTransportMatrixBuilder.jl/actions/workflows/CI.yml?query=branch%3Amain)

The purpose of this package is to build transport matrices from standard CMIP model output as part of the Transport Matrix Intercomparison Project (TMIP).

By rearranging the 3D grid of the ocean into a vector, the divergence of the flow of any tracer can be conveniently expressed in matrix form.
That is, for the vector $\boldsymbol{x}$ representing the 3D tracer concentrations $\chi$, the divergence of the flow is linear in $\chi$ and can be represented by $\mathbf{T}\boldsymbol{x}$ as a transport matrix $\mathbf{T}$ acting on $\boldsymbol{x}$.

These matrices are useful in a number of contexts, e.g., for avoiding spin ups, optimization, or novel diagnostics[^John_et_al_2020][^Pasquier_etal_2023][^Pasquier_etal_2024a][^Pasquier_etal_2024b].
The motivation for writing this package is to facilitate the use of novel diagnostics across CMIP models.
The original intended application is for validating marine Carbon Dioxide Removal (mCDR) by computing the timescales and pathways for water in the deep ocean to reemerge to the surface[^DeVries_etal_2012][^Siegel_etal_2021].

> [!WARNING]
> This is work in progress. Breaking changes expected.

## Example use

The intended use would be simply to feed the ocean transport output from a given CMIP model to the functions in this package to build the desired transport matrix.
After the required NetCDF files have been created (see some preliminary but working Python code for creating such files in [TMIP-code/notebooks](https://github.com/TMIP-code/notebooks))

```julia
using NetCDF
using YAXArrays
using OceanTransportMatrixBuilder

inputdir = "/Users/benoitpasquier/Data/TMIP/data/ACCESS-ESM1-5/historical/r1i1p1f1/Jan1990-Dec1999" # <- this is the path on my machine

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
(; T) = transportmatrix(; u, mlotst, modelgrid, indices,
    ρ = 1025.0, # mean seawater density (kg/m^3)
    κH = 500.0, # horizontal diffusivity (m^2/s)
    κVML = 0.1, # mixed-layer vertical diffusivity (m^2/s)
    κVdeep = 1e-5, # background vertical diffusivity (m^2/s)
)
```

That's it! You've got yourself the transport matrix of your dreams!

> [!TIP]
> The `test/` directory should contain up-to-date examples for building the transport matrix and doing some simple calculations.

## Acknowledgements

BP is funded through CSIRO's CarbonLock Future Science Platform.


[^John_et_al_2020]: [John et al. (2020)](10.1016/j.chemgeo.2019.119403) AWESOME OCIM: A simple, flexible, and powerful tool for modeling elemental cycling in the oceans.
[^Pasquier_etal_2023]: [Pasquier et al. (2023)](10.5194/bg-20-2985-2023) Optimal parameters for the ocean's nutrient, carbon, and oxygen cycles compensate for circulation biases but replumb the biological pump.
[^Pasquier_etal_2024a]: [Pasquier et al. (2024a)](10.5194/bg-21-3373-2024) The biological and preformed carbon pumps in perpetually slower and warmer oceans.
[^Pasquier_etal_2024b]: [Pasquier et al. (2024b)](10.1029/2024JC021043) Deoxygenation and Its Drivers Analyzed in Steady State for Perpetually Slower and Warmer Oceans.
[^DeVries_etal_2012]: [DeVries et al. (2012)](10.1029/2012GL051963) The sequestration efficiency of the biological pump.
[^Siegel_etal_2021]: [Siegel et al. (2021)](10.1088/1748-9326/ac0be0) Assessing the sequestration time scales of some ocean-based carbon dioxide reduction strategies.
