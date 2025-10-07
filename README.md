# OceanTransportMatrixBuilder

*A Julia package to build ocean transport matrices from CMIP model output.*

[![Build Status](https://github.com/TMIP-code/OceanTransportMatrixBuilder.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/TMIP-code/OceanTransportMatrixBuilder.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13864392.svg)](https://doi.org/10.5281/zenodo.13864392)


> [!WARNING]
> This is work in progress. Breaking changes expected.

The purpose of this package is to build transport matrices from standard CMIP model output as part of the Transport Matrix Intercomparison Project (TMIP).

By rearranging the 3D grid of the ocean into a vector, the divergence of the flow of any tracer can be conveniently expressed in matrix form.
That is, for the vector ***x*** representing the 3D tracer concentrations *x*(***r***) at locations ***r***, the divergence of the flow is linear in *x* and can be represented by **T** ***x*** as a transport matrix **T** acting on ***x***.
The code of this package is an extension of the work of Matt Chamberlain, who built transport matrices from the ACCESS1.3 model[^Chamberlain_etal_2019], which has been successfully used in multiple following projects[^Holzer_etal_2020][^Pasquier_etal_2023][^Pasquier_etal_2024a][^Pasquier_etal_2024b].

The main application driving this project is for the validation of marine Carbon Dioxide Removal (mCDR) by computing the timescales and pathways for water in the deep ocean to reemerge to the surface[^DeVries_etal_2012][^Siegel_etal_2021].
A preprint, titled "The sequestration efficiency of the deep ocean", has been submitted to GRL[^Pasquier_etal_2025].
However, these matrices are useful in a number of contexts, e.g., for avoiding spin ups, optimization, or novel diagnostics[^John_et_al_2020].
The motivation for sharing this package is thus to facilitate the use of novel diagnostics across CMIP models.

## Example use

The intended use would be simply to feed the ocean transport output from a given CMIP model to the functions in this package to build the desired transport matrix.
After the required NetCDF files have been created (see some preliminary but working Python code for creating such files in [TMIP-code/notebooks](https://github.com/TMIP-code/notebooks))

```julia
using NetCDF
using YAXArrays
using OceanTransportMatrixBuilder

inputdir = "/Users/benoitpasquier/Data/TMIP/data/ACCESS-ESM1-5/historical/r1i1p1f1/Jan1990-Dec1999"
# The path above is for my local machine. See test/onlinebuild.jl for an example using data from the cloud.

# Load datasets
umo_ds = open_dataset(joinpath(inputdir, "umo.nc"))
vmo_ds = open_dataset(joinpath(inputdir, "vmo.nc"))
mlotst_ds = open_dataset(joinpath(inputdir, "mlotst.nc"))
volcello_ds = open_dataset(joinpath(inputdir, "volcello.nc"))
areacello_ds = open_dataset(joinpath(inputdir, "areacello.nc"))

# Load the required variables
# Note: mass transport is expected to live on the (u,v) coordinates of an Arakawa C-grid
umo = umo_ds.umo # mass transport across "east" cell face
vmo = vmo_ds.vmo # mass transport across "north" cell face
mlotst = mlotst_ds.mlotst # MLD
areacello = areacello_ds.areacello # top area of cells
volcello = volcello_ds.volcello # volume of cells
lon = volcello_ds.lon # longitude of cell centers
lat = volcello_ds.lat # latitude of cell centers
lev = volcello_ds.lev # depth of cell centers
lon_vertices = volcello_ds.lon_verticies # cell vertices
lat_vertices = volcello_ds.lat_verticies # cell vertices
# "verticies" name is an xmip bug: https://github.com/jbusecke/xMIP/issues/369
# (my local data was preprocessed with xmip)

# Make the required data from grid geometry
gridmetrics = makegridmetrics(; areacello, volcello, lon, lat, lev, lon_vertices, lat_vertices)

# Make the indices for going back and forth between 3D and 1D
indices = makeindices(gridmetrics.v3D)

# Make arrays of the flux on each face for each grid cell
ϕ = facefluxesfrommasstransport(; umo, vmo, gridmetrics, indices)

# Some parameter values
ρ = 1035.0    # density (kg/m^3)
κH = 500.0    # horizontal diffusivity (m^2/s)
κVML = 0.1    # mixed-layer vertical diffusivity (m^2/s)
κVdeep = 1e-5 # background vertical diffusivity (m^2/s)

# Make the transport matrix (should take a few seconds)
(; T) = transportmatrix(; ϕ, mlotst, gridmetrics, indices, ρ, κH, κVML, κVdeep)
```

That's it! You've got yourself the transport matrix of your dreams!

> [!WARNING]
> This does not work for all CMIP models! See below for a list of models that pass the rudimentary tests.



## List of models tested

```
ACCESS1-0
ACCESS1-3
ACCESS-ESM1-5
ACCESS-CM2
```

> [!TIP]
> The `test/` directory contains up-to-date examples for building the transport matrix and doing some simple calculations.

## Citation



This code is © Benoît Pasquier (2025) and contributors, and it is made available under the MIT license enclosed with the software.

Over and above the legal restrictions imposed by this license, if you use this software for an academic publication then you are obliged to provide proper attribution.
This can be to this code directly,

> Benoît Pasquier (2025). **OceanTransportMatrixBuilder.jl: A Julia package to build ocean transport matrices from CMIP model output**. *Zenodo*. doi: [10.5281/zenodo.13864392](https://doi.org/10.5281/zenodo.13864392).

or to the paper that it was developed for:

> Pasquier, B., Matear, R. J., Chamberlain, M. A., Ziehn, T., Hutchinson, D. K., Primeau, F. W., et al. (2025). **The sequestration efficiency of the deep ocean**. *Geophysical Research Letters*, 52, e2025GL116799. doi: [10.1029/2025GL116799](https://doi.org/10.1029/2025GL116799).


but, ideally, both.
You can also find the citation(s) in BibTeX format in the `CITATION.bib` file.


## Acknowledgements

This package improves on published work[^Chamberlain_etal_2019].
The original Julia code was essentially translated from code by Matt Chamberlain at CSIRO.
Benoît Pasquier (@briochemc) acknowledges funding through CSIRO's CarbonLock Future Science Platform, supervised by Richard Matear at CSIRO.



[^Chamberlain_etal_2019]: [Chamberlain et al. (2019)](10.1016/j.ocemod.2019.01.005) Transport matrices from standard ocean-model output and quantifying circulation response to climate change.
[^Holzer_etal_2020]: [Holzer et al. (2020)](10.1029/2020JC016414) Climate-driven changes in the ocean's ventilation pathways and time scales diagnosed from transport matrices.
[^Pasquier_etal_2023]: [Pasquier et al. (2023)](10.5194/bg-20-2985-2023) Optimal parameters for the ocean's nutrient, carbon, and oxygen cycles compensate for circulation biases but replumb the biological pump.
[^Pasquier_etal_2024a]: [Pasquier et al. (2024a)](10.5194/bg-21-3373-2024) The biological and preformed carbon pumps in perpetually slower and warmer oceans.
[^Pasquier_etal_2024b]: [Pasquier et al. (2024b)](10.1029/2024JC021043) Deoxygenation and Its Drivers Analyzed in Steady State for Perpetually Slower and Warmer Oceans.
[^John_et_al_2020]: [John et al. (2020)](10.1016/j.chemgeo.2019.119403) AWESOME OCIM: A simple, flexible, and powerful tool for modeling elemental cycling in the oceans.
[^DeVries_etal_2012]: [DeVries et al. (2012)](10.1029/2012GL051963) The sequestration efficiency of the biological pump.
[^Siegel_etal_2021]: [Siegel et al. (2021)](10.1088/1748-9326/ac0be0) Assessing the sequestration time scales of some ocean-based carbon dioxide reduction strategies.
[^Pasquier_etal_2025]: [Pasquier et al. (2025)](https://doi.org/10.22541/essoar.174646362.20144320/v1) The sequestration efficiency of the deep ocean.