# ISMIP7

Ice sheet simulations for [ISMIP7](https://www.ismip6.org/ismip7/) using
[icepack2](https://github.com/icepack/icepack2) with Firedrake and tlm_adjoint.

## Structure

```
ismip7/
├── icepack2_tools/         # Reusable utilities (mesh, forcing, regrid, eikonal, grounding zone)
├── antarctica/             # Antarctic continent simulations
│   ├── scripts/            # Pipeline scripts
│   ├── data/               # Downloaded datasets (see data/README.md)
│   ├── mesh/               # Generated meshes
│   ├── results/            # Inversion and sensitivity results
│   └── figs/               # Figures
```

## Dependencies

- [Firedrake](https://firedrakeproject.org) 2026.4 (brings PETSc, MUMPS, mpi4py)
- [icepack2](https://github.com/icepack/icepack2) and [icepack](https://github.com/icepack/icepack)
- [icepack_tools](https://github.com/hoffmaao/icepack_tools) (mesh adaptation and the level-set front)
- [tlm_adjoint](https://github.com/jrmaddison/tlm_adjoint) (inversions)
- `xarray netCDF4 scipy rasterio pyproj shapely gmsh matplotlib`, `geopandas`
  for mesh generation, plus `earthaccess` and `globus-sdk` for downloads
- [isschecker](https://github.com/ismip/ISM_SimulationChecker) for submissions
  (needs Python >= 3.11, so it wants its own venv)

Not every part needs every one:
`antarctica/README.md` section 0 is a table of what to install and download for
each of inverting, running forward, adapting the mesh and submitting, and
section 0.5 covers running on a cluster.
