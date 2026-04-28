# ISMIP7

Ice sheet simulations for [ISMIP7](https://www.ismip6.org/ismip7/) using
[icepack2](https://github.com/icepack/icepack2) with Firedrake and tlm_adjoint.

## Structure

```
ismip7/
├── icepack2_tools/         # Reusable utilities (eikonal, grounding zone, plotting)
├── antarctica/             # Antarctic continent simulations
│   ├── scripts/            # Pipeline scripts
│   ├── data/               # Downloaded datasets (see data/README.md)
│   ├── mesh/               # Generated meshes
│   ├── results/            # Inversion and sensitivity results
│   └── figs/               # Figures
├── greenland/              # Greenland ice sheet simulations
│   ├── scripts/            # Files to get data, preprocess, run model
│   ├── data/               # Downloaded datasets and mesh outlines
│   ├── mesh/               # Generated meshes
│   ├── results/            # Results
│   └── figs/               # Figures
```

## Dependencies

- [Firedrake](https://firedrakeproject.org)
- [icepack2](https://github.com/icepack/icepack2) (Antarctica only)
- [tlm_adjoint](https://github.com/jrmaddison/tlm_adjoint) (Antarctica only)
- gmsh, rasterio, icepack, geopandas, earthaccess
- cartopy, rioxarray (Greenland only)
