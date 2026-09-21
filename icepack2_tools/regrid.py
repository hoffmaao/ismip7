r"""The ISMIP7 AIS 8 km output grid.

Only the grid is defined here. The regridding itself is the conservative
supermesh remap in ``antarctica/scripts/write_ismip7_output.py``, which takes
each variable's units, names and fill policy from the bundled variable
request. A point-sampling regridder (``VertexOnlyMesh``) used to live in this
module with its own hand-written variable table; nothing called it, point
sampling is not what the protocol asks of the mass fluxes (discussions #13
and #22), and it stamped its files with a regridding method the submission
does not use, so it is gone.

The numbers are those of the official grid description,
``isschecker/data/gdfs/gdf_ISMIP7_AIS_08000m.txt`` in
``ismip/ISM_SimulationChecker`` (xsize = ysize = 761, xfirst = yfirst =
-3040000, xinc = yinc = 8000, EPSG:3031), and ``tests/test_ismip7_output.py``
holds them to it.
"""

import numpy as np

ISMIP7_NX = 761
ISMIP7_NY = 761
ISMIP7_X0 = -3_040_000.0
ISMIP7_X1 = 3_040_000.0
ISMIP7_Y0 = -3_040_000.0
ISMIP7_Y1 = 3_040_000.0
ISMIP7_DX = 8000.0


def ismip7_grid_coords():
    r"""Return ISMIP7 AIS grid cell center coordinates."""
    x = np.linspace(ISMIP7_X0, ISMIP7_X1, ISMIP7_NX)
    y = np.linspace(ISMIP7_Y0, ISMIP7_Y1, ISMIP7_NY)
    return x, y
