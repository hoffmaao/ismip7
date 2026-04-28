# Greenland workflow

## Get data
run ```download_data.py```. The data directory can be selected by the user but defaults to ../data

## Make meshes
run ```mesh_greenland.py```. We have not selected a final mesh, and currently there are several options to choose from. The main thing is to specify an outline file--there are two in this repository, and download_data gets one from PROMICE.

All meshes vary in resolution from min-lc to 10xmin-lc, using velocity and strain rate to adapt.

Use ```--min-lc``` to control the finest resolution considered. Use ```--num-levels``` to control how many other meshes (in steps of 2 ** i) are made.

## Load data
For the initial work it is convenient to cache some data. Do this using ```load_data.py```. That script can plot, but doing so is slow.

## Test diagnostic solves
Run ```diagnostic_solves.py``` using varying numbers of cores. This plays with solver parameters (mostly damping) which can then be used to get inversions running more smoothly.

This script tries to be clever about what tests are run depending on the number of cores.

## Test inversions
Run ```inversion_tests.py```.
