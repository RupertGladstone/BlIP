
# BlIP
Blatter-Pattyn ISMIP Plus setup for Elmer/Ice

## Solver input files

SSA.sif

BP.sif


## Mesh generation

There's a python script to create a geo file for GMSH to turn into a mesh.
See options:
python meshgen/generate_geo.py --help

There's an example pipeline here:
meshgen/build_mesh.sh
meshgen/build_mesh.sh mesh/coarse_test 5 --dx-refined 2000 --dx-background 4000 --transition 40000


## VTU to Netcdf

We can use a script to convert pvtu files to gridded netcdf (linear interpolation happens internally).
This can make it easier to access a spun up state when starting on a new mesh.
Example:
python3 tools/pvtu_to_netcdf.py spinup/VTUoutputs/ssa_t0002.pvtu  --variables h  --dx 4005  --output spinup/VTUoutputs/h_4km.nc
