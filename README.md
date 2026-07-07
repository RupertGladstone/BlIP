
# BlIP

Blatter-Pattyn ISMIP Plus setup for Elmer/Ice

Notes for viewing this readme: \
The markdown should automatically be interpreted when viewing this file on the github website. \
At a Linux terminal glow is recommended.
```bash
glow -w 0 README.md
```

## Solver input files

#### SSA.sif

Transient 2D simulation using Shallow Shelf Approximation (SSA).

### SSA_r.sif

"r" for "restart", though technically it isn't a restart. It uses the griddatareader to read in a spunup geometry from a previous simulation. Apart from that it is essentially the same as SSA.sif

### BP.sif


## Mesh generation

There's a python script to create a geo file for GMSH to turn into a mesh.
See options:
python meshgen/generate_geo.py --help

There's an example pipeline here:
```bash
meshgen/build_mesh.sh
meshgen/build_mesh.sh mesh/coarse_test 5 --dx-refined 2000 --dx-background 4000 --transition 40000
```


## VTU to Netcdf

We can use a script to convert pvtu files to gridded netcdf (linear interpolation happens internally).
This can make it easier to access a spun up state when starting on a new mesh.

Example:
```bash
python3 tools/pvtu_to_netcdf.py spinup/VTUoutputs/ssa_t0002.pvtu  --variables h  --dx 4005  --output spinup/VTUoutputs/h_4km.nc
```




## Spin up

Personal notes relating to Rupert's laptop using lmod and conda; just replace the conda and module commands with whatever is relevant for your environment.

I ran SSA on the mesh generated like this:
meshgen/build_mesh.sh mesh/coarse_test 5 --dx-refined 2000 --dx-background 4000 --transition 40000

(The module dependency automatically loads the conda environment)
module load elmer/devel
```bash
cp spinup
mpirun -np 5 ElmerSolver SSA.sif
```


Then process the geometry into netcdf file:
```bash
conda activate outflow 
python3 tools/pvtu_to_netcdf.py spinup/VTUoutputs/ssa_t0301.pvtu  --variables h zs zb groundedmask  --dx 1000  --output spinup/VTUoutputs/CoarseSpunupGeom.nc
```

meshgen/build_mesh.sh mesh/medium_test 5 --dx-refined 1000 --dx-background 2000 --transition 50000

and restart:
```bash
module load elmer/devel
mpirun -np 5 ElmerSolver SSA_r.sif
```

mpirun -np 5 ElmerSolver SSA.sif > output.txt

Watching thickness NRM evolve over time:
```bash
tail -f output | grep -P 'Time|(?=.*NRM)(?=.*thickness).*' 
```

grep -P 'Time|(?=.*NRM)(?=.*thickness).*' output.txt       
