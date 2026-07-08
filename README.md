
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

#### SSA_r.sif

"r" for "restart", though technically it isn't a restart. It uses the griddatareader to read in a spunup geometry from a previous simulation. Apart from that it is essentially the same as SSA.sif

#### BP.sif


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


## Example workflow

Personal notes relating to Rupert's laptop using lmod and conda; just replace the conda and module commands with whatever is relevant for your environment.

I generated my coarse and medium resolution meshes like this:
```bash
meshgen/build_mesh.sh mesh/coarse_test 5 --dx-refined 2000 --dx-background 4000 --transition 40000
meshgen/build_mesh.sh mesh/medium_test 5 --dx-refined 1000 --dx-background 2000 --transition 50000
```

Then I ran a spinup on the coarse mesh (after copying the mesh files from the meshgen to spinup directory).
(The module dependency automatically loads the conda environment)
```bash
module load elmer/devel
cp spinup
mpirun -np 5 ElmerSolver SSA.sif > output.txt
```

Then process the spunup geometry into a netcdf file:
```bash
conda activate outflow 
python3 tools/pvtu_to_netcdf.py spinup/VTUoutputs/ssa_t0301.pvtu  --variables h zs zb groundedmask  --dx 1000  --output spinup/VTUoutputs/CoarseSpunupGeom.nc
```

Then restart with the medium resolution mesh:
```bash
module load elmer/devel
mpirun -np 5 ElmerSolver SSA_r.sif > output.txt
```

Watching thickness NRM evolve over time:
```bash
tail -f output.txt | grep -P 'Time|(?=.*NRM)(?=.*thickness).*' 
```

Or to check if it has already finished:
```bash
grep -P 'Time|(?=.*NRM)(?=.*thickness).*' output.txt
```

Then process the spunup geometry into a netcdf file (same as before but worth doing again because of the improved resolution):
```bash
conda activate outflow 
python3 tools/pvtu_to_netcdf.py spinup/VTUoutputs/ssa_r_t1001.pvtu  --variables h zs zb groundedmask ssavelocity  --dx 1000  --output spinup/MediumSpunupGeom.nc
```

```bash
--- Elmer sif snippet to read this back in (GridDataReader) ---
Solver N
  Equation = "GridDataReader"
  Procedure = "GridDataReader" "GridDataReader"
  Filename = File "MediumSpunupGeom.nc"
  X Dim Name = String "x"
  Y Dim Name = String "y"
  Variable 1 = String "h"
  Target Variable 1 = String "<elmer variable for h>"
  Variable 2 = String "zs"
  Target Variable 2 = String "<elmer variable for zs>"
  Variable 3 = String "zb"
  Target Variable 3 = String "<elmer variable for zb>"
  Variable 4 = String "groundedmask"
  Target Variable 4 = String "<elmer variable for groundedmask>"
  Variable 5 = String "ssavelocity_1"
  Target Variable 5 = String "<elmer variable for ssavelocity_1>"
  Variable 6 = String "ssavelocity_2"
  Target Variable 6 = String "<elmer variable for ssavelocity_2>"
End
```

Please make your own judgments on achieving steady state and suitable resolution. I can't guarantee that run lengths in these sif files are sufficient to achieve steady state.

I consider MediumSpunupGeom.nc to be a viable spun up SSA geometry for starting MISMIP+ standard experiments. Can also be used as a starting point for spinning up a BP run.

Using the spun up geometry to kick off a BP run:
```bash
module load elmer/devel
elmerf90 ../HydrostaticNSVec.F90 -o HydrostaticNSVec_local.so
mpirun -np 5 ElmerSolver BP_r.sif > output.txt
```
