#!/usr/bin/env python3
"""Regrid an Elmer VTU/PVTU output onto a structured grid and write NetCDF,
for use as an initial condition on a different mesh (e.g. one with different
resolution or centerline refinement) via Elmer/Ice's GridDataReader solver.

Accepts either a .pvtu (parallel, reads the referenced .vtu pieces
automatically) or a plain serial .vtu.

Two modes:
  --list-variables    Print available point-data fields and the detected
                       x/y bounding box, then exit without writing anything.
  --output PATH        Regrid the selected --variables (default: all) onto a
                       structured grid and write them to a NetCDF file.

Grid construction: the x/y extent is taken from the mesh's own bounding box.
If --dx (and optionally --dy) is not given, the domain is split into
--ncells (default 100) cells per axis, dividing it exactly. If --dx/--dy is
given and does not divide the extent exactly, the number of cells is rounded
up and the leftover is split evenly as padding on both sides, so the output
grid fully encompasses the mesh (--dy defaults to --dx, i.e. square cells,
if only --dx is given).

Points on the padded/output grid that fall outside the source mesh (VTK's
"vtkValidPointMask" == 0) are filled by nearest-neighbor extrapolation from
the nearest valid grid point, rather than left at VTK's raw (meaningless,
typically 0) value.

Vector fields (e.g. Elmer's 2-DOF SSAVelocity, which VTK output pads to 3
components with an all-zero third) are split into separate scalar NetCDF
variables per component, named "<field>_1", "<field>_2", ... matching the
convention used to reference vector components in a sif (e.g. "SSAVelocity
1"). An all-zero trailing third component is dropped.

BP/FS output is a 3D (extruded) volumetric mesh, not the inherently-2D mesh
SSA writes, so there is no single well-defined x/y surface to regrid unless
one is extracted first. If the input mesh contains volumetric cells, a
single boundary is extracted (by Elmer's "GeometryIds" cell-data field,
written when "Save Geometry Ids = Logical True" is set on the output
solver) before regridding. Use --surface-id to choose which one; if not
given on a volumetric mesh it defaults to 104 (Elmer's convention for
Boundary Condition 4, which is "lower_surface" in BP_r.sif/FS_r.sif).
Elmer's own convention (fem/src/modules/ResultOutputSolve/VtuOutputSolver.F90):
a boundary's GeometryIds value is 100 + <its Boundary Condition number in
the sif> (e.g. Boundary Condition 5 -> 105), unless "BC Id Offset" was set
to something other than the default 100 in the sif's output solver.
"""
import argparse
import math

import numpy as np
import pyvista as pv
from netCDF4 import Dataset
from scipy.interpolate import NearestNDInterpolator

# Elmer's default BC Id Offset (see VtuOutputSolver.F90): a boundary's
# GeometryIds value is 100 + <Boundary Condition number>. 104 = Boundary
# Condition 4 = "lower_surface" in BP_r.sif/FS_r.sif.
DEFAULT_SURFACE_ID = 104

# VTK cell type codes for volumetric (3D) cells likely to appear in an Elmer
# extruded mesh (linear and quadratic tetrahedra, hexahedra, wedges,
# pyramids). Used only to decide whether surface extraction is needed.
_VOLUMETRIC_CELL_TYPES = {10, 12, 13, 14, 24, 25, 26, 27, 29, 30, 31, 32}


def read_mesh(path):
    return pv.read(path)


def sanitize_name(name):
    return name.strip().replace(" ", "_")


def is_volumetric(mesh):
    """True if the mesh contains any 3D (volumetric) cells."""
    return any(int(t) in _VOLUMETRIC_CELL_TYPES for t in mesh.celltypes)


def extract_surface_by_id(mesh, geometry_id):
    """Return the subset of mesh whose "GeometryIds" cell-data field equals
    geometry_id (an Elmer boundary, e.g. 104 for Boundary Condition 4)."""
    if "GeometryIds" not in mesh.cell_data:
        raise KeyError(
            "Mesh has no 'GeometryIds' cell data - re-run the Elmer output "
            "solver with 'Save Geometry Ids = Logical True', or pass "
            "--surface-id only when that is available.")
    ids = mesh.cell_data["GeometryIds"]
    mask = ids == geometry_id
    if not np.any(mask):
        available = sorted(set(int(i) for i in ids))
        raise ValueError(
            f"No cells found with GeometryIds == {geometry_id}. "
            f"Available GeometryIds in this file: {available}")
    surface = mesh.extract_cells(mask)
    print(f"extracted surface GeometryIds={geometry_id}: "
          f"{surface.n_points} points, {surface.n_cells} cells")
    return surface


def list_variables(mesh):
    xmin, xmax, ymin, ymax, _, _ = mesh.bounds
    print(f"bounding box: x in [{xmin}, {xmax}], y in [{ymin}, {ymax}]")
    print(f"n points: {mesh.n_points}")
    print("point data fields:")
    for name in mesh.point_data.keys():
        arr = mesh.point_data[name]
        ncomp = arr.shape[1] if arr.ndim > 1 else 1
        print(f"  {name!r}: {ncomp} component(s), dtype={arr.dtype}")


def build_axis(vmin, vmax, dx, ncells_default):
    extent = vmax - vmin
    if dx is None:
        n_cells = ncells_default
        dx = extent / n_cells
        pad = 0.0
    else:
        n_cells = max(1, math.ceil(extent / dx))
        pad = n_cells * dx - extent
    v0 = vmin - pad / 2.0
    n_points = n_cells + 1
    values = v0 + dx * np.arange(n_points)
    return values


def split_components(name, arr):
    """Yield (output_name, 1D array) pairs for a point-data field, splitting
    vector fields into scalar components and dropping an all-zero padded
    third component."""
    if arr.ndim == 1:
        yield sanitize_name(name), arr
        return
    ncomp = arr.shape[1]
    if ncomp == 3 and np.all(arr[:, 2] == 0):
        ncomp = 2
    for i in range(ncomp):
        yield f"{sanitize_name(name)}_{i+1}", arr[:, i]


def regrid(mesh, variables, dx, dy, ncells):
    xmin, xmax, ymin, ymax, _, _ = mesh.bounds
    x = build_axis(xmin, xmax, dx, ncells)
    y = build_axis(ymin, ymax, dy if dy is not None else dx, ncells)

    nx, ny = len(x), len(y)
    grid = pv.ImageData(
        dimensions=(nx, ny, 1),
        spacing=(x[1] - x[0], y[1] - y[0], 1.0),
        origin=(x[0], y[0], 0.0),
    )

    # .sample() does a real 3D point-in-cell test, so a query grid fixed at
    # z=0 mostly misses a mesh whose actual z (e.g. a BP/FS surface's zb/zs,
    # often hundreds of metres away from 0) isn't ~0 everywhere - almost every
    # point then falls back to the nearest-neighbor fill below, which is only
    # meant for a few padding points. We only want x/y interpolation here (z
    # is *data*, e.g. one of the fields being regridded, not geometry that
    # should matter for it), so sample against a flattened copy instead.
    mesh_flat = mesh.copy()
    mesh_flat.points[:, 2] = 0.0
    result = grid.sample(mesh_flat)
    valid = result["vtkValidPointMask"].astype(bool)
    grid_xy = result.points[:, :2]

    fields = {}
    for name in variables:
        if name not in mesh.point_data:
            raise KeyError(f"{name!r} not found in point data; available: "
                            f"{list(mesh.point_data.keys())}")
        for out_name, values in split_components(name, result[name]):
            if not np.all(valid):
                filler = NearestNDInterpolator(grid_xy[valid], values[valid])
                values = values.copy()
                values[~valid] = filler(grid_xy[~valid])
            fields[out_name] = values.reshape(ny, nx)
    return x, y, fields


def write_netcdf(path, x, y, fields, x_name, y_name):
    with Dataset(path, "w") as nc:
        nc.createDimension(x_name, len(x))
        nc.createDimension(y_name, len(y))
        xv = nc.createVariable(x_name, "f8", (x_name,))
        yv = nc.createVariable(y_name, "f8", (y_name,))
        xv[:] = x
        yv[:] = y
        for name, values in fields.items():
            var = nc.createVariable(name, "f8", (y_name, x_name))
            var[:, :] = values


def print_sif_snippet(fields, x_name, y_name, netcdf_path):
    print("\n--- Elmer sif snippet to read this back in (GridDataReader) ---")
    print('Solver N')
    print('  Equation = "GridDataReader"')
    print('  Procedure = "GridDataReader" "GridDataReader"')
    print(f'  Filename = File "{netcdf_path}"')
    print(f'  X Dim Name = String "{x_name}"')
    print(f'  Y Dim Name = String "{y_name}"')
    for i, name in enumerate(fields, start=1):
        print(f'  Variable {i} = String "{name}"')
        print(f'  Target Variable {i} = String "<elmer variable for {name}>"')
    print('End\n')


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input", help="input .pvtu or .vtu file")
    p.add_argument("--list-variables", action="store_true",
                    help="list available point-data fields and the bounding "
                         "box, then exit without writing anything")
    p.add_argument("--variables", nargs="+", default=None,
                    help="point-data fields to regrid (default: all)")
    p.add_argument("--output", "-o", help="output NetCDF file path")
    p.add_argument("--dx", type=float, default=None,
                    help="grid cell size in x, m (default: split extent into "
                         "--ncells cells)")
    p.add_argument("--dy", type=float, default=None,
                    help="grid cell size in y, m (default: same as --dx if "
                         "given, else split extent into --ncells cells)")
    p.add_argument("--ncells", type=int, default=100,
                    help="default cells per axis when --dx/--dy are not "
                         "given (default 100)")
    p.add_argument("--x-name", default="x", help="output x dimension/variable name")
    p.add_argument("--y-name", default="y", help="output y dimension/variable name")
    p.add_argument("--surface-id", type=int, default=None,
                    help="Elmer GeometryIds value of the boundary to extract "
                         "before regridding (100 + Boundary Condition number, "
                         "e.g. 104 for Boundary Condition 4). Only relevant "
                         "for volumetric (BP/FS) output; ignored (with a "
                         "warning) if given for an already-2D mesh. Defaults "
                         f"to {DEFAULT_SURFACE_ID} on a volumetric mesh if "
                         "not given.")
    return p.parse_args()


def main():
    args = parse_args()
    mesh = read_mesh(args.input)

    volumetric = is_volumetric(mesh)
    if args.surface_id is not None:
        if not volumetric:
            print(f"warning: --surface-id {args.surface_id} given but input "
                  "mesh has no volumetric cells (already a 2D surface) - "
                  "ignoring and using the whole mesh")
        else:
            mesh = extract_surface_by_id(mesh, args.surface_id)
    elif volumetric:
        print(f"input is a volumetric (3D) mesh; no --surface-id given, "
              f"defaulting to {DEFAULT_SURFACE_ID} (see --help)")
        mesh = extract_surface_by_id(mesh, DEFAULT_SURFACE_ID)

    if args.list_variables:
        list_variables(mesh)
        return

    if not args.output:
        raise SystemExit("--output is required unless --list-variables is given")

    variables = args.variables if args.variables else list(mesh.point_data.keys())
    x, y, fields = regrid(mesh, variables, args.dx, args.dy, args.ncells)
    write_netcdf(args.output, x, y, fields, args.x_name, args.y_name)

    xmin, xmax, ymin, ymax, _, _ = mesh.bounds
    print(f"source bounds: x in [{xmin}, {xmax}], y in [{ymin}, {ymax}]")
    print(f"output grid: {len(x)} x {len(y)} points, "
          f"x in [{x[0]}, {x[-1]}], y in [{y[0]}, {y[-1]}]")
    print(f"wrote {args.output}")
    print_sif_snippet(fields, args.x_name, args.y_name, args.output)


if __name__ == "__main__":
    main()
