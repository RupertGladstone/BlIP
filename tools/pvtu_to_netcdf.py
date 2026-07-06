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
"""
import argparse
import math

import numpy as np
import pyvista as pv
from netCDF4 import Dataset
from scipy.interpolate import NearestNDInterpolator


def read_mesh(path):
    return pv.read(path)


def sanitize_name(name):
    return name.strip().replace(" ", "_")


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
    result = grid.sample(mesh)
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
    return p.parse_args()


def main():
    args = parse_args()
    mesh = read_mesh(args.input)

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
