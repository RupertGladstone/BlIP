#!/usr/bin/env python3
"""Regrid an Elmer VTU/PVTU output onto a structured grid and write NetCDF.

This is a thin shell -- the real implementation moved to the shared
elmerice_pytools package (github/local repo: ~/Source/elmerice-pytools),
also used by Fimbul, so both projects depend on one tested copy rather than
each carrying their own. See elmerice_pytools.cli.pvtu_to_netcdf for the
full docstring/options, or run this script with --help.

Requires: pip install -e ~/Source/elmerice-pytools
"""
from elmerice_pytools.cli.pvtu_to_netcdf import main

if __name__ == "__main__":
    main()
