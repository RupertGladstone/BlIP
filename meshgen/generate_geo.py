#!/usr/bin/env python3
"""Generate a Gmsh .geo file for the MISMIP+ domain.

This is a thin shell -- the real implementation moved to the shared
elmerice_pytools package (github/local repo: ~/Source/elmerice-pytools),
alongside Fimbul's own (unrelated -- real streamline-traced, not synthetic)
.geo generator, so both projects' Elmer/VTU tooling lives in one tested
place. See elmerice_pytools.cli.mismip_plus_geo for the full docstring/
options, or run this script with --help.

Requires: pip install -e ~/Source/elmerice-pytools
"""
from elmerice_pytools.cli.mismip_plus_geo import main

if __name__ == "__main__":
    main()
