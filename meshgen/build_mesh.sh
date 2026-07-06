#!/usr/bin/env bash
# Build and partition a MISMIP+ mesh: generate .geo -> mesh with gmsh ->
# convert to Elmer format with ElmerGrid -> partition with ElmerGrid/Metis.
# Fails immediately (set -e) if any step fails.
#
# Usage: build_mesh.sh <output_dir> <n_partitions> [extra generate_geo.py args...]
#
# Example (coarse test mesh, standard domain, no centerline refinement,
# 3km refined / 8km background, 5 partitions):
#   ./build_mesh.sh mesh/coarse_test 5 --dx-refined 3000 --dx-background 8000
set -euo pipefail

if [ $# -lt 2 ]; then
    echo "Usage: $0 <output_dir> <n_partitions> [extra generate_geo.py args...]" >&2
    exit 1
fi

OUTPUT_DIR=$1
N_PARTITIONS=$2
shift 2

for cmd in python3 gmsh ElmerGrid; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "error: required command '$cmd' not found on PATH" >&2
        exit 1
    fi
done

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
GEO_FILE="$OUTPUT_DIR/mismip.geo"
MSH_FILE="$OUTPUT_DIR/mismip.msh"
MESH_DIR="$OUTPUT_DIR/mesh"

mkdir -p "$OUTPUT_DIR"

echo "==> Generating .geo file"
python3 "$SCRIPT_DIR/generate_geo.py" --output "$GEO_FILE" "$@"

echo "==> Meshing with gmsh"
gmsh -2 "$GEO_FILE" -o "$MSH_FILE" -format msh2

echo "==> Converting to Elmer mesh format"
ElmerGrid 14 2 "$MSH_FILE" -out "$MESH_DIR"

echo "==> Partitioning into $N_PARTITIONS parts"
ElmerGrid 2 2 "$MESH_DIR" -metiskway "$N_PARTITIONS"

echo "==> Done."
echo "    Mesh:          $MESH_DIR"
echo "    Partitioning:  $OUTPUT_DIR/partitioning.$N_PARTITIONS"
