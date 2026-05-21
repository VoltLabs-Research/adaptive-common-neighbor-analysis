# CommonNeighborAnalysis

`CommonNeighborAnalysis` classifies atoms by local crystal environment and exports the reconstructed state consumed by DXA-compatible downstream tools.

## One-Command Install

```bash
curl -sSL https://raw.githubusercontent.com/VoltLabs-Research/CoreToolkit/main/scripts/install-plugin.sh | bash -s -- CommonNeighborAnalysis
```

## Build from source

Requires [Conan 2.x](https://docs.conan.io/2/installation.html), CMake 3.20+, and a C++23 compiler (GCC 14+ or Clang 17+).

### Prerequisites

The following Conan packages must be available in your local cache:

- `coretoolkit/1.0.0` (from the `CoreToolkit` repository)
- `structure-identification/1.0.0` (from the `StructureIdentification` repository)

For each dependency, clone its repository and create the package:

```bash
conan create <path-to-dependency-repo> --build=missing -o "hwloc/*:shared=True"
```

### Build

From the root of this repository:

```bash
conan install . -of build --build=missing -o "hwloc/*:shared=True"
cmake --preset conan-release
cmake --build build/build/Release -j
```

### Run

```bash
./build/build/Release/common-neighbor-analysis --help
```

### Package as Conan recipe

To make this plugin available as a Conan package for other projects:

```bash
conan create . --build=missing -o "hwloc/*:shared=True"
```

## CLI

Usage:

```bash
common-neighbor-analysis <lammps_file> [output_base] [options]
```

### Arguments

| Argument | Required | Description | Default |
| --- | --- | --- | --- |
| `<lammps_file>` | Yes | Input LAMMPS dump file. | |
| `[output_base]` | No | Base path for output files. | derived from input |
| `--crystal_structure <type>` | No | Input crystal structure: `FCC`, `BCC`, `HCP`, `CUBIC_DIAMOND`, `HEX_DIAMOND`. | `FCC` |
| `--dissolve_small_clusters` | No | Mark small clusters as `OTHER` after clustering. | `false` |
| `--help` | No | Print CLI help. | |
