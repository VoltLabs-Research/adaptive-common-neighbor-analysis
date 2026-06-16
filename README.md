# Adaptive Common Neighbor Analysis

Runs adaptive common neighbor analysis (ACNA/CNA), exports per-atom structure type, and generates the cluster-graph artifacts consumed by OpenDXA.

## Install

```bash
vpm install @voltlabs/adaptive-common-neighbor-analysis
```

## CLI

```bash
common-neighbor-analysis <input_dump> [output_base] [options]
```

| Argument | Required | Default | Description |
|---|---|---|---|
| `<input_dump>` | yes | — | Input LAMMPS dump. |
| `[output_base]` | no | derived from input | Base path for output files. |
| `--crystal_structure <type>` | no | `FCC` | Reference crystal structure: `BCC`, `FCC`, `HCP`, `CUBIC_DIAMOND`, `HEX_DIAMOND`. |
| `--dissolve_small_clusters` | no | `false` | Mark small clusters as `OTHER` after clustering. |

## Exports

| Output file | Exposure | Exporter → artifact |
|---|---|---|
| `{output_base}_atoms.parquet` | Structure Identification | AtomisticExporter → glb |
| `{output_base}_cna_analysis.parquet` | CNA Analysis | — |
| `{output_base}_atoms.parquet` | Structure Counts Chart | ChartExporter → chart-png |
| `{output_base}_clusters.table` | Clusters Table | — |
| `{output_base}_cluster_transitions.table` | Clusters Transitions | — |
| `{output_base}_neighbor_lattice.parquet` | Neighbor Lattice | — |

---

Full input contract and examples: https://docs.voltcloud.dev/docs/plugins/adaptive-common-neighbor-analysis
