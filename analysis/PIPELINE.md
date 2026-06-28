# CYAxiverse Parallel Pipeline

## Architecture

```
HDF5 Database (970 geoms, 50,440 axions)
         │
  create_batches.py ────► batches.json (24 batches × ~2000 axions)
         │
  ┌──────┴──────┐
  │  Run in     │  run_batch.jl (parallel via GitHub Actions / Jenkins / Make)
  │  parallel   │  Each batch: dataset gen + DM ratio + λ_self
  └──────┬──────┘
         │
  ┌──────┴──────┐
  │  Merge      │  merge_batches.jl — combines all batch outputs
  └──────┬──────┘
         │
  ┌──────┴──────┐
  │  Global     │  global_analysis.jl — stats + plots
  │  Analysis   │
  └─────────────┘
```

## What was created

| File | Purpose |
|------|---------|
| `tools/create_batches.py` | Compute optimal h11 splits for equal batch sizes |
| `scripts/run_batch.jl` | **The parallel unit**: process one h11 range end-to-end |
| `scripts/merge_batches.jl` | Combine all batch outputs into unified datasets |
| `scripts/global_analysis.jl` | Stats + plots + summary on merged data |
| `.github/workflows/pipeline.yml` | GitHub Actions: parallel matrix + merge + analysis |
| `Jenkinsfile` | Jenkins declarative pipeline (same stages, parallel) |

## Batch splits (target ~2000 axions)

```
$ python tools/create_batches.py --target 2000

 Batch   h11 range   geoms    axions
------------------------------------
     1     4-22       190      2470
     2    23-31        90      2430
     3    32-38        70      2450
     4    39-44        60      2490
     5    45-49        50      2350
     6    50-53        40      2060
     7    54-57        40      2220
     8    58-61        40      2380
     9    62-64        30      1890
    10    65-67        30      1980
    11    68-70        30      2070
    12    71-73        30      2160
    13    74-76        30      2250
    14    77-79        30      2340
    15    80-82        30      2430
    16    83-84        20      1670
    17    85-86        20      1710
    18    87-88        20      1750
    19    89-90        20      1790
    20    91-92        20      1830
    21    93-94        20      1870
    22    95-96        20      1910
    23    97-98        20      1950
    24    99-100       20      1990
------------------------------------
 Total                 970     50440
```

## Data storage

**Runners are ephemeral** — all output must go to persistent storage.

### Option A: Network mount (recommended)
Mount an NFS/SMB share or S3 bucket to `/mnt/results/` on every runner:
```
/mnt/results/
├── batches/                          # per-batch outputs
│   ├── batch_004_022/
│   │   ├── axion_dataset.jls
│   │   ├── geometry_dataset.jls
│   │   ├── dm_density_ratio_results.jls
│   │   ├── quartic_lambda_results.jls
│   │   ├── batch_info.json
│   │   └── batch.log
│   ├── batch_023_031/
│   └── ...
├── axion_dataset.jls                 # merged (all batches)
├── geometry_dataset.jls              # merged
├── dm_density_ratio_results.jls      # merged
├── quartic_lambda_results.jls        # merged
├── correlation_matrix.jls
├── pca_results.jls
├── pipeline_summary.txt
├── pipeline_merge.log
├── plots/
│   ├── hist_mass.png
│   ├── hist_decay.png
│   ├── hist_abundance.png
│   ├── scatter_mass_abundance.png
│   ├── scatter_decay_abundance.png
│   ├── scatter_mass_decay.png
│   ├── hist_light_axions.png
│   ├── hist_geometry_abundance.png
│   └── scatter_geom_DMratio_h11.png
└── logs/
    └── global_analysis.log
```

### Option B: GitHub Actions Artifacts
Jobs upload per-batch `.jls` files as artifacts; merge job downloads them.  
Limitations: 10 GB storage, 90-day retention, upload/download time.

## Running

### GitHub Actions (self-hosted runner)
```bash
# Install runner on the Mac
mkdir actions-runner && cd actions-runner
curl -o actions-runner-osx-arm64.tar.gz -L \
  https://github.com/actions/runner/releases/latest/download/actions-runner-osx-arm64.tar.gz
tar xzf actions-runner-osx-arm64.tar.gz
./config.sh --url https://github.com/YOUR_USER/YOUR_REPO --labels self-hosted,macos,docker
./run.sh
```
Then trigger from GitHub UI: Actions → CYAxiverse Parallel Pipeline → Run workflow.

### Jenkins
Configure the Jenkins pipeline to point at this repo. The `Jenkinsfile` uses
`parallel` to run all batches concurrently.

### Manual / Makefile
```bash
# Generate batch list
python tools/create_batches.py --target 2000 --format make > batches.mk

# Run 4 at a time (in separate terminals, or with GNU parallel)
for lo hi in 4-22 23-31 32-38 39-44; do
  docker exec cyax_pluto julia --project=/opt/CYAxiverse.jl/ \
    scripts/run_batch.jl --h11-min $lo --h11-max $hi --output /mnt/results &
done
wait

# Merge
docker exec cyax_pluto julia --project=/opt/CYAxiverse.jl/ \
  scripts/merge_batches.jl --input /mnt/results --output /mnt/results

# Global analysis
docker exec cyax_pluto julia --project=/opt/CYAxiverse.jl/ \
  scripts/global_analysis.jl /mnt/results
```

## Cost estimate

| Batch | Wall time (per runner) | Runners | Total wall time |
|-------|----------------------|---------|-----------------|
| Dataset gen (CYAxiverse) | ~10 min compile + ~40s work | 24 | ~11 min |
| DM ratio | ~2s | — | included |
| λ_self (HDF5) | ~5s | — | included |
| Merge | ~2s | 1 | ~2s |
| Global analysis | ~2 min (CairoMakie compile) | 1 | ~2 min |

**Total**: ~13 min wall time to process all 50,440 axions (vs ~8+ hours sequential).
