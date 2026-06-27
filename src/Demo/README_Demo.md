# Demo: Danish Transmission Networks and NPI Effects Pipeline

This folder provides a self-contained, runnable version of the full analysis pipeline described in the paper, using randomly generated artificial data in place of the restricted Danish registry data. Because the data are artificial, outputs do not reproduce the paper's results; the purpose is to illustrate how the code is structured and how each stage can be adapted.

---

## Folder structure

```
Demo/
├── data/                          # Pre-generated input data (see below)
├── python/
│   ├── transmission_functions.py  # Helper module (imported by both notebooks)
│   ├── demo_pipeline.ipynb        # Main analysis notebook — run this first
│   └── generate_demo_data.ipynb   # Optional: regenerate or customise the input data
├── R/
│   ├── Rt_estimation_and_NPI_effects.R     # Rt estimation and NPI effect plots
│   ├── partial_sequencing_adjustment.R     # Sequencing-proportion correction helper
│   ├── plotting_fns.R                      # Plotting helper functions
│   ├── stan_utility_fns.R                  # Stan utility functions
│   ├── zero_inflated_neg_bin.stan          # Stan model: zero-inflated negative binomial
│   ├── zero_inflated_poisson.stan          # Stan model: zero-inflated Poisson
│   ├── data_ps_ml_ms.csv                   # Transmission setting proportions
│   ├── mcmc_object_zibb_total_transmission_re_summary.csv  # Pre-run MCMC summary (NPI model)
│   └── stan_output/               # Pre-computed Stan model outputs
├── trees_output/                  # Sampled transmission trees (written by Step 2)
└── clusters_output/               # Transmission clusters (written by Step 3)
```

---

## Prerequisites

### Python

Python 3.11 is recommended. Install the required packages with:

```bash
pip install numpy pandas tqdm biopython regex hammingdist polars pyarrow scipy \
            networkx seaborn matplotlib joblib numpyro jax arviz netCDF4
```

> **Note:** `numpyro` and `jax` are required for the Negative Binomial model in Step 4.
> On Apple Silicon Macs, `jax` should be installed via `pip install jax-metal` or following
> the [JAX installation guide](https://jax.readthedocs.io/en/latest/installation.html).

### R

R 4.x is recommended. The following packages are required:

```r
install.packages(c("tidyverse", "tidybayes", "rstan", "patchwork",
                   "zoo", "MASS", "lme4", "jtools", "mgcv", "ggsci","ggtext","ggh4x"))
```

---

## Input data

All input files are pre-generated and located in `data/`. They are produced by `generate_demo_data.ipynb` (see below if you want to regenerate them).

| File | Description |
|---|---|
| `demo_sequences.fasta` | 1,000 artificial SARS-CoV-2-like sequences (length 5 nt) |
| `demo_ids.csv` | Sequence IDs corresponding to the FASTA file |
| `demo_sample_dates.csv` | Simulated test dates for each individual |
| `all_attributes.csv` | Individual-level attributes: age, address, school, workplace, vaccination status, NPI covariate codes |
| `demo_family_edgelist.csv` | Simulated family (non-household) relationship links |
| `school_types.csv` | Mapping of school IDs to school type |
| `hamming.csv` | Pre-computed pairwise Hamming distance matrix |
| `adjacency_demo.npz` | Pre-computed weighted adjacency matrix (plausible transmission network) |
| `dates_demo.npz` | Pre-computed pairwise test-date difference matrix |
| `rt_iar.csv` | Simulated reproduction number (Rt) and infection attack rate (IAR) time series |
| `serial_interval.csv` | Serial interval distribution (gamma fit) |
| `proportion_of_pcrpositives_sequenced.csv` | Simulated sequencing proportion over time |
| `node_attrs.csv` | Node-level tree output with NPI covariate codes, used as input to the R analysis |

---

## Running the pipeline

### Step 1 — Open `python/demo_pipeline.ipynb`

Launch Jupyter from the `python/` directory so that `transmission_functions.py` is on the import path:

```bash
cd python
jupyter notebook demo_pipeline.ipynb
```

The notebook is divided into four numbered sections. Run them in order from top to bottom.

---

#### Section 1: Build the plausible transmission network

Reads the artificial sequences and test dates, then constructs a weighted directed adjacency matrix. Each entry encodes the probability that two individuals form a transmission pair, based on the genetic distance (Hamming distance) between their sequences and the time elapsed between their test dates. The matrix is saved to `../data/adjacency_demo.npz`.

Key parameters (editable at the top of the section):

| Parameter | Default | Description |
|---|---|---|
| `weights` | `1` | Use weighted edges (`1`) or unweighted (`0`) |
| `sensitivity` | `0` | Shift the generation interval by −2 days for sensitivity analysis |

---

#### Section 2: Sample transmission trees

Reads the adjacency matrix and samples `n_trees = 100` plausible transmission trees. Two sampling strategies are applied:

- **Random trees** — each individual's infector is selected at random from plausible candidates, weighted by edge probability.
- **Priority-setting trees** — edges within a shared social setting (household, school, workplace, family) are up-weighted before sampling, as in the main analysis.

Each sampled tree is saved as `../trees_output/nodelist_{i}.csv`. The section also plots the proportion of infections attributed to each setting and the age-structured infector–infectee matrix.

---

#### Section 3: Generate transmission clusters

Uses the setting-attributed directed graph to identify setting-specific transmission chains. Individuals connected by edges within the same social setting are linked into a "settings network", whose weakly connected components define transmission clusters. Cluster data are written to `../clusters_output/`.

---

#### Section 4: Run the Negative Binomial model

Fits a Zero-Inflated Beta-Binomial (ZIBB) model using `numpyro`/`JAX` to estimate the effects of NPIs, vaccination status, age group, and SARS-CoV-2 variant on individual-level reproduction numbers (out-degree in the sampled trees). The demo runs this on a single tree (`nodelist_0.csv`) with shortened MCMC chains for speed.

Outputs (written to `python/`):

- `mcmc_object_zibb_total_transmission_re_results.nc` — full MCMC posterior (NetCDF format, readable with `arviz`)
- `mcmc_object_zibb_total_transmission_re_summary.csv` — posterior summary statistics

Pre-run versions of these outputs are already included so that the R section can be run independently.

---

### Step 2 — Run the R analysis (optional)

Set your working directory to `R/` and open `Rt_estimation_and_NPI_effects.R`. This script:

1. Reads `../data/node_attrs.csv` (the aggregated tree output) and fits a weekly zero-inflated negative binomial model in Stan to estimate the reproduction number (R) and overdispersion parameter (k) over time.
2. Reads `mcmc_object_zibb_total_transmission_re_summary.csv` and produces forest plots of NPI and risk-factor effects.

> **Note:** The Stan models can take 10–20 minutes to run. Pre-computed outputs are provided in `stan_output/` so that the plotting sections can be run immediately without re-fitting.

---

## Regenerating the input data (optional)

If you want to modify the artificial data (e.g. change the number of individuals, sequence length, or epidemic duration), open `python/generate_demo_data.ipynb` and edit the parameters at the top of each cell. Running the full notebook will overwrite the files in `../data/`.

Key parameters:

| Parameter | Default | Description |
|---|---|---|
| `n` | `1000` | Number of simulated individuals |
| `length` | `5` | Sequence length (nucleotides) |
| `total_days_epidemic` | `31` | Duration of the simulated epidemic |
| `gmean`, `gsd` | `4.87`, `1.98` | Mean and SD of the generation time distribution (days) |

---

## Expected run times

| Step | Estimated time |
|---|---|
| Section 1 (network construction) | < 1 minute |
| Section 2 (tree sampling, n=100) | 2–5 minutes |
| Section 3 (cluster detection) | 1–3 minutes |
| Section 4 (MCMC, demo settings) | 5–15 minutes |
| R: Stan overdispersion model | 10–20 minutes |

---

## Notes

- **Artificial data only.** The input data are randomly generated and bear no relation to Danish patient records. Outputs from this demo therefore do not reproduce any result in the paper.
- **Relative paths.** Both notebooks use paths relative to `python/` (e.g. `../data/`, `../trees_output/`). Always launch Jupyter from the `python/` directory, or adjust the paths if running from elsewhere.
- **`hamming.cpp`** in the Demo root is a placeholder and is not used by any notebook.
- **`partial_sequencing_adjustment.R`** contains a function that references hardcoded paths to the restricted Danish data environment; this function is not called anywhere in the demo and can be ignored.
