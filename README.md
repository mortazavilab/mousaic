# cis/trans regulatory inference explorer

This app is an interactive viewer for exploring how gene regulation changes across tissues, cell types, and mouse strains.


## Features
- Compare regulatory composition across tissues, subtypes, and strains.
- Explore transitions in regulatory assignment using Sankey diagrams.
- Drill into one gene at a time and visualize expression by allele (P1, H1, H2, P2).
- Export selected gene expression tables for downstream analysis.

Enables  fast comparative analysis and hypothesis generation without writing custom code.

## Quick start (Docker)

Fastest path:

```bash
docker pull ghcr.io/mortazavilab/mousaic:latest
docker volume create cistrans_data
docker run --rm -p 8501:8501 \
  -e DATA_ROOT=/app/data \
  -v cistrans_data:/app/data \
  ghcr.io/mortazavilab/mousaic:latest
```

Open http://localhost:8501.

On first startup, the app downloads required data and extracts `gene_count_data/` into the volume. On later startups with the same volume, download/extract is skipped.

## Run with Docker (GHCR)

Pull image:

```bash
docker pull ghcr.io/mortazavilab/mousaic:latest
```

Create a named volume for persistent data (one-time):

```bash
docker volume create cistrans_data
```

Run container (auto-downloads data on first startup):

```bash
docker run --rm -p 8501:8501 \
  -e DATA_ROOT=/app/data \
  -v cistrans_data:/app/data \
  ghcr.io/mortazavilab/mousaic:latest
```

If port 8501 is already in use, map a different host port:

```bash
docker run --rm -p 8502:8501 \
  -e DATA_ROOT=/app/data \
  -v cistrans_data:/app/data \
  ghcr.io/mortazavilab/mousaic:latest
```

Then open http://localhost:8501 (or http://localhost:8502 if using the alternate port mapping).

Notes:

- Data files are not baked into the image.
- On first startup, the container downloads and extracts data into `/app/data`.
- On subsequent startups with the same `cistrans_data` volume, download/extract is skipped.
- If startup download fails (network or URL issue), the app shows an error and exits cleanly.

Optional: persist data in a host directory instead of a named volume.

```bash
docker run --rm -p 8501:8501 \
  -e DATA_ROOT=/app/data \
  -v "$PWD/data:/app/data" \
  ghcr.io/mortazavilab/mousaic:latest
```

### Optional: Build locally

```bash
docker build -t cistrans-viewer .
```

```bash
docker run --rm -p 8501:8501 \
  -e DATA_ROOT=/app/data \
  -v cistrans_data:/app/data \
  cistrans-viewer
```

## Run locally

```bash
pip install -r requirements.txt
streamlit run app.py
```

Optional: store downloaded data outside the project root.

```bash
DATA_ROOT=./data streamlit run app.py
```

Open the URL shown (for example, http://localhost:8501).

## Data

### Automatic first-run download

On startup, the app checks for required input data and auto-downloads missing artifacts from Zenodo:

1. `cis_trans_results_table.csv`
2. `gene_count_data.tar.gz` (automatically extracted to `gene_count_data/`)

The source URLs are stored in `data_sources.txt`:

- https://zenodo.org/records/19340139/files/cis_trans_results_table.csv
- https://zenodo.org/records/19340139/files/gene_count_data.tar.gz

This bootstrap step is idempotent:

- First run: downloads missing files and extracts `gene_count_data/`.
- Later runs: skips download/extract when data already exists.

By default, data is stored under the project root. You can override this with `DATA_ROOT`.

Expected layout under your data root:

```text
cistrans_paper_viewer/  (or DATA_ROOT if overridden)
  cis_trans_results_table.csv
  gene_count_data/
    subtype/
      <tissue>/
        <strain>_xgener_input_dataframe_FILTERED.csv
        <strain>_xgener_input_metadata_FILTERED.csv
```

### Cis/trans results

- **File:** `cis_trans_results_table.csv` in the data root.
- **Key columns:** `gene`, `strain`, `subtype`, `tissue`, `subtype_tis`, `Parlog2FC`, `Hyblog2FC`, `fdr_cis`, `fdr_trans`, `reg_assignment`, `cis_prop_reordered` / `cis_prop_reordered_fixed`, optional `colors`.

### Gene expression (allele/genotype) viewer

The **Gene expression** view uses a separate directory of per-tissue, per-strain count matrices and metadata:

**Base path:** `gene_count_data/subtype`

**Layout:**

```
gene_count_data/subtype/
  <tissue>/
    <strain>_xgener_input_dataframe_FILTERED.csv   # cells × genes (index column "Unnamed: 0" or first column)
    <strain>_xgener_input_metadata_FILTERED.csv   # cells × obs (index aligned with counts)
```

- **Count matrix:** cells as rows, genes as columns; index = cell IDs.
- **Metadata:** must include `subtype` (cell type) and `Allele` (values in `P1`, `P2`, `H1`, `H2`).
- **Plot:** For a chosen gene, tissue, and cell type, the app loads each strain’s AnnData (counts + metadata), filters to that subtype, and draws one panel per strain with boxplots by allele (P1, H1, H2, P2). Colors are strain-specific (P1/H1 = B6J, P2/H2 = founder for that strain). Only strains present in the selected tissue directory are shown.

The Zenodo `gene_count_data.tar.gz` files use `_FILTERED` suffixes by default in this app configuration. If you use custom files, update `BASE_PATH`, `DATA_SUFFIX`, and `META_SUFFIX` in `utils.py` to match.

## Views

- **Cell type view:** Select tissue(s) and cell type; scatter (cis_prop vs Parlog2FC) + stacked reg_assignment bar per strain.
- **Cell type proportions - strain:** One stacked bar plot of reg_assignment proportions across strains for a single (tissue, cell type).
- **Tissue-wide composition:** One stacked bar plot per strain for a tissue (all cell types, proportions).
- **Sankey view:** Flows of reg_assignment between two conditions (tissues, strains, or subtypes); overlapping genes only.
- **Gene expression:** Gene, tissue, and cell type; boxplots by allele (P1, H1, H2, P2) per strain, with optional download as PNG.

## Performance

- Cis/trans table and pre-aggregated bar/Sankey tables are cached (`@st.cache_data`).
- Gene count data: only the selected tissue is scanned; each strain’s AnnData is loaded on demand and cached by `(tissue, strain)`.
