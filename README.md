# cis/trans regulatory inference explorer

Streamlit dashboard for exploring cis/trans regulatory inference results and gene expression (allele/genotype) across tissues, cell types, and strains.

## Data

### Download required input data

This app requires two input artifacts from Zenodo:

1. `cis_trans_results_table.csv`
   - https://zenodo.org/records/19340139/files/cis_trans_results_table.csv
2. `gene_count_data.tar.gz`
   - https://zenodo.org/records/19340139/files/gene_count_data.tar.gz

From the project root, download both files:

```bash
curl -L -o cis_trans_results_table.csv \
  https://zenodo.org/records/19340139/files/cis_trans_results_table.csv

curl -L -o gene_count_data.tar.gz \
  https://zenodo.org/records/19340139/files/gene_count_data.tar.gz
```

Extract the gene count archive in the project root:

```bash
tar -xzf gene_count_data.tar.gz
```

After extraction, your top-level layout should include:

```text
cistrans_paper_viewer/
  app.py
  utils.py
  cis_trans_results_table.csv
  gene_count_data/
    subtype/
      <tissue>/
        <strain>_xgener_input_dataframe_FILTERED.csv
        <strain>_xgener_input_metadata_FILTERED.csv
```

### Cis/trans results

- **File:** `cis_trans_results_table.csv` in the project root.
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

The Zenodo `gene_count_data.tar.gz` files use `_FILTERED` suffixes. If your local `utils.py` uses different values for `BASE_PATH`, `DATA_SUFFIX`, or `META_SUFFIX`, update those constants to match your extracted file names.

## Run

```bash
pip install -r requirements.txt
streamlit run app.py
```

Open the URL shown (e.g. http://localhost:8501).

## Run with Docker (GHCR)

Pull image:

```bash
docker pull ghcr.io/mortazavilab/mousaic:latest
```


Run container (mount local data into `/app`):

```bash
docker run --rm -p 8501:8501 \
  -v "$PWD/cis_trans_results_table.csv:/app/cis_trans_results_table.csv" \
  -v "$PWD/gene_count_data:/app/gene_count_data" \
  ghcr.io/rlweber23/mousaic:latest
```

If port 8501 is already in use, map a different host port:

```bash
docker run --rm -p 8502:8501 \
  -v "$PWD/cis_trans_results_table.csv:/app/cis_trans_results_table.csv" \
  -v "$PWD/gene_count_data:/app/gene_count_data" \
  ghcr.io/rlweber23/mousaic:latest
```

Then open http://localhost:8501 (or http://localhost:8502 if using the alternate port mapping).

Notes:

- Data files are not baked into the image; mount `cis_trans_results_table.csv` and `gene_count_data/` at runtime.
- The app expects data at `/app/cis_trans_results_table.csv` and `/app/gene_count_data/`.

### Optional: Build Locally

```bash
docker build -t cistrans-viewer .
```

```bash
docker run --rm -p 8501:8501 \
  -v "$PWD/cis_trans_results_table.csv:/app/cis_trans_results_table.csv" \
  -v "$PWD/gene_count_data:/app/gene_count_data" \
  cistrans-viewer
```

## Views

- **Cell type view:** Select tissue(s) and cell type; scatter (cis_prop vs Parlog2FC) + stacked reg_assignment bar per strain.
- **Cell type proportions - strain:** One stacked bar plot of reg_assignment proportions across strains for a single (tissue, cell type).
- **Tissue-wide composition:** One stacked bar plot per strain for a tissue (all cell types, proportions).
- **Sankey view:** Flows of reg_assignment between two conditions (tissues, strains, or subtypes); overlapping genes only.
- **Gene expression:** Gene, tissue, and cell type; boxplots by allele (P1, H1, H2, P2) per strain, with optional download as PNG.

## Performance

- Cis/trans table and pre-aggregated bar/Sankey tables are cached (`@st.cache_data`).
- Gene count data: only the selected tissue is scanned; each strain’s AnnData is loaded on demand and cached by `(tissue, strain)`.
