# cis/trans regulatory inference explorer

Streamlit dashboard for exploring cis/trans regulatory inference results and gene expression (allele/genotype) across tissues, cell types, and strains.

## Data

### Cis/trans results

- **File:** `cis_trans_results_table.csv` in the project root.
- **Key columns:** `gene`, `strain`, `subtype`, `tissue`, `subtype_tis`, `Parlog2FC`, `Hyblog2FC`, `fdr_cis`, `fdr_trans`, `reg_assignment`, `cis_prop_reordered` / `cis_prop_reordered_fixed`, optional `colors`.

### Gene expression (allele/genotype) viewer

The **Gene expression** view uses a separate directory of per-tissue, per-strain count matrices and metadata:

**Base path:** `gene_count_data/subtype/filtered`

**Layout:**

```
gene_count_data/subtype/filtered/
  <tissue>/
    <strain>_xgener_input_dataframe_FILTERED.csv   # cells × genes (index column "Unnamed: 0" or first column)
    <strain>_xgener_input_metadata_FILTERED.csv   # cells × obs (index aligned with counts)
```

- **Count matrix:** cells as rows, genes as columns; index = cell IDs.
- **Metadata:** must include `subtype` (cell type) and `Allele` (values in `P1`, `P2`, `H1`, `H2`).
- **Plot:** For a chosen gene, tissue, and cell type, the app loads each strain’s AnnData (counts + metadata), filters to that subtype, and draws one panel per strain with boxplots by allele (P1, H1, H2, P2). Colors are strain-specific (P1/H1 = B6J, P2/H2 = founder for that strain). Only strains present in the selected tissue directory are shown.

If your files use different suffixes (e.g. `_xgener_input_dataframe.csv` without `_FILTERED`), set `BASE_PATH`, `DATA_SUFFIX`, and `META_SUFFIX` in `utils.py` to match.

## Run

```bash
pip install -r requirements.txt
streamlit run app.py
```

Open the URL shown (e.g. http://localhost:8501).

## Views

- **Cell type view:** Select tissue(s) and cell type; scatter (cis_prop vs Parlog2FC) + stacked reg_assignment bar per strain.
- **Cell type proportions - strain:** One stacked bar plot of reg_assignment proportions across strains for a single (tissue, cell type).
- **Tissue-wide composition:** One stacked bar plot per strain for a tissue (all cell types, proportions).
- **Sankey view:** Flows of reg_assignment between two conditions (tissues, strains, or subtypes); overlapping genes only.
- **Gene expression:** Gene, tissue, and cell type; boxplots by allele (P1, H1, H2, P2) per strain, with optional download as PNG.

## Performance

- Cis/trans table and pre-aggregated bar/Sankey tables are cached (`@st.cache_data`).
- Gene count data: only the selected tissue is scanned; each strain’s AnnData is loaded on demand and cached by `(tissue, strain)`.
