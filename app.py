"""
Cis/trans regulatory inference explorer + Gene expression (allele/genotype) viewer.
"""
import io
import os
import textwrap
import matplotlib.pyplot as plt
import logging
import numpy as np
import pandas as pd

import streamlit as st
from streamlit_plotly_events import plotly_events
from data_bootstrap import ensure_data_ready, get_data_root

from utils import (
    BASE_PATH,
    REG_ORDER,
    SANKEY_ORDER,
    REG_COLORS,
    build_sankey_condition_index,
    get_unique_sorted,
    load_adata,
    load_results_table,
    list_available_tissues_strains,
    make_celltype_strain_figure,
    make_tissue_composition_figure,
    plot_celltype_scatter_and_reg_proportions,
    plot_celltype_interactive_scatter_plotly,
    plot_celltype_proportions_matplotlib,
    plot_gene_across_strains,
    precompute_bar_aggregates,
    subset_for_celltype_view,
    build_sankey_figure,
    build_sankey_transition_table,
    display_sankey_gene_modal,
    get_condition_labels,
    get_gene_list_for_tissue,
    get_celltype_strain_data_table,
    get_tissue_composition_data_table,
    ALLELE_ORDER,
    STRAINS_DISPLAY_ORDER,
)

DATA_PATH = str(get_data_root() / "cis_trans_results_table.csv")
PAGES = [
    "Home",
    "Cell type view",
    "Cell type proportions - strain",
    "Tissue-wide composition",
    "Sankey view",
    "Gene expression",
]


# Configure logging
log_level = os.getenv("LOG_LEVEL", "WARNING").upper()
logging.basicConfig(
    level=getattr(logging, log_level, logging.WARNING),
    format="%(asctime)s - %(levelname)s - %(message)s",
    force=True,
)


def init_session_state() -> None:
    defaults = {"page": "Home", "last_page": "Home"}
    for key, val in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = val


def render_top_navigation() -> str:
    """Render top tab-like navigation and return selected page."""
    st.title("cis/trans regulatory inference explorer")
    st.caption("Interactive dashboard for cis/trans results and gene expression (allele/genotype) viewer.")

    cols = st.columns(len(PAGES))
    for idx, page_name in enumerate(PAGES):
        is_active = st.session_state.get("page") == page_name
        label = f"[{page_name}]" if is_active else page_name
        if cols[idx].button(label, key=f"top_nav_{idx}", use_container_width=True):
            st.session_state["page"] = page_name
            st.rerun()

    st.markdown("---")
    return st.session_state.get("page", "Home")


def page_landing():
    """Landing page with overview and navigation."""
    st.markdown("# Results Viewer")
    st.caption("Choose a view below to explore cis/trans regulatory variation.")
    
    st.markdown("---")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.markdown("### Cell Type View")
        if st.button("Open Cell Type View", use_container_width=True, key="btn_cell_view"):
            st.session_state["page"] = "Cell type view"
            st.rerun()
    
    with col2:
        st.markdown("### Cell Type Proportions")
        if st.button("Open Cell Type Proportions", use_container_width=True, key="btn_proportions"):
            st.session_state["page"] = "Cell type proportions - strain"
            st.rerun()
    
    with col3:
        st.markdown("### Tissue Composition")
        if st.button("Open Tissue Composition", use_container_width=True, key="btn_tissue"):
            st.session_state["page"] = "Tissue-wide composition"
            st.rerun()
    
    st.markdown("---")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("### Sankey Diagram")
        if st.button("Open Sankey View", use_container_width=True, key="btn_sankey"):
            st.session_state["page"] = "Sankey view"
            st.rerun()
    
    with col2:
        st.markdown("### Gene Expression")
        if st.button("Open Gene Expression", use_container_width=True, key="btn_gene"):
            st.session_state["page"] = "Gene expression"
            st.rerun()


def render_sidebar_header(page: str) -> None:
    st.sidebar.title("Controls")
    st.sidebar.caption(f"Current view: {page}")
    st.sidebar.markdown("---")


def page_celltype_view(df, counts_df):
    st.subheader("Cell type view")
    st.caption("Regulatory assignments by tissue, cell type, and strain.")
    tissue_values = get_unique_sorted(df, "tissue")
    if not tissue_values:
        st.warning("No `tissue` values found in the data.")
        return

    st.sidebar.markdown("**Cell type view**")
    default_tissue = "CortexHippocampus"
    default_tissues = [default_tissue] if default_tissue in tissue_values else []
    tissue_choice = st.sidebar.multiselect("Tissues", options=tissue_values, default=default_tissues, key="ctv_tissues")
    if not tissue_choice:
        st.info("Select at least one tissue to see cell types.")
        return
    subtype_sets = []
    for tis in tissue_choice:
        subtypes_in_tis = df[df["tissue"] == tis]["subtype"].dropna().unique().tolist()
        subtype_sets.append(set(subtypes_in_tis))
    if not subtype_sets:
        st.warning("No subtypes found for the selected tissues.")
        return
    common_subtypes = set.intersection(*subtype_sets)
    if not common_subtypes:
        st.warning("No cell types are shared across the selected tissues. Try selecting fewer tissues.")
        return
    subtype_values = sorted(map(str, common_subtypes))
    default_subtype = "glutamatergic_neuron"
    default_index = subtype_values.index(default_subtype) if default_subtype in subtype_values else 0
    subtype = st.sidebar.selectbox("Cell type (subtype)", options=subtype_values, index=default_index, key="ctv_subtype")
    strain_values = get_unique_sorted(df, "strain")
    strain_mode = st.sidebar.radio("Strain selection", options=["All strains", "Select strain"], key="ctv_strain_mode")
    if strain_mode == "Select strain":
        strain_values_sel = st.sidebar.multiselect("Strains", options=strain_values, default=strain_values, key="ctv_strains")
    else:
        strain_values_sel = strain_values

    for tissue in tissue_choice:
        st.markdown(f"### Tissue: `{tissue}`")
        for strain in strain_values_sel:
            st.markdown(f"#### Strain: `{strain}`")
            subset = subset_for_celltype_view(
                df=df, subtype=subtype, tissue=tissue, strain=strain,
                genes=None, reg_include=None, fdr_cis_max=None, fdr_trans_max=None,
            )
            if subset.empty:
                st.info(f"No data for strain `{strain}`.")
                continue
            
            # Display interactive scatter plot and proportions side-by-side
            col1, col2 = st.columns([1.2, 0.4])
            
            with col1:
                # Interactive scatter plot with hover info
                fig_scatter = plot_celltype_interactive_scatter_plotly(subset)
                st.plotly_chart(fig_scatter, use_container_width=True, config={"responsive": True})
            
            with col2:
                # Proportions bar chart
                fig_bar = plot_celltype_proportions_matplotlib(subset)
                st.pyplot(fig_bar, clear_figure=True)

    st.markdown("---")
    st.markdown("### Explore table")

    table_df = df[
        (df["tissue"].isin(tissue_choice))
        & (df["strain"].isin(strain_values_sel))
        & (df["subtype"] == subtype)
    ].copy()

    if table_df.empty:
        st.info("No rows available for the current selection.")
        return

    search_gene = st.text_input("Filter gene", value="", key="ctv_gene_search")
    reg_options = get_unique_sorted(table_df, "reg_assignment")
    selected_reg = st.multiselect(
        "Filter regulatory assignment",
        options=reg_options,
        default=reg_options,
        key="ctv_reg_filter",
    )

    if search_gene:
        table_df = table_df[table_df["gene"].astype(str).str.contains(search_gene, case=False, na=False)]
    if selected_reg:
        table_df = table_df[table_df["reg_assignment"].isin(selected_reg)]
    else:
        table_df = table_df.iloc[0:0]

    display_cols = [
        "gene",
        "strain",
        "tissue",
        "subtype",
        "reg_assignment",
        "Parlog2FC",
        "Hyblog2FC",
        "fdr_cis",
        "fdr_trans",
    ]
    available_cols = [c for c in display_cols if c in table_df.columns]
    display_table = table_df[available_cols].drop_duplicates().sort_values(["tissue", "strain", "gene"]).reset_index(drop=True)

    if display_table.empty:
        st.info("No rows match your table filters.")
        return

    st.caption(f"{len(display_table)} rows")
    st.dataframe(display_table, use_container_width=True, hide_index=True)
    st.download_button(
        "Download table (CSV)",
        data=display_table.to_csv(index=False),
        file_name=f"cell_type_view_{subtype}.csv",
        mime="text/csv",
        key="ctv_download_table",
    )


def page_celltype_strain(df, props_df):
    st.subheader("Cell type proportions across strains")
    st.caption("Compare regulatory-assignment proportions across strains.")
    tissue_values = get_unique_sorted(df, "tissue")
    if not tissue_values:
        st.warning("No `tissue` values found in the data.")
        return

    st.sidebar.markdown("**Cell type proportions**")
    default_tissue = "CortexHippocampus"
    default_index = tissue_values.index(default_tissue) if default_tissue in tissue_values else 0
    tissue = st.sidebar.selectbox("Tissue", options=tissue_values, index=default_index, key="cts_tissue")
    subtype_values = get_unique_sorted(df[df["tissue"] == tissue], "subtype")
    if not subtype_values:
        st.info("No cell types found for this tissue.")
        return
    default_subtype = "glutamatergic_neuron"
    default_sub_index = subtype_values.index(default_subtype) if default_subtype in subtype_values else 0
    subtype = st.sidebar.selectbox("Cell type (subtype)", options=subtype_values, index=default_sub_index, key="cts_subtype")
    
    # Sort option for plot
    sort_by = st.sidebar.selectbox(
        "Sort x-axis by",
        options=["n", "conserved_prop", "cis_prop", "trans_prop"],
        format_func=lambda x: {"n": "Total n", "conserved_prop": "Conserved proportion", "cis_prop": "Cis proportion", "trans_prop": "Trans proportion"}[x],
        key="celltype_strain_sort"
    )
    
    # Display matplotlib chart with all strains
    fig = make_celltype_strain_figure(props_df=props_df, tissue=tissue, subtype=subtype, sort_by=sort_by)
    if fig is None:
        st.info("No data for this tissue × cell type across strains.")
        return
    
    st.pyplot(fig, clear_figure=True)
    
    # Gene table viewer with optional strain and gene filtering
    st.markdown("---")
    st.markdown("### View genes")
    
    # Get all strains for this tissue/subtype
    strain_values = get_unique_sorted(df[(df["tissue"] == tissue) & (df["subtype"] == subtype)], "strain")
    
    # Create two-column layout for filters
    # Option to filter by strain or show all
    show_all_strains = st.sidebar.checkbox("Show all strains", value=True, key="celltype_show_all_strains")
    if show_all_strains:
        selected_strains = strain_values
    else:
        selected_strains = st.sidebar.multiselect("Select strains", options=strain_values, default=strain_values[:1] if strain_values else [], key="celltype_strain_filter")

    # Gene search
    search_gene = st.sidebar.text_input("Search gene by name (optional)", value="", key="celltype_gene_search", placeholder="Leave blank to show all")
    
    # Get and filter the gene table
    if selected_strains:
        gene_table = df[(df["tissue"] == tissue) & (df["subtype"] == subtype) & (df["strain"].isin(selected_strains))]
        
        # Apply gene search filter
        if search_gene:
            gene_table = gene_table[gene_table["gene"].str.contains(search_gene, case=False, na=False)]
        
        # Select and deduplicate columns
        display_cols = ["gene", "strain", "reg_assignment", "Parlog2FC", "Hyblog2FC", "fdr_cis", "fdr_trans"]
        available_cols = [col for col in display_cols if col in gene_table.columns]
        
        filtered_table = gene_table[available_cols].drop_duplicates(subset=["gene", "strain"]).sort_values(["strain", "gene"]).reset_index(drop=True)
        
        if not filtered_table.empty:
            st.markdown(f"**Genes** ({len(filtered_table)} results)")
            
            # Display table with proper formatting
            display_table = filtered_table.copy()
            for col in ["Parlog2FC", "Hyblog2FC", "fdr_cis", "fdr_trans"]:
                if col in display_table.columns:
                    display_table[col] = display_table[col].apply(lambda x: f"{x:.4f}" if pd.notna(x) else "NA")
            
            st.dataframe(display_table, use_container_width=True, hide_index=True)
            st.download_button(
                "Download table (CSV)",
                data=display_table.to_csv(index=False),
                file_name=f"cell_type_proportions_{tissue}_{subtype}.csv",
                mime="text/csv",
                key="celltype_strain_download",
            )
        else:
            st.info("No genes match your filters.")
    else:
        st.info("Please select at least one strain.")



def page_tissue_wide(df, props_df):
    st.subheader("Tissue-wide composition view")
    st.caption("Compare regulatory-assignment proportions across cell types.")
    tissue_values = get_unique_sorted(df, "tissue")
    if not tissue_values:
        st.warning("No `tissue` values found in the data.")
        return

    st.sidebar.markdown("**Tissue-wide composition**")
    default_tissue = "CortexHippocampus"
    default_index = tissue_values.index(default_tissue) if default_tissue in tissue_values else 0
    tissue = st.sidebar.selectbox("Tissue", options=tissue_values, index=default_index, key="tw_tissue")
    
    strain_values = get_unique_sorted(df[df["tissue"] == tissue], "strain")
    if not strain_values:
        st.info("No strains found for this tissue.")
        return
    
    # Display matplotlib chart for each strain
    for strain in strain_values:
        st.markdown(f"### Strain: `{strain}`")
        fig = make_tissue_composition_figure(props_df=props_df, tissue=tissue, strain=strain)
        if fig is None:
            st.info(f"No data for tissue `{tissue}` × strain `{strain}`.")
            continue
        st.pyplot(fig, clear_figure=True)
    
    # Cell type and strain selector to view gene data
    st.markdown("---")
    st.markdown("### View genes for a cell type")
    
    subtype_values = get_unique_sorted(df[df["tissue"] == tissue], "subtype")
    if not subtype_values:
        st.info("No cell types found for this tissue.")
        return
    selected_subtype = st.sidebar.selectbox("Select cell type", options=subtype_values, key="tissue_wide_subtype_select")
    selected_strain = st.sidebar.selectbox("Select strain", options=strain_values, key="tissue_wide_strain_select")
    
    # Get and display the gene table
    gene_table = get_tissue_composition_data_table(df, tissue=tissue, subtype=selected_subtype, strain=selected_strain)
    
    if not gene_table.empty:
        st.markdown(f"**Genes in {tissue} - {selected_subtype} - {selected_strain}** ({len(gene_table)} genes)")
        
        # Add search functionality
        search_gene = st.sidebar.text_input("Search gene by name", value="", key="tissue_wide_gene_search")
        
        if search_gene:
            filtered_table = gene_table[gene_table["gene"].str.contains(search_gene, case=False, na=False)]
        else:
            filtered_table = gene_table
        
        # Display table with proper formatting
        display_table = filtered_table.copy()
        for col in ["Parlog2FC", "Hyblog2FC", "fdr_cis", "fdr_trans"]:
            if col in display_table.columns:
                display_table[col] = display_table[col].apply(lambda x: f"{x:.4f}" if pd.notna(x) else "NA")
        
        st.dataframe(display_table, use_container_width=True, hide_index=True)
        st.download_button(
            "Download table (CSV)",
            data=display_table.to_csv(index=False),
            file_name=f"tissue_composition_{tissue}_{selected_subtype}_{selected_strain}.csv",
            mime="text/csv",
            key="tissue_wide_download",
        )
    else:
        st.info(f"No genes found for {selected_subtype} in {selected_strain}.")


def page_sankey(df, condition_index):
    """Sankey diagram view for exploring gene transitions across conditions."""
    st.subheader("Sankey view")
    st.caption("Visualize transitions between regulatory assignments.")
    
    # Get unique values for all selection controls
    tissue_values = get_unique_sorted(df, "tissue")
    strain_values = get_unique_sorted(df, "strain")
    subtype_values = get_unique_sorted(df, "subtype")
    if not (tissue_values and strain_values and subtype_values):
        st.warning("Not enough metadata columns to construct a Sankey view.")
        return
    
    # Set defaults
    default_strain = "CAST_EiJ"
    default_strain_idx = strain_values.index(default_strain) if default_strain in strain_values else 0
    default_tissue = "CortexHippocampus"
    default_tissue_idx = tissue_values.index(default_tissue) if default_tissue in tissue_values else 0
    default_right_tissue = "DiencephalonPituitary"
    default_right_idx = tissue_values.index(default_right_tissue) if default_right_tissue in tissue_values else default_tissue_idx
    default_subtype = "glutamatergic_neuron"
    default_sub_idx = subtype_values.index(default_subtype) if default_subtype in subtype_values else 0

    # ========== CONTROLS SECTION ==========
    st.sidebar.markdown("**Sankey configuration**")

    mode = st.sidebar.radio(
        "Select comparison mode",
        options=[
            "Same subtype across two tissues (within a strain)",
            "Same subtype across two strains (within a tissue)",
            "Two subtypes within same tissue and strain",
        ],
        key="sankey_mode",
    )
    
    # Determine left and right conditions based on selected mode
    if mode.startswith("Same subtype across two tissues"):
        strain = st.sidebar.selectbox("Strain", options=strain_values, index=default_strain_idx, key="sankey_strain")
        left_tissue = st.sidebar.selectbox("Left tissue", options=tissue_values, index=default_tissue_idx, key="sankey_left_tissue")
        right_tissue = st.sidebar.selectbox("Right tissue", options=tissue_values, index=default_right_idx, key="sankey_right_tissue")
        subtype = st.sidebar.selectbox("Subtype", options=subtype_values, index=default_sub_idx, key="sankey_subtype_a")
            
        left = (strain, left_tissue, subtype)
        right = (strain, right_tissue, subtype)
        title = f"{subtype}: {strain} across {left_tissue} → {right_tissue}"
            
    elif mode.startswith("Same subtype across two strains"):
        tissue = st.sidebar.selectbox("Tissue", options=tissue_values, index=default_tissue_idx, key="sankey_tissue_fixed")
        subtype = st.sidebar.selectbox("Subtype", options=subtype_values, index=default_sub_idx, key="sankey_subtype_b")
        left_strain = st.sidebar.selectbox("Left strain", options=strain_values, index=default_strain_idx, key="sankey_left_strain")
        right_strain = st.sidebar.selectbox("Right strain", options=strain_values, key="sankey_right_strain")
            
        left = (left_strain, tissue, subtype)
        right = (right_strain, tissue, subtype)
        title = f"{subtype}, {tissue}: {left_strain} → {right_strain}"
            
    else:  # Two subtypes within same tissue and strain
        strain = st.sidebar.selectbox("Strain", options=strain_values, index=default_strain_idx, key="sankey_strain2")
        tissue = st.sidebar.selectbox("Tissue", options=tissue_values, index=default_tissue_idx, key="sankey_tissue2")
        left_subtype = st.sidebar.selectbox("Left subtype", options=subtype_values, key="sankey_left_subtype")
        right_subtype = st.sidebar.selectbox("Right subtype", options=subtype_values, key="sankey_right_subtype")
            
        left = (strain, tissue, left_subtype)
        right = (strain, tissue, right_subtype)
        title = f"{tissue}, {strain}: {left_subtype} → {right_subtype}"

    # ========== BUILD SANKEY FIGURE ==========
    result = build_sankey_figure(condition_index, left, right, title)
    if result is None:
        st.warning("No data available for the selected conditions.")
        return
    
    fig, total_genes, shared_genes, gene_transitions = result
    left_label, right_label = get_condition_labels(left, right)
    
    # ========== DISPLAY SECTION ==========
    st.markdown("---")
    col1, col2 = st.columns([2, 1])
    with col1:
        st.markdown("## Sankey Diagram")
    with col2:
        st.caption(f"Total genes: {total_genes} | Shared: {shared_genes}")
    
    st.plotly_chart(fig, use_container_width=True, key="sankey_diagram")
    
    # ========== GENE TRANSITION SELECTOR ==========
    st.markdown("---")
    st.markdown("## Select a transition to view genes")
    
    # Build dropdown options
    transition_options = sorted([f"{l} → {r} ({len(genes)} genes)" 
                                for (l, r), genes in gene_transitions.items()])
    
    if not transition_options:
        st.error("No transitions found in this comparison.")
        return
    
    # Use a simple key for the dropdown to persist selection
    if "sankey_transition_dropdown" not in st.session_state:
        st.session_state["sankey_transition_dropdown"] = transition_options[0]
    
    selected_option = st.selectbox(
        "Transition",
        options=transition_options,
        key="sankey_transition_dropdown",
        index=0
    )
    
    # Parse the selected option to extract the transition
    if selected_option:
        # Format: "left_state → right_state (n genes)"
        parts = selected_option.split(" → ")
        left_state = parts[0].strip()
        right_part = parts[1].strip()
        right_state = right_part.split(" (")[0].strip()
        selected_transition = (left_state, right_state)
        
        # ========== DISPLAY GENE TABLE ==========
        if selected_transition in gene_transitions:
            selected_genes = gene_transitions[selected_transition]
            
            st.markdown("---")
            st.markdown(f"### Genes: {left_state} → {right_state}")
            st.caption(f"{len(selected_genes)} genes")
            
            # Build the table inline (simple and fast)
            data = []
            for gene in sorted(selected_genes):
                left_reg = condition_index.get(left, pd.Series()).get(gene, "N/A")
                right_reg = condition_index.get(right, pd.Series()).get(gene, "N/A")
                data.append({
                    "Gene": gene,
                    left_label: left_reg,
                    right_label: right_reg,
                })
            
            gene_table = pd.DataFrame(data)
            
            # Filter option
            search_query = st.text_input(
                "Filter genes by name",
                value="",
                placeholder="Type to search...",
                key=f"sankey_search_{left_state}_{right_state}"
            )
            
            # Apply filter
            if search_query:
                display_df = gene_table[
                    gene_table["Gene"].str.contains(search_query, case=False, na=False)
                ]
            else:
                display_df = gene_table
            
            # Display table
            st.dataframe(display_df, use_container_width=True, hide_index=True)
            
            # Download button
            csv = display_df.to_csv(index=False)
            st.download_button(
                label="Download as CSV",
                data=csv,
                file_name=f"genes_{left_state}_to_{right_state}.csv",
                mime="text/csv",
                key=f"sankey_download_{left_state}_{right_state}"
            )


def page_gene_expression():
    st.subheader("Gene expression (allele/genotype) viewer")
    st.caption("Plot gene expression across strains with allele-specific detail.")
    logging.debug("Loading available tissues and strains.")
    st.sidebar.markdown("**Gene expression**")
    tissues, by_tissue = list_available_tissues_strains()
    if not tissues:
        st.warning(f"No tissue directories found under `{BASE_PATH}`. Check that gene count data is present.")
        logging.error("No tissue directories found under the base path.")
        return

    default_tissue = "CortexHippocampus"
    default_tissue_idx = tissues.index(default_tissue) if default_tissue in tissues else 0
    tissue = st.sidebar.selectbox("Tissue", options=tissues, index=default_tissue_idx, key="ge_tissue")
    logging.debug(f"Selected tissue: {tissue}")

    strains = by_tissue.get(tissue, [])
    if not strains:
        st.info(f"No strain data found for tissue `{tissue}`.")
        logging.warning(f"No strain data found for tissue {tissue}.")
        return

    # Get subtypes from first available strain for this tissue
    adata_sample = None
    strain = None  # Initialize strain variable
    for s in strains:
        try:
            adata_sample = load_adata(tissue, s)
            if adata_sample is not None:
                strain = s  # Assign the strain variable when data is successfully loaded
                break
        except Exception as e:
            logging.error(f"Error loading AnnData for strain {s}: {e}")

    if adata_sample is None:
        st.info("Could not load sample metadata to list cell types.")
        logging.error("Could not load sample metadata to list cell types.")
        return

    subtype_options = sorted(adata_sample.obs["subtype"].astype(str).unique().tolist())
    default_subtype = "glutamatergic_neuron"
    default_sub_idx = subtype_options.index(default_subtype) if default_subtype in subtype_options else 0
    subtype = st.sidebar.selectbox("Cell type (subtype)", options=subtype_options, index=default_sub_idx, key="ge_subtype")
    logging.debug(f"Selected subtype: {subtype}")

    all_genes = get_gene_list_for_tissue(tissue)
    logging.debug(f"Loaded {len(all_genes)} genes for tissue {tissue}.")
    
    # Initialize session state for gene selection and plotting
    if "ge_selected_gene" not in st.session_state:
        st.session_state["ge_selected_gene"] = "Dcc"
    if "ge_should_plot" not in st.session_state:
        st.session_state["ge_should_plot"] = True  # Plot on first load
    
    # ========== GENE SELECTION SECTION ==========
    st.sidebar.markdown("### Gene Selection")
    
    # Autocomplete dropdown for gene selection
    try:
        default_index = all_genes.index(st.session_state["ge_selected_gene"])
    except ValueError:
        default_index = 0
    
    gene = st.sidebar.selectbox(
        "Select gene (type to autocomplete):",
        options=all_genes,
        index=default_index,
        key="ge_gene_select",
    )
    st.session_state["ge_selected_gene"] = gene
    
    logging.debug(f"Selected gene: {gene}")

    st.sidebar.markdown("---")
    st.sidebar.markdown("**Plot Options**")
    layer = st.sidebar.radio("Expression layer", options=["CPM", "raw_counts"], index=0, key="ge_layer")
    show_points = st.sidebar.checkbox("Show points", value=True, key="ge_show_points")
    jitter = st.sidebar.slider("Jitter", min_value=0.0, max_value=0.5, value=0.15, step=0.05, key="ge_jitter")
    y_scale = "linear"
    y_min = 0.0
    downsample_per_allele = None

    # ========== PLOT BUTTON SECTION ==========
    st.sidebar.markdown("---")
    plot_clicked = st.sidebar.button("Plot", key="ge_plot_button", use_container_width=True)

    if st.sidebar.button("Reset", key="ge_reset_button", use_container_width=True):
        st.session_state["ge_selected_gene"] = "Dcc"
        st.session_state["ge_gene_search"] = ""
        st.session_state["ge_should_plot"] = True
        st.rerun()
    
    # Plot on button click or on first load with Dcc
    should_plot = plot_clicked or st.session_state.get("ge_should_plot", False)
    
    if should_plot and gene and gene != "":
        st.session_state["ge_should_plot"] = False  # Only auto-plot once
        try:
            with st.spinner(f"Loading plot for {gene}..."):
                fig, warnings = plot_gene_across_strains(
                    tissue=tissue,
                    subtype=subtype,
                    gene=gene,
                    layer=layer,
                    show_points=show_points,
                    jitter=jitter,
                )
            for w in warnings:
                st.warning(w)
                logging.warning(w)
            if fig is None:
                st.error("Could not generate plot. Check that the gene exists in the data and that the cell type is present.")
                logging.error("Plot generation failed. Check gene and cell type data.")
            else:
                st.pyplot(fig, clear_figure=True)

                buf = io.BytesIO()
                fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
                buf.seek(0)
                st.download_button("📥 Download figure (PNG)", data=buf, file_name=f"gene_expression_{tissue}_{subtype}_{gene}.png", mime="image/png", key="ge_download")
                plt.close(fig)
        except Exception as e:
            st.error(f"An error occurred while generating the plot: {e}")
            logging.error(f"An error occurred while generating the plot: {e}")

    st.markdown("---")
    st.markdown("### Explore table")

    results_df = load_results_table(DATA_PATH)
    strain_order = [s for s in STRAINS_DISPLAY_ORDER if s in strains]
    rows = []

    for s in strain_order:
        adata_s = load_adata(tissue, s)
        if adata_s is None or gene not in adata_s.var_names:
            continue
        if "subtype" not in adata_s.obs.columns or "Allele" not in adata_s.obs.columns:
            continue

        sub = adata_s.obs.loc[(adata_s.obs["subtype"] == subtype) & (adata_s.obs["Allele"].isin(ALLELE_ORDER))].copy()
        if sub.empty:
            continue

        gene_idx = adata_s.var_names.get_loc(gene)
        layer_data = adata_s.layers["raw_counts"] if layer == "raw_counts" else adata_s.layers[layer]
        expr = np.ravel(layer_data[np.asarray(adata_s.obs_names.get_indexer(sub.index)), gene_idx]).astype(float)
        sub["_expr"] = expr

        reg_row = results_df[
            (results_df["gene"] == gene)
            & (results_df["tissue"] == tissue)
            & (results_df["subtype"] == subtype)
            & (results_df["strain"] == s)
        ]
        reg_assignment = reg_row["reg_assignment"].iloc[0] if not reg_row.empty else "not_detected"

        for sample_id, row in sub.iterrows():
            rows.append(
                {
                    "sample_id": str(sample_id),
                    "strain": s,
                    "allele": str(row["Allele"]),
                    "gene": gene,
                    "subtype": subtype,
                    "layer": layer,
                    "expression": float(row["_expr"]),
                    "reg_assignment": reg_assignment,
                }
            )

    if not rows:
        st.info("No expression rows available for the current selection.")
    else:
        table = pd.DataFrame(rows)
        table["strain"] = pd.Categorical(table["strain"], categories=strain_order, ordered=True)
        table["allele"] = pd.Categorical(table["allele"], categories=ALLELE_ORDER, ordered=True)
        table = table.sort_values(["strain", "allele"]).reset_index(drop=True)

        st.caption(f"{len(table)} rows")
        st.dataframe(table, use_container_width=True, hide_index=True)
        st.download_button(
            "Download table (CSV)",
            data=table.to_csv(index=False),
            file_name=f"gene_expression_{tissue}_{subtype}_{gene}_{layer}.csv",
            mime="text/csv",
            key="ge_download_table",
        )

    # Add debug logs to verify data loading and plotting
    logging.debug(f"Available genes for tissue {tissue}: {len(all_genes)} total")
    logging.debug(f"Loaded AnnData for tissue {tissue}, strain {strain}: {adata_sample}")
    logging.debug(f"AnnData layers: {list(adata_sample.layers.keys()) if adata_sample else 'None'}")
    logging.debug(f"AnnData obs columns: {list(adata_sample.obs.columns) if adata_sample else 'None'}")
    logging.debug(f"Plotting data for gene {gene}, tissue {tissue}, subtype {subtype}")


def main():
    st.set_page_config(page_title="cis/trans regulatory inference explorer", layout="wide")
    init_session_state()

    with st.spinner("Checking required input data. This may take a moment..."):
        bootstrap_result = ensure_data_ready()
    if bootstrap_result.error:
        st.error(textwrap.dedent(f"""
            Failed to prepare required input data.
            Reason: {bootstrap_result.error}

            Check network access and verify data source URLs in `data_sources.txt`.
            You can also set `DATA_ROOT` to a writable directory and retry.
        """))
        return
    if bootstrap_result.downloaded:
        downloaded_items = ", ".join(bootstrap_result.downloaded)
        st.info(f"First-run data setup completed. Downloaded: {downloaded_items}")

    try:
        df = load_results_table(DATA_PATH)
    except FileNotFoundError:
        st.error(textwrap.dedent("""
            Could not find `cis_trans_results_table.csv` after startup data checks.
            Verify that `DATA_ROOT` points to a readable data directory and refresh.
        """))
        return
    if df.empty:
        st.error("Loaded CSV is empty.")
        return

    counts_df, props_df = precompute_bar_aggregates(df)
    condition_index = build_sankey_condition_index(df)
    page = render_top_navigation()

    # Re-entering Gene expression should auto-plot the current/default selection.
    prev_page = st.session_state.get("last_page", "Home")
    if page == "Gene expression" and prev_page != "Gene expression":
        st.session_state.setdefault("ge_selected_gene", "Dcc")
        st.session_state["ge_should_plot"] = True

    st.session_state["last_page"] = page

    render_sidebar_header(page)
    st.session_state["common_filters"] = {}

    if page == "Home":
        page_landing()
    elif page == "Cell type view":
        page_celltype_view(df, counts_df)
    elif page == "Cell type proportions - strain":
        page_celltype_strain(df, props_df)
    elif page == "Tissue-wide composition":
        page_tissue_wide(df, props_df)
    elif page == "Sankey view":
        page_sankey(df, condition_index)
    elif page == "Gene expression":
        page_gene_expression()


if __name__ == "__main__":
    main()
