"""
Utilities for cis/trans dashboard and gene expression (allele/genotype) viewer.
"""
import os
from typing import Dict, List, Literal, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from data_bootstrap import get_data_root
from matplotlib.patches import Rectangle

# ---------------------------------------------------------------------------
# Cis/trans constants
# ---------------------------------------------------------------------------
REG_ORDER: List[str] = ["conserved", "cis", "trans", "cisxtrans", "cis+trans"]
SANKEY_ORDER: List[str] = ["conserved", "cis", "cisxtrans", "trans"]

REG_COLORS: Dict[str, str] = {
    "conserved": "#D3D3D3",
    "cis": "#FF4500",
    "trans": "#4169E1",
    "cis+trans": "#87CEEB",
    "cisxtrans": "#228B22",
    "not_detected": "#999999",
}

TISSUE_COMPOSITION_COLORS: Dict[str, str] = {
    "conserved": "#D3D3D3",
    "cis": "#FF4500",
    "trans": "#4169E1",
    "cis+trans": "#87CEEB",
    "cisxtrans": "#228B22",
}

# ---------------------------------------------------------------------------
# Gene expression viewer constants
# ---------------------------------------------------------------------------
DATA_ROOT = str(get_data_root())
BASE_PATH = os.path.join(DATA_ROOT, "gene_count_data", "subtype")
DATA_SUFFIX = "_xgener_input_dataframe_FILTERED.csv"
META_SUFFIX = "_xgener_input_metadata_FILTERED.csv"

STRAINS_DISPLAY_ORDER = [
    "A_J",
    "NOD_ShiLtJ",
    "129S1_SvImJ",
    "NZO_HlLtJ",
    "WSB_EiJ",
    "PWK_PhJ",
    "CAST_EiJ",
]

FOUNDER_SHORTNAME: Dict[str, str] = {
    "129S1_SvImJ": "129S1J",
    "A_J": "AJ",
    "CAST_EiJ": "CASTJ",
    "NOD_ShiLtJ": "NODJ",
    "NZO_HlLtJ": "NZOJ",
    "PWK_PhJ": "PWKJ",
    "WSB_EiJ": "WSBJ",
}

GENO_DICT: Dict[str, str] = {
    "129S1J": "#DA9CC1",
    "B6J": "#C0BFBF",
    "AJ": "#F4C245",
    "CASTJ": "#55AF5B",
    "NODJ": "#4F6EAF",
    "NZOJ": "#52A5DB",
    "PWKJ": "#D83026",
    "WSBJ": "#683C91",
}

ALLELE_ORDER = ["P1", "H1", "H2", "P2"]


# ---------------------------------------------------------------------------
# Cis/trans data loading and caching
# ---------------------------------------------------------------------------
@st.cache_data(show_spinner=False)
def load_results_table(csv_path: str) -> pd.DataFrame:
    dtypes = {
        "gene": "category",
        "strain": "category",
        "subtype": "category",
        "tissue": "category",
        "subtype_tis": "category",
        "reg_assignment": "category",
        "colors": "category",
    }
    df = pd.read_csv(csv_path, dtype=dtypes, low_memory=False)
    for col in [
        "Parlog2FC",
        "Hyblog2FC",
        "fdr_cis",
        "fdr_trans",
        "cis_prop_reordered_fixed",
        "cis_prop_reordered",
    ]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


@st.cache_data(show_spinner=False)
def precompute_bar_aggregates(
    df: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if "reg_assignment" not in df.columns:
        empty = pd.DataFrame(
            columns=["strain", "tissue", "subtype", "reg_assignment", "n"]
        )
        return empty, empty.assign(proportion=0.0)
    group_cols = ["strain", "tissue", "subtype", "reg_assignment"]
    counts = (
        df.groupby(group_cols, observed=True)
        .size()
        .reset_index(name="n")
        .astype({"n": "int64"})
    )
    totals = (
        counts.groupby(["strain", "tissue", "subtype"], observed=True)["n"]
        .sum()
        .rename("total_n")
        .reset_index()
    )
    props = counts.merge(totals, on=["strain", "tissue", "subtype"], how="left")
    props["proportion"] = props["n"] / props["total_n"].replace(0, np.nan)
    return counts, props


@st.cache_data(show_spinner=False)
def build_sankey_condition_index(
    df: pd.DataFrame,
) -> Dict[Tuple[str, str, str], pd.Series]:
    required_cols = {"gene", "strain", "tissue", "subtype", "reg_assignment"}
    if not required_cols.issubset(df.columns):
        return {}
    index: Dict[Tuple[str, str, str], pd.Series] = {}
    grouped = df.groupby(["strain", "tissue", "subtype"], observed=True)
    for key, sub in grouped:
        ser = (
            sub[["gene", "reg_assignment"]]
            .dropna(subset=["gene"])
            .drop_duplicates(subset=["gene"])
            .set_index("gene")["reg_assignment"]
        )
        index[(str(key[0]), str(key[1]), str(key[2]))] = ser
    return index


def _apply_basic_filters(
    df: pd.DataFrame,
    tissue: Optional[str] = None,
    strain: Optional[str] = None,
    subtype: Optional[str] = None,
    genes: Optional[List[str]] = None,
    reg_include: Optional[List[str]] = None,
    fdr_cis_max: Optional[float] = None,
    fdr_trans_max: Optional[float] = None,
) -> pd.DataFrame:
    mask = pd.Series(True, index=df.index)
    if tissue and "tissue" in df.columns:
        mask &= df["tissue"] == tissue
    if strain and "strain" in df.columns:
        mask &= df["strain"] == strain
    if subtype and "subtype" in df.columns:
        mask &= df["subtype"] == subtype
    if genes and "gene" in df.columns:
        mask &= df["gene"].isin(genes)
    if reg_include is not None and "reg_assignment" in df.columns:
        mask &= df["reg_assignment"].isin(reg_include)
    if fdr_cis_max is not None and "fdr_cis" in df.columns:
        mask &= df["fdr_cis"] <= fdr_cis_max
    if fdr_trans_max is not None and "fdr_trans" in df.columns:
        mask &= df["fdr_trans"] <= fdr_trans_max
    return df[mask]


def plot_celltype_scatter_and_reg_proportions(
    df: pd.DataFrame,
    *,
    x_col: str = "cis_prop_reordered",
    y_col: str = "Parlog2FC",
    color_col: str = "colors",
    reg_col: str = "reg_assignment",
) -> plt.Figure:
    d = df.copy()
    if d.empty:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, "No data for this selection", ha="center", va="center", fontsize=12)
        ax.axis("off")
        return fig

    def infer_unique(col: str) -> str:
        if col not in d.columns:
            raise ValueError(f"Column '{col}' not found in dataframe.")
        vals = d[col].dropna().unique()
        if len(vals) != 1:
            raise ValueError(f"Expected exactly one unique value in '{col}', found {len(vals)}.")
        return str(vals[0])

    strain = infer_unique("strain")
    tissue = infer_unique("tissue")
    subtype = infer_unique("subtype")
    title = f"{tissue} {subtype} - {strain}"
    bar_xtick = subtype

    total_count = d[reg_col].notna().sum()
    plot_df = (
        d[reg_col]
        .value_counts(normalize=True)
        .reindex(REG_ORDER)
        .fillna(0)
        .to_frame("sample")
    )

    fig, (ax_sc, ax_bar) = plt.subplots(
        1, 2, figsize=(6.6, 5.8), gridspec_kw={"width_ratios": (3.3, 1.1)}
    )

    band = (-0.5, 0.5)
    vlines = (-1, -0.5, 0, 0.5, 1)
    vline_colors = ("forestgreen", "royalblue", "lightblue", "orangered", "forestgreen")

    ax_sc.axhspan(band[0], band[1], color="lightgray", alpha=0.4, zorder=0)
    ax_sc.grid(True, which="major", color="#d0d0d0", linewidth=0.8)
    ax_sc.grid(True, which="minor", color="#eeeeee", linewidth=0.5)
    ax_sc.minorticks_on()

    if color_col in d.columns and d[color_col].notna().any():
        point_colors = d[color_col]
    elif reg_col in d.columns:
        point_colors = d[reg_col].map(REG_COLORS)
    else:
        point_colors = "#000000"

    ax_sc.scatter(
        d[x_col], d[y_col], c=point_colors, s=20, alpha=0.9, linewidths=0, zorder=1
    )
    ax_sc.axhline(0, color="black", linewidth=1.5, zorder=4)
    for xv, col in zip(vlines, vline_colors):
        ax_sc.axvline(xv, color=col, linewidth=3, alpha=0.9, zorder=2)
    ax_sc.set_xlim(-1.1, 1.1)
    ax_sc.set_ylim(-11, 11)
    ax_sc.set_xlabel("Proportion cis", fontsize=15)
    ax_sc.set_ylabel(r"$R_P$", fontsize=16)
    ax_sc.tick_params(axis="both", labelsize=12)
    ax_sc.set_xticks([-1, -0.5, 0, 0.5, 1])
    ax_sc.set_xticklabels([0.5, 0, 0.5, 1, 0.5], fontsize=12)
    ax_sc.set_axisbelow(True)
    for spine in ["top", "right"]:
        ax_sc.spines[spine].set_visible(False)

    reg_order_for_bar = REG_ORDER
    bar_colors = [REG_COLORS[c] for c in reg_order_for_bar]
    plot_df.T.plot(kind="bar", stacked=True, width=0.75, color=bar_colors, ax=ax_bar, legend=False)
    stack_label_min = 0.03
    cumulative = 0.0
    for cat in reg_order_for_bar:
        value = float(plot_df.loc[cat, "sample"])
        if value > stack_label_min:
            ax_bar.text(0, cumulative + value / 2, f"{value:.2f}", ha="center", va="center", fontsize=11)
        cumulative += value
    ax_bar.text(0, 1.05, f"n={total_count}", ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax_bar.set_ylim(0, 1.05)
    ax_bar.set_xlabel("")
    ax_bar.set_xticklabels([bar_xtick], rotation=0, fontsize=12)
    ax_bar.tick_params(axis="y", labelsize=12)
    for spine in ["top", "right"]:
        ax_bar.spines[spine].set_visible(False)
    fig.suptitle(title, fontsize=21, y=0.99)
    handles = [Rectangle((0, 0), 1, 1, color=REG_COLORS[c]) for c in reg_order_for_bar]
    fig.legend(
        handles, reg_order_for_bar, title="Reg Assignment", title_fontsize=11, fontsize=11,
        ncol=len(reg_order_for_bar), loc="upper center", bbox_to_anchor=(0.5, -0.12), frameon=False,
    )
    plt.tight_layout(rect=[0, 0.05, 1, 0.94])
    return fig


def plot_celltype_interactive_scatter_plotly(
    df: pd.DataFrame,
    *,
    x_col: str = "cis_prop_reordered",
    y_col: str = "Parlog2FC",
    color_col: str = "colors",
    reg_col: str = "reg_assignment",
    gene_col: str = "gene",
) -> go.Figure:
    """Create interactive Plotly scatter plot with hover info (gene, reg_assignment, prop_cis)."""
    d = df.copy()
    
    if d.empty:
        fig = go.Figure()
        fig.add_annotation(
            text="No data for this selection",
            showarrow=False,
            font=dict(size=14)
        )
        fig.update_layout(
            plot_bgcolor="white",
            paper_bgcolor="white"
        )
        return fig
    
    # Ensure gene names are strings
    if gene_col in d.columns:
        # Convert categorical to string (NaNs become "nan"), then replace with "unknown"
        d[gene_col] = d[gene_col].astype(str).replace("nan", "unknown")
    else:
        d[gene_col] = "unknown"
    
    # Get color mapping
    if color_col in d.columns and d[color_col].notna().any():
        point_colors = d[color_col]
    elif reg_col in d.columns:
        point_colors = d[reg_col].map(REG_COLORS)
    else:
        point_colors = "#000000"
    
    # Create hover text with gene, reg_assignment, prop_cis, and log2FC values
    hover_text = []
    for idx, row in d.iterrows():
        gene = row.get(gene_col, "unknown")
        reg = row.get(reg_col, "unknown")
        cis_prop = row.get(x_col, np.nan)
        parlog2fc = row.get(y_col, np.nan)
        hyblog2fc = row.get("Hyblog2FC", np.nan)
        
        hover_info = f"<b>{gene}</b><br>"
        hover_info += f"Reg Assignment: {reg}<br>"
        hover_info += f"Prop Cis: {cis_prop:.3f}<br>"
        hover_info += f"Parlog2FC (R_P): {parlog2fc:.3f}<br>"
        if not np.isnan(hyblog2fc):
            hover_info += f"Hyblog2FC: {hyblog2fc:.3f}"
        hover_text.append(hover_info)
    
    d["hover"] = hover_text
    
    # Create Plotly scatter
    fig = go.Figure()
    
    # Add reference bands and lines (visual guides)
    # Band from -0.5 to 0.5
    fig.add_shape(
        type="rect",
        x0=-1.1, x1=1.1, y0=-0.5, y1=0.5,
        fillcolor="lightgray", opacity=0.2, line_width=0, layer="below"
    )
    
    # Vertical reference lines
    vlines = [(-1, "forestgreen"), (-0.5, "royalblue"), (0, "lightblue"), (0.5, "orangered"), (1, "forestgreen")]
    for xline, color in vlines:
        fig.add_vline(x=xline, line_color=color, line_width=3, opacity=0.9)
    
    # Horizontal reference line at y=0
    fig.add_hline(y=0, line_color="black", line_width=2)
    
    # Add scatter points
    fig.add_trace(go.Scatter(
        x=d[x_col],
        y=d[y_col],
        mode="markers",
        marker=dict(
            size=8,
            color=point_colors.values if isinstance(point_colors, pd.Series) else point_colors,
            opacity=0.9,
            line_width=0
        ),
        text=d["hover"],
        hovertemplate="%{text}<extra></extra>",
        name=""
    ))
    
    # Update layout
    fig.update_layout(
        title="",
        xaxis_title="Proportion cis",
        yaxis_title="R_P",
        xaxis=dict(
            range=[-1.1, 1.1],
            tickvals=[-1, -0.5, 0, 0.5, 1],
            ticktext=["0.5", "0", "0.5", "1", "0.5"],
            showgrid=True,
            gridwidth=1,
            gridcolor="#d0d0d0",
            tickcolor="black",
            tickfont=dict(color="black"),
            titlefont=dict(color="black")
        ),
        yaxis=dict(
            range=[-11, 11],
            showgrid=True,
            gridwidth=1,
            gridcolor="#d0d0d0",
            tickcolor="black",
            tickfont=dict(color="black"),
            titlefont=dict(color="black")
        ),
        hovermode="closest",
        showlegend=False,
        template="plotly_white",
        plot_bgcolor="white",
        paper_bgcolor="white",
        height=600,
        width=700,
        margin=dict(l=80, r=50, t=50, b=80)
    )
    
    return fig


def plot_celltype_proportions_matplotlib(
    df: pd.DataFrame,
    *,
    reg_col: str = "reg_assignment",
) -> plt.Figure:
    """Create matplotlib bar chart showing regulatory assignment proportions."""
    d = df.copy()
    
    if d.empty:
        fig, ax = plt.subplots(figsize=(1.5, 6), facecolor='white')
        ax.set_facecolor('white')
        ax.text(0.5, 0.5, "No data", ha="center", va="center", fontsize=12)
        ax.axis("off")
        return fig
    
    def infer_unique(col: str) -> str:
        if col not in d.columns:
            return "unknown"
        vals = d[col].dropna().unique()
        return str(vals[0]) if len(vals) > 0 else "unknown"
    
    subtype = infer_unique("subtype")
    total_count = d[reg_col].notna().sum()
    
    plot_df = (
        d[reg_col]
        .value_counts(normalize=True)
        .reindex(REG_ORDER)
        .fillna(0)
        .to_frame("sample")
    )
    
    fig, ax = plt.subplots(figsize=(2.2, 6), facecolor='white')
    ax.set_facecolor('white')
    
    reg_order_for_bar = REG_ORDER
    bar_colors = [REG_COLORS[c] for c in reg_order_for_bar]
    plot_df.T.plot(kind="bar", stacked=True, width=0.75, color=bar_colors, ax=ax, legend=False)
    
    stack_label_min = 0.03
    cumulative = 0.0
    for cat in reg_order_for_bar:
        value = float(plot_df.loc[cat, "sample"])
        if value > stack_label_min:
            ax.text(0, cumulative + value / 2, f"{value:.2f}", ha="center", va="center", fontsize=11)
        cumulative += value
    
    ax.text(0, 1.05, f"n={total_count}", ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("")
    ax.set_xticklabels([subtype], rotation=0, fontsize=11)
    ax.tick_params(axis="y", labelsize=11)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    
    # Create legend to the left of the bar plot
    handles = [Rectangle((0, 0), 1, 1, color=REG_COLORS[c]) for c in reg_order_for_bar]
    fig.legend(
        handles, reg_order_for_bar, title="Reg Assignment", title_fontsize=11, fontsize=10,
        ncol=1, loc="center left", bbox_to_anchor=(-0.6, 0.5), frameon=False,
    )
    plt.tight_layout()
    return fig


def make_celltype_strain_interactive_plotly(
    props_df: pd.DataFrame,
    tissue: str,
    subtype: str,
    sort_by: Literal["n", "conserved_prop", "cis_prop", "trans_prop"] = "n",
) -> Optional[go.Figure]:
    """Create interactive Plotly bar chart for cell type across strains with fixed strain order."""
    subset = props_df[
        (props_df["tissue"] == tissue) & (props_df["subtype"] == subtype)
    ]
    if subset.empty:
        return None
    
    # Prepare data
    strain_counts = subset.groupby("strain")["total_n"].first().rename("n").sort_index()
    df_props = (
        subset.pivot_table(
            index="strain", columns="reg_assignment", values="proportion", fill_value=0.0
        )
        .reindex(columns=REG_ORDER, fill_value=0.0)
    )
    df_props = df_props.merge(strain_counts.to_frame(), left_index=True, right_index=True)
    
    # Reorder strains to STRAINS_DISPLAY_ORDER
    df_props = df_props.reindex(STRAINS_DISPLAY_ORDER).dropna(how="all")
    
    # Create Plotly stacked bar chart
    fig = go.Figure()
    
    # Add bars for each regulation assignment
    for i, reg_cat in enumerate(REG_ORDER):
        fig.add_trace(go.Bar(
            x=df_props.index,
            y=df_props[reg_cat],
            name=reg_cat,
            marker_color=REG_COLORS.get(reg_cat, "#cccccc"),
            text=[f"{v:.2f}" if v > 0.05 else "" for v in df_props[reg_cat]],
            textposition="inside",
            hovertemplate=f"<b>{reg_cat}</b><br>Strain: %{{x}}<br>Proportion: %{{y:.3f}}<extra></extra>"
        ))
    
    # Update layout
    fig.update_layout(
        barmode="stack",
        xaxis_title="Strain",
        yaxis_title="Proportion",
        title=f"{tissue} – {subtype}",
        hovermode="x unified",
        template="plotly_white",
        plot_bgcolor="white",
        paper_bgcolor="white",
        height=600,
        width=900,
        showlegend=True,
        legend=dict(title="Reg Assignment", orientation="v", yanchor="top", y=0.99, xanchor="right", x=0.99),
        xaxis=dict(tickangle=-45),
        margin=dict(l=80, r=150, t=100, b=100),
        font=dict(color="black", size=12)
    )
    
    # Add count labels above bars
    for i, strain in enumerate(df_props.index):
        n_count = int(df_props.loc[strain, "n"])
        fig.add_annotation(
            x=strain, y=1.05,
            text=f"n={n_count}",
            showarrow=False,
            font=dict(size=10, color="black"),
            xanchor="center", yanchor="bottom"
        )
    
    return fig


def get_celltype_strain_data_table(
    df: pd.DataFrame,
    tissue: str,
    subtype: str,
    strain: str,
) -> pd.DataFrame:
    """Get a table of genes for a specific tissue/subtype/strain combo with key statistics."""
    subset = df[
        (df["tissue"] == tissue) & 
        (df["subtype"] == subtype) & 
        (df["strain"] == strain)
    ]
    
    if subset.empty:
        return pd.DataFrame()
    
    # Select relevant columns for display
    display_cols = ["gene", "reg_assignment", "Parlog2FC", "Hyblog2FC", "fdr_cis", "fdr_trans"]
    available_cols = [col for col in display_cols if col in subset.columns]
    
    result = subset[available_cols].drop_duplicates(subset=["gene"]).copy()
    result = result.sort_values("gene").reset_index(drop=True)
    
    return result


def make_tissue_composition_interactive_plotly(
    props_df: pd.DataFrame,
    tissue: str,
    sort_by: Literal["n", "conserved_prop", "cis_prop", "trans_prop"] = "n",
) -> Optional[go.Figure]:
    """Create interactive Plotly bar chart for tissue-wide composition across cell types."""
    subset = props_df[props_df["tissue"] == tissue]
    
    if subset.empty:
        return None
    
    # Prepare data
    subtype_counts = subset.groupby("subtype")["total_n"].first().rename("n").sort_index()
    df_props = (
        subset.pivot_table(
            index="subtype", columns="reg_assignment", values="proportion", fill_value=0.0
        )
        .reindex(columns=REG_ORDER, fill_value=0.0)
    )
    df_props = df_props.merge(subtype_counts.to_frame(), left_index=True, right_index=True)
    
    # Apply sorting
    if sort_by == "n":
        df_props = df_props.sort_values("n", ascending=False)
    elif sort_by == "conserved_prop":
        df_props = df_props.sort_values("conserved", ascending=False)
    elif sort_by == "cis_prop":
        df_props = df_props.sort_values("cis", ascending=False)
    elif sort_by == "trans_prop":
        df_props = df_props.sort_values("trans", ascending=False)
    
    # Create Plotly stacked bar chart
    fig = go.Figure()
    
    # Add bars for each regulation assignment
    for reg_cat in REG_ORDER:
        fig.add_trace(go.Bar(
            x=df_props.index,
            y=df_props[reg_cat],
            name=reg_cat,
            marker_color=REG_COLORS.get(reg_cat, "#cccccc"),
            text=[f"{v:.2f}" if v > 0.05 else "" for v in df_props[reg_cat]],
            textposition="inside",
            hovertemplate=f"<b>{reg_cat}</b><br>Cell Type: %{{x}}<br>Proportion: %{{y:.3f}}<extra></extra>"
        ))
    
    # Update layout
    fig.update_layout(
        barmode="stack",
        xaxis_title="Cell Type",
        yaxis_title="Proportion",
        title=f"{tissue} – Cell Type Composition",
        hovermode="x unified",
        template="plotly_white",
        plot_bgcolor="white",
        paper_bgcolor="white",
        height=600,
        width=1000,
        showlegend=True,
        legend=dict(title="Reg Assignment", orientation="v", yanchor="top", y=0.99, xanchor="right", x=0.99),
        xaxis=dict(tickangle=-45),
        margin=dict(l=80, r=150, t=100, b=120),
        font=dict(color="black", size=12)
    )
    
    # Add count labels above bars
    for i, subtype in enumerate(df_props.index):
        n_count = int(df_props.loc[subtype, "n"])
        fig.add_annotation(
            x=subtype, y=1.05,
            text=f"n={n_count}",
            showarrow=False,
            font=dict(size=10, color="black"),
            xanchor="center", yanchor="bottom"
        )
    
    return fig


def get_tissue_composition_data_table(
    df: pd.DataFrame,
    tissue: str,
    subtype: str,
    strain: Optional[str] = None,
) -> pd.DataFrame:
    """Get a table of genes for tissue composition filtering by tissue, subtype, and optionally strain."""
    subset = df[
        (df["tissue"] == tissue) & 
        (df["subtype"] == subtype)
    ]
    
    if strain:
        subset = subset[subset["strain"] == strain]
    
    if subset.empty:
        return pd.DataFrame()
    
    # Select relevant columns for display
    display_cols = ["gene", "strain", "reg_assignment", "Parlog2FC", "Hyblog2FC", "fdr_cis", "fdr_trans"]
    available_cols = [col for col in display_cols if col in subset.columns]
    
    result = subset[available_cols].drop_duplicates(subset=["gene", "strain"]).copy()
    result = result.sort_values(["strain", "gene"]).reset_index(drop=True)
    
    return result


def make_celltype_strain_figure(
    props_df: pd.DataFrame,
    tissue: str,
    subtype: str,
    sort_by: Literal["n", "conserved_prop", "cis_prop", "trans_prop"] = "n",
) -> Optional[plt.Figure]:
    subset = props_df[
        (props_df["tissue"] == tissue) & (props_df["subtype"] == subtype)
    ]
    if subset.empty:
        return None
    strain_counts = subset.groupby("strain")["total_n"].first().rename("n").sort_index()
    df_props = (
        subset.pivot_table(
            index="strain", columns="reg_assignment", values="proportion", fill_value=0.0
        )
        .reindex(columns=REG_ORDER, fill_value=0.0)
    )
    df_props = df_props.merge(strain_counts.to_frame(), left_index=True, right_index=True)
    
    # Apply sort_by logic if requested
    if sort_by == "n":
        df_props = df_props.sort_values("n", ascending=False)
    elif sort_by == "conserved_prop":
        df_props = df_props.sort_values("conserved", ascending=False)
    elif sort_by == "cis_prop":
        df_props = df_props.sort_values("cis", ascending=False)
    elif sort_by == "trans_prop":
        df_props = df_props.sort_values("trans", ascending=False)
    
    # Reorder by STRAINS_DISPLAY_ORDER to ensure consistent display order
    strain_order_available = [s for s in STRAINS_DISPLAY_ORDER if s in df_props.index]
    df_props = df_props.reindex(strain_order_available)
    
    df_n = df_props["n"]
    df_vals = df_props[REG_ORDER]
    n_bars = len(df_vals)
    fig_w = max(6.0, 0.7 * n_bars)
    fig, ax = plt.subplots(figsize=(fig_w, 6.0))
    x = np.arange(n_bars)
    bottom = np.zeros(n_bars)
    for cat in REG_ORDER:
        vals = df_vals[cat].values
        ax.bar(x, vals, bottom=bottom, color=REG_COLORS.get(cat, "#cccccc"), width=0.8, edgecolor="none", label=cat)
        for i, v in enumerate(vals):
            if v > 0.05:
                ax.text(x[i], bottom[i] + v / 2, f"{v:.2f}", ha="center", va="center", fontsize=10, color="black")
        bottom += vals
    for i, strain in enumerate(df_vals.index):
        ax.text(x[i], 1.03, f"{int(df_n.loc[strain])}", ha="center", va="bottom", fontsize=11, color="black")
    ax.set_title(f"{tissue} – {subtype}", fontsize=18)
    ax.set_xlabel("Strain", fontsize=14)
    ax.set_ylabel("Proportion", fontsize=14)
    ax.set_ylim(0, 1.15)
    ax.set_xticks(x)
    ax.set_xticklabels(df_vals.index, rotation=45, ha="right", fontsize=12)
    ax.tick_params(axis="y", labelsize=12, colors="black")
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.legend(title="reg_assignment", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=10, title_fontsize=11, frameon=False)
    fig.tight_layout()
    return fig


def make_tissue_composition_figure(
    props_df: pd.DataFrame,
    tissue: str,
    strain: str,
    sort_by: Literal["n", "conserved_prop", "cis_prop", "trans_prop"] = "n",
) -> Optional[plt.Figure]:
    subset = props_df[
        (props_df["tissue"] == tissue) & (props_df["strain"] == strain)
    ]
    if subset.empty:
        return None
    subtype_counts = subset.groupby("subtype")["total_n"].first().rename("n").sort_index()
    df_props = (
        subset.pivot_table(
            index="subtype", columns="reg_assignment", values="proportion", fill_value=0.0
        )
        .reindex(columns=REG_ORDER, fill_value=0.0)
    )
    df_props = df_props.merge(subtype_counts.to_frame(), left_index=True, right_index=True)
    if sort_by == "n":
        df_props = df_props.sort_values("n", ascending=False)
    elif sort_by == "conserved_prop":
        df_props = df_props.sort_values("conserved", ascending=False)
    elif sort_by == "cis_prop":
        df_props = df_props.sort_values("cis", ascending=False)
    elif sort_by == "trans_prop":
        df_props = df_props.sort_values("trans", ascending=False)
    df_n = df_props["n"]
    df_vals = df_props[REG_ORDER]
    n_bars = len(df_vals)
    fig_w = max(6.0, 0.7 * n_bars)
    fig, ax = plt.subplots(figsize=(fig_w, 6.0))
    x = np.arange(n_bars)
    bottom = np.zeros(n_bars)
    for cat in REG_ORDER:
        vals = df_vals[cat].values
        ax.bar(x, vals, bottom=bottom, color=REG_COLORS.get(cat, "#cccccc"), width=0.8, edgecolor="none", label=cat)
        for i, v in enumerate(vals):
            if v > 0.05:
                ax.text(x[i], bottom[i] + v / 2, f"{v:.2f}", ha="center", va="center", fontsize=10, color="black")
        bottom += vals
    for i, subtype in enumerate(df_vals.index):
        ax.text(x[i], 1.03, f"{int(df_n.loc[subtype])}", ha="center", va="bottom", fontsize=11, color="black")
    ax.set_title(f"{tissue} – {strain}", fontsize=18)
    ax.set_xlabel("Cell type (subtype)", fontsize=14)
    ax.set_ylabel("Proportion", fontsize=14)
    ax.set_ylim(0, 1.15)
    ax.set_xticks(x)
    ax.set_xticklabels(df_vals.index, rotation=45, ha="right", fontsize=12)
    ax.tick_params(axis="y", labelsize=12, colors="black")
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.legend(title="reg_assignment", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=10, title_fontsize=11, frameon=False)
    fig.tight_layout()
    return fig



def _prepare_sankey_series(
    condition_index: Dict[Tuple[str, str, str], pd.Series],
    strain: str, tissue: str, subtype: str,
) -> Optional[pd.Series]:
    key = (strain, tissue, subtype)
    ser = condition_index.get(key)
    if ser is None or ser.empty:
        return None
    ser = ser.where(ser != "cis+trans").dropna()
    return ser


def get_condition_labels(
    left: Tuple[str, str, str],
    right: Tuple[str, str, str],
) -> Tuple[str, str]:
    """Generate meaningful column labels based on which condition differs.
    
    Args:
        left: (strain, tissue, subtype) tuple
        right: (strain, tissue, subtype) tuple
    
    Returns:
        (left_label, right_label) tuple
    """
    if left[0] != right[0]:  # Strain differs
        return f"reg_assignment ({left[0]})", f"reg_assignment ({right[0]})"
    elif left[1] != right[1]:  # Tissue differs
        return f"reg_assignment ({left[1]})", f"reg_assignment ({right[1]})"
    else:  # Subtype differs
        return f"reg_assignment ({left[2]})", f"reg_assignment ({right[2]})"


@st.cache_data(show_spinner=False)
@st.cache_data(show_spinner=False)
def build_sankey_transition_table(
    genes_tuple: Tuple[str, ...],
    condition_index_key_left: Tuple[str, str, str],
    condition_index_key_right: Tuple[str, str, str],
    condition_index_left_data: Tuple[Tuple[str, str], ...],
    condition_index_right_data: Tuple[Tuple[str, str], ...],
    left_label: str,
    right_label: str,
) -> pd.DataFrame:
    """Build the gene table dataframe for a transition (cached for performance).
    
    Args:
        genes_tuple: Tuple of gene names
        condition_index_key_left: Left condition tuple (strain, tissue, subtype)
        condition_index_key_right: Right condition tuple (strain, tissue, subtype)
        condition_index_left_data: Tuple of (gene, reg_assignment) pairs for left
        condition_index_right_data: Tuple of (gene, reg_assignment) pairs for right
        left_label: Column label for left condition
        right_label: Column label for right condition
    
    Returns:
        DataFrame with gene, left_label, right_label columns
    """
    # Convert tuples back to dicts
    left_dict = dict(condition_index_left_data)
    right_dict = dict(condition_index_right_data)
    
    data = []
    for gene in genes_tuple:
        left_reg = left_dict.get(gene, "N/A")
        right_reg = right_dict.get(gene, "N/A")
        data.append({
            "gene": gene,
            left_label: left_reg,
            right_label: right_reg,
        })
    
    return pd.DataFrame(data)


def display_sankey_gene_modal(
    df: pd.DataFrame,
    transition: Tuple[str, str],
    genes: List[str],
    condition_index: Dict[Tuple[str, str, str], pd.Series],
    left_condition: Tuple[str, str, str],
    right_condition: Tuple[str, str, str],
    left_label: str = "From reg_assignment",
    right_label: str = "To reg_assignment",
) -> None:
    """Display a modal with genes for a specific Sankey transition."""
    if not genes:
        st.warning("No genes found for this transition.")
        return
    
    left_state, right_state = transition
    st.markdown(f"### Genes: **{left_state}** → **{right_state}** ({len(genes)} genes)")
    
    # Build simplified dataframe with only gene and reg_assignments
    data = []
    for gene in genes:
        left_reg = condition_index.get(left_condition, pd.Series()).get(gene, "N/A")
        right_reg = condition_index.get(right_condition, pd.Series()).get(gene, "N/A")
        data.append({
            "gene": gene,
            left_label: left_reg,
            right_label: right_reg,
        })
    
    simplified_df = pd.DataFrame(data)
    
    # Display search box for gene filtering with dynamic key to reset on transition change
    transition_key = f"{left_state}_to_{right_state}"
    search_query = st.text_input("Filter genes by name", value="", key=f"sankey_search_{transition_key}")
    
    # Apply gene name filter
    display_df = simplified_df.copy()
    if search_query:
        display_df = display_df[display_df['gene'].str.contains(search_query, case=False, na=False)]
    
    # Display the dataframe
    st.dataframe(display_df, use_container_width=True, hide_index=True)
    
    # Download button with dynamic key
    csv = display_df.to_csv(index=False)
    st.download_button(
        label="Download as CSV",
        data=csv,
        file_name=f"genes_{left_state}_to_{right_state}.csv",
        mime="text/csv",
        key=f"sankey_download_{transition_key}"
    )


def build_sankey_figure(
    condition_index: Dict[Tuple[str, str, str], pd.Series],
    left: Tuple[str, str, str],
    right: Tuple[str, str, str],
    title: str,
) -> Optional[Tuple[go.Figure, int, int, Dict[Tuple[str, str], List[str]]]]:
    left_strain, left_tissue, left_subtype = left
    right_strain, right_tissue, right_subtype = right
    left_ser = _prepare_sankey_series(condition_index, left_strain, left_tissue, left_subtype)
    right_ser = _prepare_sankey_series(condition_index, right_strain, right_tissue, right_subtype)
    if left_ser is None or right_ser is None:
        return None
    shared_genes = sorted(set(left_ser.index) & set(right_ser.index))
    if not shared_genes:
        return None
    left_vals = left_ser.reindex(shared_genes)
    right_vals = right_ser.reindex(shared_genes)
    pair_counts: Dict[Tuple[str, str], int] = {}
    gene_transitions: Dict[Tuple[str, str], List[str]] = {}
    for gene, (l, r) in zip(shared_genes, zip(left_vals, right_vals)):
        if (l not in SANKEY_ORDER) or (r not in SANKEY_ORDER):
            continue
        pair_counts[(l, r)] = pair_counts.get((l, r), 0) + 1
        if (l, r) not in gene_transitions:
            gene_transitions[(l, r)] = []
        gene_transitions[(l, r)].append(gene)
    if not pair_counts:
        return None
    left_nodes = [f"{lbl} (L)" for lbl in SANKEY_ORDER]
    right_nodes = [f"{lbl} (R)" for lbl in SANKEY_ORDER]
    labels = left_nodes + right_nodes
    color_map = [REG_COLORS.get(lbl.split()[0], "#cccccc") for lbl in SANKEY_ORDER]
    node_colors = color_map + color_map
    label_to_index = {}
    for i, lbl in enumerate(SANKEY_ORDER):
        label_to_index[("L", lbl)] = i
        label_to_index[("R", lbl)] = i + len(SANKEY_ORDER)
    sources, targets, values, link_colors = [], [], [], []
    for (l, r), count in pair_counts.items():
        sources.append(label_to_index[("L", l)])
        targets.append(label_to_index[("R", r)])
        values.append(count)
        link_colors.append(REG_COLORS.get(l, "#cccccc"))
    fig = go.Figure(
        data=[
            go.Sankey(
                node=dict(pad=15, thickness=15, line=dict(color="black", width=0.5), label=labels, color=node_colors),
                link=dict(source=sources, target=targets, value=values, color=link_colors),
            )
        ]
    )
    fig.update_layout(title=title, font_size=12, height=650, margin=dict(l=40, r=40, t=60, b=40))
    total_genes = len(set(left_ser.index) | set(right_ser.index))
    shared_n = len(shared_genes)
    return fig, total_genes, shared_n, gene_transitions


def get_unique_sorted(df: pd.DataFrame, col: str) -> List[str]:
    if col not in df.columns:
        return []
    return sorted(map(str, df[col].dropna().unique().tolist()))


def gene_search_matches(df: pd.DataFrame, query: str) -> List[str]:
    if "gene" not in df.columns or not query:
        return []
    q = query.strip().lower()
    genes = df["gene"].astype(str).unique()
    return sorted([g for g in genes if q in g.lower()])[:100]


def subset_for_celltype_view(
    df: pd.DataFrame,
    subtype: str,
    tissue: Optional[str],
    strain: Optional[str],
    genes: Optional[List[str]],
    reg_include: Optional[List[str]],
    fdr_cis_max: Optional[float],
    fdr_trans_max: Optional[float],
) -> pd.DataFrame:
    return _apply_basic_filters(
        df=df, tissue=tissue, strain=strain, subtype=subtype,
        genes=genes, reg_include=reg_include, fdr_cis_max=fdr_cis_max, fdr_trans_max=fdr_trans_max,
    )


# ---------------------------------------------------------------------------
# Gene expression viewer: data discovery and loading
# ---------------------------------------------------------------------------
@st.cache_data(show_spinner=False)
def list_available_tissues_strains(_base_path: str = BASE_PATH) -> Tuple[List[str], Dict[str, List[str]]]:
    """Return (tissues_sorted, {tissue: [strain, ...]}). Uses suffixes to find strain names."""
    tissues: List[str] = []
    by_tissue: Dict[str, List[str]] = {}
    if not os.path.isdir(_base_path):
        return [], {}
    for name in sorted(os.listdir(_base_path)):
        tissue_dir = os.path.join(_base_path, name)
        if not os.path.isdir(tissue_dir):
            continue
        tissues.append(name)
        strains = []
        for f in os.listdir(tissue_dir):
            if f.endswith(DATA_SUFFIX):
                # e.g. CAST_EiJ_xgener_input_dataframe_FILTERED.csv -> CAST_EiJ
                stem = f[: -len(DATA_SUFFIX)]
                strains.append(stem)
        by_tissue[name] = sorted(strains)
    return tissues, by_tissue


def _normalize_gene_names(genes: pd.Index) -> pd.Index:
    return pd.Index([str(g).strip().strip('"') for g in genes])


@st.cache_data(show_spinner=False)
def load_adata(tissue: str, strain: str, _base_path: str = BASE_PATH) -> Optional["AnnData"]:
    import anndata as ad
    tissue_dir = os.path.join(_base_path, tissue)
    data_path = os.path.join(tissue_dir, strain + DATA_SUFFIX)
    meta_path = os.path.join(tissue_dir, strain + META_SUFFIX)
    if not os.path.isfile(data_path) or not os.path.isfile(meta_path):
        return None
    counts = pd.read_csv(data_path, index_col=0)
    meta = pd.read_csv(meta_path, index_col=0)
    common = counts.index.intersection(meta.index)
    if len(common) == 0:
        return None
    counts = counts.loc[common]
    meta = meta.loc[common]
    X = counts.values.astype(np.float32)
    adata = ad.AnnData(X=X, obs=meta, var=pd.DataFrame(index=counts.columns))
    adata.var_names = _normalize_gene_names(adata.var_names)
    adata.layers["raw_counts"] = X.copy()
    totals = X.sum(axis=1, keepdims=True)
    totals[totals == 0] = 1
    adata.layers["CPM"] = (X / totals) * 1e6
    return adata


def _allele_palette_for_strain(strain: str) -> Dict[str, str]:
    P1_COLOR = GENO_DICT["B6J"]
    short = FOUNDER_SHORTNAME.get(strain, strain)
    P2_COLOR = GENO_DICT.get(short, "#888888")
    return {"P1": P1_COLOR, "H1": P1_COLOR, "H2": P2_COLOR, "P2": P2_COLOR}


def add_colorbars_below_ticks_for_strain(
    ax,
    order=("P1", "H1", "H2", "P2"),
    *,
    strain: str = "",
    P1_COLOR: str = "#C0BFBF",
    P2_COLOR: str = "#888888",
    add_genotype_text: bool = True,
    bar_height: float = 0.08,
    bar_spacing: float = 0.03,
    y_offset: float = 0.25,
    # text controls
    genotype_text_fs: int = 7,
    genotype_text_weight: str = "bold",
    geno_text_y_frac: float = -0.5,
):
    """
    Add two color bars below the boxplot:
    TOP bar: Allele-of-origin (P1 color for P1/H1, P2 color for P2/H2)
    BOTTOM bar: Genotype (solid blocks for P1/P2, diagonal split for F1 hybrids)
    
    Uses axis-relative coordinates (transAxes) so bars are always visible within each subplot.
    Positions drawn below the xaxis but within the axis bounds.
    
    Args:
        ax: matplotlib axis
        order: tuple of alleles in order (e.g., ("P1", "H1", "H2", "P2"))
        strain: strain name (e.g., "CAST_EiJ") for getting parent names
        P1_COLOR: color for P1 allele (B6J)
        P2_COLOR: color for P2 allele (strain-specific)
        add_genotype_text: whether to add genotype labels
        bar_height: height of each bar in axis coordinates (default: 0.08)
        bar_spacing: spacing between bars in axis coordinates (default: 0.03)
        y_offset: vertical offset from bottom of plot (default: 0.25)
        genotype_text_fs: font size for genotype labels
        genotype_text_weight: font weight for genotype labels
        geno_text_y_frac: vertical position fraction within genotype bar
    """
    from matplotlib.patches import Polygon
    
    # Get parent names
    p1_name = "B6J"
    p2_name = FOUNDER_SHORTNAME.get(strain, strain) if strain else "P2"
    
    n_ticks = len(order)
    
    # Position bars using axis coordinates (0-1 scale)
    # y_offset controls how far below the plot bottom the bars are positioned
    y_allele = -y_offset
    y_geno = y_allele - bar_height - bar_spacing
    
    # x positions in axis coordinates (0 = left edge, 1 = right edge)
    x_left = 0
    x_right = 1
    width_per_tick = (x_right - x_left) / n_ticks
    
    # y position for genotype text (within genotype bar)
    geno_text_y = y_geno + bar_height * geno_text_y_frac
    
    # track hybrid span for centered F1 label
    hybrid_x0 = None
    hybrid_x1 = None
    
    for i, tick_label in enumerate(order):
        x_start = x_left + i * width_per_tick
        x0 = x_start
        x1 = x_start + width_per_tick
        
        # -------------------------
        # TOP BAR: Allele-of-origin
        # -------------------------
        ao_color = P1_COLOR if tick_label in ("P1", "H1") else P2_COLOR
        rect_allele = plt.Rectangle(
            (x_start, y_allele), width_per_tick, bar_height,
            transform=ax.transAxes,
            facecolor=ao_color, edgecolor="black", linewidth=0.5, clip_on=False
        )
        ax.add_patch(rect_allele)
        
        # -------------------------
        # BOTTOM BAR: Genotype
        # -------------------------
        if tick_label == "P1":
            rect_geno = plt.Rectangle(
                (x_start, y_geno), width_per_tick, bar_height,
                transform=ax.transAxes,
                facecolor=P1_COLOR, edgecolor="black", linewidth=0.5, clip_on=False
            )
            ax.add_patch(rect_geno)
            if add_genotype_text:
                ax.text(
                    (x0 + x1) / 2, geno_text_y,
                    p1_name,
                    transform=ax.transAxes,
                    ha="center", va="center",
                    fontsize=genotype_text_fs,
                    fontweight=genotype_text_weight,
                )
        
        elif tick_label == "P2":
            rect_geno = plt.Rectangle(
                (x_start, y_geno), width_per_tick, bar_height,
                transform=ax.transAxes,
                facecolor=P2_COLOR, edgecolor="black", linewidth=0.5, clip_on=False
            )
            ax.add_patch(rect_geno)
            if add_genotype_text:
                ax.text(
                    (x0 + x1) / 2, geno_text_y,
                    p2_name,
                    transform=ax.transAxes,
                    ha="center", va="center",
                    fontsize=genotype_text_fs,
                    fontweight=genotype_text_weight,
                )
        
        else:
            # Hybrid (H1 or H2) - record span for single centered F1 label
            hybrid_x0 = x0 if hybrid_x0 is None else min(hybrid_x0, x0)
            hybrid_x1 = x1 if hybrid_x1 is None else max(hybrid_x1, x1)
            
            # Diagonal split: upper-left P1, lower-right P2
            y0 = y_geno
            y1 = y_geno + bar_height
            
            # Upper-left triangle (P1)
            tri1 = Polygon(
                [(x0, y1), (x0, y0), (x1, y1)],
                closed=True,
                transform=ax.transAxes,
                facecolor=P1_COLOR,
                edgecolor="none",
                linewidth=0,
                clip_on=False
            )
            ax.add_patch(tri1)
            
            # Lower-right triangle (P2)
            tri2 = Polygon(
                [(x1, y0), (x1, y1), (x0, y0)],
                closed=True,
                transform=ax.transAxes,
                facecolor=P2_COLOR,
                edgecolor="none",
                linewidth=0,
                clip_on=False
            )
            ax.add_patch(tri2)
            
            # Border around hybrid box
            rect_border = plt.Rectangle(
                (x0, y0), width_per_tick, bar_height,
                transform=ax.transAxes,
                facecolor="none",
                edgecolor="black",
                linewidth=0.5,
                clip_on=False
            )
            ax.add_patch(rect_border)
    
    # -------------------------
    # ONE centered F1 label (across H1+H2)
    # -------------------------
    if add_genotype_text and (hybrid_x0 is not None) and (hybrid_x1 is not None):
        ax.text(
            (hybrid_x0 + hybrid_x1) / 2,
            geno_text_y,
            "F1",
            transform=ax.transAxes,
            ha="center", va="center",
            fontsize=genotype_text_fs,
            fontweight=genotype_text_weight,
        )


def plot_gene_boxplot_for_strain(
    adata: "AnnData",
    gene: str,
    subtype: str,
    strain: str,
    *,
    ax: Optional[plt.Axes] = None,
    layer: str = "CPM",
    show_points: bool = True,
    point_size: float = 8,
    jitter: float = 0.15,
    y_scale: Literal["linear", "log1p"] = "linear",
    y_min: float = 0,
    downsample_per_allele: Optional[int] = None,
    reg_assignment: Optional[str] = None,
) -> Optional[Tuple[plt.Figure, plt.Axes, int]]:
    import anndata as ad
    import logging
    
    if gene not in adata.var_names:
        logging.warning(f"[{strain}] Gene '{gene}' not found in var_names. Available: {len(adata.var_names)} genes")
        return None
    
    obs = adata.obs
    if "subtype" not in obs.columns or "Allele" not in obs.columns:
        logging.warning(f"[{strain}] Missing columns. Has subtype: {'subtype' in obs.columns}, Has Allele: {'Allele' in obs.columns}. Available: {obs.columns.tolist()}")
        return None
    
    sub = obs.loc[obs["subtype"] == subtype]
    if sub.empty:
        logging.warning(f"[{strain}] No cells found for subtype '{subtype}'. Available subtypes: {obs['subtype'].unique().tolist()}")
        return None
    
    logging.debug(f"[{strain}] {len(sub)} cells with subtype '{subtype}'")
    
    sub = sub.copy()
    try:
        # Get gene index in var_names
        gene_idx = adata.var_names.get_loc(gene)
        
        # Extract expression for the filtered cells and gene
        if layer == "raw_counts":
            raw = adata.layers["raw_counts"][np.asarray(adata.obs_names.get_indexer(sub.index)), gene_idx]
        else:
            raw = adata.layers[layer][np.asarray(adata.obs_names.get_indexer(sub.index)), gene_idx]
        
        sub["_expr"] = np.ravel(raw)
        logging.debug(f"[{strain}] Expression extracted: min={sub['_expr'].min():.2f}, max={sub['_expr'].max():.2f}, mean={sub['_expr'].mean():.2f}")
    except Exception as e:
        logging.error(f"[{strain}] Failed to extract expression for gene '{gene}': {e}")
        return None
    
    sub = sub[sub["Allele"].isin(ALLELE_ORDER)]
    if sub.empty:
        logging.warning(f"[{strain}] No cells with valid alleles. Expected one of {ALLELE_ORDER}. Available: {obs['Allele'].unique().tolist()}")
        return None
    
    logging.debug(f"[{strain}] {len(sub)} cells after allele filtering. Alleles: {sub['Allele'].unique().tolist()}")
    if downsample_per_allele is not None and downsample_per_allele > 0:
        sub = sub.groupby("Allele", group_keys=False).apply(
            lambda g: g.sample(n=min(len(g), downsample_per_allele), random_state=42)
        ).reset_index(drop=True)
    palette = _allele_palette_for_strain(strain)
    order = [a for a in ALLELE_ORDER if a in sub["Allele"].unique()]
    if not order:
        return None
    n_cells = len(sub)
    short = FOUNDER_SHORTNAME.get(strain, strain)
    if ax is None:
        fig, ax = plt.subplots(figsize=(4, 2.5))
    else:
        fig = ax.figure
    positions = np.arange(len(order))
    boxes_data = [sub.loc[sub["Allele"] == a, "_expr"].values for a in order]
    if y_scale == "log1p":
        boxes_data = [np.log1p(np.maximum(b, y_min)) for b in boxes_data]
    else:
        boxes_data = [np.maximum(b, y_min) for b in boxes_data]
    bp = ax.boxplot(boxes_data, positions=positions, widths=0.5, patch_artist=True, showfliers=False)
    for i, (a, patch) in enumerate(zip(order, bp["boxes"])):
        patch.set_facecolor(palette.get(a, "#cccccc"))
    if show_points:
        for i, a in enumerate(order):
            y = boxes_data[i]
            x_jitter = np.random.uniform(-jitter, jitter, size=len(y)) + positions[i]
            ax.scatter(x_jitter, y, s=point_size, c="black", alpha=1.0, edgecolors="none")
    ax.set_xticks(positions)
    ax.set_xticklabels(order)
    ax.set_ylabel("CPM" if layer == "CPM" else "Counts")
    # Ensure y-axis starts at 0
    ax.set_ylim(bottom=0)
    title = f"{strain} trio"
    if reg_assignment:
        title += f" | {reg_assignment}"
    ax.set_title(title)
    
    # Add color bars below the plot
    P1_COLOR = palette.get("P1", "#cccccc")
    P2_COLOR = palette.get("P2", "#cccccc")
    add_colorbars_below_ticks_for_strain(
        ax,
        order=tuple(order),
        strain=strain,
        P1_COLOR=P1_COLOR,
        P2_COLOR=P2_COLOR,
        add_genotype_text=True,
        bar_height=0.08,
        bar_spacing=0.03,
        y_offset=0.25,
    )
    
    return fig, ax, n_cells


def plot_gene_across_strains(
    tissue: str,
    subtype: str,
    gene: str,
    *,
    layer: str = "CPM",
    show_points: bool = True,
    point_size: float = 8,
    jitter: float = 0.15,
) -> Tuple[Optional[plt.Figure], List[str]]:  
    """Returns (figure with 2-row layout showing all 7 strains, list of warning messages)."""
    import matplotlib.gridspec as gridspec
    
    # Load results table to check for gene-tissue-subtype combinations
    results_path = os.path.join(DATA_ROOT, "cis_trans_results_table.csv")
    results_df = load_results_table(results_path)
    
    # Check if gene-tissue-subtype combination exists in results
    gene_in_results = results_df[
        (results_df["gene"] == gene) & 
        (results_df["tissue"] == tissue) & 
        (results_df["subtype"] == subtype)
    ]
    
    if gene_in_results.empty:
        # Gene not detected in cis_trans_results_table for this tissue-subtype
        return None, [f"Gene '{gene}' for subtype '{subtype}' not found in cis_trans results for tissue '{tissue}."]
    
    _, by_tissue = list_available_tissues_strains()
    strains_in_tissue = by_tissue.get(tissue, [])
    strains_to_plot = [s for s in STRAINS_DISPLAY_ORDER if s in strains_in_tissue]
    warnings: List[str] = []
    
    # Load data for all strains
    strain_data = {}
    for strain in strains_to_plot:
        adata = load_adata(tissue, strain)
        if adata is None:
            strain_data[strain] = None
            warnings.append(f"Could not load data for {strain}.")
        else:
            strain_data[strain] = adata
    
    # Create 2-row figure with GridSpec (4 cols on top, 3 cols on bottom)
    # Extra height and bottom margin to accommodate color bars below
    fig = plt.figure(figsize=(12, 8))
    gs = gridspec.GridSpec(2, 4, figure=fig, hspace=0.55, wspace=0.3, bottom=0.15, top=0.91)
    
    # Top row: 4 strains
    axes_top = [fig.add_subplot(gs[0, i]) for i in range(4)]
    # Bottom row: 3 strains
    axes_bottom = [fig.add_subplot(gs[1, i]) for i in range(3)]
    all_axes = axes_top + axes_bottom
    
    # Map strains to axes positions
    has_data = False
    
    for idx, strain in enumerate(strains_to_plot):
        if idx >= len(all_axes):
            break
        ax = all_axes[idx]
        adata = strain_data[strain]
        
        if adata is None:
            ax.text(0.5, 0.5, f"{strain}\nnot detected", ha="center", va="center", fontsize=10)
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.axis("off")
        elif gene not in adata.var_names:
            ax.text(0.5, 0.5, f"{strain}\nnot detected", ha="center", va="center", fontsize=10)
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.axis("off")
        elif subtype not in adata.obs["subtype"].astype(str).values:
            ax.text(0.5, 0.5, f"{strain}\nnot detected", ha="center", va="center", fontsize=10)
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.axis("off")
        else:
            # Check if this gene-tissue-subtype-strain combo exists in results table
            strain_reg = results_df[
                (results_df["gene"] == gene) & 
                (results_df["tissue"] == tissue) & 
                (results_df["subtype"] == subtype) &
                (results_df["strain"] == strain)
            ]
            
            if strain_reg.empty:
                # Gene not found in results table for this strain
                ax.text(0.5, 0.5, f"{strain}\nnot detected", ha="center", va="center", fontsize=10)
                ax.set_xlim(0, 1)
                ax.set_ylim(0, 1)
                ax.axis("off")
            else:
                reg_assignment = strain_reg["reg_assignment"].iloc[0]
                plot_gene_boxplot_for_strain(
                    adata, gene, subtype, strain,
                    ax=ax,
                    layer=layer, show_points=show_points, point_size=point_size, jitter=jitter,
                    y_scale="linear", y_min=0.0, downsample_per_allele=None,
                    reg_assignment=reg_assignment,
                )
                has_data = True
    
    # Hide unused subplots (if fewer than 7 strains available)
    for idx in range(len(strains_to_plot), len(all_axes)):
        all_axes[idx].axis("off")
    
    if not has_data:
        return None, warnings
    
    fig.suptitle(f"{subtype} | {gene}", fontsize=14, fontweight="bold")
    return fig, warnings

# ---------------------------------------------------------------------------
# Gene expression viewer: gene list for tissue
# ---------------------------------------------------------------------------   
@st.cache_data(show_spinner=False)
def get_gene_list_for_tissue(tissue: str, _base_path: str = BASE_PATH) -> List[str]:
    """Return sorted list of gene names from the first available strain in a tissue (fast loading)."""
    _, by_tissue = list_available_tissues_strains()
    strains = by_tissue.get(tissue, [])
    
    # Load from just the first available strain (fast - 1 file instead of 7)
    for strain in strains:
        adata = load_adata(tissue, strain)
        if adata is not None:
            return sorted(adata.var_names.tolist())
    
    return []
