#!/usr/bin/env python3
"""
Streamlit GUI for multisource lollipop cleavage-site plots.

Place this file in the same folder as `run_analysis_lollipop_any_protein.py`, then run:

    streamlit run streamlit_lollipop_app.py

The app supports two workflows:
  1. Extract a target protein from the three source Excel datasets.
  2. Plot user-supplied custom lollipop points.
"""

from __future__ import annotations

import io
import tempfile
import zipfile
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple

import pandas as pd
import streamlit as st

from run_analysis_lollipop_any_protein import (
    build_cluster_summary,
    build_lollipop_points,
    infer_protein_length,
    make_output_prefix,
    make_target_spec,
    plot_cluster_matrix,
    plot_lollipop_multisource_hotspots,
    read_hunter_arabidopsis,
    read_metacaspase,
    read_semitryptome,
    safe_numeric,
    split_alias_string,
)


# -----------------------------------------------------------------------------
# Page configuration
# -----------------------------------------------------------------------------

st.set_page_config(
    page_title="Protein lollipop plotter",
    page_icon="🧬",
    layout="wide",
)

DATASET_ORDER = ["Semitryptome", "HUNTER MJ", "Metacaspase matrix"]


# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

def save_uploaded_file(uploaded_file, destination: Path) -> Path:
    """Write a Streamlit UploadedFile object to disk and return the file path."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(uploaded_file.getbuffer())
    return destination


def csv_bytes(df: pd.DataFrame) -> bytes:
    """Return a UTF-8 CSV representation of a dataframe."""
    return df.to_csv(index=False).encode("utf-8")


def make_download_zip(files: Dict[str, bytes]) -> bytes:
    """Create an in-memory ZIP archive from a mapping of archive names to bytes."""
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, mode="w", compression=zipfile.ZIP_DEFLATED) as zf:
        for archive_name, content in files.items():
            zf.writestr(archive_name, content)
    buffer.seek(0)
    return buffer.getvalue()


def read_optional_file_bytes(path: Path) -> Optional[bytes]:
    return path.read_bytes() if path.exists() else None


def normalize_alias_text(alias_text: str) -> list[str]:
    aliases: list[str] = []
    for token in split_alias_string(alias_text):
        aliases.append(token)
    return aliases


def render_download_panel(result_key: str) -> None:
    """Render generated plot, tables, and download buttons from session state."""
    result = st.session_state.get(result_key)
    if not result:
        return

    st.divider()
    st.subheader("Generated output")

    metrics = result.get("metrics", {})
    if metrics:
        cols = st.columns(len(metrics))
        for col, (label, value) in zip(cols, metrics.items()):
            col.metric(label, value)

    if result.get("main_png"):
        st.image(result["main_png"], caption="Lollipop plot", use_container_width=True)

    dcols = st.columns(4)
    with dcols[0]:
        st.download_button(
            "Download PNG",
            data=result["main_png"],
            file_name=result["main_png_name"],
            mime="image/png",
            use_container_width=True,
        )
    with dcols[1]:
        st.download_button(
            "Download PDF",
            data=result["main_pdf"],
            file_name=result["main_pdf_name"],
            mime="application/pdf",
            use_container_width=True,
        )
    with dcols[2]:
        st.download_button(
            "Download points CSV",
            data=result["points_csv"],
            file_name=result["points_csv_name"],
            mime="text/csv",
            use_container_width=True,
        )
    with dcols[3]:
        st.download_button(
            "Download all ZIP",
            data=result["zip_bytes"],
            file_name=result["zip_name"],
            mime="application/zip",
            use_container_width=True,
        )

    with st.expander("Preview lollipop points table", expanded=False):
        st.dataframe(result["points_df"], use_container_width=True)

    with st.expander("Preview cluster summary table", expanded=False):
        st.dataframe(result["clusters_df"], use_container_width=True)

    if result.get("matrix_png"):
        with st.expander("Show cluster evidence matrix", expanded=False):
            st.image(result["matrix_png"], caption="Cluster evidence matrix", use_container_width=True)


def build_result_payload(
    prefix: str,
    out_dir: Path,
    points: pd.DataFrame,
    clusters: pd.DataFrame,
    main_png: Path,
    main_pdf: Path,
    matrix_png: Optional[Path] = None,
    matrix_pdf: Optional[Path] = None,
    extra_tables: Optional[Dict[str, pd.DataFrame]] = None,
    metrics: Optional[Dict[str, int]] = None,
) -> Dict[str, object]:
    """Collect plot/table bytes in a Streamlit-friendly result dictionary."""
    points_csv_name = f"{prefix}_lollipop_points.csv"
    clusters_csv_name = f"{prefix}_lollipop_clusters.csv"

    zip_files: Dict[str, bytes] = {
        f"plots/{main_png.name}": main_png.read_bytes(),
        f"plots/{main_pdf.name}": main_pdf.read_bytes(),
        f"tables/{points_csv_name}": csv_bytes(points),
        f"tables/{clusters_csv_name}": csv_bytes(clusters),
    }

    matrix_png_bytes = read_optional_file_bytes(matrix_png) if matrix_png else None
    matrix_pdf_bytes = read_optional_file_bytes(matrix_pdf) if matrix_pdf else None
    if matrix_png_bytes and matrix_png:
        zip_files[f"plots/{matrix_png.name}"] = matrix_png_bytes
    if matrix_pdf_bytes and matrix_pdf:
        zip_files[f"plots/{matrix_pdf.name}"] = matrix_pdf_bytes

    for table_name, table in (extra_tables or {}).items():
        zip_files[f"tables/{table_name}"] = csv_bytes(table)

    return {
        "main_png": main_png.read_bytes(),
        "main_pdf": main_pdf.read_bytes(),
        "matrix_png": matrix_png_bytes,
        "points_csv": csv_bytes(points),
        "clusters_csv": csv_bytes(clusters),
        "zip_bytes": make_download_zip(zip_files),
        "points_df": points,
        "clusters_df": clusters,
        "metrics": metrics or {},
        "main_png_name": main_png.name,
        "main_pdf_name": main_pdf.name,
        "points_csv_name": points_csv_name,
        "clusters_csv_name": clusters_csv_name,
        "zip_name": f"{prefix}_lollipop_output.zip",
    }


def coerce_custom_points(df: pd.DataFrame) -> pd.DataFrame:
    """Validate and standardize manually supplied lollipop points."""
    required = {"dataset", "cleavage_site", "weight"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Custom table is missing required columns: {', '.join(missing)}")

    out = df.copy()
    out["dataset"] = out["dataset"].astype(str).str.strip()
    out["cleavage_site"] = safe_numeric(out["cleavage_site"])
    out["weight"] = safe_numeric(out["weight"]).fillna(1)

    out = out.dropna(subset=["cleavage_site"])
    out = out[out["dataset"].isin(DATASET_ORDER)]

    if out.empty:
        raise ValueError(
            "No valid custom rows remained after filtering. Dataset must be one of: "
            + ", ".join(DATASET_ORDER)
        )

    out["cleavage_site"] = out["cleavage_site"].astype(int)
    if "peptide" not in out.columns:
        out["peptide"] = ""
    if "start" not in out.columns:
        out["start"] = out["cleavage_site"]
    if "end" not in out.columns:
        out["end"] = out["cleavage_site"]
    if "n_observations" not in out.columns:
        out["n_observations"] = 1

    return out[["dataset", "cleavage_site", "weight", "peptide", "start", "end", "n_observations"]]


# -----------------------------------------------------------------------------
# Header and workflow selection
# -----------------------------------------------------------------------------

st.title("Protein lollipop cleavage-site plotter")
st.write(
    "Use this interface to generate multisource lollipop plots for a selected protein. "
    "You can either upload the three original Excel datasets and extract a target automatically, "
    "or provide custom residue-level points manually."
)

mode = st.sidebar.radio(
    "Input mode",
    ["Extract from Excel datasets", "Custom lollipop points"],
)

st.sidebar.markdown("---")
st.sidebar.caption(
    "Expected dataset labels: Semitryptome, HUNTER MJ, Metacaspase matrix. "
    "Dot size is scaled within each dataset track."
)


# -----------------------------------------------------------------------------
# Workflow 1: source Excel extraction
# -----------------------------------------------------------------------------

if mode == "Extract from Excel datasets":
    st.header("Extract target protein from source Excel files")

    with st.form("excel_form"):
        st.subheader("1. Upload input files")
        c1, c2, c3 = st.columns(3)
        with c1:
            semi_file = st.file_uploader(
                "Arabidopsis semitryptome Excel",
                type=["xlsx"],
                key="semi_file",
            )
        with c2:
            hunter_file = st.file_uploader(
                "HUNTER shoot MJ Excel",
                type=["xlsx"],
                key="hunter_file",
            )
        with c3:
            meta_file = st.file_uploader(
                "ProteinBased P10/P10 metacaspase Excel",
                type=["xlsx"],
                key="meta_file",
            )

        st.subheader("2. Define target protein")
        target_col, uniprot_col, alias_col = st.columns([1.1, 1.0, 2.0])
        with target_col:
            target_query = st.text_input(
                "Primary target ID or name",
                value="AT3G09260.1",
                help="Usually the TAIR protein ID for Arabidopsis, e.g. AT3G09260.1.",
            )
        with uniprot_col:
            uniprot_id = st.text_input(
                "UniProt ID",
                value="Q9SR37",
                help="Optional UniProt accession for the same protein, e.g. Q9SR37.",
            )
        with alias_col:
            alias_text = st.text_input(
                "Aliases / gene names",
                value="PYK10; BGLU23; LEB; PSR3.1",
                help="Separate terms with semicolons, commas, or pipes.",
            )

        st.subheader("3. Plot options")
        o1, o2, o3, o4 = st.columns(4)
        with o1:
            protein_length = st.number_input(
                "Protein length",
                min_value=0,
                value=525,
                step=1,
                help="Use 0 to infer length from the extracted data.",
            )
        with o2:
            cluster_gap = st.number_input(
                "Cluster gap, aa",
                min_value=1,
                value=10,
                step=1,
            )
        with o3:
            shade_support = st.number_input(
                "Shade clusters supported by ≥ datasets",
                min_value=1,
                max_value=3,
                value=2,
                step=1,
            )
        with o4:
            x_tick_gap = st.number_input(
                "X-axis tick gap",
                min_value=5,
                value=20,
                step=5,
            )

        m1, m2 = st.columns(2)
        with m1:
            plot_every_row = st.checkbox(
                "Plot every extracted row instead of aggregating exact duplicate sites",
                value=False,
            )
        with m2:
            include_isoforms = st.checkbox(
                "Match metacaspase Isoforms/Description fields",
                value=False,
                help="Default matching uses the primary Protein column only.",
            )

        submitted = st.form_submit_button("Generate lollipop plot", use_container_width=True)

    if submitted:
        if not all([semi_file, hunter_file, meta_file]):
            st.error("Please upload all three Excel files before generating the plot.")
        else:
            try:
                with st.spinner("Extracting target rows and generating plots..."):
                    run_dir = Path(tempfile.mkdtemp(prefix="lollipop_streamlit_"))
                    input_dir = run_dir / "input"
                    out_dir = run_dir / "output"
                    plots_dir = out_dir / "plots"
                    tables_dir = out_dir / "tables"
                    plots_dir.mkdir(parents=True, exist_ok=True)
                    tables_dir.mkdir(parents=True, exist_ok=True)

                    semi_path = save_uploaded_file(semi_file, input_dir / semi_file.name)
                    hunter_path = save_uploaded_file(hunter_file, input_dir / hunter_file.name)
                    meta_path = save_uploaded_file(meta_file, input_dir / meta_file.name)

                    target = make_target_spec(
                        query=target_query,
                        uniprot_id=uniprot_id,
                        aliases=normalize_alias_text(alias_text),
                        protein_length=int(protein_length) if protein_length else None,
                    )
                    prefix = target.output_prefix

                    semi = read_semitryptome(semi_path, target)
                    hunter = read_hunter_arabidopsis(hunter_path, target)
                    meta = read_metacaspase(meta_path, target, match_isoforms=include_isoforms)

                    points = build_lollipop_points(
                        semi,
                        hunter,
                        meta,
                        aggregate_sites=not plot_every_row,
                    )
                    if points.empty:
                        raise ValueError(
                            "No lollipop points were produced. Check the target ID/aliases and the uploaded files."
                        )

                    clusters = build_cluster_summary(points, gap=int(cluster_gap))
                    plot_len = int(protein_length) if protein_length else infer_protein_length(semi, hunter, meta)

                    main_png = plots_dir / f"{prefix}_lollipop_multisource_hotspots.png"
                    main_pdf = plots_dir / f"{prefix}_lollipop_multisource_hotspots.pdf"
                    matrix_png = plots_dir / f"{prefix}_cluster_evidence_matrix.png"
                    matrix_pdf = plots_dir / f"{prefix}_cluster_evidence_matrix.pdf"

                    plot_lollipop_multisource_hotspots(
                        points=points,
                        cluster_summary=clusters,
                        protein_len=plot_len,
                        out_png=main_png,
                        out_pdf=main_pdf,
                        shade_min_dataset_support=int(shade_support),
                        x_tick_gap=int(x_tick_gap),
                        title=f"{target.display_name}: cleavage-site evidence",
                        x_label=f"{target.query} residue coordinate",
                    )
                    plot_cluster_matrix(
                        cluster_summary=clusters,
                        out_png=matrix_png,
                        out_pdf=matrix_pdf,
                        title=f"{target.query} cluster-level evidence matrix",
                    )

                    st.session_state["excel_result"] = build_result_payload(
                        prefix=prefix,
                        out_dir=out_dir,
                        points=points,
                        clusters=clusters,
                        main_png=main_png,
                        main_pdf=main_pdf,
                        matrix_png=matrix_png,
                        matrix_pdf=matrix_pdf,
                        extra_tables={
                            f"{prefix}_semitryptome.csv": semi,
                            f"{prefix}_hunter_mj.csv": hunter,
                            f"{prefix}_metacaspase.csv": meta,
                        },
                        metrics={
                            "Semitryptome rows": len(semi),
                            "HUNTER rows": len(hunter),
                            "Metacaspase rows": len(meta),
                            "Lollipop points": len(points),
                            "Clusters": len(clusters),
                            "Protein length": plot_len,
                        },
                    )
                st.success("Plot generated successfully.")
            except Exception as exc:
                st.error(str(exc))

    render_download_panel("excel_result")


# -----------------------------------------------------------------------------
# Workflow 2: custom point plotting
# -----------------------------------------------------------------------------

else:
    st.header("Plot custom lollipop points")
    st.write(
        "Provide residue coordinates directly. This mode is useful when sites were manually curated, "
        "or when you want to plot a protein without using the source Excel extraction step."
    )

    default_custom = pd.DataFrame(
        {
            "dataset": ["Semitryptome", "HUNTER MJ", "Metacaspase matrix"],
            "cleavage_site": [23, 25, 284],
            "weight": [5.0, 2.3, 0.61],
            "peptide": ["", "", "DSQDGASIDR"],
        }
    )

    with st.form("custom_form"):
        st.subheader("1. Protein and plotting options")
        c1, c2, c3, c4 = st.columns(4)
        with c1:
            custom_target = st.text_input("Protein label", value="AT3G09260.1")
        with c2:
            custom_display = st.text_input("Plot title label", value="AT3G09260.1 / PYK10")
        with c3:
            custom_length = st.number_input("Protein length", min_value=1, value=525, step=1)
        with c4:
            custom_prefix = st.text_input("Output prefix", value="custom_lollipop")

        c5, c6, c7 = st.columns(3)
        with c5:
            custom_cluster_gap = st.number_input("Cluster gap, aa", min_value=1, value=10, step=1)
        with c6:
            custom_shade_support = st.number_input(
                "Shade clusters supported by ≥ datasets",
                min_value=1,
                max_value=3,
                value=2,
                step=1,
            )
        with c7:
            custom_x_tick_gap = st.number_input("X-axis tick gap", min_value=5, value=20, step=5)

        st.subheader("2. Enter custom sites")
        input_style = st.radio(
            "How do you want to provide points?",
            ["Editable table", "Paste CSV"],
            horizontal=True,
        )

        if input_style == "Editable table":
            custom_df = st.data_editor(
                default_custom,
                num_rows="dynamic",
                use_container_width=True,
                column_config={
                    "dataset": st.column_config.SelectboxColumn(
                        "dataset",
                        options=DATASET_ORDER,
                        required=True,
                    ),
                    "cleavage_site": st.column_config.NumberColumn("cleavage_site", min_value=1, step=1),
                    "weight": st.column_config.NumberColumn("weight", min_value=0.0, step=0.1),
                },
            )
            pasted_csv = ""
        else:
            pasted_csv = st.text_area(
                "Paste CSV with columns: dataset, cleavage_site, weight, optional peptide/start/end",
                value=default_custom.to_csv(index=False),
                height=180,
            )
            custom_df = pd.DataFrame()

        custom_submit = st.form_submit_button("Generate custom plot", use_container_width=True)

    if custom_submit:
        try:
            with st.spinner("Generating custom lollipop plot..."):
                if input_style == "Paste CSV":
                    raw_df = pd.read_csv(io.StringIO(pasted_csv))
                else:
                    raw_df = custom_df

                points = coerce_custom_points(raw_df)
                clusters = build_cluster_summary(points, gap=int(custom_cluster_gap))

                run_dir = Path(tempfile.mkdtemp(prefix="custom_lollipop_streamlit_"))
                plots_dir = run_dir / "plots"
                plots_dir.mkdir(parents=True, exist_ok=True)

                prefix = make_output_prefix(custom_prefix or custom_target)
                main_png = plots_dir / f"{prefix}_lollipop_multisource_hotspots.png"
                main_pdf = plots_dir / f"{prefix}_lollipop_multisource_hotspots.pdf"
                matrix_png = plots_dir / f"{prefix}_cluster_evidence_matrix.png"
                matrix_pdf = plots_dir / f"{prefix}_cluster_evidence_matrix.pdf"

                plot_lollipop_multisource_hotspots(
                    points=points,
                    cluster_summary=clusters,
                    protein_len=int(custom_length),
                    out_png=main_png,
                    out_pdf=main_pdf,
                    shade_min_dataset_support=int(custom_shade_support),
                    x_tick_gap=int(custom_x_tick_gap),
                    title=f"{custom_display}: cleavage-site evidence",
                    x_label=f"{custom_target} residue coordinate",
                )
                plot_cluster_matrix(
                    cluster_summary=clusters,
                    out_png=matrix_png,
                    out_pdf=matrix_pdf,
                    title=f"{custom_target} cluster-level evidence matrix",
                )

                st.session_state["custom_result"] = build_result_payload(
                    prefix=prefix,
                    out_dir=run_dir,
                    points=points,
                    clusters=clusters,
                    main_png=main_png,
                    main_pdf=main_pdf,
                    matrix_png=matrix_png,
                    matrix_pdf=matrix_pdf,
                    metrics={
                        "Lollipop points": len(points),
                        "Clusters": len(clusters),
                        "Protein length": int(custom_length),
                    },
                )
            st.success("Custom plot generated successfully.")
        except Exception as exc:
            st.error(str(exc))

    render_download_panel("custom_result")
