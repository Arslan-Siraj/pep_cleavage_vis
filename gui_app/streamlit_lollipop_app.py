"""
Streamlit GUI for Arabidopsis multisource cleavage-site plotting and maize
local-context conservation screening.
The app provides two analysis pages:
  1. Arabidopsis multisource cleavage-site plot from the three Arabidopsis Excel datasets.
  2. Maize local-context conservation screen comparing Arabidopsis and maize cleavage-site windows.
"""

from __future__ import annotations

import io
import re
import tempfile
import zipfile
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

from run_analysis import (
    build_cluster_summary,
    build_exact_site_summary,
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
    page_title="Protein cleavage-site plotter",
    page_icon="🧬",
    layout="wide",
)

DATASET_ORDER = ["Semitryptome", "HUNTER MJ", "Metacaspase matrix"]


# -----------------------------------------------------------------------------
# Shared helper functions
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




# -----------------------------------------------------------------------------
# Maize local cleavage-window conservation helper functions
# -----------------------------------------------------------------------------

def read_fasta(path: str | Path) -> Dict[str, str]:
    """Minimal FASTA parser with tolerant Arabidopsis ID aliasing."""
    seqs: Dict[str, str] = {}
    current_ids: list[str] = []
    chunks: list[str] = []

    def add_aliases(ids: list[str], header: str) -> list[str]:
        out: list[str] = []
        for value in ids:
            value = str(value).strip()
            if not value:
                continue
            out.append(value)
            m = re.match(r"^(AT[1-5CM]G\d{5})\.\d+$", value, flags=re.I)
            if m:
                out.append(m.group(1))
        for token in re.findall(r"AT[1-5CM]G\d{5}(?:\.\d+)?", header, flags=re.I):
            out.append(token)
            m = re.match(r"^(AT[1-5CM]G\d{5})\.\d+$", token, flags=re.I)
            if m:
                out.append(m.group(1))
        sym_match = re.search(r"Symbols:\s*([^|]+)", header, flags=re.I)
        if sym_match:
            for sym in re.split(r"[,;\s]+", sym_match.group(1)):
                sym = sym.strip()
                if sym:
                    out.append(sym)
        return list(dict.fromkeys(out))

    def flush() -> None:
        nonlocal chunks, current_ids
        if current_ids and chunks:
            seq = "".join(chunks).upper().replace("*", "")
            seq = re.sub(r"[^A-Z]", "", seq)
            for key in current_ids:
                if key and key not in seqs:
                    seqs[key] = seq
                key_upper = key.upper()
                if key_upper and key_upper not in seqs:
                    seqs[key_upper] = seq
        current_ids = []
        chunks = []

    with open(path, "r", encoding="utf-8-sig") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header = line[1:].strip()
                first = header.split()[0]
                ids = [first]
                if "|" in first:
                    ids.extend([x.strip() for x in first.split("|") if x.strip()])
                current_ids = add_aliases(ids, header)
            else:
                chunks.append(re.sub(r"[^A-Za-z*]", "", line))
        flush()
    return seqs


def choose_col(df: pd.DataFrame, candidates: Iterable[str]) -> Optional[str]:
    lower = {str(c).lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    return None


def normalize_site(value) -> Optional[int]:
    try:
        if pd.isna(value):
            return None
        return int(round(float(value)))
    except Exception:
        return None


def reconstruct_window(seq: str, cleavage_site: int, window: int) -> tuple[str, str, str]:
    """
    Reconstruct cleavage-centered window from protein sequence.

    cleavage_site is treated as the first residue of the observed neo-N terminus,
    using 1-based coordinates. The cleavage bond is therefore between
    cleavage_site - 1 and cleavage_site.
    """
    if not seq or cleavage_site is None:
        return "", "", ""
    idx = int(cleavage_site) - 1
    if idx < 0 or idx >= len(seq):
        return "", "", ""
    before = seq[max(0, idx - window):idx]
    after = seq[idx:min(len(seq), idx + window)]
    return before, after, before + after


def split_maize_window(row: pd.Series, window: int) -> tuple[str, str, str]:
    """Use maize cleavage_sequence_10 where possible; otherwise combine before/after columns."""
    raw = row.get("cleavage_sequence_10", "")
    raw = "" if pd.isna(raw) else re.sub(r"[^A-Za-z]", "", str(raw)).upper()
    if len(raw) >= 2:
        if len(raw) >= 2 * window:
            before = raw[:window]
            after = raw[window:window * 2]
        else:
            mid = len(raw) // 2
            before = raw[:mid]
            after = raw[mid:]
        return before, after, before + after

    before = row.get("ten_residues_before", row.get("five_residues_before", ""))
    after = row.get("ten_residues_after", row.get("five_residues_after", ""))
    before = "" if pd.isna(before) else re.sub(r"[^A-Za-z]", "", str(before)).upper()
    after = "" if pd.isna(after) else re.sub(r"[^A-Za-z]", "", str(after)).upper()
    return before[-window:], after[:window], before[-window:] + after[:window]


def sequence_identity(a: str, b: str) -> float:
    """Position-wise identity over the shared length."""
    a = re.sub(r"[^A-Za-z]", "", str(a)).upper()
    b = re.sub(r"[^A-Za-z]", "", str(b)).upper()
    n = min(len(a), len(b))
    if n == 0:
        return 0.0
    return sum(aa == bb for aa, bb in zip(a[:n], b[:n])) / n


def centered_core(before: str, after: str, core: int) -> str:
    """Return residues nearest the cleavage bond: core before + core after."""
    before = str(before).upper()
    after = str(after).upper()
    return before[-core:] + after[:core]


def load_arabidopsis_points_from_df(df: pd.DataFrame, fasta: Dict[str, str], protein_query: str, window: int) -> pd.DataFrame:
    site_col = choose_col(df, ["cleavage_site", "first_aa_position", "site", "position"])
    dataset_col = choose_col(df, ["dataset", "source"])
    peptide_col = choose_col(df, ["peptide", "sequence"])
    window_col = choose_col(df, ["cleavage_sequence_10", "window", "local_window", "sequence_window"])

    if site_col is None:
        raise ValueError("Arabidopsis points table needs a cleavage-site column, e.g. cleavage_site or first_aa_position.")

    protein_col = choose_col(df, ["protein", "protein_id", "target", "target_id", "accession"])
    if protein_col is not None:
        mask = df[protein_col].astype(str).str.contains(re.escape(protein_query), case=False, na=False)
        if mask.any():
            df = df[mask].copy()

    seq = fasta.get(protein_query) or fasta.get(protein_query.upper())
    locus = re.sub(r"\.\d+$", "", protein_query, flags=re.I)
    if seq is None and locus != protein_query:
        seq = fasta.get(locus) or fasta.get(locus.upper())

    rows = []
    skipped = []
    seq_len = len(seq) if seq else None

    for i, row in df.iterrows():
        site = normalize_site(row[site_col])
        if site is None:
            skipped.append({"row": int(i), "reason": "missing/non-numeric cleavage site", "site": row.get(site_col, "")})
            continue

        before = after = combined = ""
        if window_col is not None and not pd.isna(row.get(window_col)):
            raw = re.sub(r"[^A-Za-z]", "", str(row[window_col])).upper()
            if len(raw) >= 2:
                if len(raw) >= 2 * window:
                    before, after = raw[:window], raw[window:window * 2]
                else:
                    mid = len(raw) // 2
                    before, after = raw[:mid], raw[mid:]
                combined = before + after

        if not combined and seq:
            before, after, combined = reconstruct_window(seq, site, window)

        if not combined:
            if seq:
                reason = f"site {site} outside FASTA sequence length {len(seq)}"
            else:
                available = ", ".join(list(fasta.keys())[:12]) if fasta else "no FASTA records loaded"
                reason = f"protein {protein_query!r} not found in FASTA; available IDs include: {available}"
            skipped.append({"row": int(i), "reason": reason, "site": site})
            continue

        rows.append({
            "species": "Arabidopsis",
            "protein": protein_query,
            "dataset": row.get(dataset_col, "Arabidopsis") if dataset_col else "Arabidopsis",
            "cleavage_site": site,
            "peptide": row.get(peptide_col, "") if peptide_col else "",
            "before": before,
            "after": after,
            "window_sequence": combined,
        })

    if not rows:
        site_values = pd.to_numeric(df[site_col], errors="coerce").dropna().astype(int).head(20).tolist()
        raise ValueError(
            "Could not build any Arabidopsis sequence windows.\n"
            f"Protein requested: {protein_query}\n"
            f"FASTA sequence length found: {seq_len}\n"
            f"First cleavage-site values in points file: {site_values}\n"
            f"Examples of skipped rows: {skipped[:5]}\n\n"
            "Most likely cause: the Arabidopsis points CSV is not for this protein, or the cleavage coordinates "
            "are outside the FASTA sequence length. Re-run the Arabidopsis multisource plotting workflow for this protein, "
            "then use that lollipop-points CSV in the conservation screen."
        )

    if skipped:
        st.warning(f"Skipped {len(skipped)} Arabidopsis point rows that could not be windowed. First examples: {skipped[:3]}")

    return pd.DataFrame(rows).drop_duplicates()


def load_maize_windows(path: str | Path, maize_query: str, window: int) -> pd.DataFrame:
    xl = pd.ExcelFile(path)
    frames = []
    for sheet in xl.sheet_names:
        df = pd.read_excel(path, sheet_name=sheet)
        if not {"gene", "protein"}.issubset(set(df.columns)):
            continue
        q = str(maize_query)
        mask = (
            df["gene"].astype(str).str.contains(re.escape(q), case=False, na=False)
            | df["protein"].astype(str).str.contains(re.escape(q), case=False, na=False)
        )
        if mask.any():
            temp = df[mask].copy()
            temp["source_sheet"] = sheet
            frames.append(temp)
    if not frames:
        raise ValueError(f"No maize rows matched query: {maize_query}")

    df = pd.concat(frames, ignore_index=True)
    rows = []
    for _, row in df.iterrows():
        site = normalize_site(row.get("first_aa_position"))
        preceding = normalize_site(row.get("preceding_aa_position"))
        before, after, combined = split_maize_window(row, window)
        rows.append({
            "species": "Maize",
            "protein": row.get("protein", ""),
            "gene": row.get("gene", ""),
            "source_sheet": row.get("source_sheet", ""),
            "cleavage_site": site,
            "preceding_site": preceding,
            "peptide": row.get("peptide", ""),
            "before": before,
            "after": after,
            "window_sequence": combined,
            "logFC": row.get("hunter_limma_logFC", pd.NA),
            "adj_p": row.get("hunter_limma_adj_p_val", pd.NA),
            "regulation": row.get("hunter_limma_regulation", ""),
        })
    out = pd.DataFrame(rows)
    return out.dropna(subset=["cleavage_site"]).drop_duplicates()


def compare_windows(ath: pd.DataFrame, maize: pd.DataFrame, core: int, min_core_identity: float, min_window_identity: float) -> pd.DataFrame:
    hits = []
    for _, a in ath.iterrows():
        a_core = centered_core(a["before"], a["after"], core)
        for _, m in maize.iterrows():
            m_core = centered_core(m["before"], m["after"], core)
            core_id = sequence_identity(a_core, m_core)
            window_id = sequence_identity(a["window_sequence"], m["window_sequence"])
            if core_id >= min_core_identity and window_id >= min_window_identity:
                hits.append({
                    "arabidopsis_protein": a["protein"],
                    "arabidopsis_dataset": a["dataset"],
                    "arabidopsis_site": a["cleavage_site"],
                    "arabidopsis_peptide": a.get("peptide", ""),
                    "arabidopsis_context": f'{a["before"]}|{a["after"]}',
                    "maize_gene": m.get("gene", ""),
                    "maize_protein": m.get("protein", ""),
                    "maize_site": m["cleavage_site"],
                    "maize_peptide": m.get("peptide", ""),
                    "maize_context": f'{m["before"]}|{m["after"]}',
                    "core_identity": round(core_id, 3),
                    "window_identity": round(window_id, 3),
                    "arabidopsis_core": a_core,
                    "maize_core": m_core,
                    "maize_logFC": m.get("logFC", pd.NA),
                    "maize_adj_p": m.get("adj_p", pd.NA),
                    "maize_regulation": m.get("regulation", ""),
                })
    if not hits:
        return pd.DataFrame(columns=[
            "arabidopsis_protein", "arabidopsis_dataset", "arabidopsis_site", "arabidopsis_peptide",
            "arabidopsis_context", "maize_gene", "maize_protein", "maize_site", "maize_peptide",
            "maize_context", "core_identity", "window_identity", "arabidopsis_core", "maize_core",
            "maize_logFC", "maize_adj_p", "maize_regulation"
        ])
    return pd.DataFrame(hits).sort_values(["core_identity", "window_identity"], ascending=False)


def make_context_figure(ath: pd.DataFrame, maize: pd.DataFrame, hits: pd.DataFrame, title: str):
    import plotly.graph_objects as go

    fig = go.Figure()
    y_map = {"Arabidopsis": 1, "Maize": 0}
    for df, species in [(ath, "Arabidopsis"), (maize, "Maize")]:
        fig.add_trace(go.Scatter(
            x=df["cleavage_site"],
            y=[y_map[species]] * len(df),
            mode="markers",
            name=species,
            marker={"size": 14},
            customdata=df[["peptide", "window_sequence"]].fillna("").to_numpy(),
            hovertemplate=(
                f"<b>{species}</b><br>Site=%{{x}}<br>Peptide=%{{customdata[0]}}"
                "<br>Window=%{customdata[1]}<extra></extra>"
            ),
        ))

    for _, h in hits.head(50).iterrows():
        fig.add_shape(
            type="line",
            x0=h["arabidopsis_site"], y0=y_map["Arabidopsis"],
            x1=h["maize_site"], y1=y_map["Maize"],
            line={"width": 1, "dash": "dot"},
        )

    fig.update_layout(
        title=title,
        xaxis_title="Residue coordinate within each species protein",
        yaxis={"tickmode": "array", "tickvals": [0, 1], "ticktext": ["Maize", "Arabidopsis"]},
        hovermode="closest",
        template="plotly_white",
        height=440,
    )
    return fig


def build_task5_payload(prefix: str, ath: pd.DataFrame, maize: pd.DataFrame, hits: pd.DataFrame, fig_html: str) -> Dict[str, object]:
    files = {
        f"tables/{prefix}_arabidopsis_windows.csv": csv_bytes(ath),
        f"tables/{prefix}_maize_windows.csv": csv_bytes(maize),
        f"tables/{prefix}_candidate_local_context_matches.csv": csv_bytes(hits),
        f"plots/{prefix}_paired_context_plot.html": fig_html.encode("utf-8"),
    }
    return {
        "ath": ath,
        "maize": maize,
        "hits": hits,
        "fig_html": fig_html,
        "ath_csv": csv_bytes(ath),
        "maize_csv": csv_bytes(maize),
        "hits_csv": csv_bytes(hits),
        "html_bytes": fig_html.encode("utf-8"),
        "zip_bytes": make_download_zip(files),
        "prefix": prefix,
    }


def render_task5_panel(result_key: str) -> None:
    result = st.session_state.get(result_key)
    if not result:
        return

    st.divider()
    st.subheader("Conservation-screen output")
    c1, c2, c3 = st.columns(3)
    c1.metric("Arabidopsis windows", len(result["ath"]))
    c2.metric("Maize windows", len(result["maize"]))
    c3.metric("Candidate matches", len(result["hits"]))

    st.markdown("#### Candidate local-context matches")
    if result["hits"].empty:
        st.info("No candidate local-context matches passed the selected thresholds. Try lowering the identity cutoffs or testing another deeply conserved protein.")
    else:
        st.dataframe(result["hits"], use_container_width=True)

    st.markdown("#### Paired local-context plot")
    st.caption("Coordinates remain species-specific. Dotted links indicate local sequence-context matches, not raw coordinate equivalence.")
    components.html(result["fig_html"], height=500, scrolling=True)

    d1, d2, d3, d4 = st.columns(4)
    prefix = result["prefix"]
    with d1:
        st.download_button("Download candidate matches", result["hits_csv"], f"{prefix}_candidate_local_context_matches.csv", "text/csv", use_container_width=True)
    with d2:
        st.download_button("Download Arabidopsis windows", result["ath_csv"], f"{prefix}_arabidopsis_windows.csv", "text/csv", use_container_width=True)
    with d3:
        st.download_button("Download plot HTML", result["html_bytes"], f"{prefix}_paired_context_plot.html", "text/html", use_container_width=True)
    with d4:
        st.download_button("Download all ZIP", result["zip_bytes"], f"{prefix}_conservation_output.zip", "application/zip", use_container_width=True)

    with st.expander("Preview maize windows", expanded=False):
        st.dataframe(result["maize"], use_container_width=True)
    with st.expander("Preview Arabidopsis windows", expanded=False):
        st.dataframe(result["ath"], use_container_width=True)


# -----------------------------------------------------------------------------
# Header and workflow selection
# -----------------------------------------------------------------------------

st.title("Protein cleavage-site plotter")
st.write(
    "This interface first builds a reproducible Arabidopsis multisource cleavage-site plot from the "
    "three Arabidopsis datasets. A separate conservation screen can then compare Arabidopsis cleavage-site "
    "sequence windows with maize cleavage-site windows from the maize HUNTER dataset."
)

workflow = st.sidebar.radio(
    "Analysis page",
    [
        "Arabidopsis multisource cleavage-site plot",
        "Maize local-context conservation screen",
    ],
)

st.sidebar.markdown("---")
st.sidebar.caption(
    "Use the Arabidopsis page to extract and plot cleavage evidence. Use the maize page for candidate "
    "cross-species conservation based on local amino-acid sequence context."
)


# -----------------------------------------------------------------------------
# Arabidopsis multisource source Excel extraction
# -----------------------------------------------------------------------------

if workflow == "Arabidopsis multisource cleavage-site plot":
    st.header("Arabidopsis multisource cleavage-site plot")
    st.markdown(
        "This page extracts one Arabidopsis protein from the three source Excel datasets and places all "
        "cleavage or neo-N-terminal peptide evidence on one protein-coordinate axis. The resulting plot is "
        "the main Arabidopsis evidence summary and also produces the lollipop-points table that can be reused "
        "in the maize conservation screen."
    )
    with st.expander("What this page needs and produces", expanded=True):
        st.markdown(
            """
            **Inputs**
            - Arabidopsis semitryptome Excel file.
            - Arabidopsis HUNTER shoot methyl-jasmonate Excel file.
            - Protein-based P10/P10 metacaspase Excel file.
            - A target Arabidopsis protein identifier, with UniProt ID and aliases.

            **Outputs**
            - A lollipop plot with one track per Arabidopsis dataset.
            - A cluster evidence matrix summarizing nearby cleavage regions.
            - CSV tables for extracted source rows, combined lollipop points, clusters, and exact residue-level overlaps.

            **Interpretation**
            Dot position is the residue coordinate. Dot size is scaled within each dataset: PSM for semitryptome, maximum absolute logFC for HUNTER, and maximum metacaspase score for the matrix. Shaded regions mark nearby-site clusters supported by multiple datasets.
            """
        )

    with st.form("excel_form"):
        st.subheader("1. Upload input files")
        c1, c2, c3 = st.columns(3)
        with c1:
            semi_file = st.file_uploader(
                "Arabidopsis semitryptome Excel",
                type=["xlsx"],
                key="semi_file",
                help="Excel table containing semitryptome cleavage evidence. Matching rows are identified using the target protein ID, UniProt accession, gene symbol, or aliases.",
            )
        with c2:
            hunter_file = st.file_uploader(
                "HUNTER shoot MJ Excel",
                type=["xlsx"],
                key="hunter_file",
                help="Arabidopsis methyl-jasmonate HUNTER dataset. The app uses peptide start coordinates and summarizes the strongest absolute logFC across time points.",
            )
        with c3:
            meta_file = st.file_uploader(
                "ProteinBased P10/P10 metacaspase Excel",
                type=["xlsx"],
                key="meta_file",
                help="Metacaspase matrix. Rows are matched strictly by the primary Protein column to avoid pulling homologous proteins by alias.",
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
                help="UniProt accession for the same protein, e.g. Q9SR37.",
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
                help="Nearby cleavage coordinates are assigned to the same processing region when adjacent sites are separated by no more than this number of amino acids.",
            )
        with o3:
            shade_support = st.number_input(
                "Supported by ≥ datasets",
                min_value=1,
                max_value=3,
                value=2,
                step=1,
                help="Background shading is added only for clusters that contain evidence from at least this many Arabidopsis datasets.",
            )
        with o4:
            x_tick_gap = st.number_input(
                "X-axis tick gap",
                min_value=5,
                value=20,
                step=5,
                help="Spacing between residue-coordinate tick marks on the x-axis.",
            )

        m1, m2 = st.columns(2)
        with m1:
            plot_every_row = st.checkbox(
                "Plot every extracted row instead of aggregating exact duplicate sites",
                value=False,
                help="When off, exact duplicate coordinates from the same dataset are collapsed into one dot to avoid overplotting. When on, every extracted row is plotted.",
            )
        with m2:
            st.info(
                "Metacaspase extraction is strict: only rows whose primary Protein "
                "column equals the Primary target ID are retained."
            )
            include_isoforms = False

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
                    meta = read_metacaspase(meta_path, target, match_isoforms=False)

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
                    exact_sites = build_exact_site_summary(points)
                    exact_three_dataset_sites = exact_sites[
                        exact_sites["is_exact_three_dataset_site"] == True
                    ].copy()
                    plot_len = int(protein_length) if protein_length else infer_protein_length(semi, hunter, meta)

                    main_png = plots_dir / f"{prefix}_lollipop_multisource_hotspots.png"
                    main_pdf = plots_dir / f"{prefix}_lollipop_multisource_hotspots.pdf"
                    matrix_png = plots_dir / f"{prefix}_cluster_evidence_matrix.png"
                    matrix_pdf = plots_dir / f"{prefix}_cluster_evidence_matrix.pdf"

                    plot_lollipop_multisource_hotspots(
                        points=points,
                        cluster_summary=clusters,
                        exact_site_summary=exact_sites,
                        highlight_exact_dataset_support=3,
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
                            f"{prefix}_exact_site_summary.csv": exact_sites,
                            f"{prefix}_exact_three_dataset_sites.csv": exact_three_dataset_sites,
                        },
                        metrics={
                            "Semitryptome rows": len(semi),
                            "HUNTER rows": len(hunter),
                            "Metacaspase rows": len(meta),
                            "Lollipop points": len(points),
                            "Clusters": len(clusters),
                            "Exact 3-dataset sites": len(exact_three_dataset_sites),
                            "Protein length": plot_len,
                        },
                    )
                st.success("Plot generated successfully.")
            except Exception as exc:
                st.error(str(exc))
                st.exception(exc)

    render_download_panel("excel_result")


# -----------------------------------------------------------------------------
# Maize local-context conservation screen
# -----------------------------------------------------------------------------

else:
    st.header("Maize local-context conservation screen")
    st.write(
        "This workflow compares Arabidopsis and maize cleavage sites by local amino-acid context, "
        "not by raw residue coordinate. It is intended for conserved proteins such as PsbQ1, ribosomal proteins, "
        "or histone/nucleosome-associated proteins."
    )

    st.info(
        "Recommended use: first generate an Arabidopsis plot for a conserved Arabidopsis protein, for example AT4G21280.1 / PsbQ1. "
        "Then run this conservation screen using the Arabidopsis lollipop-points table, the Arabidopsis protein FASTA, and the maize Excel file."
    )

    with st.expander("What this conservation screen does", expanded=True):
        st.markdown(
            """
            This analysis asks whether Arabidopsis and maize cleavage evidence occurs in similar local amino-acid sequence environments. It deliberately avoids direct raw-coordinate matching, because an Arabidopsis residue number and a maize residue number are not automatically equivalent.

            The app reconstructs an Arabidopsis cleavage window from the protein FASTA and compares it with the maize `cleavage_sequence_10` window. Candidate matches are retained when both the cleavage-proximal core and the full local window pass the selected identity thresholds.
            """
        )

    with st.expander("Required inputs and why they are needed", expanded=True):
        st.markdown(
            """
            - **Arabidopsis lollipop-points CSV**: contains the Arabidopsis cleavage-site coordinates extracted from the three Arabidopsis datasets. You can either use the current app-session result or upload the `<prefix>_lollipop_points.csv` file.
            - **Arabidopsis protein FASTA**: contains the amino-acid sequence for the same Arabidopsis protein. This is used to extract residues before and after each cleavage site.
            - **Maize HUNTER Excel file**: contains maize cleavage evidence and the local `cleavage_sequence_10` windows.
            - **Arabidopsis and maize identifiers**: specify which Arabidopsis protein and maize gene/protein should be compared, for example `AT4G21280.1` and `PSBQ1`.
            """
        )

    with st.expander("How to interpret the output", expanded=False):
        st.markdown(
            """
            - **Arabidopsis windows** are cleavage-centered sequence contexts reconstructed from the uploaded FASTA.
            - **Maize windows** are cleavage-centered sequence contexts read from the maize Excel file.
            - **Candidate matches** are Arabidopsis–maize window pairs that pass both similarity cutoffs.
            - **Dotted lines in the plot** connect candidate local-context matches. They do not mean the residue coordinates are identical between species.
            """
        )

    with st.expander("How to download the Arabidopsis protein FASTA", expanded=False):
        st.markdown(
            """
            1. Open the TAIR sequence bulk download page.
            2. Select **Araport11 proteins** or **Araport11 peptide/protein sequences** as the dataset. Do not select transcripts, because transcripts are nucleotide sequences.
            3. Paste the Arabidopsis protein ID, for example `AT4G21280.1`.
            4. Select **FASTA** output.
            5. Download and save the file as something like `AT4G21280_protein.fasta`.

            The FASTA header should contain the target ID, for example:

            `>AT4G21280.1 | Symbols: PSBQ-1, PSBQA, PSBQ | ...`
            """
        )

    points_source = st.radio(
        "Arabidopsis cleavage-site input",
        ["Use Arabidopsis plot result from this app session", "Upload lollipop-points CSV"],
        horizontal=True,
        help="Use the app-session result if you just generated the Arabidopsis plot above. Upload the CSV if you generated the points table earlier or outside this app.",
    )

    with st.form("task5_form"):
        st.subheader("1. Upload or select required inputs")
        c1, c2, c3 = st.columns(3)
        with c1:
            if points_source == "Upload lollipop-points CSV":
                ath_points_file = st.file_uploader(
                    "Arabidopsis lollipop-points CSV",
                    type=["csv"],
                    help="This is the Arabidopsis plot output table named <prefix>_lollipop_points.csv.",
                    key="task5_ath_points_file",
                )
            else:
                ath_points_file = None
                st.caption("The app will use the current `excel_result` lollipop points table if available.")
        with c2:
            ath_fasta_file = st.file_uploader(
                "Arabidopsis protein FASTA",
                type=["fasta", "fa", "faa", "txt"],
                help="Required unless the points CSV already contains a cleavage_sequence_10/window column.",
                key="task5_fasta_file",
            )
        with c3:
            maize_xlsx_file = st.file_uploader(
                "Maize HUNTER Excel",
                type=["xlsx"],
                help="Example: SHMJ007_maize_shoot_rep_exclusion.xlsx.",
                key="task5_maize_file",
            )

        st.subheader("2. Define Arabidopsis and maize proteins")
        p1, p2, p3 = st.columns(3)
        with p1:
            ath_protein = st.text_input("Arabidopsis protein ID", value="AT4G21280.1", help="Protein identifier present in the FASTA header and corresponding to the Arabidopsis lollipop-points table.")
        with p2:
            maize_query = st.text_input("Maize gene/protein query", value="PSBQ1", help="Gene symbol or protein accession used to find matching rows in the maize Excel file.")
        with p3:
            task5_prefix = st.text_input("Output prefix", value="psbq1_local_context", help="Prefix used for the downloadable CSV, HTML plot, and ZIP files.")

        st.subheader("3. Similarity thresholds")
        s1, s2, s3, s4 = st.columns(4)
        with s1:
            window = st.number_input("Window size each side", min_value=3, max_value=30, value=10, step=1, help="Number of amino acids taken before and after each cleavage position. A value of 10 creates a 20-aa context window.")
        with s2:
            core = st.number_input("Core size each side", min_value=1, max_value=15, value=5, step=1, help="Number of amino acids immediately next to the cleavage bond used for the stricter core similarity score.")
        with s3:
            min_core_identity = st.number_input("Minimum core identity", min_value=0.0, max_value=1.0, value=0.40, step=0.05, help="Minimum fraction of identical amino acids in the cleavage-proximal core. For a 10-aa core, 0.40 means about 4 matching positions.")
        with s4:
            min_window_identity = st.number_input("Minimum full-window identity", min_value=0.0, max_value=1.0, value=0.25, step=0.05, help="Minimum fraction of identical amino acids across the complete before+after window. For a 20-aa window, 0.25 means about 5 matching positions.")

        st.markdown("**Threshold interpretation:** lower cutoffs return more exploratory candidates; higher cutoffs return fewer but more sequence-similar candidates. The default values are intended for screening deeply conserved proteins, not for proving conservation alone.")

        task5_submit = st.form_submit_button("Run local-context comparison", use_container_width=True)

    if task5_submit:
        try:
            with st.spinner("Comparing Arabidopsis and maize cleavage-site windows..."):
                if points_source == "Use Arabidopsis plot result from this app session":
                    current = st.session_state.get("excel_result")
                    if not current or "points_df" not in current:
                        raise ValueError("No Arabidopsis plot result is available in this app session. Generate the Arabidopsis plot first, or choose 'Upload lollipop-points CSV'.")
                    ath_points = current["points_df"].copy()
                else:
                    if ath_points_file is None:
                        raise ValueError("Please upload the Arabidopsis lollipop-points CSV.")
                    ath_points = pd.read_csv(ath_points_file)

                if ath_fasta_file is None:
                    has_window_col = choose_col(ath_points, ["cleavage_sequence_10", "window", "local_window", "sequence_window"])
                    if has_window_col is None:
                        raise ValueError("Please upload the Arabidopsis protein FASTA, or provide a points CSV containing a local window column.")
                    fasta = {}
                    fasta_path = None
                else:
                    run_dir = Path(tempfile.mkdtemp(prefix="task5_local_context_"))
                    fasta_path = save_uploaded_file(ath_fasta_file, run_dir / ath_fasta_file.name)
                    fasta = read_fasta(fasta_path)

                if maize_xlsx_file is None:
                    raise ValueError("Please upload the maize HUNTER Excel file.")
                if 'run_dir' not in locals():
                    run_dir = Path(tempfile.mkdtemp(prefix="task5_local_context_"))
                maize_path = save_uploaded_file(maize_xlsx_file, run_dir / maize_xlsx_file.name)

                ath = load_arabidopsis_points_from_df(ath_points, fasta, ath_protein, int(window))
                maize = load_maize_windows(maize_path, maize_query, int(window))
                hits = compare_windows(ath, maize, int(core), float(min_core_identity), float(min_window_identity))

                fig = make_context_figure(
                    ath,
                    maize,
                    hits,
                    f"Local cleavage-context comparison: {ath_protein} vs {maize_query}",
                )
                fig_html = fig.to_html(include_plotlyjs="cdn", full_html=True)
                prefix = make_output_prefix(task5_prefix or f"{ath_protein}_{maize_query}_local_context")
                st.session_state["task5_result"] = build_task5_payload(prefix, ath, maize, hits, fig_html)

            st.success("Local-context comparison completed.")
        except Exception as exc:
            st.error(str(exc))
            st.exception(exc)

    render_task5_panel("task5_result")
