from __future__ import annotations
import argparse
import html
import re
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple
import pandas as pd


MAIZE_REQUIRED_COLUMNS = {
    "protein", "gene", "peptide", "first_aa_position", "preceding_aa_position", "cleavage_sequence_10"
}


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
            # Also index the locus without transcript suffix, e.g. AT4G21280 for AT4G21280.1.
            m = re.match(r"^(AT[1-5CM]G\d{5})\.\d+$", value, flags=re.I)
            if m:
                out.append(m.group(1))
        # Capture all ATxGxxxxx(.isoform) tokens anywhere in the header.
        for token in re.findall(r"AT[1-5CM]G\d{5}(?:\.\d+)?", header, flags=re.I):
            out.append(token)
            m = re.match(r"^(AT[1-5CM]G\d{5})\.\d+$", token, flags=re.I)
            if m:
                out.append(m.group(1))
        # Capture symbol aliases after "Symbols:" when present.
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
    lower = {c.lower(): c for c in df.columns}
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

    cleavage_site is treated as the first residue of the observed neo-N terminus, using 1-based coordinates.
    The cleavage bond is therefore between cleavage_site - 1 and cleavage_site.
    Returns before, after, combined before+after.
    """
    if not seq or cleavage_site is None:
        return "", "", ""
    idx = int(cleavage_site) - 1  # 0-based first residue after cleavage
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
        # SHMJ007 cleavage_sequence_10 is usually ten residues before + ten residues after.
        mid = min(window, len(raw) // 2)
        # For the expected 20 aa string this gives exactly 10 + 10.
        if len(raw) >= 2 * window:
            before = raw[:window]
            after = raw[window:window * 2]
        else:
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


def load_arabidopsis_points(path: str | Path, fasta: Dict[str, str], protein_query: str, window: int) -> pd.DataFrame:
    df = pd.read_csv(path)
    site_col = choose_col(df, ["cleavage_site", "first_aa_position", "site", "position"])
    dataset_col = choose_col(df, ["dataset", "source"])
    peptide_col = choose_col(df, ["peptide", "sequence"])
    window_col = choose_col(df, ["cleavage_sequence_10", "window", "local_window", "sequence_window"])

    if site_col is None:
        raise ValueError("Arabidopsis points table needs a cleavage-site column, e.g. cleavage_site or first_aa_position.")

    # Keep rows matching the requested protein when a suitable protein column exists.
    protein_col = choose_col(df, ["protein", "protein_id", "target", "target_id", "accession"])
    if protein_col is not None:
        mask = df[protein_col].astype(str).str.contains(re.escape(protein_query), case=False, na=False)
        if mask.any():
            df = df[mask].copy()

    # Tolerant FASTA lookup: exact, uppercase, and locus without isoform suffix.
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
            f"Input points file: {path}\n"
            f"Protein requested: {protein_query}\n"
            f"FASTA sequence length found: {seq_len}\n"
            f"First cleavage-site values in points file: {site_values}\n"
            f"Examples of skipped rows: {skipped[:5]}\n\n"
            "Most likely cause: the Arabidopsis points CSV is not for this protein, or the cleavage coordinates "
            "are outside the FASTA sequence length. Re-run the Arabidopsis lollipop workflow for AT4G21280.1 and "
            "use its output/tables/AT4G21280_lollipop_points.csv."
        )

    if skipped:
        print(f"WARNING: skipped {len(skipped)} Arabidopsis point rows that could not be windowed.")
        for item in skipped[:5]:
            print(f"  - row {item['row']}: {item['reason']}")

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


def make_plot(ath: pd.DataFrame, maize: pd.DataFrame, hits: pd.DataFrame, out_html: Path, title: str) -> None:
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

    # Draw simple visual links between matched local contexts. Coordinates remain species-specific; links are evidence links.
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
        height=420,
    )
    fig.write_html(out_html, include_plotlyjs="cdn")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--arabidopsis-points", required=True, help="Arabidopsis lollipop-points CSV for a conserved target protein.")
    ap.add_argument("--arabidopsis-fasta", default=None, help="Protein FASTA containing the Arabidopsis target. Required unless the points CSV already has windows.")
    ap.add_argument("--arabidopsis-protein", required=True, help="Arabidopsis target ID, e.g. AT4G21280.1.")
    ap.add_argument("--maize-xlsx", required=True, help="SHMJ007 maize Excel file.")
    ap.add_argument("--maize-query", required=True, help="Maize gene/protein query, e.g. PSBQ1 or Q41048.")
    ap.add_argument("--window", type=int, default=10, help="Residues before and after the cleavage site.")
    ap.add_argument("--core", type=int, default=5, help="Residues nearest the cleavage bond used for core similarity.")
    ap.add_argument("--min-core-identity", type=float, default=0.40)
    ap.add_argument("--min-window-identity", type=float, default=0.25)
    ap.add_argument("--out-prefix", default="point5_local_context")
    ap.add_argument("--plot", action="store_true")
    args = ap.parse_args()

    fasta = read_fasta(args.arabidopsis_fasta) if args.arabidopsis_fasta else {}
    ath = load_arabidopsis_points(args.arabidopsis_points, fasta, args.arabidopsis_protein, args.window)
    maize = load_maize_windows(args.maize_xlsx, args.maize_query, args.window)
    hits = compare_windows(ath, maize, args.core, args.min_core_identity, args.min_window_identity)

    prefix = Path(args.out_prefix)
    ath_out = prefix.with_name(prefix.name + "_arabidopsis_windows.csv")
    maize_out = prefix.with_name(prefix.name + "_maize_windows.csv")
    hits_out = prefix.with_name(prefix.name + "_candidate_local_context_matches.csv")

    ath.to_csv(ath_out, index=False)
    maize.to_csv(maize_out, index=False)
    hits.to_csv(hits_out, index=False)

    print(f"Arabidopsis windows: {len(ath)} -> {ath_out}")
    print(f"Maize windows: {len(maize)} -> {maize_out}")
    print(f"Candidate local-context matches: {len(hits)} -> {hits_out}")

    if args.plot:
        html_out = prefix.with_name(prefix.name + "_paired_context_plot.html")
        make_plot(ath, maize, hits, html_out, f"Local cleavage-context comparison: {args.arabidopsis_protein} vs {args.maize_query}")
        print(f"Plot: {html_out}")


if __name__ == "__main__":
    main()
