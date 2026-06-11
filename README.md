# Arabidopsis cleavage-site visualization

This Streamlit app visualizes cleavage-site evidence for a selected **Arabidopsis protein** using three Arabidopsis datasets, with an optional maize local-context comparison.

## Inputs

- `Arabidopsis_semitryptome_23112020.xlsx`
- `HUNTER_data_shoot_MJ.xlsx`
- `ProteinBased_www_P10'P10.xlsx`
- Optional: `SHMJ007_maize_shoot_rep_exclusion.xlsx` for the maize comparison
- Optional: Arabidopsis protein FASTA for local sequence-window comparison

## What the app does

The app searches for the selected target protein using the provided protein ID, UniProt ID, and aliases. It extracts cleavage or peptide-start evidence and plots the results on a common protein-coordinate lollipop plot.

Evidence tracks:

- **Semitryptome:** peptide-spectrum-match support
- **HUNTER MJ:** methyl-jasmonate-responsive N-terminal peptides
- **Metacaspase matrix:** predicted or scored metacaspase substrate evidence

The maize screen compares Arabidopsis and maize cleavage sites by local amino-acid sequence context, not by raw residue number.

## Web app

[Open the Streamlit app](https://arslan-siraj-pep-cleavage--gui-appstreamlit-lollipop-app-ukwaqf.streamlit.app/)

## Outputs

- Lollipop plot as PNG/PDF
- Cleavage-site table as CSV
- Cluster and exact-site summaries as CSV
- Optional maize local-context match table and HTML plot

## Output of workflow
[output folder](https://github.com/Arslan-Siraj/pep_cleavage_vis/output)