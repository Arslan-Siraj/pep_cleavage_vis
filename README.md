# Integrated visualization of Arabidopsis PYK10 cleavage-site evidence

## Objective

This analysis visualizes cleavage-site and peptide-coverage evidence for **Arabidopsis AT3G09260.1**, also known as **LEB/BGLU23/PYK10/PSR3.1**, by integrating three Arabidopsis datasets and adding a separate exploratory fourth analysis against maize methyl-jasmonate N-terminomics.

## Input datasets

1. `Arabidopsis_semitryptome_23112020.xlsx`: semitryptic peptide evidence and candidate cleavage sites across Arabidopsis proteins.
2. `HUNTER_data_shoot_MJ.xlsx`: Arabidopsis shoot methyl-jasmonate N-terminomics; rows are N-terminally modified peptides with time-resolved log2 fold-changes.
3. `ProteinBased_www_P10'P10.xlsx`: protein-based metacaspase substrate matrix with peptide windows, peptide starts, and MC1/MC2/MC4/MC9 scoring blocks.
4. `SHMJ007_maize_shoot_rep_exclusion.xlsx`: maize shoot methyl-jasmonate N-terminomics used separately for Task 5.

## Main computational strategy

The code extracts all rows corresponding to PYK10/AT3G09260.1 from the three Arabidopsis datasets. Each row is converted to a common coordinate system using the peptide start or cleavage-site coordinate. The final plot contains three evidence tracks:

- **Semitryptome track:** peptide coverage segments; line thickness increases with peptide-spectrum-match support.
- **HUNTER MJ track:** N-terminal cleavage/start sites; marker colour and height encode the mean log2 fold-change across MJ time points.
- **Metacaspase matrix track:** peptide segments from the metacaspase substrate matrix; opacity encodes the maximum available substrate score.

access the GUI at: (StreamlitApp)[https://arslan-siraj-pep-cleavage--gui-appstreamlit-lollipop-app-ukwaqf.streamlit.app/]