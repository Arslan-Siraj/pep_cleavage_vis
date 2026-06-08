# Integrated visualization of Arabidopsis PYK10 cleavage-site evidence

## Objective

This analysis visualizes cleavage-site and peptide-coverage evidence for **Arabidopsis AT3G09260.1**, also known as **LEB/BGLU23/PYK10/PSR3.1**, by integrating three Arabidopsis datasets and adding a separate exploratory fourth analysis against maize methyl-jasmonate N-terminomics.

## Input datasets

1. `Arabidopsis_semitryptome_23112020(2).xlsx`: semitryptic peptide evidence and candidate cleavage sites across Arabidopsis proteins.
2. `HUNTER_data_shoot_MJ(1).xlsx`: Arabidopsis shoot methyl-jasmonate N-terminomics; rows are N-terminally modified peptides with time-resolved log2 fold-changes.
3. `ProteinBased_www_P10'P10 1(1).xlsx`: protein-based metacaspase substrate matrix with peptide windows, peptide starts, and MC1/MC2/MC4/MC9 scoring blocks.
4. `SHMJ007_maize_shoot_rep_exclusion(1).xlsx`: maize shoot methyl-jasmonate N-terminomics used separately for Task 5.

## Main computational strategy

The code extracts all rows corresponding to PYK10/AT3G09260.1 from the three Arabidopsis datasets. Each row is converted to a common coordinate system using the peptide start or cleavage-site coordinate. The final plot contains three evidence tracks:

- **Semitryptome track:** peptide coverage segments; line thickness increases with peptide-spectrum-match support.
- **HUNTER MJ track:** N-terminal cleavage/start sites; marker colour and height encode the mean log2 fold-change across MJ time points.
- **Metacaspase matrix track:** peptide segments from the metacaspase substrate matrix; opacity encodes the maximum available substrate score.

The inferred PYK10 length used for plotting is approximately **531 amino acids**. The final visualization also marks the predicted signal-peptide cleavage region near residue 25, which is present in the semitryptome/HUNTER annotations.

## Key extracted counts

| Dataset | PYK10 rows extracted |
|---|---:|
| Arabidopsis semitryptome | 171 |
| Arabidopsis HUNTER shoot MJ significant peptides | 26 |
| Metacaspase substrate matrix | 70 |

## Biological interpretation

The integrated plot is designed to distinguish three related but non-identical signals. The semitryptome provides broad peptide-level evidence for where semi-tryptic peptide starts occur. The HUNTER dataset provides treatment-responsive N-terminal peptides and therefore highlights cleavage/start sites that change after methyl jasmonate treatment. The metacaspase matrix provides a separate substrate-oriented perspective and helps identify whether any peptide starts overlap with evidence from MC-related experiments.

For PYK10, the visualization typically reveals a strong N-terminal/signal-peptide region near residue 25 and multiple internal regions with N-terminal peptide evidence, including MJ-responsive internal sites. Overlap or close proximity between HUNTER sites and semitryptome/metacaspase peptide windows is more informative than exact identity alone, because different enrichment, digestion, and scoring pipelines can report nearby peptide starts for the same underlying proteolytic region.

## Task 5: maize comparison

The maize comparison is provided as a separate exploratory screen rather than as a definitive orthology analysis. Because the supplied maize table does not by itself provide explicit Arabidopsis–maize orthology mappings, the script compares local cleavage windows directly. It reports Arabidopsis and maize peptide-start windows with high local amino-acid identity and strong maize MJ response. This is a reasonable first-pass approach for deeply conserved proteins, but it should be followed by explicit orthologue assignment and sequence alignment before making biological claims about evolutionary conservation.

## Output files

- `plots/final_pyk10_integrated_coverage.png`: final integrated peptide coverage and cleavage-site plot.
- `plots/final_pyk10_integrated_coverage.pdf`: vector/PDF version of the final plot.
- `plots/intermediate_site_distribution.png`: intermediate QC plot showing cleavage-site density by dataset.
- `plots/intermediate_hunter_logfc_heatmap.png`: intermediate HUNTER time-course heatmap for PYK10 peptides.
- `plots/task5_conserved_candidates.png`: exploratory Arabidopsis–maize conserved-window candidates.
- `tables/pyk10_semitryptome.csv`: extracted semitryptome rows for PYK10.
- `tables/pyk10_hunter_mj.csv`: extracted HUNTER rows for PYK10.
- `tables/pyk10_metacaspase.csv`: extracted metacaspase-matrix rows for PYK10.
- `tables/task5_conserved_candidates.csv`: candidate conserved cleavage windows between Arabidopsis and maize.

## Reproducibility

Run the workflow from the project directory with:

```bash
python scripts/run_analysis.py --input-dir /path/to/input_xlsx_files --output-dir /path/to/output_folder
```

The script does not require manual row selection. Protein identifiers, gene symbols, peptide coordinates, cleavage windows, and quantitative columns are parsed programmatically.

## Separate Task 5 implementation

In addition to the PYK10-focused `task5_conserved_candidates.csv`, I included a broader cross-species screen in `scripts/task5_conserved_screen.py`. This compares the Arabidopsis and maize HUNTER `All peptides` sheets across all proteins using local cleavage-window and peptide-sequence similarity. The ranked output is `tables/task5_conserved_hunter_cross_species_candidates.csv`, and the corresponding plot is `plots/task5_conserved_hunter_cross_species_candidates.png`. A dedicated explanation is provided in `TASK5_conserved_cleavage_sites.md`.
