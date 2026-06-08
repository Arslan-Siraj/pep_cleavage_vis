# Task 5: exploratory Arabidopsis–maize conserved cleavage-site screen

## Purpose

Task 5 asks whether cleavage sites detected in Arabidopsis are evolutionarily conserved, or at least occur in reasonable proximity to conserved sites, in maize. The supplied files do not include a formal Arabidopsis–maize orthology table or aligned protein sequences. Therefore, the implemented solution performs a conservative **sequence-window based screen** rather than claiming definitive conservation.

## Input used

- Arabidopsis: `HUNTER_data_shoot_MJ(1).xlsx`, sheet `All peptides`
- Maize: `SHMJ007_maize_shoot_rep_exclusion(1).xlsx`, sheet `All peptides`

The script `scripts/task5_conserved_screen.py` reads both tables, extracts N-terminal peptides, cleavage/start coordinates, local cleavage windows, and MJ log2 fold-change values.

## Method

For each Arabidopsis and maize N-terminal peptide, the script builds two comparable sequence features:

1. the local cleavage window, `cleavage_sequence_10`, and
2. the detected N-terminal peptide sequence.

The algorithm then searches for cross-species pairs with high local sequence identity. To keep the comparison computationally tractable, it first indexes maize sequences by 6-mers, then evaluates only candidate pairs that share at least one local 6-mer. Candidate pairs are retained when either the peptide or the local cleavage window has a local identity of at least 0.65.

This detects candidate conserved processing regions in proteins such as photosystem and other deeply conserved proteins, but it remains a hypothesis-generating screen. A definitive result would require orthologue assignment and alignment of full-length proteins.

## Main outputs

- `tables/task5_conserved_hunter_cross_species_candidates.csv`: ranked Arabidopsis–maize candidate pairs.
- `plots/task5_conserved_hunter_cross_species_candidates.png`: bar plot of the top candidate pairs.

## Interpretation

The strongest candidates include highly conserved photosystem proteins, for example Arabidopsis chloroplast-encoded PSBB paired with maize psbB at a matching peptide/cleavage region. These cases are credible first-pass candidates because both species show highly similar peptide sequences and local cleavage context. However, entries where peptide identity is high but window identity is low should be interpreted more cautiously; they may reflect conserved peptide sequence detected at nearby but not identical local cleavage contexts.

## Recommended follow-up

For a manuscript-quality evolutionary analysis, the next step would be:

1. map Arabidopsis and maize proteins to orthologue groups,
2. retrieve full protein sequences,
3. perform pairwise or multiple sequence alignment,
4. project cleavage coordinates onto the alignment, and
5. quantify whether cleavage sites fall within a small aligned-distance threshold.

The current solution provides an automated, transparent first-pass screen and generates candidate pairs for that stricter downstream analysis.
