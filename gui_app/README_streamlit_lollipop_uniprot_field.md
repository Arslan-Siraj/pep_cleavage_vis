# Streamlit lollipop GUI with separate UniProt input

This version keeps the earlier simple GUI and matching logic, but separates the UniProt accession from the alias/gene-name field.

## Files

- `streamlit_lollipop_app.py` — Streamlit GUI
- `run_analysis_lollipop_any_protein.py` — analysis and plotting functions
- `requirements_streamlit_lollipop.txt` — Python dependencies

## Recommended conda environment

```bash
conda create -n lollipop_gui python=3.10
conda activate lollipop_gui
pip install -r requirements_streamlit_lollipop.txt
```

## Run

```bash
streamlit run streamlit_lollipop_app.py
```

## Example PYK10 inputs

- Primary target ID or name: `AT3G09260.1`
- UniProt ID: `Q9SR37`
- Aliases / gene names: `PYK10; BGLU23; LEB; PSR3.1`
- Protein length: `525`

The UniProt ID is included as a search term internally, but it is no longer mixed into the alias/gene-name input box.

## Command-line example

```bash
python run_analysis_lollipop_any_protein.py \
  --input-dir data \
  --output-dir output_pyk10 \
  --target AT3G09260.1 \
  --uniprot-id Q9SR37 \
  --aliases PYK10 BGLU23 LEB PSR3.1 \
  --protein-length 525
```
