# Streamlit lollipop GUI: strict metacaspase extraction

This update keeps the Streamlit GUI workflow but changes metacaspase extraction to be strict.

## Key fix

The metacaspase matrix now retains only rows where the primary `Protein` column exactly matches the primary target ID supplied by the user.

For example, if the target is:

```text
AT3G09260.1
```

then the metacaspase table will keep only:

```text
Protein == AT3G09260.1
```

It will not keep homologous or shared mappings such as:

```text
AT1G66270.2
AT1G66280.1
AT3G16420.2
AT3G21370.1
```

even if those rows contain the target protein in `Isoforms` or `Description`.

## Files

Place these files in the same folder:

```text
streamlit_lollipop_app_strict_metacaspase.py
run_analysis_lollipop_any_protein.py
requirements_streamlit_lollipop.txt
```

The ZIP includes the backend under the import name `run_analysis_lollipop_any_protein.py`, so the GUI can run without editing the import line.

Run:

```bash
streamlit run streamlit_lollipop_app_strict_metacaspase.py
```

## Notes

The GUI still accepts UniProt ID and aliases for semitryptome and HUNTER matching, but metacaspase extraction is intentionally restricted to the primary `Protein` column.
