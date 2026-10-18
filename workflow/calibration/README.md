# Calibration (Table 1)

`run_mouse_calib.sh` documents how the counting settings for mouse were chosen: the
mouse Ensembl STAR alignments (made from the ENA reads by `../per_database/`) are
recounted with Ensembl 110/112 and several featureCounts settings, and the setting
closest to the published DEG list is kept (Ensembl 110, `-s 0`, featureCounts defaults;
85 DEGs vs 69 published). The script's paths point to our cluster storage.

The calibration itself is part of `run_all.sh` and uses only the public data:
`scripts/calibration.R` compares the mouse and pig Ensembl results with the published
ones (mouse 85 vs 69; pig 234 and 59, identical gene lists), and the dog calibration
(protein-coding genes, as in the published analysis) is computed on the 22 SRA libraries
(`tables/calibration_dog.txt`).
