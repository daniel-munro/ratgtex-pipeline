# v5 sex-concordance review

Open `index.html` for a self-contained report, or `datasets.pdf` for one figure per configured source dataset. Each figure shows the current score and cutoff by metadata sex, Xist versus Y expression, and held-out alternative predictions. `samples.tsv` provides individual values and predictions; `datasets.tsv` gives directional mismatch counts and separation; `models.tsv` records fold training sizes and standardized coefficients.

Regenerate from the project root:

```bash
snakemake --profile slurm -n -p qc_sex_concordance_report
snakemake --profile slurm qc_sex_concordance_report
```

The script requires numpy, pandas, PyYAML, matplotlib, and scikit-learn. Dataset membership comes from `config/config.yaml`, expression from every configured dataset's unnormalized Pantry BED, and labels from `geno/genotyping_log.csv` (last record per RFID, matching existing QC). Held-out datasets also exclude their rats from training in other datasets. Models use fixed regularization and features, without tuning. The report does not change the production predictor or metadata. Metadata agreement does not establish biological sex or identity; review residuals against identity QC and independent genotype evidence before adopting the model.

## Frozen production predictor

`config/sex_concordance_model.json` now supplies routine `qc_sex_concordance` predictions. It preserves the evaluated per-dataset held-out fits, including scaler parameters, coefficients, training sample counts, source-table SHA256, and scikit-learn version. There is no runtime retraining. The original HTML/PDF remains the historical comparison against the old threshold method.

To reproduce the frozen artifact from the reviewed sample table (requires scikit-learn):

```bash
OPENBLAS_NUM_THREADS=1 python reports/sex_concordance/freeze_model.py reports/sex_concordance/samples.tsv config/sex_concordance_model.json
```

Export asserts that every fitted probability reproduces the reviewed report before writing the artifact. Replacing the training table or extending dataset coverage requires another reviewed export. Frozen folds exclude the sample IDs present in the reviewed dataset snapshot; additions sharing rats with the old training data require renewed evaluation. Runtime prediction needs only numpy and pandas. In the reviewed v5 data, all 3,909 predictions are reproduced, with three metadata disagreements belonging to E21DCD3672 across three tissues.
