"""Predict sex with frozen, dataset-held-out Xist + Y expression models."""

import argparse
import io
import json
from pathlib import Path

import numpy as np
import pandas as pd


def predict(expr_path, model_path, dataset):
    model = json.loads(Path(model_path).read_text())
    parameters = model['datasets'][dataset]
    # Retain only sex markers while streaming potentially large expression BEDs.
    with open(expr_path) as handle:
        header = next(handle)
        lines = [line for line in handle
                 if line.startswith('chrY\t') or line.split('\t', 4)[3] == 'Xist']
    expr = pd.read_csv(io.StringIO(header + ''.join(lines)), sep='\t').set_index('phenotype_id')
    y = expr.loc[expr['#chr'] == 'chrY'].drop(columns=['#chr', 'start', 'end']).astype(float)
    xist = expr.loc[['Xist']].drop(columns=['#chr', 'start', 'end']).astype(float)
    if y.empty or len(xist) != 1:
        raise ValueError('Expression BED must contain chrY genes and exactly one Xist row')
    for values in (y, xist):
        if not np.isfinite(values.to_numpy()).all() or (values < 0).any().any():
            raise ValueError('Sex-marker TPM values must be finite and nonnegative')
    result = pd.DataFrame({'TPM': y.clip(upper=model['max_tpm']).sum(), 'Xist_TPM': xist.iloc[0]})
    features = np.log1p(result[['TPM', 'Xist_TPM']].to_numpy())
    standardized = (features - np.array(parameters['mean'])) / np.array(parameters['scale'])
    log_odds = standardized @ np.array(parameters['coefficients']) + parameters['intercept']
    result['p_male'] = np.exp(-np.logaddexp(0, -log_odds))
    result['pred'] = np.where(result.p_male >= model['male_probability_cutoff'], 'M', 'F')
    return result, len(y), model


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('expr', help='Unnormalized expression BED with TPM values')
    parser.add_argument('meta', help='Metadata CSV with rfid and sex columns')
    parser.add_argument('out', help='Sex-concordance text report')
    parser.add_argument('--dataset', required=True, help='Source dataset identifying the frozen held-out model')
    parser.add_argument('--model', default='config/sex_concordance_model.json', help='Frozen model JSON')
    args = parser.parse_args()
    expr, n_y_genes, model = predict(args.expr, args.model, args.dataset)
    meta = pd.read_csv(args.meta, dtype=str).drop_duplicates('rfid', keep='last').set_index('rfid')
    expr = expr.join(meta[['sex']].rename(columns={'sex': 'label'}), how='left')
    labeled = expr.loc[expr.label.isin(['F', 'M'])]
    with open(args.out, 'w') as out:
        out.write(f"Model: {model['method']}; held-out dataset: {args.dataset}.\n")
        out.write(f"Model file: {args.model}; training SHA256: {model['training_sha256']}.\n")
        out.write(f"chrY genes: {n_y_genes}; per-gene TPM cap: {model['max_tpm']}; male probability cutoff: {model['male_probability_cutoff']}.\n")
        out.write(f'F/M sex labels present for {len(labeled)} of {len(expr)} samples.\n')
        out.write('\nExpression and male probability by metadata sex:\n')
        out.write(labeled.groupby('label')[['TPM', 'Xist_TPM', 'p_male']].agg(['count', 'min', 'median', 'max']).to_string() + '\n\n')
        out.write(f'{(labeled.pred == labeled.label).sum()} matches between expression-predicted sex and metadata labels.\n')
        mismatches = labeled.loc[labeled.pred != labeled.label]
        if len(mismatches):
            out.write('Samples with mismatch:\n' + mismatches.to_string() + '\n')
        out.write('\nAll results (TPM = capped chrY sum; p_male = Xist + Y model probability):\n')
        out.write(expr.to_string() + '\n')


if __name__ == '__main__':
    main()
