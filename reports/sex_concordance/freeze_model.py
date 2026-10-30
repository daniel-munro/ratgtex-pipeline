"""Freeze the evaluated dataset-held-out Xist + Y models for routine QC."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import sklearn
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('samples', help='Reviewed report samples.tsv')
    parser.add_argument('output', help='Frozen JSON model')
    args = parser.parse_args()
    data = pd.read_csv(args.samples, sep='\t', dtype={'sample': str})
    models = {}
    for dataset in data.dataset.unique():
        test = data.dataset.eq(dataset)
        train = ~test & ~data['sample'].isin(data.loc[test, 'sample']) & data.label.isin(['F', 'M'])
        model = make_pipeline(StandardScaler(), LogisticRegression(C=1, class_weight='balanced', max_iter=1000))
        model.fit(np.log1p(data.loc[train, ['score', 'Xist']]), data.loc[train, 'label'])
        probability = model.predict_proba(np.log1p(data.loc[test, ['score', 'Xist']]))[:, 1]
        np.testing.assert_allclose(probability, data.loc[test, 'xy_logistic_p_male'], rtol=1e-12, atol=1e-12)
        models[dataset] = {
            'mean': model[0].mean_.tolist(), 'scale': model[0].scale_.tolist(),
            'coefficients': model[-1].coef_[0].tolist(), 'intercept': float(model[-1].intercept_[0]),
            'training_rows': int(train.sum()), 'training_rats': int(data.loc[train, 'sample'].nunique()),
        }
    result = {
        'method': 'dataset-held-out Xist + capped chrY logistic regression',
        'features': ['log1p(capped_chrY_TPM_sum)', 'log1p(Xist_TPM)'],
        'max_tpm': 5.0, 'male_probability_cutoff': 0.5,
        'training_source': args.samples,
        'training_sha256': hashlib.sha256(Path(args.samples).read_bytes()).hexdigest(),
        'sklearn_version': sklearn.__version__,
        'training': 'StandardScaler; LogisticRegression(C=1, class_weight=balanced, max_iter=1000). Exclude held-out dataset and all its sample IDs from training.',
        'datasets': models,
    }
    Path(args.output).write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
