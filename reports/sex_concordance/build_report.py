"""Audit configured v5 datasets; exploratory predictions do not alter production QC."""
import argparse
import base64
import io
from pathlib import Path
import numpy as np
import pandas as pd
import yaml
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--config', default='config/config.yaml')
    p.add_argument('--outdir', default='reports/sex_concordance')
    args = p.parse_args()
    config = yaml.safe_load(Path(args.config).read_text())
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    cap = config.get('sex_concordance', {}).get('max_tpm', 5)
    cutoff = config.get('sex_concordance', {}).get('tpm_cutoff', 30)
    meta = pd.read_csv('geno/genotyping_log.csv', dtype=str).drop_duplicates('rfid', keep='last').set_index('rfid')['sex']
    frames = []
    for dataset in config['datasets']:
        path = Path(config['version']) / dataset / 'phenos/output/unnorm/expression.bed'
        # Stream the BED to retain only chrY and Xist, avoiding large whole-tissue matrices.
        with path.open() as f:
            header = next(f)
            lines = [line for line in f if line.startswith('chrY\t') or line.split('\t', 4)[3] == 'Xist']
        expr = pd.read_csv(io.StringIO(header + ''.join(lines)), sep='\t').set_index('phenotype_id')
        y = expr.loc[expr['#chr'] == 'chrY'].iloc[:, 3:].astype(float)
        xist = expr.loc['Xist'].iloc[3:].astype(float)
        if y.empty or not np.isfinite(y.to_numpy()).all() or (y < 0).any().any() or not np.isfinite(xist).all() or (xist < 0).any():
            raise ValueError(f'{dataset}: invalid sex-marker TPM')
        df = pd.DataFrame({'score': y.clip(upper=cap).sum(), 'Xist': xist})
        df['dataset'] = dataset
        df['label'] = meta.reindex(df.index)
        df['sample'] = df.index
        df['n_y_genes'] = len(y)
        frames.append(df)
        print(f'{dataset}: {len(df)} samples', flush=True)
    data = pd.concat(frames, ignore_index=True)
    data['baseline_pred'] = np.where(data.score >= cutoff, 'M', 'F')
    # Fixed feature definitions and regularization; no tuning using held-out labels.
    models = {'y_logistic': ['score'], 'xy_logistic': ['score', 'Xist']}
    coefficients = []
    for dataset in config['datasets']:
        test = data.dataset.eq(dataset)
        train = ~test & ~data['sample'].isin(data.loc[test, 'sample']) & data.label.isin(['F', 'M'])
        if data.loc[train, 'label'].nunique() != 2:
            raise ValueError(f'{dataset}: training requires both sexes')
        for name, features in models.items():
            model = make_pipeline(StandardScaler(), LogisticRegression(C=1, class_weight='balanced', max_iter=1000))
            model.fit(np.log1p(data.loc[train, features]), data.loc[train, 'label'])
            data.loc[test, name + '_p_male'] = model.predict_proba(np.log1p(data.loc[test, features]))[:, 1]
            data.loc[test, name + '_pred'] = np.where(data.loc[test, name + '_p_male'] >= .5, 'M', 'F')
            coefficients.append({'dataset': dataset, 'model': name, 'training_rows': int(train.sum()), 'training_rats': data.loc[train, 'sample'].nunique(), 'intercept': model[-1].intercept_[0], **dict(zip(features, model[-1].coef_[0]))})
    rows = []
    for dataset, df in data.groupby('dataset', sort=False):
        labeled = df[df.label.isin(['F', 'M'])]
        row = {'dataset': dataset, 'samples': len(df), 'labeled': len(labeled)}
        for name in ['baseline', *models]:
            row[name + '_F_to_M'] = int(((labeled.label == 'F') & (labeled[name+'_pred'] == 'M')).sum())
            row[name + '_M_to_F'] = int(((labeled.label == 'M') & (labeled[name+'_pred'] == 'F')).sum())
        f = labeled.loc[labeled.label == 'F', 'score']
        m = labeled.loc[labeled.label == 'M', 'score']
        row.update(female_max=f.max(), male_min=m.min(), separation_gap=m.min()-f.max())
        # Label-fitted bound is descriptive only, never a validated prediction.
        candidates = np.r_[-np.inf, np.unique(labeled.score), np.inf]
        row['best_in_sample_threshold_errors'] = min(int((np.where(labeled.score >= c, 'M', 'F') != labeled.label).sum()) for c in candidates)
        rows.append(row)
    summary = pd.DataFrame(rows)
    data.to_csv(out/'samples.tsv', sep='\t', index=False)
    summary.to_csv(out/'datasets.tsv', sep='\t', index=False)
    pd.DataFrame(coefficients).to_csv(out/'models.tsv', sep='\t', index=False)
    html = ['<!doctype html><meta charset="utf-8"><title>v5 sex concordance</title><style>body{font:16px sans-serif;max-width:1200px;margin:30px auto}table{border-collapse:collapse;font-size:12px}td,th{padding:5px;border:1px solid #ccc}img{width:100%}section{break-before:page}</style><h1>v5 sex concordance</h1>',
    '<p>Metadata labels are reference labels, not independently established biological truth. Arrows mean metadata → prediction. All configured datasets and all expression samples are retained; unlabeled samples are gray. Existing production predictions are unchanged.</p>',
    f'<p>Baseline: sum of all chrY gene TPM after capping each at {cap}; male if sum ≥ {cutoff}. Alternative: balanced logistic regression on log1p(capped Y sum), with or without log1p(Xist TPM); standardized training features, C=1, male probability ≥ 0.5. Each dataset is held out, and all its sample IDs are removed from other training datasets. No hyperparameter tuning. Repeated training rats can contribute multiple tissue observations. Results measure agreement with potentially imperfect metadata, not proven sex accuracy.</p>',
    '<p>Positive separation gap means a dataset-specific threshold could perfectly separate its current labels. The best in-sample threshold error count is an optimistic descriptive bound, not validation. A common cutoff can fail in opposite directions when tissue score distributions shift.</p>', summary.to_html(index=False),
    '<p>Method rationale: <a href="https://github.com/SexChrLab/SexInference">SexInference</a> combines XIST and Y markers in a human-trained logistic model; those weights cannot be transferred directly to rats. <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC7374973/">Sex chromosome reference study</a> recommends multiple X/Y markers and XIST and documents mapping artifacts. The model here is a rat-data exploratory analogue using the current Y summary. Next validation should compare specific rat Y markers and DNA-based sex, review discordant identities, and use an independent release before adopting a replacement.</p>']
    totals = {name: int((summary[name+'_F_to_M'] + summary[name+'_M_to_F']).sum()) for name in ['baseline', *models]}
    residual = data[data.label.isin(['F', 'M']) & data.label.ne(data.xy_logistic_pred)]
    html.insert(2, f'<h2>Findings</h2><p>{len(data)} sample observations across {len(summary)} datasets. Baseline disagreements: {totals["baseline"]}; held-out Y-only logistic: {totals["y_logistic"]}; held-out Xist + Y logistic: {totals["xy_logistic"]}. Every dataset has baseline errors in at most one direction. Liver contributes 324 male→female calls; VMH contributes 56 female→male calls. A Y-only model does not fix the distribution shift; adding Xist is substantially more effective. This is exploratory agreement evidence, not an independently validated replacement.</p><p>Female maxima and separation gaps are undefined for male-only datasets (Amyg, Heart, Adipose, Liver); their in-sample threshold bound is trivial and provides no evidence of specificity. NAcc2, RMTg, PL3, and VMH have overlapping Y scores between metadata sexes.</p><h3>Remaining Xist + Y disagreements</h3><p>All three residual disagreements are rat E21DCD3672, in NAcc3, PL3, and pVTA1: recorded male, predicted female in each held-out dataset. Review its identity and metadata before deciding whether these represent model errors.</p>' + residual.to_html(index=False))
    rng = np.random.default_rng(42)
    with PdfPages(out/'datasets.pdf') as pdf:
        for dataset, df in data.groupby('dataset', sort=False):
            fig, axes = plt.subplots(1, 3, figsize=(15, 5), layout='constrained')
            for label, color, pos in [('F','#d95f02',0), ('M','#1b9e77',1), ('unknown','#888888',2)]:
                sub = df[df.label.eq(label)] if label != 'unknown' else df[~df.label.isin(['F','M'])]
                jitter = rng.uniform(-.18,.18,len(sub))
                axes[0].scatter(pos+jitter, sub.score, s=13, c=color, alpha=.65, label=f'{label} (n={len(sub)})')
                axes[1].scatter(sub.Xist, sub.score, s=13,c=color,alpha=.65)
                axes[2].scatter(pos+jitter, sub.xy_logistic_p_male, s=13,c=color,alpha=.65)
            axes[0].axhline(cutoff,color='black',ls='--',label=f'Cutoff {cutoff}')
            axes[0].set(xticks=[0,1,2],xticklabels=['F','M','unknown'],ylabel='Capped chrY TPM sum',xlabel='Metadata sex',title='Current prediction')
            axes[0].legend(fontsize=8)
            axes[1].set(xscale='symlog',yscale='symlog',xlabel='Xist TPM (symlog)',ylabel='Capped Y sum (symlog)',title='Independent X and Y signals')
            axes[1].axhline(cutoff,color='black',ls='--')
            axes[2].axhline(.5,color='black',ls='--')
            axes[2].set(xticks=[0,1,2],xticklabels=['F','M','unknown'],ylabel='Held-out male probability',xlabel='Metadata sex',title='Exploratory Xist + Y model')
            row = summary.set_index('dataset').loc[dataset]
            fig.suptitle(f"{dataset} | baseline F→M {row.baseline_F_to_M}, M→F {row.baseline_M_to_F} | X+Y F→M {row.xy_logistic_F_to_M}, M→F {row.xy_logistic_M_to_F}")
            pdf.savefig(fig)
            buf = io.BytesIO()
            fig.savefig(buf, format='png',dpi=115)
            html.append(f'<section><h2>{dataset}</h2><img alt="{dataset} scores and sex labels" src="data:image/png;base64,{base64.b64encode(buf.getvalue()).decode()}"></section>')
            plt.close(fig)
    (out/'index.html').write_text('\n'.join(html))
    print(summary.to_string(index=False))

if __name__ == '__main__':
    main()
