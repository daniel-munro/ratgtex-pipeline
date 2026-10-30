"""Check the frozen predictor's feature order, boundary, and invalid-input behavior."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

spec = importlib.util.spec_from_file_location('sex_concordance', Path(__file__).resolve().parents[1] / 'scripts/qc/sex_concordance.py')
qc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(qc)


class SexPredictionTests(unittest.TestCase):
    def test_features_and_boundary(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bed = root / 'expression.bed'
            bed.write_text('#chr\tstart\tend\tphenotype_id\tA\tB\nchrY\t0\t1\tY1\t100\t1\nchrX\t0\t1\tXist\t0\t9\n')
            model = root / 'model.json'
            model.write_text(json.dumps({'max_tpm': 5, 'male_probability_cutoff': .5, 'datasets': {'test': {'mean': [0, 0], 'scale': [1, 1], 'coefficients': [1, -1], 'intercept': -float(np.log1p(5))}}}))
            result, n_y, _ = qc.predict(bed, model, 'test')
            self.assertEqual(n_y, 1)
            self.assertEqual(result.loc['A', 'TPM'], 5)
            self.assertAlmostEqual(result.loc['A', 'p_male'], .5)
            self.assertEqual(result.pred.tolist(), ['M', 'F'])
            with self.assertRaises(KeyError):
                qc.predict(bed, model, 'unknown')
            bed.write_text(bed.read_text().replace('\t100\t1', '\t-1\t1'))
            with self.assertRaises(ValueError):
                qc.predict(bed, model, 'test')


if __name__ == '__main__':
    unittest.main()
