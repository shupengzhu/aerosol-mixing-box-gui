"""Result rounding must not invent NaN from finite simulation outputs."""
import json
import unittest

import numpy as np

from scrambox.server import _round


class ResultRoundingTest(unittest.TestCase):
    def test_subnormal_finite_values_remain_finite_and_serializable(self):
        raw = np.array([[0.0, 6.7422242e-310, 1.937e-321, -5e-324],
                        [1e-303, 1e-304, 1e-305, 1e308]])
        result = _round(raw)
        self.assertEqual(np.asarray(result).shape, raw.shape)
        self.assertTrue(np.all(np.isfinite(result)))
        json.dumps(result, allow_nan=False)
        self.assertEqual(result[0][1], raw[0, 1])
        self.assertEqual(result[0][2], raw[0, 2])
        self.assertEqual(result[0][3], raw[0, 3])

    def test_ordinary_values_keep_five_significant_digits(self):
        result = _round([123.456789, 0.00123456789, -98.7654321, 0.0])
        self.assertEqual(result, [123.46, 0.0012346, -98.765, 0.0])

    def test_raw_nonfinite_is_not_silently_relabelled_as_valid(self):
        for invalid in (np.nan, np.inf, -np.inf):
            with self.subTest(invalid=invalid):
                with self.assertRaises(ValueError):
                    _round([invalid])


if __name__ == "__main__":
    unittest.main()
