import os
from unittest import mock, TestCase

from scipy import stats

import dynast.estimation.p_c as p_c
import dynast.estimation.p_e as p_e
import dynast.preprocessing.aggregation as aggregation

from .. import mixins


class TestPC(mixins.TestMixin, TestCase):

    def test_read_p_c(self):
        p_c.read_p_c(self.umi_p_c_path, group_by=['barcode'])

        result = p_c.read_p_c(self.p_c_int_path, group_by=['group'])
        self.assertEqual({'1': 0.5, '2': 0.75}, result)

    def test_binomial_pmf(self):
        self.assertAlmostEqual(p_c.binomial_pmf(5, 10, 0.5), stats.binom.pmf(5, 10, 0.5))
        self.assertAlmostEqual(p_c.binomial_pmf(10, 10, 0.5), stats.binom.pmf(10, 10, 0.5))
        self.assertAlmostEqual(p_c.binomial_pmf(0, 10, 0.5), stats.binom.pmf(0, 10, 0.5))

    def test_estimate_p_c(self):
        p_c_path = os.path.join(self.temp_dir, 'p_c.csv')
        df_aggregates = aggregation.read_aggregates(self.umi_aggregates_path)
        with mock.patch('dynast.estimation.p_c.utils.as_completed_with_progress', mixins.tqdm_mock):
            self.assertEqual(
                p_c_path,
                p_c.estimate_p_c(
                    df_aggregates,
                    p_e.read_p_e(self.umi_p_e_path, group_by=['barcode']),
                    p_c_path,
                    group_by=['barcode'],
                    threshold=0,
                    n_threads=2
                )
            )
            self.assertTrue(mixins.files_equal(self.umi_p_c_path, p_c_path))

    def test_estimate_p_c_nasc(self):
        p_c_path = os.path.join(self.temp_dir, 'p_c.csv')
        df_aggregates = aggregation.read_aggregates(self.nasc_aggregates_path)
        with mock.patch('dynast.estimation.p_c.utils.as_completed_with_progress', mixins.tqdm_mock):
            self.assertEqual(
                p_c_path,
                p_c.estimate_p_c(
                    df_aggregates,
                    p_e.read_p_e(self.nasc_p_e_path, group_by=['barcode']),
                    p_c_path,
                    group_by=['barcode'],
                    threshold=0,
                    n_threads=2,
                    nasc=True,
                )
            )
            self.assertTrue(mixins.files_equal(self.nasc_p_c_path, p_c_path))
