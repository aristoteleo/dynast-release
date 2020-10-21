from unittest import mock, TestCase

import numpy as np
from scipy import sparse, stats

import dynast.estimation.p_c as p_c
from tests import mixins


class TestPC(mixins.TestMixin, TestCase):

    def test_read_p_c(self):
        pass

    def test_binomial_pmf(self):
        self.assertEqual(stats.binom.pmf(0, 0, 0), p_c.binomial_pmf(0, 0, 0))
        self.assertEqual(stats.binom.pmf(0, 1, 0.5), p_c.binomial_pmf(0, 1, 0.5))
        self.assertEqual(stats.binom.pmf(1, 1, 0), p_c.binomial_pmf(1, 1, 0))

    def test_filter_aggregates_part(self):
        with mock.patch('dynast.estimation.p_c.binomial_pmf') as binomial_pmf:
            binomial_pmf.side_effect = [1] * 16
            values = np.array([
                [0, 0, 3],
                [1, 3, 5],
                [3, 3, 2],
            ])
            p_e = 0.5
            sp, mask = p_c.filter_aggregates_part(values, p_e)
            self.assertTrue((np.array([
                [3, 0, 0, 0],
                [0, 0, 0, 5],
                [0, 0, 0, 0],
                [0, 0, 0, 2],
            ]) == sp.toarray()).all())
            self.assertEqual([(0, 3), (1, 3), (2, 3)], mask)

    def test_filter_aggregates(self):
        pass

    def test_expectation_maximization_part(self):
        sp = sparse.lil_matrix(np.array([
            [3, 0, 0, 0],
            [0, 0, 0, 5],
            [0, 0, 0, 0],
            [0, 0, 0, 2],
        ]))
        mask = [(0, 3), (1, 3), (2, 3)]
        new_sp, estimated = p_c.expectation_maximization_part(sp, mask)
        self.assertTrue(
            np.allclose(
                np.array([
                    [3, 0, 0, 0.08767744],
                    [0, 0, 0, 0.74597309],
                    [0, 0, 0, 2.11561777],
                    [0, 0, 0, 2],
                ]), new_sp.toarray()
            )
        )
        self.assertAlmostEqual(0.7393152527850702, estimated)

    def test_expectation_maximization(self):
        pass

    def test_estimate_p_c(self):
        pass
