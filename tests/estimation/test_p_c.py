import os
from unittest import mock, TestCase

import numpy as np
import pandas as pd
from scipy import sparse, stats

import dynast.estimation.p_c as p_c
from tests import mixins


class TestPC(mixins.TestMixin, TestCase):

    def test_read_p_c_by_barcode(self):
        self.assertEqual({
            'AAGATGTCCAGG': 0.2,
            'AGCGCGAATTAT': 0.0
        }, p_c.read_p_c(self.p_c_barcode_path, group_by=['barcode']))

    def test_read_p_c_by_barcode_gene(self):
        self.assertEqual({('AGCGCGAATTAT', 'ENSG00000186092'): 0.0,
                          ('AAGATGTCCAGG', 'ENSG00000279928'): 0.2},
                         p_c.read_p_c(self.p_c_barcode_gene_path, group_by=['barcode', 'GX']))

    def test_read_p_c_by_none(self):
        self.assertEqual(0.125, p_c.read_p_c(self.p_c_none_path, group_by=None))

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

    def test_filter_aggregates_by_barcode(self):
        df = pd.read_csv(self.aggregate_paths['TC'])
        p_e = {
            'AAGATGTCCAGG': 0.2,
            'AGCGCGAATTAT': 0.3,
        }
        with mock.patch('dynast.estimation.p_c.tqdm'), \
            mock.patch('dynast.estimation.p_c.ProcessPoolExecutor') as executor, \
            mock.patch('dynast.estimation.p_c.as_completed', mixins.tqdm_mock):

            mock1 = mock.MagicMock()
            mock2 = mock.MagicMock()
            executor.return_value.__enter__.return_value.submit.side_effect = [mock1, mock2]

            self.assertEqual({
                'AAGATGTCCAGG': mock1.result.return_value,
                'AGCGCGAATTAT': mock2.result.return_value,
            }, p_c.filter_aggregates(df, p_e, group_by=['barcode']))

            executor.assert_called_once()
            executor.return_value.__enter__.return_value.submit.assert_has_calls([
                mock.call(p_c.filter_aggregates_part, mock.ANY, 0.2),
                mock.call(p_c.filter_aggregates_part, mock.ANY, 0.3),
            ])

    def test_filter_aggregates_by_gene(self):
        df = pd.read_csv(self.aggregate_paths['TC'])
        p_e = {
            'ENSG00000279928': 0.2,
            'ENSG00000186092': 0.3,
        }
        with mock.patch('dynast.estimation.p_c.tqdm'), \
            mock.patch('dynast.estimation.p_c.ProcessPoolExecutor') as executor, \
            mock.patch('dynast.estimation.p_c.as_completed', mixins.tqdm_mock):

            mock1 = mock.MagicMock()
            mock2 = mock.MagicMock()
            executor.return_value.__enter__.return_value.submit.side_effect = [mock1, mock2]

            self.assertEqual({
                'ENSG00000279928': mock2.result.return_value,
                'ENSG00000186092': mock1.result.return_value,
            }, p_c.filter_aggregates(df, p_e, group_by=['GX']))

            executor.assert_called_once()
            executor.return_value.__enter__.return_value.submit.assert_has_calls([
                mock.call(p_c.filter_aggregates_part, mock.ANY, 0.3),
                mock.call(p_c.filter_aggregates_part, mock.ANY, 0.2),
            ])

    def test_filter_aggregates_by_none(self):
        df = pd.read_csv(self.aggregate_paths['TC'])
        p_e = 0.5
        with mock.patch('dynast.estimation.p_c.filter_aggregates_part') as filter_aggregates_part, \
            mock.patch('dynast.estimation.p_c.tqdm'), \
            mock.patch('dynast.estimation.p_c.ProcessPoolExecutor') as executor, \
            mock.patch('dynast.estimation.p_c.as_completed', mixins.tqdm_mock):
            self.assertEqual(filter_aggregates_part.return_value, p_c.filter_aggregates(df, p_e, group_by=None))
            filter_aggregates_part.assert_called_once_with(mock.ANY, 0.5)
            executor.assert_not_called()

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

    def test_expectation_maximization_by_barcode(self):
        sp = sparse.lil_matrix(np.array([
            [3, 0, 0, 0],
            [0, 0, 0, 5],
            [0, 0, 0, 0],
            [0, 0, 0, 2],
        ]))
        mask = [(0, 3), (1, 3), (2, 3)]
        filtered = {'BC': (sp, mask)}
        with mock.patch('dynast.estimation.p_c.ProcessPoolExecutor') as executor, \
            mock.patch('dynast.estimation.p_c.as_completed', mixins.tqdm_mock), \
            mock.patch('dynast.estimation.p_c.tqdm'):
            m = mock.MagicMock()
            executor.return_value.__enter__.return_value.submit.return_value = m
            self.assertEqual({'BC': m.result.return_value},
                             p_c.expectation_maximization(filtered, group_by=['barcode']))
            executor.return_value.__enter__.return_value.submit.assert_called_once_with(
                p_c.expectation_maximization_part, sp, mask
            )

    def test_expectation_maximization_by_gene(self):
        sp = sparse.lil_matrix(np.array([
            [3, 0, 0, 0],
            [0, 0, 0, 5],
            [0, 0, 0, 0],
            [0, 0, 0, 2],
        ]))
        mask = [(0, 3), (1, 3), (2, 3)]
        filtered = {'GENE': (sp, mask)}
        with mock.patch('dynast.estimation.p_c.ProcessPoolExecutor') as executor, \
            mock.patch('dynast.estimation.p_c.as_completed', mixins.tqdm_mock), \
            mock.patch('dynast.estimation.p_c.tqdm'):
            m = mock.MagicMock()
            executor.return_value.__enter__.return_value.submit.return_value = m
            self.assertEqual({'GENE': m.result.return_value}, p_c.expectation_maximization(filtered, group_by=['gene']))
            executor.return_value.__enter__.return_value.submit.assert_called_once_with(
                p_c.expectation_maximization_part, sp, mask
            )

    def test_expectation_maximization_by_none(self):
        sp = sparse.lil_matrix(np.array([
            [3, 0, 0, 0],
            [0, 0, 0, 5],
            [0, 0, 0, 0],
            [0, 0, 0, 2],
        ]))
        mask = [(0, 3), (1, 3), (2, 3)]
        filtered = (sp, mask)
        with mock.patch('dynast.estimation.p_c.expectation_maximization_part') as expectation_maximization_part, \
            mock.patch('dynast.estimation.p_c.ProcessPoolExecutor') as executor, \
            mock.patch('dynast.estimation.p_c.as_completed', mixins.tqdm_mock), \
            mock.patch('dynast.estimation.p_c.tqdm'):
            self.assertEqual(
                expectation_maximization_part.return_value, p_c.expectation_maximization(filtered, group_by=None)
            )
            expectation_maximization_part.assert_called_once_with(sp, mask)
            executor.assert_not_called()

    def test_estimate_p_c_by_barcode(self):
        df = pd.read_csv(self.aggregate_paths['GC'])
        p_e = dict(pd.read_csv(self.p_e_barcode_path).set_index('barcode')['p_e'])
        p_c_path = os.path.join(self.temp_dir, 'p_c.csv')
        aggregate_path = os.path.join(self.temp_dir, 'aggregate.csv')
        with mock.patch('dynast.estimation.p_c.tqdm'):
            result = p_c.estimate_p_c(
                df, p_e, p_c_path, aggregate_path, group_by=['barcode'], value_columns=['GC', 'G', 'count']
            )
            self.assertEqual(({'AGCGCGAATTAT': 0.0, 'AAGATGTCCAGG': 0.2}, p_c_path, aggregate_path), result)
            self.assertTrue(mixins.files_equal(self.p_c_barcode_path, p_c_path))
            self.assertTrue(mixins.files_equal(self.aggregate_barcode_path, aggregate_path))

    def test_estimate_p_c_by_gene(self):
        df = pd.read_csv(self.aggregate_paths['GC'])
        p_e = dict(pd.read_csv(self.p_e_gene_path).set_index('GX')['p_e'])
        p_c_path = os.path.join(self.temp_dir, 'p_c.csv')
        aggregate_path = os.path.join(self.temp_dir, 'aggregate.csv')
        with mock.patch('dynast.estimation.p_c.tqdm'):
            result = p_c.estimate_p_c(
                df, p_e, p_c_path, aggregate_path, group_by=['GX'], value_columns=['GC', 'G', 'count']
            )
            self.assertEqual(({'ENSG00000186092': 0.0, 'ENSG00000279928': 0.2}, p_c_path, aggregate_path), result)
            self.assertTrue(mixins.files_equal(self.p_c_gene_path, p_c_path))
            self.assertTrue(mixins.files_equal(self.aggregate_gene_path, aggregate_path))

    def test_estimate_p_c_by_barcode_gene(self):
        df = pd.read_csv(self.aggregate_paths['GC'])
        p_e = dict(pd.read_csv(self.p_e_barcode_gene_path).set_index(['barcode', 'GX'])['p_e'])
        p_c_path = os.path.join(self.temp_dir, 'p_c.csv')
        aggregate_path = os.path.join(self.temp_dir, 'aggregate.csv')
        with mock.patch('dynast.estimation.p_c.tqdm'):
            result = p_c.estimate_p_c(
                df, p_e, p_c_path, aggregate_path, group_by=['barcode', 'GX'], value_columns=['GC', 'G', 'count']
            )
            self.assertEqual(({('AGCGCGAATTAT', 'ENSG00000186092'): 0.0,
                               ('AAGATGTCCAGG', 'ENSG00000279928'): 0.2}, p_c_path, aggregate_path), result)
            self.assertTrue(mixins.files_equal(self.p_c_barcode_gene_path, p_c_path))
            self.assertTrue(mixins.files_equal(self.aggregate_barcode_gene_path, aggregate_path))

    def test_estimate_p_c_by_none(self):
        df = pd.read_csv(self.aggregate_paths['GC'])
        with open(self.p_e_none_path, 'r') as f:
            p_e = float(f.read())
        p_c_path = os.path.join(self.temp_dir, 'p_c.csv')
        aggregate_path = os.path.join(self.temp_dir, 'aggregate.csv')
        with mock.patch('dynast.estimation.p_c.tqdm'):
            result = p_c.estimate_p_c(
                df, p_e, p_c_path, aggregate_path, group_by=None, value_columns=['GC', 'G', 'count']
            )
            self.assertEqual((0.125, p_c_path, aggregate_path), result)
            self.assertTrue(mixins.files_equal(self.p_c_none_path, p_c_path))
            self.assertTrue(mixins.files_equal(self.aggregate_none_path, aggregate_path))
