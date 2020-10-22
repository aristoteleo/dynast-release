import os
from unittest import mock, TestCase

import pandas as pd

import dynast.estimation.pi as pi
from tests import mixins


class TestPi(mixins.TestMixin, TestCase):

    def test_fit_stan(self):
        values = [[1, 2, 1], [4, 5, 2]]
        p_e = mock.MagicMock()
        p_c = mock.MagicMock()
        model = mock.MagicMock()
        model.sampling.return_value.extract.return_value = {'alpha': [1, 2], 'beta': [3, 4]}
        n_chains = 5
        self.assertEqual((1.5, 3.5), pi.fit_stan(values, p_e, p_c, model, n_chains=n_chains))
        model.sampling.assert_called_once_with(
            data={
                'N': 3,
                'contents': [2, 5, 5],
                'conversions': [1, 4, 4],
                'p_c': p_c,
                'p_e': p_e
            },
            n_jobs=mock.ANY,
            iter=mock.ANY,
            chains=mock.ANY,
            init=[{
                'log_alpha': 0.0,
                'log_beta': 0.0,
                'pi_g': 0.5
            }] * n_chains,
            control=mock.ANY
        )
        model.sampling.return_value.extract.assert_called_once_with(('alpha', 'beta'))

    def test_estimate_pi_by_barcode_by_barcode_gene(self):
        df = pd.read_csv(self.aggregate_paths['GC'])
        p_e = dict(pd.read_csv(self.p_e_barcode_path).set_index('barcode')['p_e'])
        p_c = dict(pd.read_csv(self.p_c_barcode_path).set_index('barcode')['p_c'])
        filter_dict = {'barcode': ['AAGATGTCCAGG']}
        pi_path = os.path.join(self.temp_dir, 'pi.csv')
        with mock.patch('dynast.estimation.pi.ProcessPoolExecutor') as executor, \
            mock.patch('dynast.estimation.pi.tqdm'), \
            mock.patch('dynast.estimation.pi.pystan'), \
            mock.patch('dynast.estimation.pi.as_completed', mixins.tqdm_mock):
            m = mock.MagicMock()
            m.result.return_value = (2, 2)
            executor.return_value.__enter__.return_value.submit.return_value = m

            self.assertEqual(
                pi_path,
                pi.estimate_pi(
                    df,
                    p_e,
                    p_c,
                    pi_path,
                    filter_dict=filter_dict,
                    p_group_by=['barcode'],
                    group_by=['barcode', 'GX'],
                    value_columns=['GC', 'G', 'count']
                )
            )
            self.assertTrue(mixins.files_equal(self.pi_barcode_barcode_gene_path, pi_path))
