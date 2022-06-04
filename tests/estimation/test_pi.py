import os
import pickle
from unittest import TestCase, mock

from scipy import stats

import dynast.estimation.p_c as p_c
import dynast.estimation.p_e as p_e
import dynast.estimation.pi as pi
import dynast.preprocessing.aggregation as aggregation

from .. import mixins


class TestPi(mixins.TestMixin, TestCase):

    def test_read_pi(self):
        pi.read_pi(self.umi_pi_path, group_by=['barcode', 'GX'])

        result = pi.read_pi(self.pi_int_path, group_by=['group'])
        self.assertEqual({'1': 0.5, '2': 0.75}, result)

    def test_beta_mean(self):
        self.assertAlmostEqual(stats.beta.mean(1, 1), pi.beta_mean(1, 1))
        self.assertAlmostEqual(stats.beta.mean(1, 2), pi.beta_mean(1, 2))

    def test_beta_mode(self):
        self.assertAlmostEqual(0, pi.beta_mode(1, 2))

    def test_guess_beta_parameters(self):
        alpha, beta = pi.guess_beta_parameters(0.8, strength=5)
        self.assertEqual(0.8, stats.beta.mean(alpha, beta))

    def test_estimate_pi(self):
        pi_path = os.path.join(self.temp_dir, 'pi.csv')
        with mock.patch('dynast.estimation.pi.pystan.StanModel') as StanModel, \
            mock.patch('dynast.estimation.pi.utils.as_completed_with_progress', mixins.tqdm_mock):
            model = mock.MagicMock()
            StanModel.return_value = model
            model.sampling.return_value.extract.return_value = {'alpha': [2], 'beta': [2], 'pi_g': [0.5]}

            try:
                self.assertEqual(
                    pi_path,
                    pi.estimate_pi(
                        aggregation.read_aggregates(self.umi_aggregates_path),
                        p_e.read_p_e(self.umi_p_e_path, group_by=['barcode']),
                        p_c.read_p_c(self.umi_p_c_path, group_by=['barcode']),
                        pi_path,
                        group_by=['barcode', 'GX'],
                        p_group_by=['barcode'],
                        n_threads=1,
                        threshold=1,
                        seed=None,
                    )
                )
                with open(pi_path, 'r') as f:
                    self.assertTrue(
                        f.read().
                        startswith('barcode,GX,guess,alpha,beta,pi\nAAACCCAACGTA,ENSG00000172009,0.99,2.0,2.0,0.5\n')
                    )
            except pickle.PicklingError:
                pass
