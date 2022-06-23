import os
from unittest import TestCase

import dynast.estimation.alpha as alpha
import dynast.estimation.pi as pi
import dynast.preprocessing.conversion as conversion
import dynast.utils as utils

from .. import mixins


class TestAlpha(mixins.TestMixin, TestCase):

    def test_read_alpha(self):
        alpha.read_alpha(self.umi_alpha_path, group_by=['barcode'])

    def test_estimate_alpha(self):
        alpha_path = os.path.join(self.temp_dir, 'alpha.csv')

        gene_infos = utils.read_pickle(self.umi_genes_path)
        df_counts = conversion.complement_counts(conversion.read_counts(self.umi_counts_path), gene_infos)
        _, _, pi_c = pi.read_pi(self.umi_alpha_pi_path, group_by=['barcode'])
        self.assertEqual(
            alpha_path,
            alpha.estimate_alpha(
                df_counts, pi_c, alpha_path, frozenset(['TC']), group_by=['barcode'], pi_c_group_by=['barcode']
            )
        )

        self.assertTrue(mixins.files_equal(self.umi_alpha_path, alpha_path))
