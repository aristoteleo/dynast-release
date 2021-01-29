import os
from unittest import TestCase

import dynast.estimation.p_e as p_e
import dynast.preprocessing.aggregation as aggregation
import dynast.preprocessing.conversion as conversion
import dynast.utils as utils

from .. import mixins


class TestPE(mixins.TestMixin, TestCase):

    def test_read_p_e(self):
        p_e.read_p_e(self.umi_p_e_path, group_by=['barcode'])

    def test_estimate_p_e_control(self):
        p_e_path = os.path.join(self.temp_dir, 'p_e.csv')
        self.assertEqual(
            p_e_path,
            p_e.estimate_p_e_control(
                conversion.complement_counts(
                    conversion.read_counts(self.control_counts_path),
                    utils.read_pickle(self.control_genes_path),
                ),
                p_e_path,
                conversion='TC'
            )
        )
        self.assertAlmostEqual(p_e.read_p_e(self.control_p_e_path), p_e.read_p_e(p_e_path))

    def test_estimate_p_e(self):
        p_e_path = os.path.join(self.temp_dir, 'p_e.csv')
        self.assertEqual(
            p_e_path,
            p_e.estimate_p_e(
                conversion.complement_counts(
                    conversion.read_counts(self.umi_counts_path),
                    utils.read_pickle(self.umi_genes_path),
                ),
                p_e_path,
                conversion='TC',
                group_by=['barcode']
            )
        )
        p_e_test = p_e.read_p_e(p_e_path, group_by=['barcode'])
        for barcode, value in p_e.read_p_e(self.umi_p_e_path, group_by=['barcode']).items():
            self.assertAlmostEqual(value, p_e_test[barcode])

    def test_estimate_p_e_nasc(self):
        p_e_path = os.path.join(self.temp_dir, 'p_e.csv')
        self.assertEqual(
            p_e_path,
            p_e.estimate_p_e_nasc(
                aggregation.read_rates(self.nasc_rates_path), p_e_path, conversion='TC', group_by=['barcode']
            )
        )
        p_e_test = p_e.read_p_e(p_e_path, group_by=['barcode'])
        for barcode, value in p_e.read_p_e(self.nasc_p_e_path, group_by=['barcode']).items():
            self.assertAlmostEqual(value, p_e_test[barcode])
