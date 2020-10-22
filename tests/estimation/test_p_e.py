import os
from unittest import mock, TestCase

import pandas as pd

import dynast.estimation.p_e as p_e
from tests import mixins


class TestPE(mixins.TestMixin, TestCase):

    def test_read_p_e(self):
        estimate = p_e.read_p_e(self.p_e_barcode_path, group_by='barcode')
        self.assertEqual(estimate['AAGATGTCCAGG'], 0.0625)

    def test_estimate_p_e_by_barcode(self):
        df = pd.read_csv(self.rates_barcode_path)
        p_e_path = os.path.join(self.temp_dir, 'p_e.csv')
        self.assertEqual((mock.ANY, p_e_path), p_e.estimate_p_e(df, p_e_path, group_by=['barcode']))
        self.assertTrue(mixins.files_equal(self.p_e_barcode_path, p_e_path))

    def test_estimate_p_e_by_gene(self):
        df = pd.read_csv(self.rates_gene_path)
        p_e_path = os.path.join(self.temp_dir, 'p_e.csv')
        self.assertEqual((mock.ANY, p_e_path), p_e.estimate_p_e(df, p_e_path, group_by=['GX']))
        self.assertTrue(mixins.files_equal(self.p_e_gene_path, p_e_path))

    def test_estimate_p_e_by_barcode_gene(self):
        df = pd.read_csv(self.rates_barcode_gene_path)
        p_e_path = os.path.join(self.temp_dir, 'p_e.csv')
        self.assertEqual((mock.ANY, p_e_path), p_e.estimate_p_e(df, p_e_path, group_by=['barcode', 'GX']))
        self.assertTrue(mixins.files_equal(self.p_e_barcode_gene_path, p_e_path))

    def test_estimate_p_e_by_none(self):
        df = pd.read_csv(self.rates_none_path)
        p_e_path = os.path.join(self.temp_dir, 'p_e.csv')
        self.assertEqual((mock.ANY, p_e_path), p_e.estimate_p_e(df, p_e_path, group_by=None))
        self.assertTrue(mixins.files_equal(self.p_e_none_path, p_e_path))
