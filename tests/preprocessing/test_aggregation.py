import os
from unittest import mock, TestCase

import pandas as pd

import dynast.preprocessing.aggregation as aggregation
from tests import mixins


class TestAggregation(mixins.TestMixin, TestCase):

    def test_read_rates(self):
        df = aggregation.read_rates(self.rates_barcode_path)
        self.assertEqual(['barcode'] + self.conversions, list(df.columns))

    def test_read_aggregates(self):
        df = aggregation.read_aggregates(self.aggregate_paths['TC'])
        self.assertEqual(['barcode', 'GX', 'TC', 'T', 'count'], list(df.columns))

    def test_calculate_mutation_rates_by_barcode(self):
        rates_path = os.path.join(self.temp_dir, 'rates.csv')
        df_counts = pd.read_csv(self.counts_path)

        self.assertEqual(rates_path, aggregation.calculate_mutation_rates(df_counts, rates_path, group_by='barcode'))
        self.assertTrue(mixins.files_equal(self.rates_barcode_path, rates_path))

    def test_calculate_mutation_rates_by_gene(self):
        rates_path = os.path.join(self.temp_dir, 'rates.csv')
        df_counts = pd.read_csv(self.counts_path)

        self.assertEqual(rates_path, aggregation.calculate_mutation_rates(df_counts, rates_path, group_by='GX'))
        self.assertTrue(mixins.files_equal(self.rates_gene_path, rates_path))

    def test_calculate_mutation_rates_by_none(self):
        rates_path = os.path.join(self.temp_dir, 'rates.csv')
        df_counts = pd.read_csv(self.counts_path)

        self.assertEqual(rates_path, aggregation.calculate_mutation_rates(df_counts, rates_path, group_by=None))
        self.assertTrue(mixins.files_equal(self.rates_none_path, rates_path))

    def test_aggregate_counts(self):
        df_counts = pd.read_csv(self.counts_path)
        df_genes = pd.read_csv(self.genes_path)

        with mock.patch('dynast.preprocessing.aggregation.tqdm', mixins.tqdm_mock):
            result = aggregation.aggregate_counts(df_counts, df_genes, self.temp_dir)
        self.assertEqual({
            conversion: os.path.join(self.temp_dir, f'{conversion}.csv')
            for conversion in self.conversions
        }, result)

        for conversion in self.conversions:
            self.assertTrue(mixins.files_equal(self.aggregate_paths[conversion], result[conversion]))
