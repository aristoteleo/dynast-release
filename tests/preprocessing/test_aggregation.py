import os
from unittest import TestCase

import pandas as pd

import dynast.preprocessing.aggregation as aggregation
import dynast.preprocessing.conversion as conversion
import dynast.utils as utils

from .. import mixins


class TestAggregation(mixins.TestMixin, TestCase):

    def test_read_rates(self):
        aggregation.read_rates(self.umi_rates_path)

    def test_read_aggregates(self):
        aggregation.read_aggregates(self.umi_aggregates_path)

    def test_merge_aggregates(self):
        df1 = pd.read_csv(self.umi_aggregates_path)
        df2 = pd.read_csv(self.umi_aggregates_path)
        aggregation.merge_aggregates(df1, df2)

    def test_calculate_mutation_rates(self):
        rates_path = os.path.join(self.temp_dir, 'rates.csv')
        self.assertEqual(
            rates_path,
            aggregation.calculate_mutation_rates(
                conversion.complement_counts(pd.read_csv(self.umi_counts_path), utils.read_pickle(self.umi_genes_path)),
                rates_path,
                group_by='barcode'
            )
        )
        self.assertTrue(mixins.dataframes_equal(self.umi_rates_path, rates_path, sort_by='barcode'))

    def test_aggregate_counts(self):
        df = conversion.complement_counts(pd.read_csv(self.umi_counts_path), utils.read_pickle(self.umi_genes_path))
        aggregates_path = os.path.join(self.temp_dir, 'aggregates.csv')
        self.assertEqual(aggregates_path, aggregation.aggregate_counts(df, aggregates_path, conversions=['TC']))
        self.assertTrue(mixins.dataframes_equal(self.umi_aggregates_path, aggregates_path, sort_by=['barcode', 'GX']))
