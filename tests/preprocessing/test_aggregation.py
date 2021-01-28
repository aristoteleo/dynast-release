import os
from unittest import TestCase

import pandas as pd

import dynast.preprocessing.aggregation as aggregation
import dynast.preprocessing.conversion as conversion
import dynast.utils as utils
from tests import mixins


class TestAggregation(mixins.TestMixin, TestCase):

    def test_read_rates(self):
        aggregation.read_rates(self.umi_rates_path)

    def test_read_aggregates(self):
        aggregation.read_aggregates(self.umi_aggregates_paths['transcriptome']['TC'])

    def test_merge_aggregates(self):
        df1 = pd.read_csv(self.umi_aggregates_paths['spliced']['TC'])
        df2 = pd.read_csv(self.umi_aggregates_paths['unspliced']['TC'])
        aggregation.merge_aggregates(df1, df2, conversion='TC')

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
        self.assertTrue(mixins.files_equal(self.umi_rates_path, rates_path))

    def test_aggregate_counts(self):
        df = conversion.complement_counts(pd.read_csv(self.umi_counts_path), utils.read_pickle(self.umi_genes_path))
        self.assertEqual({c: os.path.join(self.temp_dir, f'{c}.csv')
                          for c in self.conversions},
                         aggregation.aggregate_counts(df[df['transcriptome']], self.temp_dir))
        for c in self.conversions:
            self.assertTrue(
                mixins.files_equal(
                    self.umi_aggregates_paths['transcriptome'][c], os.path.join(self.temp_dir, f'{c}.csv')
                )
            )
