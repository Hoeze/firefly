import pandas as pd
import numpy as np

import pyspark.sql.functions as f

import firefly
import firefly.spark_functions.identifiers
from firefly.spark_functions import displayHead, flatten


def test_read_vcf(example_variants_sdf):
    print(displayHead(flatten(example_variants_sdf)).T)


def test_split_normalize(split_normalized_variants_sdf):
    print(displayHead(flatten(split_normalized_variants_sdf)).T)


def test_variant_identifiers_biallelic(split_normalized_variants_sdf):
    variant_ids_sdf = split_normalized_variants_sdf.select(
        firefly.spark_functions.identifiers.variant_SPDI(
            contig=f.col("contigName"),
            start=f.col("start"),
            ref=f.col("referenceAllele"),
            alt=f.col("alternateAlleles")[0],
        ).alias("variant_ids")
    )

    variant_ids_pdf = variant_ids_sdf.toPandas()
    expected_variant_ids_pdf = pd.read_csv("resources/glow_example_SPDI_biallelic.csv")

    assert np.all(np.sort(expected_variant_ids_pdf.values.flatten()) == np.sort(variant_ids_pdf.values.flatten()))


def test_variant_identifiers_multiallelic(example_variants_sdf):
    variant_ids_sdf = example_variants_sdf.select(
        firefly.spark_functions.identifiers.variant_SPDI(
            contig=f.col("contigName"),
            start=f.col("start"),
            ref=f.col("referenceAllele"),
            alt=f.col("alternateAlleles"),
            join_alternate_alleles=True
        ).alias("variant_ids")
    )

    variant_ids_pdf = variant_ids_sdf.toPandas()
    expected_variant_ids_pdf = pd.read_csv("resources/glow_example_SPDI_multiallelic.csv")

    assert np.all(np.sort(expected_variant_ids_pdf.values.flatten()) == np.sort(variant_ids_pdf.values.flatten()))
