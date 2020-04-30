import pytest

import pyspark

import glow

EXAMPLE_VCF_PATH = "resources/glow_example.vcf"
EXAMPLE_REF_GENOME = "resources/hg38_chr20.fa"


@pytest.fixture
def glow_session(spark_session) -> pyspark.sql.SparkSession:
    glow.register(spark_session)

    return spark_session


@pytest.fixture
def example_variants_sdf(glow_session: pyspark.sql.SparkSession):
    df = (
        glow_session
            .read
            .format('vcf')
            .option("includeSampleIds", True)
            # .option("splitToBiallelic", True)
            .load(EXAMPLE_VCF_PATH)
    )

    return df


@pytest.fixture
def split_normalized_variants_sdf(glow_session, example_variants_sdf):
    split_variants_df = glow.transform(
        "split_multiallelics",
        example_variants_sdf
    )

    normalized_variants_df = glow.transform(
        "normalize_variants",
        split_variants_df,
        reference_genome_path=EXAMPLE_REF_GENOME,
    ).repartition(2)

    return normalized_variants_df
