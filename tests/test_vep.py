import pytest

import yaml

import firefly.spark_functions.identifiers

with open("pytest_env.yml", "r") as fd:
    VEP_TEST_CFG = yaml.safe_load(fd)["vep_test_opts"]

import firefly
from firefly.spark_functions import displayHead

_GLOW_DEBUG_InVcfHeader = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO"""


@pytest.fixture
def vep_transformer():
    init_args = VEP_TEST_CFG["init_args"]
    predefined_vcf_annotation_args = VEP_TEST_CFG["predefined_vcf_annotation"]
    cadd_dir = VEP_TEST_CFG["cadd_dir"]

    vep_transformer = firefly.vep.VEPtransformer(**init_args)
    vep_transformer.add_sift()
    vep_transformer.add_polyphen()
    vep_transformer.add_predefined_vcf_annotation(**predefined_vcf_annotation_args)
    vep_transformer.add_cadd_plugin(cadd_dir=cadd_dir)
    vep_transformer.add_colocated_variants()

    return vep_transformer


def test_VEP_transformer_compile_cmd(vep_transformer):
    expected_cmd = VEP_TEST_CFG["expected_cmd"]

    assert vep_transformer.compile_cmd() == expected_cmd


def test_VEP_transformer_schema(vep_transformer):
    expected_transcript_consequences = {
        "cadd_phred",
        "cadd_raw",
        "sift_prediction",
        "sift_score",
        "polyphen_prediction",
        "polyphen_score",
    }
    transcript_consequence_fields = {
        *vep_transformer.get_output_struct_field("transcript_consequences").elementType.fieldNames()
    }
    assert all(f in transcript_consequence_fields for f in expected_transcript_consequences)

    expected_colocated_variant_fields = {'allele_string', 'end', 'id', 'seq_region_name', 'start', 'strand'}
    colocated_variants_fields = {
        *vep_transformer.get_output_struct_field("colocated_variants").elementType.fieldNames()
    }
    assert all(f in expected_colocated_variant_fields for f in colocated_variants_fields)

    expected_custom_annotations_fields = {'GnomADg_v3'}
    custom_annotations_fields = {
        *vep_transformer.get_output_struct_field("custom_annotations").fieldNames()
    }
    assert all(f in expected_custom_annotations_fields for f in custom_annotations_fields)


def test_run_VEP_formatting(
        glow_session,
        split_normalized_variants_sdf,
        example_variants_sdf,
        vep_transformer: firefly.vep.VEPtransformer
):
    # df = vep_transformer.apply_on(split_normalized_variants_sparkdf)
    import json
    import pyspark.sql.functions as f

    import glow
    input_df = split_normalized_variants_sdf.select([
        "contigName",
        "start",
        "end",
        f.array(firefly.spark_functions.identifiers.variant_SPDI(
            contig="contigName",
            start="start",
            ref="referenceAllele",
            alt="alternateAlleles",
            join_alternate_alleles=True,
        )).alias("names"),
        "referenceAllele",
        "alternateAlleles",
    ])

    vep_transformed_df = glow.transform(
        "pipe",
        df=input_df,
        cmd=json.dumps(["cat"]),
        inputFormatter='vcf',
        inVcfHeader="infer",
        outputFormatter='vcf',
    )

    vep_transformed_pdf = vep_transformed_df.toPandas()
    expexted_vep_transformed_pdf = split_normalized_variants_sdf.toPandas()

    assert vep_transformed_pdf.shape[0] == expexted_vep_transformed_pdf.shape[0]


def test_run_VEP(glow_session, split_normalized_variants_sdf, vep_transformer: firefly.vep.VEPtransformer):
    df = vep_transformer.transform(split_normalized_variants_sdf)

    print(displayHead(df).T)
