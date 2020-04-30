# ---
# jupyter:
#   jupytext:
#     cell_metadata_json: true
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python [conda env:3-2018.12-glow-test]
#     language: python
#     name: conda-env-3-2018.12-glow-test-py
# ---

# +
import pyspark
from pyspark.sql import SparkSession
import glow
import os
import pandas as pd
#sc = pyspark.SparkContext()

import pyspark.sql.types as t
import pyspark.sql.functions as f
# -

os.environ.get("ASDF", "")

import json

# +
from pyspark.sql import functions as F
from pyspark.sql.types import DataType, StructType, ArrayType
from pyspark.sql import DataFrame
import re


def __rename_nested_field__(in_field: DataType, fieldname_normaliser):
	if isinstance(in_field, ArrayType):
		dtype = ArrayType(__rename_nested_field__(in_field.elementType, fieldname_normaliser), in_field.containsNull)
	elif isinstance(in_field, StructType):
		dtype = StructType()
		for field in in_field.fields:
			dtype.add(fieldname_normaliser(field.name), __rename_nested_field__(field.dataType, fieldname_normaliser))
	else:
		dtype = in_field
	return dtype


def __normalise_fieldname__(raw: str):
	return re.sub('[^A-Za-z0-9_]+', '.', raw.strip())


def __get_fields_info__(dtype: DataType, name: str = ""):
	ret = []
	if isinstance(dtype, StructType):
		for field in dtype.fields:
			for child in __get_fields_info__(field.dataType, field.name):
				wrapped_child = ["{prefix}{suffix}".format(
					prefix=("" if name == "" else "`{}`.".format(name)), suffix=child[0])] + child[1:]
				ret.append(wrapped_child)
	elif isinstance(dtype, ArrayType) and (
			isinstance(dtype.elementType, ArrayType) or isinstance(dtype.elementType, StructType)):
		for child in __get_fields_info__(dtype.elementType):
			wrapped_child = ["`{}`".format(name)] + child
			ret.append(wrapped_child)
	else:
		return [["`{}`".format(name)]]
	return ret


def normalise_fields_names(df: DataFrame, fieldname_normaliser=__normalise_fieldname__):
	return df.select([
		F.col("`{}`".format(field.name)).cast(__rename_nested_field__(field.dataType, fieldname_normaliser))
			.alias(fieldname_normaliser(field.name)) for field in df.schema.fields
	])


def flatten(df: DataFrame, fieldname_normaliser=__normalise_fieldname__):
	cols = []
	for child in __get_fields_info__(df.schema):
		if len(child) > 2:
			ex = "x.{}".format(child[-1])
			for seg in child[-2:0:-1]:
				if seg != '``':
					ex = "transform(x.{outer}, x -> {inner})".format(outer=seg, inner=ex)
			ex = "transform({outer}, x -> {inner})".format(outer=child[0], inner=ex)
		else:
			ex = ".".join(child)
		cols.append(F.expr(ex).alias(fieldname_normaliser("_".join(child).replace('`', ''))))
	return df.select(cols)


# -

## define custom display Head to simplify looking at spark DataFrames.
def displayHead(df, nrows = 5):
    return df.limit(nrows).toPandas()


# +
MEM = os.popen("ulimit -m").read()
if MEM.startswith("unlimited"):
    print("Memory not constrained, using all available memory...")
    import psutil
    MEM = psutil.virtual_memory().available / 1024
MEM = int(MEM)

N_CPU = int(os.popen("nproc").read())

print("memory: %dk" % MEM)
print("number of cores: %d" % N_CPU)

# +
# MEM = MEM // 2
# N_CPU = N_CPU // 2
# -

os.environ['PYSPARK_SUBMIT_ARGS'] = " ".join([
    "--packages io.projectglow:glow_2.11:0.3.1-SNAPSHOT,io.delta:delta-core_2.11/0.5.0",
    "--driver-memory %sk" % MEM,
    "--executor-memory %sk" % MEM,
    "pyspark-shell"
])
os.environ['PYSPARK_SUBMIT_ARGS']

spark = (
    SparkSession.builder
    .appName('abc')
    #.config("spark.local.dir", os.environ.get("TMP"))
    .config("spark.local.dir", "/data/cephrbgssd/scratch/hoelzlwi") # TODO: change this back to real scratch partition
    .config("spark.sql.execution.arrow.enabled", "true")
    .config("spark.sql.shuffle.partitions", "2001")
#     .config("spark.network.timeout", "1800s")
#     .config("spark.executor.heartbeatInterval", "600s")
    .config("spark.databricks.io.cache.enabled", "true")
    .config("spark.databricks.io.cache.maxDiskUsage", "500G")
#     .config("spark.driver.maxResultSize", "48G")
    .config("spark.driver.maxResultSize", "%dk" % MEM)
    .getOrCreate()
)
glow.register(spark)
spark

INPUT_PATH  = '/s/project/variantDatabase/gtex/v8/hg38/normalized.delta/'
OUTPUT_PATH = '/s/project/variantDatabase/gtex/v8/hg38/vep.delta/'

# +
VEP_CACHE_DIR = "/opt/modules/i12g/conda-ensembl-vep/99/cachedir"
VEP_PLUGIN_DIR = "/opt/modules/i12g/conda-ensembl-vep/99/cachedir/Plugins"

# GNOMAD_VCF_FILE = "/s/raw/ensembl/pub/data_files/homo_sapiens/GRCh37/variation_genotype/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz"
# GNOMAD_VCF_FILE = "/s/raw/ensembl/pub/data_files/homo_sapiens/GRCh38/variation_genotype/gnomad.genomes.r2.0.1.sites.GRCh38.noVEP.vcf.gz"
# GNOMAD_FIELDS = [
#     "AF",
#     "AF_AFR",
#     "AF_AMR",
#     "AF_ASJ",
#     "AF_EAS",
#     "AF_FIN",
#     "AF_NFE",
#     "AF_OTH",
# ]
GNOMAD_VCF_FILE = "/s/raw/gnomad/3.0/hg38/gnomad.genomes.r3.0.sites.vcf.gz"
GNOMAD_AF_FIELDS = [
    "AF",
    "AF_asj_female",
    "AF_eas_female",
    "AF_afr_male",
    "AF_female",
    "AF_fin_male",
    "AF_oth_female",
    "AF_ami",
    "AF_oth",
    "AF_male",
    "AF_ami_female",
    "AF_afr",
    "AF_eas_male",
    "AF_sas",
    "AF_nfe_female",
    "AF_asj_male",
    "AF_raw",
    "AF_oth_male",
    "AF_nfe_male",
    "AF_asj",
    "AF_amr_male",
    "AF_amr_female",
    "AF_sas_female",
    "AF_fin",
    "AF_afr_female",
    "AF_sas_male",
    "AF_amr",
    "AF_nfe",
    "AF_eas",
    "AF_ami_male",
    "AF_fin_female"
]

# CADD_DIR = "/s/genomes/human/hg19/CADD/v1.3"
CADD_DIR = "/s/raw/cadd/v1.5"

ASSEMBLY = "GRCh38"
# -

df = (
    spark
    .read
    .format('delta')
    .load(INPUT_PATH)
)

input_df = df.select(["contigName", "start", "end", "referenceAllele", "alternateAlleles"])

input_df.printSchema()

displayHead(input_df)

cmd = [
#    "/data/ouga04b/ag_gagneur/home/hoelzlwi/.local/bin/vep.sh",
   "vep",
#     "--no_check_variants_order",
    "--dir_cache", VEP_CACHE_DIR,
    "--dir_plugin", VEP_PLUGIN_DIR,
    "--plugin", f"CADD,{CADD_DIR}/whole_genome_SNVs.tsv.gz,{CADD_DIR}/InDels.tsv.gz",
    "--format", "vcf",
    "--output_file", "STDOUT",
    #   "--fasta", "/mnt/dbnucleus/dbgenomics/grch37_merged_vep_96/data/human_g1k_v37.fa",
    "--allele_number",
    "--assembly", ASSEMBLY,
    "--port", "3337",
    "--sift", "b",
    "--polyphen", "s",
    "--total_length", # Give cDNA, CDS and protein positions as Position/Length
    "--numbers", # Adds affected exon and intron numbering to to output. Format is Number/Total
    "--symbol", # Adds gene symbol
#     "--hgvs",
    "--ccds", # Adds the CCDS transcript identifer (where available) to the output
    "--xref_refseq", # Output aligned RefSeq mRNA identifier for transcript
    "--uniprot", # Adds best match accessions for translated protein products from three UniProt-related databases (SWISSPROT, TREMBL and UniParc) to the output
    "--af", # Add the global allele frequency (AF) from 1000 Genomes Phase 3 data for any known co-located variant to the output
    "--af_gnomad", # Include allele frequency from Genome Aggregation Database (gnomAD) exome populations
    "--max_af", # Report the highest allele frequency observed in any population from 1000 genomes, ESP or gnomAD
    "--pubmed", # Report Pubmed IDs for publications that cite existing variant
    "--canonical", # Adds a flag indicating if the transcript is the canonical transcript for the gene
    "--biotype", # Adds the biotype of the transcript or regulatory feature
    "--no_stats",
    "--cache",
    "--custom", f"{GNOMAD_VCF_FILE},gnomADg,vcf,exact,0," + ",".join(GNOMAD_AF_FIELDS),
    "--offline",
    "--json",
    "--merged"
]
" ".join(cmd)

# + {"active": ""}
# # just for debugging purposes
#
# output_df = glow.transform(
#     "pipe",
#     input_df,
#     cmd=json.dumps(["cat"]),
#     inputFormatter='vcf',
#     inVcfHeader='infer',
#     outputFormatter='text',
# )
# -

output_df = glow.transform(
    "pipe",
    input_df, 
    cmd=json.dumps(cmd), 
    inputFormatter='vcf', 
    inVcfHeader='infer', 
    outputFormatter='text',
)

print(displayHead(output_df).text[0])

vep_schema = t.StructType([
    t.StructField("seq_region_name", t.StringType()),
    t.StructField("start", t.IntegerType()),
    t.StructField("end", t.IntegerType()),
    t.StructField("id", t.StringType()),
    t.StructField("strand", t.IntegerType()),
    t.StructField("assembly_name", t.StringType()),
    t.StructField("allele_string", t.StringType()),
    t.StructField("ancestral", t.StringType()),
    t.StructField("variant_class", t.StringType()),
    t.StructField("custom_annotations", t.StructType([
        t.StructField("gnomADg", t.ArrayType(t.StructType([
            t.StructField("name", t.StringType()),
            t.StructField("allele", t.StringType()),
            t.StructField("fields", t.StructType([
                t.StructField(f, t.FloatType()) for f in GNOMAD_AF_FIELDS
            ])),
        ])))
    ])),
#     t.StructField("colocated_variants", t.ArrayType(t.StructType([
#         t.StructField("aa_allele", t.StringType()),
#         t.StructField("aa_maf", t.FloatType()),
#         t.StructField("afr_allele", t.StringType()),
#         t.StructField("afr_maf", t.FloatType()),
#         t.StructField("allele_string", t.StringType()),
#         t.StructField("amr_allele", t.StringType()),
#         t.StructField("amr_maf", t.FloatType()),
#         t.StructField("clin_sig", t.ArrayType(t.StringType())),
#         t.StructField("end", t.IntegerType()),
#         t.StructField("eas_allele", t.StringType()),
#         t.StructField("eas_maf", t.FloatType()),
#         t.StructField("ea_allele", t.StringType()),
#         t.StructField("ea_maf", t.FloatType()),
#         t.StructField("eur_allele", t.StringType()),
#         t.StructField("eur_maf", t.FloatType()),
#         t.StructField("exac_adj_allele", t.StringType()),
#         t.StructField("exac_adj_maf", t.FloatType()),
#         t.StructField("exac_allele", t.StringType()),
#         t.StructField("exac_afr_allele", t.StringType()),
#         t.StructField("exac_afr_maf", t.FloatType()),
#         t.StructField("exac_amr_allele", t.StringType()),
#         t.StructField("exac_amr_maf", t.FloatType()),
#         t.StructField("exac_eas_allele", t.StringType()),
#         t.StructField("exac_eas_maf", t.FloatType()),
#         t.StructField("exac_fin_allele", t.StringType()),
#         t.StructField("exac_fin_maf", t.FloatType()),
#         t.StructField("exac_maf", t.FloatType()),
#         t.StructField("exac_nfe_allele", t.StringType()),
#         t.StructField("exac_nfe_maf", t.FloatType()),
#         t.StructField("exac_oth_allele", t.StringType()),
#         t.StructField("exac_oth_maf", t.FloatType()),
#         t.StructField("exac_sas_allele", t.StringType()),
#         t.StructField("exac_sas_maf", t.FloatType()),
#         t.StructField("id", t.StringType()),
#         t.StructField("minor_allele", t.StringType()),
#         t.StructField("minor_allele_freq", t.FloatType()),
#         t.StructField("phenotype_or_disease", t.IntegerType()),
#         t.StructField("pubmed", t.ArrayType(t.IntegerType())),
#         t.StructField("sas_allele", t.StringType()),
#         t.StructField("sas_maf", t.FloatType()),
#         t.StructField("somatic", t.IntegerType()),
#         t.StructField("start", t.IntegerType()),
#         t.StructField("strand", t.IntegerType())
#     ]))),
    t.StructField("context", t.StringType()),
    t.StructField("input", t.StringType()),
    t.StructField("intergenic_consequences", t.ArrayType(t.StructType([
            t.StructField("allele_num", t.IntegerType()),
            t.StructField("consequence_terms", t.ArrayType(t.StringType())),
            t.StructField("impact", t.StringType()),
            t.StructField("minimised", t.IntegerType()),
            t.StructField("variant_allele", t.StringType())
    ]))),
    t.StructField("most_severe_consequence", t.StringType()),
    t.StructField("motif_feature_consequences", t.ArrayType(t.StructType([
            t.StructField("allele_num", t.IntegerType()),
            t.StructField("consequence_terms", t.ArrayType(t.StringType())),
            t.StructField("high_inf_pos", t.StringType()),
            t.StructField("impact", t.StringType()),
            t.StructField("minimised", t.IntegerType()),
            t.StructField("motif_feature_id", t.StringType()),
            t.StructField("motif_name", t.StringType()),
            t.StructField("motif_pos", t.IntegerType()),
            t.StructField("motif_score_change", t.FloatType()),
            t.StructField("strand", t.IntegerType()),
            t.StructField("variant_allele", t.StringType())
    ]))),
    t.StructField("regulatory_feature_consequences", t.ArrayType(t.StructType([
            t.StructField("allele_num", t.IntegerType()),
            t.StructField("biotype", t.StringType()),
            t.StructField("consequence_terms", t.ArrayType(t.StringType())),
            t.StructField("impact", t.StringType()),
            t.StructField("minimised", t.IntegerType()),
            t.StructField("regulatory_feature_id", t.StringType()),
            t.StructField("variant_allele", t.StringType())
    ]))),
    t.StructField("transcript_consequences", t.ArrayType(t.StructType([
            t.StructField("allele_num", t.IntegerType()),
            t.StructField("amino_acids", t.StringType()),
            t.StructField("appris", t.StringType()),
            t.StructField("biotype", t.StringType()),
            t.StructField("cadd_phred", t.FloatType()),
            t.StructField("cadd_raw", t.FloatType()),
            t.StructField("canonical", t.IntegerType()),
            t.StructField("ccds", t.StringType()),
            t.StructField("cdna_start", t.IntegerType()),
            t.StructField("cdna_end", t.IntegerType()),
            t.StructField("cds_end", t.IntegerType()),
            t.StructField("cds_start", t.IntegerType()),
            t.StructField("codons", t.StringType()),
            t.StructField("consequence_terms", t.ArrayType(t.StringType())),
            t.StructField("distance", t.IntegerType()),
            t.StructField("domains", t.ArrayType(t.StructType([
                t.StructField("db", t.StringType()),
                t.StructField("name", t.StringType())
            ]))),
            t.StructField("exon", t.StringType()),
            t.StructField("gene_id", t.StringType()),
            t.StructField("gene_pheno", t.IntegerType()),
            t.StructField("gene_symbol", t.StringType()),
            t.StructField("gene_symbol_source", t.StringType()),
            t.StructField("hgnc_id", t.StringType()),
            t.StructField("hgvsc", t.StringType()),
            t.StructField("hgvsp", t.StringType()),
            t.StructField("hgvs_offset", t.IntegerType()),
            t.StructField("impact", t.StringType()),
            t.StructField("intron", t.StringType()),
            t.StructField("lof", t.StringType()),
            t.StructField("lof_flags", t.StringType()),
            t.StructField("lof_filter", t.StringType()),
            t.StructField("lof_info", t.StringType()),
            t.StructField("minimised", t.IntegerType()),
            t.StructField("polyphen_prediction", t.StringType()),
            t.StructField("polyphen_score", t.FloatType()),
            t.StructField("protein_end", t.IntegerType()),
            t.StructField("protein_start", t.IntegerType()),
            t.StructField("protein_id", t.StringType()),
            t.StructField("refseq_transcript_ids", t.ArrayType(t.StringType())),
            t.StructField("sift_prediction", t.StringType()),
            t.StructField("sift_score", t.FloatType()),
            t.StructField("source", t.StringType()),
            t.StructField("strand", t.IntegerType()),
            t.StructField("swissprot", t.StringType()),
            t.StructField("transcript_id", t.StringType()),
            t.StructField("trembl", t.StringType()),
            t.StructField("tsl", t.IntegerType()),
            t.StructField("uniparc", t.StringType()),
            t.StructField("used_ref", t.StringType()),
            t.StructField("variant_allele", t.StringType())
    ]))),
])


# +
# vep_schema.jsonValue()

# +
# t.StructType.fromJson(vep_schema.jsonValue())

# +
#schema = spark.range(1).select(f.schema_of_json(x["text"]).alias("json")).collect()[0].json
#schema
# -

vep_df = output_df.withColumn('data', f.from_json('text', vep_schema)).select(f.col("data.*")).drop("input")
vep_df.printSchema()

vep_df.write.format("delta").save(OUTPUT_PATH)

check_vep_df = (
    spark
    .read
    .format('delta')
    .load(OUTPUT_PATH)
)

with pd.option_context('display.max_rows', None, 'display.max_columns', None):
    display(displayHead(flatten(check_vep_df)).T)


