# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.2.4
#   kernelspec:
#     display_name: Python [conda env:anaconda-florian3]
#     language: python
#     name: conda-env-anaconda-florian3-py
# ---

# +
import pyspark
from pyspark.sql import SparkSession
import glow
import os
import pandas as pd
#sc = pyspark.SparkContext()

import pyspark.sql.functions as f
# -

import pyranges as pr

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
    "--packages io.projectglow:glow_2.11:0.3.0,org.bdgenomics.adam:adam-assembly-spark2_2.12:0.29.0,io.delta:delta-core_2.11:0.5.0",
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


## define custom display Head to simplify looking at spark DataFrames.
def displayHead(df, nrows = 5):
    return df.limit(nrows).toPandas()


# GTEX_PATH = '/s/project/variantDatabase/gtex/v8/hg19/normalized.delta/'
# GTEX_PATH = '/s/project/variantDatabase/deltaLake/gtex_normalized'
GTEX_PATH = '/s/project/variantDatabase/deltaLake/gtex_normalized.bak'

gtf = pr.read_gtf("/s/project/variantDatabase/gtex/v8/hg19/reference/annotation/hg19.ensGene.gtf.gz")
gtf

transcripts = gtf.df.query("Feature == 'transcript'")
transcripts


# +
def promoter_region(start, end, strand, offset=1000):
    if strand == "+":
        return start - offset, start
    else:
        return end, end + offset

def apply_promoter_region(row):
    row.Start, row.End = promoter_region(row.Start, row.End, row.Strand)
    return row

promoter_regions = transcripts.copy()    
promoter_regions = promoter_regions.apply(apply_promoter_region, 1)
promoter_regions = promoter_regions.drop(["Score", "Frame", "exon_number", "exon_id", "Feature"], axis=1)
promoter_regions
# -

spark_promoter_regions = spark.createDataFrame(promoter_regions)


def zip_explode_cols(df: pyspark.sql.dataframe.DataFrame, cols: list, target_struct, target_colnames=None):
    if target_colnames is None:
        target_colnames = cols
    
    df = df.withColumn(target_struct, f.explode(f.arrays_zip(*cols)))
    df = df.withColumn(target_struct, f.struct(*[
        f.col(target_struct + "." + actualName).alias(targetName) 
        for targetName, actualName in zip(target_colnames, df.schema[target_struct].dataType.fieldNames())
    ]))
    
    return df


gtex_df = (
    spark
    .read
    .format('delta')
    .load(GTEX_PATH)
    .withColumn("num_alt_alleles", f.expr("genotype_states(genotypes)"))
    .drop("names")
)

gtex_df.printSchema()

# +
# filter out low-quality variants
gtex_df = gtex_df.filter(f.array_contains(f.col("filters"), "PASS") & (f.col("qual") >= 99))

# convert one-sized arrays to scalars
gtex_df = (
    gtex_df
    .withColumnRenamed("alternateAlleles", "alternateAllele").withColumn("alternateAllele", f.expr("alternateAllele[0]"))
    .withColumn("INFO_AC", f.expr("INFO_AC[0]"))
    .withColumn("INFO_AF", f.expr("INFO_AF[0]"))
    .withColumn("INFO_MLEAC", f.expr("INFO_MLEAC[0]"))
    .withColumn("INFO_MLEAF", f.expr("INFO_MLEAF[0]"))
)

# count unique elements
gtex_df = (
    gtex_df
    .withColumn(
        "zygocity",
        f.struct(
            f.expr('size(num_alt_alleles)').alias("n_samples"),
            f.expr('size(filter(num_alt_alleles, x -> x == 2))').alias("n_homo"),
            f.expr('size(filter(num_alt_alleles, x -> x == 1))').alias("n_hetero"),
            f.expr('size(filter(num_alt_alleles, x -> x == 0))').alias("n_ref"),
            f.expr('size(filter(num_alt_alleles, x -> x == -1))').alias("n_not_measured"),
            f.expr('array_distinct(filter(num_alt_alleles, x -> not x in (-1, 0, 1, 2)))').alias("other")
        )
    )
#     .withColumn("n_homozygous", f.expr('size(filter(num_alt_alleles, x -> x == 2))'))
#     .withColumn("n_heterozygous", f.expr('size(filter(num_alt_alleles, x -> x == 1))'))
#     .withColumn("n_ref", f.expr('size(filter(num_alt_alleles, x -> x == 0))'))
#     .withColumn("n_not_measured", f.expr('size(filter(num_alt_alleles, x -> x == -1))'))
#     .withColumn("n_not_measured", f.expr('size(filter(num_alt_alleles, x -> x == -1))'))
)
# -

gtex_df.printSchema()

x = displayHead(gtex_df)
x.T

gtex_singlegt = (
    zip_explode_cols(
        gtex_df,
        cols=[
            "genotypes",
            "num_alt_alleles",
        ],
        target_colnames=[
            "gt",
            "num_alt_alleles",
        ],
        target_struct="genotype"
    )
    .drop(
        "genotypes",
        "num_alt_alleles",
    )
)
gtex_singlegt.printSchema()

# +
min_heterozygous = 2
min_homozygous = 2

# filter for variants that have either an exclusive heterozygous or an exclusive homozygous sample
gtex_singlegt = gtex_singlegt.filter( 
    f"(genotype.num_alt_alleles == 1 AND zygocity.n_hetero <= {min_heterozygous})"
    + f" OR (genotype.num_alt_alleles == 2 AND zygocity.n_homo <= {min_homozygous})"
)
# -

# absolute paths:
CACHE_DIR="/s/project/rep/cache/"
RAW_DATA_DIR="/s/project/rep/raw/"
PROCESSED_DATA_DIR="/s/project/rep/processed/"
# per default, relative to PROCESSED_DATA_DIR:
MODEL_DIR="gtex/OUTRIDER"
# per default, relative to MODEL_DIR:
FIGURE_DIR="figures"

# +
MODEL_DIR=os.path.join(PROCESSED_DATA_DIR, MODEL_DIR)
FIGURE_DIR=os.path.join(MODEL_DIR, FIGURE_DIR)

if not os.path.exists(FIGURE_DIR):
    os.mkdir(FIGURE_DIR)

# +
import xarray as xr

xrds = xr.open_zarr(os.path.join(PROCESSED_DATA_DIR, "gtex/OUTRIDER/gtex_unstacked.zarr"))
xrds
# -

stacked_xrds = xrds.stack(observations=["subtissue", "individual", "gene"])
stacked_xrds

hilo_padj = stacked_xrds.hilo_padj.sel(observations= ~ stacked_xrds.missing)
hilo_padj

hilo_padj = hilo_padj.to_pandas()

hilo_padj

displayHead(vep_df)

# +
IMPACT_LEVELS = ["HIGH", "MODERATE", "LOW", "MODIFIER"]
LOF_LEVELS = ["HC", "LC"]

csq_terms = {}
for csq_type in ['intergenic_consequences','motif_feature_consequences', 'regulatory_feature_consequences', 'transcript_consequences']:
    unique_consequences = (
        vep_df
        .selectExpr(csq_type + ".consequence_terms as " + csq_type)
        .withColumn(csq_type, f.explode(csq_type))
        .withColumn(csq_type, f.explode(csq_type))
        .drop_duplicates()
    )
    # unique_consequences.printSchema()

    unique_consequences = unique_consequences.toPandas()
    csq_terms[csq_type] = unique_consequences

# -

csq_terms

exploded_df = (
    vep_df
    .withColumn("transcript_consequence", f.expr("explode(transcript_consequences)")).drop("transcript_consequences")
    .drop(
        "colocated_variants", 
        "intergenic_consequences",
        'most_severe_consequence',
        'motif_feature_consequences',
        'regulatory_feature_consequences',
    )
)

exploded_df.printSchema()

transcript_consequences = exploded_df.select(
    "*", 
    f.struct(*[
        f.array_contains("transcript_consequence.consequence_terms", name).cast("byte").alias(name)
        for name in csq_terms["transcript_consequences"].iloc[:,0].values
    ]).alias("consequences"),
    f.struct(*[
        (f.col("transcript_consequence.impact") == name).alias(name)
        for name in IMPACT_LEVELS
    ]).alias("impact_level"),
    f.struct(*[
        (f.col("transcript_consequence.lof") == name).alias(name)
        for name in LOF_LEVELS
    ]).alias("lof_level"),
)
transcript_consequences.printSchema()

# + {"active": ""}
# # unused
#
# def array_to_onehot_columns(df, input_col, output_col):
#     from pyspark.ml.feature import CountVectorizer
#     
#     vectorizer = CountVectorizer(inputCol=input_col, outputCol=output_col, binary=True).fit(df)
#     transformed_df = vectorizer.transform(df)
#     vocabulary = vectorizer.vocabulary
#     
#     udf_to_array = f.udf(lambda v: v.toArray().tolist(), 'array<double>')
#     
#     transformed_df = (
#         transformed_df
#         .withColumn(output_col, udf_to_array(output_col)) \
#         .withColumn(output_col, f.struct(*[ 
#             f.col(output_col)[i].astype('boolean').alias(vocabulary[i]) for i in range(len(vocabulary))
#         ]))
#     )
#     
#     return transformed_df
# -

gtex_with_consequences = gtex_singlegt.join(
    transcript_consequences.withColumnRenamed("seq_region_name", "contigName"),
    ["contigName", "start", "end"],
    how="inner"
)

gtex_with_consequences.printSchema()

grouped_df = gtex_with_consequences.groupby(
    "contigName", 
    "genotype.gt.sampleId", 
    "genotype.num_alt_alleles", 
    "transcript_consequence.gene_id"
)

# +
counts = grouped_df.agg(
    f.struct(*[
        # sum booleans as (0, 1) integer types
        f.sum(
            f.col("consequences." + c).cast("int")
        ).alias(c) 
        for c in transcript_consequences.schema["consequences"].dataType.fieldNames()
    ]).alias("consequences"),
    f.struct(*[
        # sum booleans as (0, 1) integer types
        f.sum(
            f.col("impact_level." + c).cast("int")
        ).alias(c) 
        for c in transcript_consequences.schema["impact_level"].dataType.fieldNames()
    ]).alias("impact_level"),
    f.struct(*[
        # sum booleans as (0, 1) integer types
        f.sum(
            f.col("lof_level." + c).cast("int")
        ).alias(c) 
        for c in transcript_consequences.schema["lof_level"].dataType.fieldNames()
    ]).alias("lof_level"),
    f.struct(*[
        f.min(f.col("transcript_consequence.sift_score")).alias("sift_min"),
        f.mean(f.col("transcript_consequence.sift_score")).alias("sift_mean"),
        f.max(f.col("transcript_consequence.sift_score")).alias("sift_max"),
        f.min(f.col("transcript_consequence.polyphen_score")).alias("polyphen_min"),
        f.mean(f.col("transcript_consequence.polyphen_score")).alias("polyphen_mean"),
        f.max(f.col("transcript_consequence.polyphen_score")).alias("polyphen_max"),
        f.min(f.col("transcript_consequence.cadd_raw")).alias("cadd_min"),
        f.mean(f.col("transcript_consequence.cadd_raw")).alias("cadd_mean"),
        f.max(f.col("transcript_consequence.cadd_raw")).alias("cadd_max"),
    ]).alias("scores"),
    
)
# -

counts.printSchema()

counts = counts.groupby("contigName", "gene_id", "sampleId").pivot("num_alt_alleles", [1, 2]).agg(
    f.struct([
        f.struct(*[
            f.sum(f.col("consequences." + c)).alias(c) 
            for c in counts.schema["consequences"].dataType.fieldNames()
        ]).alias("consequences"),
        f.struct(*[
            f.sum(f.col("impact_level." + c)).alias(c) 
            for c in counts.schema["impact_level"].dataType.fieldNames()
        ]).alias("impact_level"),
        f.struct(*[
            f.sum(f.col("lof_level." + c)).alias(c) 
            for c in counts.schema["lof_level"].dataType.fieldNames()
        ]).alias("lof_level"),
        f.struct(*[
            f.sum(f.col("scores." + c)).alias(c) 
            for c in counts.schema["scores"].dataType.fieldNames()
        ]).alias("scores"),
    ])
)

counts = counts.withColumnRenamed("1", "heterozygous").withColumnRenamed("2", "homozygous")

counts.printSchema()

counts.write.format("delta").save('/s/project/rep/processed/VEP/vep_counts2.deltalake')

written_counts = (
    spark
    .read
    .format('delta')
    .load('/s/project/rep/processed/VEP/vep_counts2.deltalake')
)
displayHead(written_counts)

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
	return re.sub('[^A-Za-z0-9_]+', '_', raw.strip())


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

flattened_counts = flatten(written_counts, lambda name: name)
flattened_counts.printSchema()

flattened_counts.write.format("delta").save('/s/project/rep/processed/VEP/vep_counts_flattened2.deltalake')

displayHead(flattened_counts)


