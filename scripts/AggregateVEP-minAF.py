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
import pyspark.sql.functions as f
import pyspark.sql.types as t

import glow
import os
# -

import json

# +
MEM = int(os.popen("ulimit -m").read())
N_CPU = int(os.popen("nproc").read())

print("memory: %d" % MEM)
print("number of cores: %d" % N_CPU)
# -

MEM = int(MEM * 0.8)

os.environ['PYSPARK_SUBMIT_ARGS'] = '--packages io.projectglow:glow_2.11:0.2.0,org.bdgenomics.adam:adam-assembly-spark2_2.12:0.29.0,io.delta:delta-core_2.11:0.4.0 --driver-memory %dk pyspark-shell' % MEM
os.environ['PYSPARK_SUBMIT_ARGS']

MAX_FAILURES=4

spark = (
    SparkSession.builder
    .appName('abc')
    #.config("spark.local.dir", os.environ.get("TMP"))
    .config("spark.local.dir", "/data/cephrbg/scratch/hoelzlwi")
    .config("spark.master", f"local[{N_CPU},{MAX_FAILURES}]")
    .config("spark.sql.execution.arrow.enabled", "true")
    #.config("spark.network.timeout", "1800s")
    #.config("spark.executor.heartbeatInterval", "600s")
    .config("spark.driver.maxResultSize", "48G")
#     .config("spark.databricks.io.cache.enabled", "true") # only enable when local storage is actually on local SSD
    .config("spark.task.maxFailures", MAX_FAILURES)
    .getOrCreate()
)
glow.register(spark)
spark

spark.sparkContext.getConf().get("spark.local.dir")

spark.sparkContext.getConf().get("spark.task.maxFailures")

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

## define custom display Head to simplify looking at spark DataFrames.
def displayHead(df, nrows = 5):
    return df.limit(nrows).toPandas()


# GTEX_PATH = '/s/project/variantDatabase/deltaLake/gtex_normalized'
GTEX_PATH = '/s/project/variantDatabase/deltaLake/gtex_normalized.bak'
GTEX_VEP_PATH = '/s/project/variantDatabase/deltaLake/GTEx_VEP'


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

# + {"active": ""}
# # TODO: remove
# gtex_df = (
#     gtex_df
#     .withColumn("call_summary_stats", f.expr("call_summary_stats(genotypes)"))
#     .withColumn("dp_summary_stats", f.expr("dp_summary_stats(genotypes)"))
#     .withColumn("gq_summary_stats", f.expr("gq_summary_stats(genotypes)"))
#     .drop("names")
# )

# +
# filter out low-quality variants
gtex_df = gtex_df.filter(f.array_contains(f.col("filters"), "PASS") | (f.col("qual") >= 99))

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
#min_heterozygous = 2
#min_homozygous = 2

## filter for variants that have either an exclusive heterozygous or an exclusive homozygous sample
#gtex_singlegt = gtex_singlegt.filter( 
#    f"(genotype.num_alt_alleles == 1 AND zygocity.n_hetero <= {min_heterozygous})"
#    + f" OR (genotype.num_alt_alleles == 2 AND zygocity.n_homo <= {min_homozygous})"
#)
# -

vep_df = (
    spark
    .read
    .format('delta')
    .load(GTEX_VEP_PATH)
    .withColumn("colocated_variants", f.col("colocated_variants").getItem(0))
    .withColumnRenamed("colocated_variants", "colocated_variant")
)

vep_df.printSchema()

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

exploded_vep = (
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

exploded_vep.printSchema()

transcript_consequences = exploded_vep.select(
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

gtex_with_consequences = gtex_singlegt.withColumn(
    "start", f.expr("start + 1"),
).join(
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
        # min of AF for each consequence
        f.max(
            f.when(f.col("consequences." + c) == True, (1 - f.col("INFO_AF"))).otherwise(0)
        ).alias(c) 
        for c in transcript_consequences.schema["consequences"].dataType.fieldNames()
    ]).alias("consequences"),
    f.struct(*[
        # min of AF for each impact level
        f.max(
            f.when(f.col("impact_level." + c) == True, (1 - f.col("INFO_AF"))).otherwise(0)
        ).alias(c) 
        for c in transcript_consequences.schema["impact_level"].dataType.fieldNames()
    ]).alias("impact_level"),
    f.struct(*[
        # min of AF for each LoF level
        f.max(
            f.when(f.col("lof_level." + c) == True, (1 - f.col("INFO_AF"))).otherwise(0)
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

counts.write.format("delta").save('/s/project/rep/processed/VEP/vep_counts_minAF.deltalake')

written_counts = (
    spark
    .read
    .format('delta')
    .load('/s/project/rep/processed/VEP/vep_counts_minAF.deltalake')
)
displayHead(written_counts)

flattened_counts = flatten(written_counts, lambda name: name)
flattened_counts.printSchema()

flattened_counts.write.format("delta").save('/s/project/rep/processed/VEP/vep_counts_minAF_flattened.deltalake')

displayHead(flattened_counts)


