# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.1
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
# import databricks.koalas as ks
#sc = pyspark.SparkContext()

import pyspark.sql.functions as f

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


def zip_explode_cols(df: pyspark.sql.dataframe.DataFrame, cols: list, target_struct, target_colnames=None):
    if target_colnames is None:
        target_colnames = cols
    
    df = df.withColumn(target_struct, f.explode(f.arrays_zip(*cols)))
    df = df.withColumn(target_struct, f.struct(*[
        f.col(target_struct + "." + actualName).alias(targetName) 
        for targetName, actualName in zip(target_colnames, df.schema[target_struct].dataType.fieldNames())
    ]))
    
    return df


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

MAX_FAILURES=4

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
    .config("spark.task.maxFailures", MAX_FAILURES)
    .config("spark.driver.extraJavaOptions=-verbose:class")
    .getOrCreate()
)
glow.register(spark)
spark

spark.sparkContext.getConf().get("spark.local.dir")

spark.sparkContext.getConf().get("spark.task.maxFailures")

input_path  = '/s/raw/gnomad/3.0/hg38/gnomad.genomes.r3.0.sites.vcf.gz'
output_path = '/s/project/variantDatabase/gnomad/3.0/hg38/normalized.delta/'
ref_genome  = '/s/project/variantDatabase/gnomad/3.0/hg38/dna.fa' # contigName = chr<#>

df = (
    spark
    .read
    .format('vcf')
    .option("includeSampleIds", True)
    .load(input_path)
)

displayHead(df)

# +
normalized_variants_df = glow.transform(
  "normalize_variants",
  df,
  reference_genome_path=ref_genome,
  mode='split_and_normalize'
)

normalized_variants_df = (
    normalized_variants_df
#    .repartition(10000, shuffle=True)
#    .persist()
#     .orderBy("contigName", "start", "end", "referenceAllele", "alternateAlleles")
#    .repartitionByRange(10000, "contigName", "start")
#    .sortWithinPartitions("contigName", "start", "end", "referenceAllele", "alternateAlleles")
#     .withColumn("call_summary_stats", f.expr("call_summary_stats(genotypes)"))
#     .withColumn("dp_summary_stats", f.expr("dp_summary_stats(genotypes)"))
#     .withColumn("gq_summary_stats", f.expr("gq_summary_stats(genotypes)"))
)
# -

displayHead(normalized_variants_df)

(
    normalized_variants_df
    .write
    .format("delta")
    .save(output_path)
)

# ## Check if data can be loaded

displayHead(
    spark
    .read
    .format("delta")
    .load(output_path)
)


