from typing import Dict, Union, List, Tuple
import re

import logging

log = logging.getLogger(__name__)

import pyspark.sql.types as t
import pyspark.sql.functions as f
import pyspark


def displayHead(df: pyspark.sql.DataFrame, nrows: int = 5):
    """
    returns the first `nrow` lines of a Spark dataframe as Pandas dataframe

    Args:
        df: The spark dataframe
        nrows: number of rows

    Returns: Pandas dataframe
    """
    return df.limit(nrows).toPandas()


def zip_explode_cols(
        df: pyspark.sql.dataframe.DataFrame,
        cols: list,
        result_name: str,
        rename_fields: Dict[str, str] = None
):
    """
    Explode multiple equally-sized arrays into one struct by zipping all arrays into one `ArrayType[StructType]`

    Args:
        df: The input Spark DataFrame
        cols: The array columns that should be zipped
        result_name: The name of the column that will contain the newly created struct
        rename_fields: dictionary mapping column names to new struct field names.
            Used to rename columns in the newly created struct.

    Returns: `df.withColumn(result_name, zip(explode(cols)))`

    """
    df = df.withColumn(result_name, f.explode(f.arrays_zip(*cols)))

    if rename_fields:  # create schema of new struct by simply renaming the top-level struct fields
        old_schema: t.StructType = df.schema[result_name].dataType

        # rename field if field ist in `old_schema.fieldNames()`
        new_field_names = [
            rename_fields[field] if field in rename_fields else field for field in old_schema.fieldNames()
        ]

        new_schema = t.StructType([
            t.StructField(name, field.dataType) for name, field in zip(new_field_names, old_schema.fields)
        ])

        df = df.withColumn(result_name, f.col(result_name).cast(new_schema))

        # # old method using withColumn and a new struct; breaks with PySpark 3.0
        # df = df.withColumn(target_struct, f.struct(*[
        #     f.col(target_struct + "." + actualName).alias(targetName)
        #     for targetName, actualName in zip(target_colnames, df.schema[target_struct].dataType.fieldNames())
        # ]))

    return df


def __rename_nested_field__(in_field: t.DataType, fieldname_normaliser):
    if isinstance(in_field, t.ArrayType):
        dtype = t.ArrayType(__rename_nested_field__(in_field.elementType, fieldname_normaliser), in_field.containsNull)
    elif isinstance(in_field, t.StructType):
        dtype = t.StructType()
        for field in in_field.fields:
            dtype.add(fieldname_normaliser(field.name), __rename_nested_field__(field.dataType, fieldname_normaliser))
    else:
        dtype = in_field
    return dtype


def __normalise_fieldname__(raw: str):
    return re.sub('[^A-Za-z0-9_]+', '.', raw.strip())


def __get_fields_info__(dtype: t.DataType, name: str = ""):
    ret = []
    if isinstance(dtype, t.StructType):
        for field in dtype.fields:
            for child in __get_fields_info__(field.dataType, field.name):
                wrapped_child = ["{prefix}{suffix}".format(
                    prefix=("" if name == "" else "`{}`.".format(name)), suffix=child[0])] + child[1:]
                ret.append(wrapped_child)
    elif isinstance(dtype, t.ArrayType) and (
            isinstance(dtype.elementType, t.ArrayType) or isinstance(dtype.elementType, t.StructType)):
        for child in __get_fields_info__(dtype.elementType):
            wrapped_child = ["`{}`".format(name)] + child
            ret.append(wrapped_child)
    else:
        return [["`{}`".format(name)]]
    return ret


def normalise_fields_names(df: pyspark.sql.DataFrame, fieldname_normaliser=__normalise_fieldname__):
    return df.select([
        f.col("`{}`".format(field.name)).cast(__rename_nested_field__(field.dataType, fieldname_normaliser))
            .alias(fieldname_normaliser(field.name)) for field in df.schema.fields
    ])


def flatten(df: pyspark.sql.DataFrame, fieldname_normaliser=__normalise_fieldname__):
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
        cols.append(f.expr(ex).alias(fieldname_normaliser("_".join(child).replace('`', ''))))
    return df.select(cols)
