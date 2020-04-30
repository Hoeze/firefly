import pyspark

import firefly


def test_spark_context(spark_session):
    print(spark_session)


def test_zip_explode_cols(spark_session):
    from pyspark.sql.functions import split

    df = (
        spark_session
            .createDataFrame([('abc, def, ghi', '1.0, 2.0, 3.0')])
            .toDF("column_1", "column_2")
            .withColumn("column_1", split("column_1", " *, *"))
            .withColumn("column_2", split("column_2", " *, *"))
    )

    df_zip_exploded = firefly.spark_functions.zip_explode_cols(
        df,
        cols=["column_1", "column_2"],
        result_name="zip_exploded",
        rename_fields={"column_1": "column_renamed_1"}
    )

    expected_schema = {
        'type': 'struct',
        'fields': [
            {
                'name': 'column_1',
                'type': {'type': 'array', 'elementType': 'string', 'containsNull': True},
                'nullable': True,
                'metadata': {}
            },
            {
                'name': 'column_2',
                'type': {'type': 'array', 'elementType': 'string', 'containsNull': True},
                'nullable': True,
                'metadata': {}
            }
        ]
    }

    print(df_zip_exploded.printSchema())
    print(df_zip_exploded.toPandas())

    assert df.schema.jsonValue() == expected_schema
