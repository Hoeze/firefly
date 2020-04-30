from pyspark.sql import functions as f


def variant_SPDI(
        contig="contigName",
        start="start",
        # end="end",
        ref="referenceAllele",
        alt="alternateAllele",
        join_alternate_alleles=False
):
    # if alt allele column is an array, concat with ','
    if join_alternate_alleles:
        alt = f.array_join(alt, ",")

    return f.concat(
        contig,
        f.lit(":"),
        start,
        f.lit(":"),
        ref,
        f.lit(":"),
        alt
    )