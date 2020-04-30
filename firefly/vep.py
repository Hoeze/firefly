from typing import Union, List, Tuple

import shlex
import json
# import collections
# import copy

import logging

import firefly.spark_functions.identifiers

log = logging.getLogger(__name__)

import pyspark
import pyspark.sql.types as t
import pyspark.sql.functions as f

import firefly.resources

# VEP_FORMAT_OPTIONS = {
#     "ensembl", "vcf", "hgvs", "id", "region", "spdi"
# }

# Dict["vep_version", Dict["json_schema"]]
VEP_BASIC_SCHEMA = firefly.resources.get("vep/basic_schema")

# Dict["vep_version", List["args"]]
VEP_BASIC_ARGS = firefly.resources.get("vep/basic_args")

# Dict["dataset", Dict["args", "dtype"]]
VEP_CUSTOM_ANNOT = firefly.resources.get("vep/predefined_custom_annot")


# def _get_field_in_list(list_of_dicts: List[Dict[str, any]], key, value):
#     for f in list_of_dicts:
#         if f.get(key) == key:
#             return f
#     return None


class VEPtransformer:
    call_args: List[str]
    output_schema: t.StructType
    genome_assembly: str

    def __init__(
            self,
            vep_command="vep",
            genome_assembly="GRCh38",
            offline=True,
            use_cache=True,
            vep_cache_dir=None,
            vep_plugin_dir=None,
            vep_version="version_99",
    ):
        self.call_args = [
            *shlex.split(vep_command),
            *VEP_BASIC_ARGS[vep_version]
        ]
        self.output_schema = t.StructType.fromJson(VEP_BASIC_SCHEMA[vep_version])
        self.genome_assembly = genome_assembly

        if vep_cache_dir:
            self.add_call_arg("--dir_cache", vep_cache_dir)
        if vep_plugin_dir:
            self.add_call_arg("--dir_plugin", vep_plugin_dir)
        if offline:
            self.add_call_arg("--offline")
        if use_cache:
            self.add_call_arg("--cache")

        self.add_call_arg("--assembly", genome_assembly)

    def add_call_arg(self, *args):
        """
        Adds arguments to the VEP call

        Args:
            *args: arguments

        Returns: self
        """
        self.call_args.extend(args)

        return self

    def add_custom_annotation(
            self,
            file_path: str,
            file_type: str,
            short_name: str,
            annotation_type: str = "overlap",
            force_report_coordinates: bool = False,
            vcf_fields: List[Tuple[str, t.DataType]] = ()
    ):
        """

        Args:
            file_path:
                The path to the file. For tabix indexed files, the VEP will check that both the file and the
                corresponding .tbi file exist. For remote files, VEP will check that the tabix index is accessible on
                startup.
            short_name:
                A name for the annotation that will appear as the key in the key=value pairs in the results.
                If not defined, this will default to the annotation filename for the first set of annotation added
                (e.g. "myPhenotypes.bed.gz" in the second example below if the short name was missing).
            file_type: any of "bed", "gff", "gtf", "vcf" or "bigwig"
            annotation_type: "exact" or "overlap".
                When using "exact" only annotations whose coordinates match exactly those of the variant will be
                reported. This would be suitable for position specific information such as conservation scores, allele
                frequencies or phenotype information. Using "overlap", any annotation that overlaps the variant by even
                1bp will be reported.
            force_report_coordinates:
                If set to "True", this forces VEP to output the coordinates of an overlapping custom feature instead of
                any found identifier (or value in the case of bigWig) field. If set to "False" (the default), VEP will
                output the identifier field if one is found. If none is found, then the coordinates are used instead.
            vcf_fields:
                List of Tuple[field name, data type];
                You can specify any info type (e.g. "AC") present in the INFO field of the custom input VCF, to add
                these as custom annotations:
                    - If using "exact" annotation type, allele-specific annotation will be retrieved.
                    - The INFO field name will be prefixed with the short name, e.g. using short name "test", the INFO
                      field "foo" will appear as "test_FOO" in the VEP output.
                    - In VCF files the custom annotations are added to the CSQ INFO field.
                    - Alleles in the input and VCF entry are trimmed in both directions in an attempt to match complex
                      or poorly formatted entries.

        Returns: self
        """

        if force_report_coordinates:
            force_report_coordinates = "1"
        else:
            force_report_coordinates = "0"

        arg = ",".join([
            file_path,
            short_name,
            file_type,
            annotation_type,
            force_report_coordinates,
            *[field_name for field_name, field_dtype in vcf_fields]
        ])

        log.debug("Adding option to VEP call: --custom %s", arg)
        self.add_call_arg("--custom", arg)

        empty_field = t.StructField('custom_annotations', t.StructType())

        self.get_output_struct_field(
            "custom_annotations",
            create_if_nonexistent=empty_field
        ).add(
            t.StructField(short_name, t.ArrayType(t.StructType([
                t.StructField("name", t.StringType()),
                t.StructField("allele", t.StringType()),
                t.StructField("fields", t.StructType([
                    t.StructField(field, field_type) for field, field_type in vcf_fields
                ])),
            ])))
        )

        return self

    def add_predefined_vcf_annotation(self, file_path, dataset="GnomADg_v3", annotation_type="overlap"):
        vcf_fields = VEP_CUSTOM_ANNOT[dataset]
        vcf_fields = [
            (field_name, t._parse_datatype_json_value(field_dtype)) for field_name, field_dtype in vcf_fields.items()
        ]

        return self.add_custom_annotation(
            file_path=file_path,
            file_type="vcf",
            short_name=dataset,
            annotation_type=annotation_type,
            vcf_fields=vcf_fields,
        )

    def add_colocated_variants(self):
        empty_field = t.StructField('colocated_variants', t.ArrayType(t.StructType()))

        f = self.get_output_struct_field("colocated_variants", create_if_nonexistent=empty_field).elementType
        f.add(t.StructField("seq_region_name", t.StringType()))
        f.add(t.StructField("strand", t.IntegerType()))
        f.add(t.StructField("start", t.LongType()))
        f.add(t.StructField("end", t.LongType()))
        f.add(t.StructField("id", t.StringType()))
        f.add(t.StructField("allele_string", t.StringType()))

        return self

    def add_sift(self, annotype="b"):
        """
	  	Species limited SIFT predicts whether an amino acid substitution affects protein function based on sequence
	  	homology and the physical properties of amino acids. VEP can output the prediction term, score or both.

        Args:
            annotype:
                - p: prediction term
                - s: score
                - b: both

        Returns: self
        """
        self.add_call_arg("--sift", annotype)

        f = self.get_output_struct_field("transcript_consequences").elementType
        f.add(t.StructField("sift_prediction", t.StringType()))
        f.add(t.StructField("sift_score", t.FloatType()))

        return self

    def add_polyphen(self, annotype="s"):
        """
        Human only. PolyPhen is a tool which predicts possible impact of an amino acid substitution on the structure and
        function of a human protein using straightforward physical and comparative considerations. VEP can output the
        prediction term, score or both. VEP uses the humVar score by default - use --humdiv to retrieve the humDiv
        score.

        Args:
            annotype:
                - p: prediction term
                - s: score
                - b: both

        Returns: self
        """
        self.add_call_arg("--polyphen", annotype)

        f = self.get_output_struct_field("transcript_consequences").elementType
        f.add(t.StructField("polyphen_prediction", t.StringType()))
        f.add(t.StructField("polyphen_score", t.FloatType()))

        return self

    def add_cadd_plugin(self, cadd_dir):
        """
        Adds CADD predictions to the transcript consequences

        Args:
            cadd_dir: directory of CADD data

        Returns: self
        """
        arg = f"CADD,{cadd_dir}/whole_genome_SNVs.tsv.gz,{cadd_dir}/InDels.tsv.gz"

        self.add_call_arg("--plugin", arg)

        f = self.get_output_struct_field("transcript_consequences").elementType
        f.add(t.StructField("cadd_phred", t.FloatType()))
        f.add(t.StructField("cadd_raw", t.FloatType()))

        return self

    def get_output_struct_field(self, key: Union[str, List[str]] = None, create_if_nonexistent: t.StructField = None):
        """
        Returns a struct field in the output schema

        Args:
            key: Key of struct; If the target field is inside a nested struct, a list of keys can be used as query.
            create_if_nonexistent: StructField that will be added in case the `key` does not exist yet

        Returns: StructField of requested key

        """
        if not key:
            return self.output_schema

        if type(key) == str:
            key = [key]

        # recursively call schema = schema[k].dataType for all values in 'key'
        schema = self.output_schema
        for k in key:
            # check if key is non-existent
            if k not in schema.fieldNames():
                if create_if_nonexistent:
                    schema.add(create_if_nonexistent)
                    return create_if_nonexistent.dataType
                else:
                    raise ValueError("Key %s does not exist in schema", k)

            # return field
            schema = schema[k].dataType

        return schema

    def compile_cmd(self) -> str:
        return ' '.join([
            shlex.quote(arg) for arg in self.call_args
        ])

    def _transform(
            self,
            input_df: f.DataFrame,
            contig: pyspark.sql.column.Column,
            start: pyspark.sql.column.Column,
            end: pyspark.sql.column.Column,
            ref: pyspark.sql.column.Column,
            alt: pyspark.sql.column.Column,
            id: pyspark.sql.column.Column,
    ):
        """
        Runs Ensembl VEP on a Spark DataFrame with VEP.
        The DataFrame needs to provide the following fields:
            - "contigName"
            - "start"
            - "end
            - "referenceAllele"
            - "alternateAlleles"

        Args:
            df: Spark DataFrame with contigNamem start, end, ref and alt
            contig: contig name column
            start: variant position column
            ref: reference allele column
            alt: array of alternate alleles column
            id: array of id's

        Returns:
            Spark DataFrame with single column `text` that contains json-formatted VEP output as string
        """
        import glow
        input_df = input_df.select([
            contig,
            start,
            end,
            id,
            ref,
            alt,
        ])

        vep_transformed_df = glow.transform(
            "pipe",
            input_df,
            cmd=json.dumps(self.call_args),
            inputFormatter='vcf',
            inVcfHeader='infer',
            outputFormatter='text',
        )

        return vep_transformed_df

    def _parse_text(self, vep_transformed_df: f.DataFrame):
        """
        Parses json-formatted VEP output string

        Args:
            vep_transformed_df: output of `self._transform()`

        Returns:
            Spark DataFrame with the schema as defined by `self.output_schema`
        """
        vep_df = (
            vep_transformed_df
                .withColumn('data', f.from_json('text', self.output_schema))
                .select(f.expr("data.*"))
                .drop("input")
        )

        return vep_df

    def transform(
            self,
            df: f.DataFrame,
            contig: str = "contigName",
            start: str = "start",
            end: str = "end",
            ref: str = "referenceAllele",
            alt: str = "alternateAlleles",
            id: str = None,
            alt_is_array=True,
            id_is_array=True,
    ):
        """
        Runs Ensembl VEP on a Spark DataFrame with VEP.
        The DataFrame needs to provide the following fields:
            - "contigName"
            - "start"
            - "end
            - "referenceAllele"
            - "alternateAlleles"

        Args:
            df: Spark DataFrame with contigNamem start, end, ref and alt
            contig: contig name column
            start: variant position column
            ref: reference allele column
            alt: (array of) alternate alleles column
            id: (array of) variant ids
            alt_is_array: Set this to True, if `alt` is already an array of alternate alleles.
            id_is_array: Set this to True, if `id` is already an array of alternate alleles or `id` is None.
        Returns:
            Spark DataFrame with the schema as defined by `self.output_schema`

        """
        alt = f.expr(alt)
        if not alt_is_array:
            alt = f.array(alt)

        if id:
            id = f.expr(id)
        else:
            id = f.array(firefly.spark_functions.identifiers.variant_SPDI(
                contig=contig,
                start=start,
                ref=ref,
                alt=alt,
                join_alternate_alleles=alt_is_array,
            ))

        if not id_is_array:
            id = f.array(id)

        vep_transformed_df = self._transform(
            df,
            contig=f.expr(contig).alias("contigName"),
            start=f.expr(start).alias("start"),
            end=f.expr(end).alias("end"),
            ref=f.expr(ref).alias("referenceAllele"),
            alt=alt.alias("alternateAlleles"),
            id=id.alias("names")
        )

        return self._parse_text(vep_transformed_df)
