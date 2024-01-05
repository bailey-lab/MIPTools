import sys
sys.path.append("/opt/src")
import mip_functions as mip

wdir=snakemake.params['wdir']
settings_file=snakemake.params['settings_file']
settings = mip.get_analysis_settings(wdir+'/'+settings_file)

vcf_file="/opt/analysis/variants.vcf.gz"
geneid_to_genename=snakemake.params.geneid_to_genename
target_aa_annotation=snakemake.params.target_aa_annotation
aggregate_nucleotides=snakemake.params.aggregate_nucleotides
aggregate_aminoacids=snakemake.params.aggregate_aminoacids
target_nt_annotation=snakemake.params.target_nt_annotation
annotate=snakemake.params.annotate
decompose_options=snakemake.params.decompose_options
annotated_vcf=snakemake.params.annotated_vcf
aggregate_none=snakemake.params.aggregate_none
output_prefix=snakemake.params.output_prefix
min_site_qual=1
min_target_site_qual=-1
min_genotype_qual=-1
min_mean_alt_qual=-1

vcf_file = vcf_file.split("/")[-1]
mip.vcf_to_tables_fb(
     vcf_file,
     settings=settings,
     settings_file=settings_file,
     annotate=annotate,
     geneid_to_genename=geneid_to_genename,
     target_aa_annotation=target_aa_annotation,
     aggregate_aminoacids=aggregate_aminoacids,
     target_nt_annotation=target_nt_annotation, 
     aggregate_nucleotides=aggregate_nucleotides, 
     decompose_options=decompose_options,
     annotated_vcf=annotated_vcf,
     aggregate_none=aggregate_none,
     min_site_qual=min_site_qual,
     min_target_site_qual=min_target_site_qual,
     min_genotype_qual=min_genotype_qual,
     min_mean_alt_qual=min_mean_alt_qual,
     output_prefix=output_prefix)
