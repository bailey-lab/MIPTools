'''
parses the original wrangler output file to produce abbreviated versions of the
same, separated into data.tsv, samples.tsv, and unique_haplotypes.csv files.
'''
import sys
sys.path.append("/opt/src")
import mip_functions as mip

wdir=snakemake.params['wdir']+'/'
settings_file=snakemake.params['settings_file']
info_files=snakemake.params['info_files']
sample_sheets=[snakemake.params['sample_sheets']]
sample_groups=snakemake.params['sample_groups']
settings = mip.get_analysis_settings(wdir+'/'+settings_file)

if len(info_files) > 1:
	mip.combine_info_files(wdir, settings_file, info_files, sample_sheets,
	settings["mipsterFile"],
	sample_sets=sample_groups)
else:
	mip.process_info_file(wdir, settings_file, info_files, sample_sheets,
	settings["mipsterFile"], sample_sets=sample_groups)
