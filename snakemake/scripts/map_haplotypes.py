import sys
sys.path.append("/opt/src")
import mip_functions as mip

wdir=snakemake.params['wdir']
settings_file=snakemake.params['settings_file']
settings = mip.get_analysis_settings(wdir+'/'+settings_file)

mip.map_haplotypes(settings)
mip.get_haplotype_counts(settings)
