'''
The code this function depends on could probably be significantly streamlined by
reading and writing settings file dictionaries with yaml instead of Ozkan's
custom functions.
Also, not sure why mipSetKey needs an empty entry added on to the existing list.
'''
import sys
sys.path.append("/opt/src")
import subprocess
import mip_functions as mip

temp_settings_file = snakemake.params['template_settings']
processor_number=snakemake.params['processor_number']
bwa_extra=snakemake.params['bwa_extra']
species=snakemake.params['species']
probe_set = snakemake.params['probe_set']
wdir=snakemake.params['wdir']
min_haplotype_barcodes=snakemake.params['min_haplotype_barcodes']
min_haplotype_samples=snakemake.params['min_haplotype_samples']
min_haplotype_sample_fraction=snakemake.params['min_haplotype_sample_fraction']
freebayes_threads=snakemake.params['freebayes_threads']

# extract the settings template
settings = mip.get_analysis_settings(temp_settings_file)
settings["bwaOptions"]=[settings['bwaOptions']]+bwa_extra
settings['species']=species
settings['freebayes_threads']=freebayes_threads
settings['processorNumber']=processor_number
settings['mipSetKey'] = [probe_set,'']
settings['minHaplotypeBarcodes']=min_haplotype_barcodes
settings['minHaplotypeSamples']=min_haplotype_samples
settings['minHaplotypeSampleFraction']=min_haplotype_sample_fraction
mip.write_analysis_settings(settings, wdir+'/settings.txt')

#this (below) is fairly silly - it just writes a list of the names of all probes to probe_sets.json
#and it doesn't even really work - includes the name of the mip arms text file
#but this json file is required for other parts of the current pipeline to function correctly.
try:
    mip.update_probe_sets("/opt/project_resources/mip_ids/mipsets.csv",
                         "/opt/project_resources/mip_ids/probe_sets.json")
except IOError:
    pass
