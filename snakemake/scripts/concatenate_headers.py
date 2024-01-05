import sys
import subprocess
import gzip
import os
import yaml
sys.path.append("/opt/src")
import mip_functions as mip

contig_vcf_gz_paths_yaml = open("/opt/analysis/contig_vcf_gz_paths.yaml",'r')
contig_vcf_gz_paths = yaml.safe_load(contig_vcf_gz_paths_yaml)

vcf_file="/opt/analysis/variants.vcf.gz"

wdir=snakemake.params['wdir']
settings_file=snakemake.params['settings_file']
options=snakemake.params['freebayes_settings']
settings = mip.get_analysis_settings(wdir+'/'+settings_file)

# concatanate contig vcfs. The number of contigs may be high, so we'll
# write the vcf paths to a file and bcftools will read from that file
cvcf_paths_file = os.path.join(wdir, "contig_vcfs", "vcf_file_list.txt")
with open(cvcf_paths_file, "w") as outfile:
    outfile.write("\n".join(contig_vcf_gz_paths) + "\n")
subprocess.run(["bcftools", "concat", "-f", cvcf_paths_file, "-Oz",
                "-o", vcf_file], check=True)
subprocess.run(["bcftools", "index", "-f", vcf_file], check=True)

# fix vcf header if --gvcf option has been used
if "--gvcf" in options:
    temp_vcf_path = os.path.join(wdir, "temp.vcf.gz")
    mip.vcf_reheader(os.path.basename(vcf_file), temp_vcf_path, wdir=wdir)
    old_vcf_path = os.path.join(wdir, "unfixed.vcf.gz")
    subprocess.run(["mv", vcf_file, old_vcf_path])
    subprocess.run(["mv", temp_vcf_path, vcf_file])
    subprocess.run(["bcftools", "index", "-f", vcf_file], check=True)
    print('did a reheader')

with open('/opt/analysis/freebayes_reheader_check.txt','w') as file:
    file.write('reheader done')

