import subprocess
sys.path.append("/opt/src")
#import mip_functions_freebayes_call_edit as mip

freebayes_command_dict = snakemake.params.freebayes_command_dict

vcf_file = snakemake.wildcards.contig
#vcf_file = snakemake.params.contig_name

command = freebayes_command_dict[vcf_file]

freebayes_status=subprocess.run(command,shell=True)
if freebayes_status.returncode==0:
	subprocess.call(f"bgzip -f /opt/analysis/contig_vcfs/{vcf_file}.vcf",shell=True)
	subprocess.call(f"bcftools index -f /opt/analysis/contig_vcfs/{vcf_file}.vcf.gz",shell=True)
