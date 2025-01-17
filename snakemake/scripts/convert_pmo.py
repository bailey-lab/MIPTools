sys.path.append("/opt/src")
sys.path.append("/opt/src/pmo_scripts")
from transformer import transform_panel_info, transform_mhap_info
import pandas as pd
import io

wrangler_table=snakemake.input.final_table
mip_arms=snakemake.input.mip_arms
gff_URL=snakemake.params.gff_URL
genome_name=snakemake.params.genome_name
taxon_ID=snakemake.params.species_ID
genome_URL=snakemake.params.genome_URL
panel_ID=snakemake.params.panel_ID
version=snakemake.params.genome_version
bioinfo_ID=snakemake.params.wrangled_name
panel_pmo=open(snakemake.output.panel_pmo, 'w')
microhaplotype_pmo=open(snakemake.output.microhaplotype_pmo, 'w')
combined_pmo=open(snakemake.output.combined_pmo, 'w')
selected_additional_fields=None

def make_panel_pmo():
	'''
	creates a json formatted PMO from the MIP arms file (the panel of MIPs from
	the project resources folder)
	'''
	genome_info={"name": genome_name,
			"taxon_id": taxon_ID,
			"url": genome_URL,
			"version": version}
	field_mapping={'target_id': 'mip_id', 'forward_primers': 'extension_arm', 'reverse_primers': 'ligation_arm'}
	df=pd.read_csv(mip_arms, sep='\t')
	transformed_df=transform_panel_info(df, panel_ID, field_mapping, genome_info, 
										selected_additional_fields)
	panel_pmo.write(transformed_df)
	return transformed_df

def make_microhaplotype_pmo():
	'''
	creates a json formatted PMO from the allInfo file (the wrangled output of
	miptools wrangler).
	'''
	field_mapping={'sampleID': 's_Sample', 'locus': 'p_targetName', 'asv': 'c_seq', 'reads': 'c_barcodeCnt'}
	df=pd.read_csv(wrangler_table, sep='\t')
	transformed_df=transform_mhap_info(
	df, bioinfo_ID, field_mapping, selected_additional_fields)
	microhaplotype_pmo.write(transformed_df)
	return transformed_df

panel_pmo=make_panel_pmo()
microhaplotype_pmo=make_microhaplotype_pmo()

combined_pmo.write(panel_pmo)
combined_pmo.write(microhaplotype_pmo)