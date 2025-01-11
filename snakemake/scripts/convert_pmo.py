sys.path.append("/opt/src")
sys.path.append("/opt/src/pmo_scripts")
from transformer import transform_panel_info
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
output_pmo=open(snakemake.output.pmo_json, 'w')
selected_additional_fields=None

genome_info={"name": genome_name,
			"taxon_id": taxon_ID,
			"url": genome_URL,
			"version": version}
field_mapping={'target_id': 'mip_id', 'forward_primers': 'extension_arm', 'reverse_primers': 'ligation_arm'}

df=pd.read_csv(mip_arms, sep='\t')
transformed_df=transform_panel_info(df, panel_ID, field_mapping, genome_info, 
	selected_additional_fields)

output_pmo.write(transformed_df)

