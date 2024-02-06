project_resources=snakemake.params.project_resources
species_resources=snakemake.params.species_resources
wrangler_directory=snakemake.params.wrangler_directory
output_directory=snakemake.params.output_directory
#miptools_directory=snakemake.params.miptools_directory
output_profile=open(snakemake.output.profile, 'w')

output_profile.write('use-singularity: True\n')
output_profile.write(f'singularity-args: "-B {project_resources}:/opt/project_resources\n')
output_profile.write(f'  -B {species_resources}:/opt/species_resources\n')
output_profile.write(f'  -B {wrangler_directory}:/opt/data\n')
output_profile.write(f'  -B {output_directory}:/opt/analysis\n')
#output_profile.write(f'  -B {miptools_directory}:/opt/src\n')
output_profile.write(f'  --app jupyter"\n')
output_profile.write('printshellcmds: True\n')
output_profile.write('cores: 16\n')
output_profile.write('keep-going: True\n')
output_profile.write('rerun-incomplete: True\n')
output_profile.write('use-conda: True\n')
output_profile.write('latency-wait: 60\n')
output_profile.write('#keep-incomplete: True\n')
output_profile.write('#restart-times: 3\n')
