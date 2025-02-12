from urllib.request import urlopen
import plotly.express as px
import json
import pandas as pd

def calculate_prevalences(metadata_file, prevalences_input_table, mutations, output_summary_table):
	def create_site_dict(metadata_file):
		'''
		Takes a metadata csv of format Sites,Sampleid 
		and creates a dictionary {Sampleid+UMI_suffix:Site}
		'''
		site_dict = {}
		for line_number, line in enumerate(open(metadata_file)):
			if line_number > 0: #discard header
				line = line.strip().split(',')
				sample = line[1]
				site = line[0]
				site_dict[sample] = site
		return site_dict

	def get_counts(prevalences_input_table, site_dict):
		'''
		Creates a dictionary of format {site:{column_number:[alt_count, cov_count]}}
		column_number is used instead of mutation_name because occasionally a mutation_name appears twice in the input_table
		'''
		count_dict = {}
		count_dict['overall'] = {}
		base, top = 0, 0
		for line_number, line in enumerate(open(prevalences_input_table, 'r')):
			if "Mutation Name" in line:
				mutations_list = line.strip().split(',')[1:]
				for mutation in mutations_list:
					if mutations_list.count(mutation) > 1:
						print(f"error: {mutation} appears in the prevalences_input_table more than once")
				mutation_dict = dict(enumerate(mutations_list))
			if line_number >= 6:
				line = line.strip().split(',')
				sample = line[0]
				sample_tallies = line[1:]
				if sample in site_dict:
					site = site_dict[sample]
					if site not in count_dict:
						count_dict[site] = {}
					for column_number, tally in enumerate(sample_tallies):
						# mutation_name = column_number
						if tally == '':
							if column_number not in count_dict[site]:
								count_dict[site][column_number] = [0, 0]
							if column_number not in count_dict['overall']:
								count_dict['overall'][column_number] = [0, 0]
						if tally == '0.0':
							if column_number not in count_dict[site]:
								count_dict[site][column_number] = [0, 0]
							count_dict[site][column_number][1] += 1
							if column_number not in count_dict['overall']:
								count_dict['overall'][column_number] = [0, 0]
							count_dict['overall'][column_number][1] += 1
						if tally == '1.0':
							if column_number not in count_dict[site]:
								count_dict[site][column_number] = [0, 0]
							count_dict[site][column_number][1] += 1
							count_dict[site][column_number][0] += 1
							if column_number not in count_dict['overall']:
								count_dict['overall'][column_number] = [0, 0]
							count_dict['overall'][column_number][1] += 1
							count_dict['overall'][column_number][0] += 1
		return count_dict, mutation_dict

	def create_output_file(count_dict, mutations, output_summary_table):
		output_file = open(output_summary_table, 'w')
		output_file.write('Sites')
		sites = list(count_dict.keys())
		sites.remove('overall')
		sites.sort()
		for column_number in mutation_dict:
			if mutation_dict[column_number] in mutations:
				output_file.write('\t'+mutation_dict[column_number])
		for site in sites:
			output_file.write('\n'+site)
			for column_number in count_dict[site]:
				if mutation_dict[column_number] in mutations:
					alt = count_dict[site][column_number][0]
					cov = count_dict[site][column_number][1]
					if alt == 0:
						prevalence = 0
					else:
						prevalence = alt/cov
					output_file.write(f"\t{prevalence} ({alt}/{cov})")
		output_file.write('\n'+'overall')
		for column_number in count_dict['overall']:
			if mutation_dict[column_number] in mutations:
				alt = count_dict['overall'][column_number][0]
				cov = count_dict['overall'][column_number][1]
				if alt == 0:
					prevalence = 0
				else:
					prevalence = alt/cov
				output_file.write(f"\t{prevalence} ({alt}/{cov})")

	site_dict = create_site_dict(metadata_file)
	count_dict, mutation_dict = get_counts(prevalences_input_table, site_dict)
	create_output_file(count_dict, mutations, output_summary_table)
	# print(count_dict)

def load_json(region_of_interest):
	with urlopen('https://raw.githubusercontent.com/simkin-bioinformatics/geojson_files/main/'+region_of_interest+'.geojson') as response:
		region_df = json.load(response)
	return region_df

def read_prevalence_table(sample_set, variants_of_interest):
	summary_table = sample_set+'_prevalence_summary.tsv'
	prevalence_summary_df = pd.read_csv(summary_table, sep='\t')
	for variant in variants_of_interest:
		prevalence_summary_df[variant]=[float(x.split()[0]) for x in prevalence_summary_df[variant]]
	return prevalence_summary_df

def display_choropleth(region_of_interest, variant, sample_set, variants_of_interest):
	region_dict = {'uganda':{'zoom':5.8, 'coordinates':{"lat": 1.3733, "lon": 32.2903}},
                   'tanzania':{'zoom':4, 'coordinates':{"lat": -6.3690, "lon": 34.8888}}
                  }
	json_file = load_json(region_of_interest)
	prevalence_summary_df = read_prevalence_table(sample_set, variants_of_interest)
	fig = px.choropleth_mapbox(prevalence_summary_df, 
                                geojson=json_file, 
                                locations='Sites', 
                                color=variant,
                                color_continuous_scale="reds",
                                mapbox_style="open-street-map",
                                featureidkey="properties.name",
                                zoom=region_dict[region_of_interest]['zoom'], 
                                center = region_dict[region_of_interest]['coordinates'],
                                labels={'Sites':'Site'},
                                range_color=(0,1),
    )
	fig.update_layout(margin={"r":0,"t":40,"l":0,"b":0})
	fig.update_layout(height=500, title=sample_set)
	return fig