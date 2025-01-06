from urllib.request import urlopen
import plotly.express as px
import json
import pandas as pd

def calculate_prevalences(metadata_file, prevalences_input_table, mutations, output_summary_table, sample_column, summarize_column):
	def create_summary_dict(metadata_file, sample_column, summarize_column):
		'''
		returns a dictionary keyed by sample that returns a list of the
		location, (from summarize_column), latitude, and longitude
		associated with each sample. Assumes that the metadata file has columns
		named 'Latitude' and 'Longitude'
		'''
		summarize_dict = {}
		for line_number, line in enumerate(open(metadata_file)):
			line = line.strip().split(',')
			if line_number==0:
				sample_column = line.index(sample_column)
				summarize_column = line.index(summarize_column)
				latitude=line.index('Latitude')
				longitude=line.index('Longitude')
			if line_number > 0: #discard header
				summarize_dict[line[sample_column]] = [line[summarize_column], float(line[latitude]), float(line[longitude])]
		return summarize_dict

	def get_counts(prevalences_input_table, summary_dict):
		'''
		Creates a dictionary of format {location:{mutation_name:[alt_count, cov_count]}}
		If a mutation occurs twice an error is reported
		'''
		count_dict, coord_dict = {}, {}
		count_dict.setdefault('overall', {})
		base, top = 0, 0
		for line_number, line in enumerate(open(prevalences_input_table, 'r')):
			if "Mutation Name" in line:
				mutations_list = line.strip().split(',')[1:]
				for mutation in mutations_list:
					if mutations_list.count(mutation) > 1:
						print(f"error: {mutation} appears in the prevalences_input_table more than once")
						exit()
			if line_number >= 6:
				line = line.strip().split(',')
				sample = line[0]
				sample_tallies = line[1:]
				if sample in summary_dict:
					location, latitude, longitude = summary_dict[sample]
					coord_dict.setdefault(location, [[],[]])
					count_dict.setdefault(location, {})
					count_dict.setdefault('overall', {})
					coord_dict[location][0].append(latitude)
					coord_dict[location][1].append(longitude)
					for column_number, tally in enumerate(sample_tallies):
						# mutation_name = column_number
						mutation=mutations_list[column_number]
						count_dict[location].setdefault(mutation, [0,0])
						count_dict['overall'].setdefault(mutation, [0,0])
						if tally:
							count_dict[location][mutation][1] += 1
							count_dict['overall'][mutation][1] += 1
							tally_float=float(tally)
							if tally_float>0:
								count_dict[location][mutation][0] += tally_float
								count_dict['overall'][mutation][0] += tally_float
		return count_dict, coord_dict

	def parse_counts(count_dict, found_mutations, location, output_line):
		for mutation in found_mutations:
			alt = count_dict[location][mutation][0]
			if int(alt)==alt:
				alt=int(alt)
			cov = count_dict[location][mutation][1]
			if alt == 0:
				prevalence = 0
			else:
				prevalence = round(alt/cov, 4)
			output_line.append(f'{prevalence} ({round(alt, 4)}/{cov})')
		return output_line

	def create_output_file(count_dict, coord_dict, mutations, output_summary_table, summarize_column):
		output_file = open(output_summary_table, 'w')

		#write the name of the column from the metadata sheet that we summarized by
		header_line=[]
		header_line.append(summarize_column)
		locations = list(count_dict.keys())
		locations.remove('overall')
		locations.sort()
		header_line.extend(['Latitude', 'Longitude'])

		#filter down to the subset of mutations of interest that were present in the AA tables
		found_mutations=[mutation for mutation in mutations if mutation in count_dict['overall']]
		for mutation in found_mutations:
			if mutation in count_dict['overall']:
				header_line.append(mutation)
		output_file.write('\t'.join(header_line)+'\n')
		for location in locations:
			output_line=[location]
			latitude=round(sum(coord_dict[location][0])/len(coord_dict[location][0]), 7)
			longitude=round(sum(coord_dict[location][1])/len(coord_dict[location][1]), 7)
			output_line.append(str(latitude))
			output_line.append(str(longitude))
			output_line=parse_counts(count_dict, found_mutations, location, output_line)
			output_file.write('\t'.join(output_line)+'\n')
		output_line=['overall', '', '']
		location='overall'
		output_line=parse_counts(count_dict, found_mutations, location, output_line)
		output_file.write('\t'.join(output_line)+'\n')
	summary_dict = create_summary_dict(metadata_file, sample_column, summarize_column)
	count_dict, coord_dict = get_counts(prevalences_input_table, summary_dict)
	create_output_file(count_dict, coord_dict, mutations, output_summary_table, summarize_column)
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