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