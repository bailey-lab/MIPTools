def calculate_prevalences(metadata_file, prevalences_input_table, UMI_suffix, mutations, output_summary_table):
	def create_site_dict(metadata_file, UMI_suffix):
		site_dict = {}
		for line in open (metadata_file):
			if 'Sites' not in line:
				line = line.strip().split(',')
				sample = line[1]+UMI_suffix
				site = line[0]
				site_dict[sample] = site
		return site_dict

	def get_counts(prevalences_input_table, site_dict):
		count_dict = {}
		count_dict['overall'] = {}
		base, top = 0, 0
		for line_number, line in enumerate(open(prevalences_input_table, 'r')):
			if "Mutation Name" in line:
				mutation_dict = dict(enumerate(line.strip().split(',')[1:]))
				# need to print an error message if two mutations have the same name
			if line_number >= 6:
				line = line.strip().split(',')
				sample = line[0]
				if sample in site_dict:
					site = site_dict[sample]
					if site not in count_dict:
						count_dict[site] = {}
					for tally_number, tally in enumerate(line[1:]):
						# mutation_name = tally_number
						if tally == '':
							if tally_number not in count_dict[site]:
								count_dict[site][tally_number] = [0, 0]
							if tally_number not in count_dict['overall']:
								count_dict['overall'][tally_number] = [0, 0]
						if tally == '0.0':
							if tally_number not in count_dict[site]:
								count_dict[site][tally_number] = [0, 1]
							count_dict[site][tally_number][1] += 1
							if tally_number not in count_dict['overall']:
								count_dict['overall'][tally_number] = [0, 1]
							count_dict['overall'][tally_number][1] += 1
						if tally == '1.0':
							if tally_number not in count_dict[site]:
								count_dict[site][tally_number] = [1, 1]
							count_dict[site][tally_number][1] += 1
							count_dict[site][tally_number][0] += 1
							if tally_number not in count_dict['overall']:
								count_dict['overall'][tally_number] = [1, 1]
							count_dict['overall'][tally_number][1] += 1
							count_dict['overall'][tally_number][0] += 1
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

	site_dict = create_site_dict(metadata_file, UMI_suffix)
	count_dict, mutation_dict = get_counts(prevalences_input_table, site_dict)
	create_output_file(count_dict, mutations, output_summary_table)
# print(count_dict)

