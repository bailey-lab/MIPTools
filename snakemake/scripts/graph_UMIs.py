import sys
sys.path.append("/opt/src")
#import mip_functions_freebayes_call_edit as mip
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px
import os

sample_summary_csv = snakemake.input.sample_summary_csv
umi_vs_probe_graph = snakemake.output.umi_vs_probe_graph

wdir=snakemake.params['wdir']

def make_graphing_list(UMI_file):
	import math
	graphing_list, rows=[],[]
	for line_number, line in enumerate(open(UMI_file)):
		line=line.strip().split(',')
		if line_number==0:
			columns=line[1:]
		if line_number>2:
			rows.append(line[0])
			int_line=list(map(int, list(map(float, line[1:]))))
			log_line=[math.log(number+1, 2) for number in int_line]
			graphing_list.append(log_line)
	return graphing_list, columns, rows

def plot_heatmap(graphing_list, x_values, y_values, x_title, y_title, count_title, output_path, width=2000, height=4000):
	import plotly.express as px
#	print(graphing_list)
	fig = px.imshow(graphing_list, aspect='auto', labels=dict(x=x_title, y=y_title,
	color=count_title), x=x_values, y=y_values)
	fig.update_xaxes(side="top")
	fig.update_layout(width=width, height=height, autosize=False)
	fig.write_html(output_path)

def plot_umi_vs_probe(sample_summary_csv, umi_vs_probe_graph):
	sample_summary = pd.read_csv(sample_summary_csv)
	fig = px.scatter(sample_summary, 
        x="UMI Count",
        y="targets_with_>=10_UMIs", 
        hover_name="Sample ID",
        title="UMI Count vs. Probe Coverage",
        hover_data="Read Count"
	)
	fig.show()
	fig.write_html(umi_vs_probe_graph)


graphing_list, x_values, y_values=make_graphing_list(wdir+'/UMI_counts.csv')
plot_heatmap(graphing_list, x_values, y_values, 'mips', 'samples', 'log2 of umi_counts+1', '/opt/analysis/umi_heatmap.html')
plot_umi_vs_probe(sample_summary_csv, umi_vs_probe_graph)