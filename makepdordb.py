#!/usr/bin/python3
"""Make P-DOR database -  a well curated genomic database for the species you choose"""

import argparse
import numpy as np
import os
import pandas as pd
pd.options.mode.chained_assignment = None
import sys
import select
from datetime import datetime

def main():
	args = parseArgs().parse_args()
	if args.species is not None:
	#create a tmp directory to store intermediate results
		if not os.path.exists(".tmp"):
			os.mkdir(".tmp")
		os.chdir(".tmp")
		species = '"' + args.species + '"'
		now = datetime.now()
		dt_string = now.strftime("_%d-%m-%Y")

### DOWNLOADING and reading last release of BV-BRC and filtering based on genome length, number of contigs and GC content
	#whole_list : last release of bvbrc genome_metadata
		whole_list = os.system("wget ftp://ftp.bvbrc.org/RELEASE_NOTES/genome_metadata")
	#sTable : filter only selected species from whole list and remove only plasmidic sequences
		sTable = filter_species("genome_metadata", args.species)
		genomelength_filter = filter_IQRbased(sTable, 'genome_length')
		contigs_filter = filter_IQRbased(genomelength_filter, 'contigs')
		gc_filter = filter_IQRbased(contigs_filter, 'gc_content')

### DOWNLOAD
	#Do you really want to start the download?
	#assuming 1 base = ~ 1 byte
		download_list = list(gc_filter['genome_id'])
		spacerequired =  humansize(gc_filter['genome_length'].sum())
		num_genomes = len(download_list)
		yes_or_no = lambda prompt: 'n' not in timeout_input(prompt + "? (Y/n)", 30, default="y")[1].lower()
		print(("You're going to download %s genomes from the BV-BRC database for a total of %s \nPlease check you have enough space on disk") %(num_genomes, spacerequired))
		answer = yes_or_no('Continue')

		#YES download!
		if answer == True:
			print("Downloading genomes...")
			name_folder = "PDORdb_" + (args.species).replace(' ', '') + dt_string
			os.chdir("..")
			if not os.path.exists(name_folder):
				 os.mkdir(name_folder)
			os.chdir(name_folder)
			download = download_genome(download_list)
		## Check all the genomes have been downloaded
			downloaded_genomes = [s.strip('.fna') for s in os.listdir()]
			check_download = is_identical(download_list, downloaded_genomes)
			if check_download == False:
				pass
			else: #If a genome(/a few genomes) is missing, try to download it again
				for genome in check_download:
					try_download = download_genome(genome)
			final_db_list = [s.strip('.fna') for s in os.listdir()]
			print(("A final number of %s genomes have been downloaded") % len(final_db_list))
			final_db_tb = pd.DataFrame(final_db_list, columns=["genome_id"], dtype=str)
			FinalTable = pd.merge(gc_filter, final_db_tb, how='inner')
			FinalTable.to_csv("../"+"_".join(species.strip('"').split(" "))+dt_string+"_list.tsv", sep='\t')


### MASH SKETCH
			os.system("ls * > sketches")
			os.system('mash sketch -l sketches')
			os.system('rm sketches')
			os.system('mv sketches.msh ../.')
			os.system('rm -r ../.tmp')
			print("...Done! Now you can run PDOR analysis")

	#Exiting
		else:
			print("Exiting...")
			os.system("rm -r ../.tmp")

#### ONLY SKETCH MODE
	if args.sketch_only is not None:
		os.chdir(args.sketch_only)
		os.system("ls * > sketches")
		os.system('mash sketch -l sketches')
		os.system('mv sketches.msh ../.')
		print("...Done! Now you can run PDOR analysis: python P-DOR.py -q [query genome folder] -sd sketches.msh -ref [reference genome] -snp_thr infl")



### DEF
def filter_species(input, species):
	completeDF = pd.read_csv(input, sep='\t', dtype=str)
	sdf = completeDF.loc[completeDF.genome_name.str.contains(species, case=False)]
	species_table = sdf[sdf['genome_status'] != "Plasmid"]
	return species_table

def filter_IQRbased(df, column):
	df.fillna(0,inplace=True)
	df[column] = df[column].astype('float')
	Q1 = df[column].quantile(0.25)
	Q3 = df[column].quantile(0.75)
	IQR = Q3 -Q1
	df_filt = df[~((df[column] < (Q1 - 2 * IQR)) |(df[column] > (Q3 + 2 * IQR)))]
	return df_filt

def humansize(nbytes):
	suffixes = ['B', 'KB', 'MB', 'GB', 'TB', 'PB']
	i = 0
	while nbytes >= 1024 and i < len(suffixes)-1:
		nbytes /= 1024.
		i += 1
	f = ('%.2f' % nbytes).rstrip('0').rstrip('.')
	return '%s %s' % (f, suffixes[i])

def download_genome(dwlist):
	for genome in dwlist:
		cmd = ("wget -qN ftp://ftp.patricbrc.org/genomes/%s/%s.fna") % (genome, genome)
		os.system(cmd)

def is_identical(list_a, list_b):
	if len(list_a) != len(list_b):
		try_new_download = []
		for i in list_a:
			if i not in list_b:
				try_new_download.append(i)
		return try_new_download
	else:
		return False

def timeout_input(prompt, timeout=10, default=""):
    print(prompt, end=': ', flush=True)
    inputs, outputs, errors = select.select([sys.stdin], [], [], timeout)
    print()
    return (0, sys.stdin.readline().strip()) if inputs else (-1, default)

### HELP
def parseArgs():
	parser = argparse.ArgumentParser(prog='makepdordb', description='makepdordb: create a well curated genomic database for P-DOR analysis')
	fullpipe = parser.add_mutually_exclusive_group(required=True) # options mutually exclusive
	fullpipe.add_argument('-s', '--species', help=' Name of the species you want to analyse  e.g. makepdordb -s "Acinetobacter baumannii" ', required=False, metavar='')
	fullpipe.add_argument('-sketch', '--sketch_only', help=' Sketching a user-defined database e.g. makepdordb -sketch mygenomes/', required=False, metavar='')
	return parser



# Call main function
if __name__ == '__main__':
	main()
