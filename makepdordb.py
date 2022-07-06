#!/usr/bin/python3
"""Make P-DOR database -  a well curated genomic database for the species you choose"""

import sys
import os
import argparse
import pandas as pd
import numpy as np
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
### OBTAINING INFO FROM PATRIC
	#whole_list  = of all the PATRIC genomes for the selected species
		whole_list =("p3-all-genomes --eq genome_name,%s | p3-get-genome-data -a genome_name -a genome_status -a genome_quality -a genome_length -a contigs | p3-match --col 4 Good > %s_list.txt") %(species,"_".join(species.strip('"').split(" ")))
		os.system(whole_list)
	#lst of column names which needs to be string
		lst_str_cols = ['genome_id', 'genome_name']
    #use dictionary comprehension to make dict of dtypes
		dict_dtypes = {x : 'str'  for x in lst_str_cols}
		whole_db = pd.read_csv("_".join(species.strip('"').split(" "))+"_list.txt", sep='\t', skiprows=1, names = ['genome_id', 'genome_name', 'genome_status', 'genome_quality', 'genome_length', 
'genome_contigs'], dtype=dict_dtypes)
		#whole_db['genome_id'] = whole_db['genome_id'].astype(str)

### FILTERING out bad genomes by size and number of contigs
		fsize_db = filter_size(whole_db)
		fcontigs_db = filter_contigs(fsize_db)
		fcontigs_db.to_csv("../"+"_".join(species.strip('"').split(" "))+dt_string+"_list.csv")
		todownload = fcontigs_db[['genome_id']]
		todownload.columns = ['genome.genome_id']
		todownload.to_csv("genomes_to_download", index=False)

### DOWNLOAD
	#Do you really want to start the download?
	#assuming 1 base = ~ 1 byte
		spacerequired =  humansize(fcontigs_db['genome_length'].sum())
		num_genomes = len(fcontigs_db)
		print(("You're going to download %s genomes from the PATRIC database for a total of %s \nPlease check you have enough space on disk") %(num_genomes, spacerequired))
		Join = input('Do you want to continue? [y/n]\n')

		#YES download!
		if Join.lower() == 'yes' or Join.lower() == 'y':
			print("Downloading genomes...")
			name_folder = "PDORdb_" + (args.species).replace(' ', '') + dt_string
			os.chdir("..")
			if not os.path.exists(name_folder):
				 os.mkdir(name_folder)
			os.chdir(name_folder)
			os.system('for i in `cat ../.tmp/genomes_to_download`; do wget -qN "ftp://ftp.patricbrc.org/genomes/$i/$i.fna";done')
			os.system('rm -r ../.tmp')
### MASH SKETCH
			os.system("ls * > sketches")
			os.system('mash sketch -l sketches')
			print("...Done! Now you can run PDOR analysis")

	#Exiting
		if Join.lower() == 'no' or Join.lower() == 'n':
			print("Exiting...")
			os.system("rm *")

#### ONLY SKETCH MODE
	if args.sketch_only is not None:
		os.chdir(args.sketch_only)
		os.system("ls * > sketches")
		os.system('mash sketch -l sketches')
		print("...Done! Now you can run PDOR analysis")



def filter_size(df):
	Q1 = df['genome_length'].quantile(0.25)
	Q3 = df['genome_length'].quantile(0.75)
	IQR = Q3 -Q1
	size_filtered = df[~((df['genome_length'] < (Q1 - 2 * IQR)) |(df['genome_length'] > (Q3 + 2 * IQR)))]
	return size_filtered

def filter_contigs(dfiltersize):
	Q1 = dfiltersize['genome_contigs'].quantile(0.25)
	Q3 = dfiltersize['genome_contigs'].quantile(0.75)
	IQR = Q3 -Q1
	contigs_filtered = dfiltersize[~((dfiltersize['genome_contigs'] > (Q3 + 2 * IQR)))]
	return contigs_filtered


def humansize(nbytes):
	suffixes = ['B', 'KB', 'MB', 'GB', 'TB', 'PB']
	i = 0
	while nbytes >= 1024 and i < len(suffixes)-1:
		nbytes /= 1024.
		i += 1
	f = ('%.2f' % nbytes).rstrip('0').rstrip('.')
	return '%s %s' % (f, suffixes[i])


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
