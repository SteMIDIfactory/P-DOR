#!/usr/bin/python3
"""Make P-DOR database -  a well curated genomic database for the species you choose"""

import argparse, os, sys, select, copy
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None
from datetime import datetime
from urllib import request

def main():
	args = parseArgs().parse_args()
	mash_threads=args.threads
	max_threads=os.popen("grep -c '^processor' /proc/cpuinfo").read().strip()
	if mash_threads>int(max_threads)/2:
		mash_threads=int(round(int(max_threads)/2,0))
		print("\n\nWARNING: The number of threads allocated for sketching genomes was automatically lowered to %i" %(mash_threads))
		print("Please bear in mind that MASH could still malfunction if other intensive processes are running on your machine.\nThis is due to a knonw MASH bug. Plase see https://github.com/marbl/Mash/issues/140\n\n")

	if args.subcommand=="download":
	#create a tmp directory to store intermediate results
		if not os.path.exists(".tmp"):
			os.mkdir(".tmp")
		os.chdir(".tmp")
		species = '"' + args.species + '"'
		now = datetime.now()
		dt_string = now.strftime("_%d-%m-%Y")
		speciesOUT="_".join(str(args.species).split())
### DOWNLOADING and reading last release of BV-BRC and filtering based on genome length, number of contigs and GC content
	#whole_list : last release of bvbrc genome_metadata


		whole_list = os.system("wget ftp://ftp.bvbrc.org/RELEASE_NOTES/genome_metadata")


	#sTable : filter only selected species from whole list and remove only plasmidic sequences
		sTable = filter_species("genome_metadata", args.species)
		genomelength_filter = filter_IQRbased(sTable, 'genome_length')
		contigs_filter = filter_IQRbased(genomelength_filter, 'contigs')
		gc_filter = filter_IQRbased(contigs_filter, 'gc_content')
		os.chdir("..")
		download_list = list(gc_filter['genome_id'])
		spacerequired =  gc_filter['genome_length'].sum()
		num_orig_genomes = len(download_list)
		if args.folder is not None:
			print("\n\nChecking genomes in folder %s" %(args.folder))
			extragenomes=[]
			Astats=os.popen("ls %s/*.fna | parallel -j %i 'assembly-stats -u {} | cut -f1,2'" %(args.folder,args.threads)).readlines()
			gc_filter=gc_filter.set_index("genome_id",drop=False)
			DICT=gc_filter.to_dict(orient="index")
			for AS in Astats:
				AS=AS.split()
				name=AS[0]
				id=AS[0].split("/")[-1].replace(".fna","")
				length=int(AS[1])
				try:
					if int(DICT[id]["genome_length"])==length:
						download_list.remove(id)
						spacerequired=spacerequired-length
				except:
					extragenomes.append(name)
			if len(extragenomes)>0:
				extrafolder=args.folder.strip().rstrip("/").split("/")
				extrafolder[-1]="Excluded_"+extrafolder[-1]
				extrafolder=os.path.abspath("/".join(extrafolder))
				os.mkdir(extrafolder)
				print("\n\n%i genomes in your folder must be excluded from the sketch dataset. Moving them to folder %s" %(len(extragenomes),extrafolder))
				for e in extragenomes:
					os.system("mv %s %s" %(e,extrafolder))


### DOWNLOAD
	#Do you really want to start the download?
	#assuming 1 base = ~ 1 byte
		spacerequired =  humansize(spacerequired)
		num_genomes = len(download_list)
		yes_or_no = lambda prompt: 'n' not in timeout_input(prompt + "? (Y/n)", 30, default="y")[1].lower()
		print(("\n\n%i genomes were selected for download from the BV-BRC database. %i are already present in your folder.\nYou're going to download %i genomes for a total of %s \nPlease double-check that you have enough space on disk") %(num_orig_genomes,(num_orig_genomes-num_genomes),num_genomes, spacerequired))
		answer = yes_or_no('Continue')

		#YES download!
		if answer == True:
			print("\n\nDownloading genomes...")

			if args.folder is None:
				name_folder = "PDORdb_" + (args.species).replace(' ', '') + dt_string
				if not os.path.exists(name_folder):
					os.mkdir(name_folder)
			else:
				name_folder=args.folder
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
			empties=os.popen("find . -maxdepth 1 -type f -size -100k -ls | grep '.fna' | cut -f2 -d'/'").readlines()
			for e in empties:
				os.system("rm %s" %(e))
			final_db_list = [s.strip('.fna') for s in os.listdir()]
			print(("\n\nThe dataset folder contains a total of %s genomes") % len(final_db_list))
			gc_filter=gc_filter.reset_index(drop=True)
			final_db_tb = pd.DataFrame(final_db_list, columns=["genome_id"], dtype=str)
			FinalTable = pd.merge(gc_filter, final_db_tb, how='inner')
			FinalTable.to_csv("../"+"_".join(species.strip('"').split(" "))+dt_string+"_list.tsv", sep='\t')


### MASH SKETCH
			os.system("ls *.fna > sketches_list")
			os.system('mash sketch -l sketches_list -p %i -o sketches.msh' %(mash_threads))
			os.system('rm sketches_list')
			os.system('mv sketches.msh ../%s.msh' %(speciesOUT))
			os.system('rm -r ../.tmp')


			print("...Done! The SD sketch file is stored as %s.msh\nThe SD genomes are stored in folder %s\nNow you can run the P-DOR analysis!" %(speciesOUT,name_folder))

	#Exiting
		else:
			print("Exiting...")

			os.system("rm -r .tmp")

#### ONLY SKETCH MODE
	elif args.subcommand=="sketch":
		os.chdir(args.folder)
		numberFNA = str((os.popen("ls *.fna | wc -l").read().strip()))
		print("You have %s FNA files in your folder\nWARNING: makepdordb uses only files .fna\n" %(numberFNA))
		os.system("ls *.fna > sketches")
		os.system('mash sketch -l sketches -p %i' %(mash_threads))
		os.system('mv sketches.msh ../.')
		print("\n...Done! The SD sketch file is stored as sketches.msh\nNow you can run the P-DOR analysis!")



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
		try:
			request.urlretrieve("ftp://ftp.bvbrc.org/genomes/%s/%s.fna" %(genome, genome),"%s.fna" %(genome))
		except:
			print("Genome %s was not found on the BV-BRC ftp site" %(genome))



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
	parser.add_argument('-t', '--threads', type=int, help='number of threads to be used with mash and assembly-stats (used when checking existing genomes)',required=False, default=1)
	subparsers = parser.add_subparsers(help='You can either sketch a collection of your genomes or download them from the BV-BRC database', dest="subcommand")
	parser_a = subparsers.add_parser('sketch', help='This option sketches the genomes contained in a user-defined folder')
	parser_a.add_argument('-f', '--folder', type=str, help='path to the folder where the genomes you want to sketch are stored',required=True)

	parser_b = subparsers.add_parser('download', help='This options downloads genomes from the BV-BRC database, filters them for quality and sketches them')
	parser_b.add_argument('-f', '--folder', type=str, help='path to the folder where you want to store genomes or a previous version of the dataset is stored. Genomes are checked and not downloaded if already present. If no folder name is provided, a new folder will be created automatically',required=False)
	parser_b.add_argument('-s', '--species', help=' Name of the species you want to analyse', required=True)

	return parser

# Call main function
if __name__ == '__main__':
	main()
