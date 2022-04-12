#!/usr/bin/python3

## LIBRARIES

import argparse
from Bio import SeqIO
import sys
import subprocess
import time
import datetime
import glob
import os
from Skynet_SNPs_library_4_0 import core_snps_2_fasta
from Skynet_SNPs_library_4_0 import core_snps_two_files
from Skynet_SNPs_library_4_0 import core_snps_list_path
from Skynet_SNPs_library_4_0 import snp_calling_mauve
from Skynet_SNPs_library_4_0 import card_resistance_genes

## FUNCTIONS


def parse_args():

	parser=argparse.ArgumentParser(description='''This is PDOR manual ''',usage='%(prog)s [options] -o <output_dir>',
		prog='pdor_v01.py',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		argument_default=argparse.SUPPRESS,
		epilog="Pdor is the Son of Kmer, the mystic god of the Carnic Mountains")

	requiredNamed = parser.add_argument_group('Input data')
	requiredNamed.add_argument('-q', help='query folder containing genomes in .fna format',metavar="<dirname>",dest="query_folder",type=str,required=True)
	requiredNamed.add_argument("-db", help="PATRIC database folder",metavar="<dirname>",dest="db_sketch", required=True)
	requiredNamed.add_argument("-ref", help="reference genome", dest="ref", metavar="<filename>",required=True)



	optional = parser.add_argument_group('Additional arguments')
#	advanced = parser.add_argument_group('Advanced_options')
	optional.add_argument("-gff", help="annotation file in gff format; if not specified (default) a dummy gff is generated",required=False,default="",nargs="?")
	optional.add_argument("-bkg_folder", help="folder containing the genomes from which the background sketch was created; if not specified (default) nearest genomes are downloaded from the PATRIC-DB",required=False,default="",nargs="?")
	optional.add_argument("-call", help="Snps calling method",required=False, choices=['purple', 'mummer'],default="mummer")
	#advanced.add_argument("-call", help="Snps calling method",required=False, choices=['Purple', 'Mummer', 'Combined'],default="Mummer")
	optional.add_argument('-n', type=int,help='Maximum nearest genomes from database',metavar="<int>",dest="near",default=20)
	#optional.add_argument('-recomb', metavar='',type=str,help="Remove recombinations [Gubbins]",default=False)
	optional.add_argument('-t', type=int,help='number of threads',metavar="<int>",dest="threads",default=10)

	parser.add_argument('-v', '--version', action='version', version='P-DOR v0.1')
	#parser.print_help()

	args = parser.parse_args()
	#if args.verbose: ### version
	return args


def check_folder(folder):
	folder_path=os.path.abspath(folder)
	if os.path.isdir(folder):
		sys.exit("\n\nERROR: the {} already exists!".format(folder_path))
		#sys.stderr.write("%s %s\n" % (foo, bar))
		#raise Exception("{} already exists".format(folder_path)) ## \nuse -f overwrite???


def Mummer_snp_call(threads,ref):
	os.chdir("Align")
	os.system('ls | parallel -j %i "nucmer -p {} %s {}; show-snps -H -C -I -r -T {}.delta | cut -f1,2,3 >{}.snp"' %(threads,ref))
	os.chdir("..")
	os.system("ls Align/*.snp >SNP_genomes.list")
	core_snps_list_path("SNP_genomes.list", 2, "SNP_positions")
	core_snps_2_fasta("SNP_genomes.list", "SNP_positions.core.list",ref,"SNP_alignment")


def Purple_snp_call(ref,gff,cpus,Results_folder_name):
	os.system("cp -r ../Purple .")
	os.system("mv Align Purple/")
	os.chdir("..")
	gffname=check_gff(ref)
	cmd="cp %s %s/Purple/" %(gffname,Results_folder_name)
	os.system(cmd)
	cmd="cp %s %s/Purple/" %(ref,Results_folder_name)
	os.system(cmd)
	os.chdir(Results_folder_name)
	os.chdir("Purple")
	localref=ref.split("/")[-1]
	cmd="perl Purple_v1.22.2.pl -r %s -g %s -f %s -cpu %i" %(localref,gffname,genome_folder,cpus)
	os.system(cmd)
	#####mv del risultato di Purple alla cartella Results
	### chdir alla cartella Results


def create_dummy_gff(ref):
	for index, record in enumerate(SeqIO.parse(ref,"fasta")):
		name=i.strip().split("/")[-1]
		if index==0:
			cmd="echo '%s	.	CDS	1	800' >%s.gff" %(record.id,name)
			os.system(cmd)
	return name+".gff"


def check_gff(ref):

	if args.gff:
		print ("using %s annotated file..." %args.gff)
		gffname=args.gff

	else:
		gffname=create_dummy_gff(ref)
		print ("generating dummy annotated file...")
	return gffname

### MAIN
args = parse_args()

timefunct=datetime.datetime.now()
Results_folder_name="Results_%s-%s-%s_%s-%s" %(str(timefunct.year),str(timefunct.month),str(timefunct.day),str(timefunct.hour),
	str(timefunct.minute))
os.mkdir(Results_folder_name)

start_time = time.time()


query_folder=parse_args().query_folder
db_sketch=parse_args().db_sketch
ref=parse_args().ref
threads=parse_args().threads
nearest=parse_args().near

ref=os.path.abspath(ref)
ABS_query_folder=os.path.abspath(query_folder)
Qglobber="%s/*.fna" %(ABS_query_folder)
genomes_query = glob.glob(Qglobber)
genomes_query= [i.strip().split("/")[-1] for i in genomes_query]


NEAREST=[]
for i in genomes_query:
	cmd="mash dist %s/%s %s -p %i | sort -gr -k5 | head -n %i"  %(query_folder,i,db_sketch,threads,nearest)
	NEAREST.append(subprocess.check_output(cmd,shell=True).read().split("\n"))

NEAREST=list(set(NEAREST))

print ("Sketches completed...")

os.chdir(Results_folder_name)

with open("query_genome_list","w") as q:
	for i in genomes_query:
		q.write("%s\n" %(i))

with open("backgroud_genome_list","w") as q:
	for i in NEAREST:
		q.write("%s\n" %(i))

os.mkdir("Background")

if args.bkg_folder:
	print ("retrieving nearest genomes from the %s folder" %args.bkg_folder)
	os.chdir("..")
	background_folder=os.path.abspath(bkg_folder)
	os.chdir(Results_folder_name)
	for N in NEAREST:
		cmd="cp %s/%s Background/" %(background_folder,N)
else:
	os.chdir("Background")
	print ("etrieving nearest genomes from the PATRIC-DB")
	for N in NEAREST:
		cmd='wget -qN "ftp://ftp.patricbrc.org/genomes/%s/%s.fna"' %(N,N)
		os.system(cmd)
	os.chdir("..")

os.mkdir("Align")
os.system("cp Background/* Align; cp ABS_query_folder/*.fna Align")

if args.call == 'mummer':
	Mummer_snp_call(threads,ref)
elif args.call == 'purple':
	Purple_snp_call(ref,gff,threads,Results_folder_name)


"""

cmd="raxmlHPC-PTHREADS -f a -x 12345 -p 12345 -# 100 -m ASC_GTRGAMMA -s SNP_alignment.core.fasta -n SNP_Phylo.nwk -T %i --no-bfgs --asc-corr=lewis" %(threads)
os.system(cmd)



############ R Plot

cmd="Rscript Snp_breaker.R"
os.system(cmd)




time=float(time.time() - start_time)
hours=int(time/3600)
minutes=((time/3600)-hours)*60

print("\n\n\n\n\n\n\n\n\n")
print("####################\n")
print ("Analysis completed in %i hours and %f minutes" %(hours,
minutes))
print("####################\n")


"""
