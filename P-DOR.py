#!/usr/bin/python3

## LIBRARIES

import argparse
from Bio import SeqIO
import sys
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

	parser=argparse.ArgumentParser(description='''This is the P-DOR help ''',usage='%(prog)s [options]',
		prog='P-DOR.py',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		argument_default=argparse.SUPPRESS,
		epilog="According to the legend, P-dor is the Son of K-mer, but it also likes SNPs")

	requiredNamed = parser.add_argument_group('Input data')
	requiredNamed.add_argument('-q', help='query folder containing genomes in .fna format',metavar="<dirname>",dest="query_folder",type=str,required=True)
	requiredNamed.add_argument("-db", help="background sketch file",metavar="<dirname>",dest="db_sketch", required=True)
	requiredNamed.add_argument("-ref", help="reference genome", dest="ref", metavar="<filename>",required=True)


	optional = parser.add_argument_group('Additional arguments')
	optional.add_argument("-meta", help="metadata file; see example file for formatting",default=None,required=False)
	optional.add_argument("-gff", help="annotation file in gff format; if not specified (default) a dummy gff is generated",required=False,default="",nargs="?")
	optional.add_argument("-bkg_folder", help="folder containing the genomes from which the background sketch was created; if not specified (default) nearest genomes are downloaded from the PATRIC-DB",required=False,default="",nargs="?")
	optional.add_argument("-call", help="Snps calling method",required=False, choices=['purple', 'mummer'],default="mummer")
	optional.add_argument("-snp_thr", type=int,help="Threshold number of SNPs to define an epidemic cluster",required=False,default=20)
	optional.add_argument('-n', type=int,help='Maximum nearest genomes from database',metavar="<int>",dest="near",default=20)
	optional.add_argument('-t', type=int,help='number of threads',metavar="<int>",dest="threads",default=10)

	parser.add_argument('-v', '--version', action='version', version='P-DOR v0.1')

	args = parser.parse_args()
	return args


def check_folder(folder):
	folder_path=os.path.abspath(folder)
	if os.path.isdir(folder):
		sys.exit("\n\nERROR: the {} already exists!".format(folder_path))
		#sys.stderr.write("%s %s\n" % (foo, bar))
		#raise Exception("{} already exists".format(folder_path)) ## \nuse -f overwrite???


def Mummer_snp_call(threads,ref):
	print("Aligning genomes with Mummer4")
	os.chdir("Align")
	os.system('ls | parallel -j %i "nucmer -p {} %s {}; show-snps -H -C -I -r -T {}.delta | cut -f1,2,3 >{}.snp"' %(threads,ref))
	os.chdir("..")
	os.system("ls Align/*.snp >SNP_genomes.list")
	core_snps_list_path("SNP_genomes.list", 2, "SNP_positions")
	core_snps_2_fasta("SNP_genomes.list", "SNP_positions.core.list",ref,"SNP_alignment")


def Purple_snp_call(ref,cpus,Results_folder_name):
	print("Aligning genomes with Purple (Mauve-based)")
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
	condaprefix=os.popen("echo $CONDA_PREFIX").read().strip()
	libpath="%s/lib/perl5/site_perl/5.22.0/" %(condaprefix)
	os.environ['PERL5LIB'] = libpath
	cmd="perl Purple_v1.22.2.pl -r %s -g %s -f Align -cpu %i" %(localref,gffname,cpus)
	os.system(cmd)
	os.system("mv Align_OUT/Results/Fasta_SNPs/Align.core_SNPs.fasta ../SNP_alignment.core.fasta")
	os.chdir("..")



def create_dummy_gff(ref):
	for index, record in enumerate(SeqIO.parse(ref,"fasta")):
		name=i.strip().split("/")[-1]
		if index==0:
			cmd="echo '%s	.	CDS	1	800' >%s.gff" %(record.id,name)
			os.system(cmd)
	return name+".gff"


def check_gff(ref):

	if args.gff:
		print ("Using %s annotation file..." %args.gff)
		gffname=args.gff

	else:
		gffname=create_dummy_gff(ref)
		print ("Generating dummy annotation file...")
	return gffname

### MAIN
args = parse_args()

timefunct=datetime.datetime.now()
Results_folder_name="Results_%s-%s-%s_%s-%s" %(str(timefunct.year),str(timefunct.month),str(timefunct.day),str(timefunct.hour),
	str(timefunct.minute))
os.mkdir(Results_folder_name)

start_time = time.time()


query_folder=args.query_folder
db_sketch=args.db_sketch
ref=args.ref
threads=args.threads
nearest=args.near
bkg_folder=args.bkg_folder

ref=os.path.abspath(ref)
ABS_query_folder=os.path.abspath(query_folder)
Qglobber="%s/*.fna" %(ABS_query_folder)
genomes_query = glob.glob(Qglobber)
genomes_query= [i.strip().split("/")[-1] for i in genomes_query]


NEAREST=[]
for i in genomes_query:
	cmd="mash dist %s/%s %s -p %i | sort -gr -k5 | head -n %i | cut -f2"  %(query_folder,i,db_sketch,threads,nearest)
	NEAREST=NEAREST+(os.popen(cmd).read().strip().strip('\n').split('\n'))

NEAREST=list(set(NEAREST))

NEAREST=[i.split("/")[-1] for i in NEAREST]

print ("Sketches completed...")

os.chdir(Results_folder_name)

with open("query_genome_list.txt","w") as q:
	for i in genomes_query:
		q.write("%s\n" %(i))

with open("backgroud_genome_list.txt","w") as q:
	for i in NEAREST:
		q.write("%s\n" %(i))

os.mkdir("Background")

if args.bkg_folder:
	print ("Retrieving nearest genomes from the %s folder" %args.bkg_folder)
	os.chdir("..")
	background_folder=os.path.abspath(bkg_folder)
	os.chdir(Results_folder_name)
	for N in NEAREST:
		cmd="cp %s/%s Background/" %(background_folder,N)
		os.system(cmd)
else:
	os.chdir("Background")
	print ("Retrieving nearest genomes from the PATRIC-DB")
	for N in NEAREST:
		N=N.replace(".fna","").replace(".fasta","")
		print(N)
		cmd='wget -qN "ftp://ftp.patricbrc.org/genomes/%s/%s.fna"' %(N,N)
		print(cmd)
		os.system(cmd)
	os.chdir("..")

os.mkdir("Align")


for i in glob.glob("Background/*"):
	name=i.strip().split("/")[-1]
	os.system("cp %s Align/DB_%s" %(i,name))

os.system("cp %s/*.fna Align" %(ABS_query_folder))



if args.call == 'mummer':
	Mummer_snp_call(threads,ref)
elif args.call == 'purple':
	Purple_snp_call(ref,threads,Results_folder_name)


print("Performing phylogeny reconstruction")

cmd="raxmlHPC-PTHREADS -f a -x 12345 -p 12345 -# 100 -m ASC_GTRGAMMA -s SNP_alignment.core.fasta -n SNP_Phylo.nwk -T %i --no-bfgs --asc-corr=lewis" %(threads)
os.system(cmd)




os.system("Rscript ../Snpbreaker.R SNP_alignment.core.fasta %i" %(args.snp_thr))
if args.meta is not None:
	os.system("Rscript ../Timeline.R ../%s" %(args.meta))
	os.system("Rscript ../Tubeplot.R ../%s" %(args.meta))

os.chdir("..")

time=float(time.time() - start_time)
hours=int(time/3600)
minutes=((time/3600)-hours)*60

print("\n\n\n\n\n\n\n\n\n")
print("####################\n")
print ("Analysis completed in %i hours and %f minutes" %(hours,minutes))
print("####################\n")
