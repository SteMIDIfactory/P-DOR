#!/usr/bin/python3

## LIBRARIES

import argparse
from Bio import SeqIO
import sys
import time
import datetime
import glob
import os
#from os.path import join
#from glob import glob
from PDOR_lib import core_snps_2_fasta
from PDOR_lib import core_snps_two_files
from PDOR_lib import core_snps_list_path
from PDOR_lib import snp_calling_mauve


## FUNCTIONS

def parse_args():

        parser=argparse.ArgumentParser(description='''This is the P-DOR help ''',usage='%(prog)s [options]',
                prog='P-DOR.py',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                argument_default=argparse.SUPPRESS,
                epilog="According to the legend, P-dor is the Son of K-mer, but it also likes SNPs")

        requiredNamed = parser.add_argument_group('Input data')
        requiredNamed.add_argument('-q', help='MA QUESTO SCRIPT LO AVETE FATTO VOI?query folder containing genomes in .fna format',metavar="<dirname>",dest="query_folder",type=str,required=True)
        requiredNamed.add_argument("-db", help="background sketch file",metavar="<dirname>",dest="db_sketch", required=True)
        requiredNamed.add_argument("-ref", help="reference genome", dest="ref", metavar="<filename>",required=True)
        
        requiredNamed.add_argument("-snp_thr", dest="snp_threshold", help="Threshold number of SNPs to define an epidemic cluster: choices are an integer number or type 'infl' to calculate it by the inflection point of SNPs distribution ",required=True)
        

        optional = parser.add_argument_group('Additional arguments')
        optional.add_argument("-meta", help="metadata file; see example file for formatting",default=None,required=False)
        optional.add_argument("-gff", help="annotation file in gff format; if not specified (default) a dummy gff is generated",required=False,default="",nargs="?")
        optional.add_argument("-bkg_folder", help="folder containing the genomes from which the background sketch was created",required=False,default="",nargs="?")
        optional.add_argument("-borders", type=int, help="length of the regions at the contig extremities from which SNPs are not called",metavar="<int>",default=20)
        optional.add_argument("-min_contig_length", type=int, help="minimum contig length to be retained for SNPs analysis",metavar="<int>",default=500)
        optional.add_argument("-snp_spacing", type=int, help="Number of bases surrounding the mutated position to call a SNP",metavar="<int>",default=10)
        optional.add_argument('-n', type=int,help='Maximum closest genomes from database',metavar="<int>",dest="near",default=20)
        optional.add_argument('-t', type=int,help='number of threads',metavar="<int>",dest="threads",default=2)

        parser.add_argument('-v', '--version', action='version', version='P-DOR v1.0')

        args = parser.parse_args()
        snp_threshold=args.snp_threshold
        if snp_threshold is not None and snp_threshold.isdigit():
              snp_threshold = int(snp_threshold)

   
   
### CHECK N 1
   
        elif snp_threshold!="infl":             

              sys.exit('\nERROR: to set the inflation point the following argument is required: -snp_thr infl\n')


        return args



args = parse_args()


def logfile():
	with open("PDOR.log","w") as p:
		for arg, value in sorted(vars(args).items()):
			p.write("%s\t%s\n" %(arg,value))




def check_res_vir(threads):
	print ("checking for resistance and virulence genes...")
	os.chdir(path_dir+"/"+Results_folder_name+"/"+"Align")
	path_res=path_dir+"/"+Results_folder_name+"/"+"Align"
	cmd=("cp %s %s") %(os.path.abspath(ref),path_res)
	os.system(cmd)
	os.system("conda env list >path_pdor")
	path_inF=open("path_pdor","r")
	for i in path_inF.readlines():
		i=i.strip().split()
		
		try:
			if i[1].strip()=="*":
				abricate_db_path=i[2].strip()+"/db/all_db"
				abs_path=i[2].strip()+"/db"
				
				cmd="mkdir -p %s" %abricate_db_path
				os.system(cmd)
				cmd="cat %s/*/*sequences | perl -pe 's/[[:^ascii:]]//g' >%s/all_db/sequences" %(abs_path,abs_path)
				os.system(cmd)
				cmd="makeblastdb -in %s/sequences -title all_db -dbtype nucl -hash_index" %abricate_db_path
				os.system(cmd)

				os.system("ls *fna | parallel -j %i 'abricate {} --minid 80 --mincov 60 --db all_db >{}.report'" %(threads))
				os.system("cat *report >summary_resistance_virulence")
				os.system("rm *report")
				rm_ref=("rm %s/%s") %(path_res,ref.split("/")[-1])
				os.system(rm_ref)
				os.system("mv summary_resistance_virulence ../")


		except IndexError:
		
			pass



				

def Mummer_snp_call(threads,ref):
	#os.rename(".fasta",".fna").replace(".fa","fna")
	print("Aligning genomes with Mummer4...\n")
	os.chdir(path_dir+"/"+Results_folder_name+"/"+"Align")
	os.system('ls *fna| parallel -j %i "nucmer --maxgap=500 -p {} %s {}; delta-filter -1 {}.delta > {}_filtered.delta; show-snps -H -C -I -r -T {}_filtered.delta | cut -f1,2,3 >{}.snp"' %(threads,ref))
	os.chdir("..")
	os.system("ls Align/*.snp >SNP_genomes.list")
	core_snps_list_path("SNP_genomes.list", snp_spacing, "SNP_positions")
	core_snps_2_fasta("SNP_genomes.list", "SNP_positions.core.list",ref,"SNP_alignment")
	if os.stat("SNP_alignment.core.fasta").st_size > 10:
		print ("Core-SNPs aligmment detected...\n")
	else:
                sys.exit("\nERROR: Core-SNPs alignment not present, check your input files...\n")




def create_dummy_gff(ref):
	for index, record in enumerate(SeqIO.parse(ref,"fasta")):
		name=i.strip().split("/")[-1]
		if index==0:
			cmd="echo '%s	.	CDS	1	800' >%s.gff" %(record.id,name)
			os.system(cmd)
	return name+".gff"





### MAIN
args = parse_args()

timefunct=datetime.datetime.now()
Results_folder_name="Results_%s-%s-%s_%s-%s-%s" %(str(timefunct.year),str(timefunct.month),str(timefunct.day),str(timefunct.hour),
	str(timefunct.minute),str(timefunct.second))
os.mkdir(Results_folder_name)
logfile()
path_dir = os.path.abspath( os.path.dirname( Results_folder_name ) ) 

os.system("mv PDOR.log %s" %(Results_folder_name))



start_time = time.time()


query_folder=args.query_folder
db_sketch=args.db_sketch
ref=args.ref
threads=args.threads
nearest=args.near
borders=args.borders
contig_length=args.min_contig_length
snp_spacing=args.snp_spacing
bkg_folder=args.bkg_folder


ref=os.path.abspath(ref)
ABS_query_folder=os.path.abspath(query_folder)


#Qglobber="%s/*.fna" %(ABS_query_folder)
#genomes_query = glob.glob(Qglobber)

exts=['*.fasta', '*.fna', '*.fa']

files = [f for ext in exts for f in glob.glob(os.path.join(ABS_query_folder, ext))]


genomes_query = [i.strip().split("/")[-1] for i in files]
#genomes_query=[i.replace(".fasta",".fna").replace(".fa",".fna") for i in genomes_query]



### CHECK N 2

print ("Check input formats...\n")

ref_file=open(ref,"r")

for id in ref_file.readlines()[0:1]:
      
      if id[0].strip()!=">":
              sys.exit('\nERROR: the reference file need to be a FASTA file!!!')



### CHECK N 3

if db_sketch.split(".")[-1]!="msh":
	
	sys.exit("\n\nERROR: check you sketch file!!!")


### CHECK N 4


for gen in genomes_query:

	if gen.endswith(".fasta") or gen.endswith(".fna") or gen.endswith(".fa"):
		n_genomes=len(genomes_query)
		print ("Detected %i genomes in the query folder\n" %n_genomes)
	
		break
	
	else:
		sys.exit("\n\nERROR: your query genome folder does not contain fasta files!!!")

	


for gen_names in genomes_query:
	ext=gen_names.strip().split(".")[-1]
	if (ext == "fa" or ext == "fasta") :
		
		#print ("%s/%s","%s/%s.fna" %(query_folder,gen_names,query_folder,gen_names.strip(ext).strip(".")))
		os.rename ("%s/%s" %(query_folder,gen_names), "%s/%s.fna" %(query_folder,gen_names.strip(ext).strip("."))) 

	
NEAREST=[]
for i in genomes_query:
	cmd="mash dist %s/%s %s -p %i | sort -gr -k5 | head -n %i | cut -f2"  %(query_folder,i,db_sketch,threads,nearest)
	NEAREST=NEAREST+(os.popen(cmd).read().strip().strip('\n').split('\n'))
	#print (cmd)


#assert 1==2

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
	print ("Retrieving nearest genomes from the PATRIC-DB...\n")
	for N in NEAREST:
		N=N.replace(".fna","").replace(".fasta","").replace(".fa","")
		print(N)
		cmd='wget -qN "ftp://ftp.patricbrc.org/genomes/%s/%s.fna"' %(N,N)
		print(cmd)
		os.system(cmd)
	os.chdir("..")
	
os.mkdir("Align")





for i in glob.glob("Background/*"):
	name=i.strip().split("/")[-1]
	#trim and select contigs
	oF=open("Align/DB_%s" %(name),"w")
	for x in SeqIO.parse(i,"fasta"):
		if len(str(x.seq))>=contig_length+2*borders:
			x.seq=x.seq[borders:-borders]
			SeqIO.write(x,oF,"fasta")
	oF.close()

for i in glob.glob("%s/*" %(ABS_query_folder)):
	name=i.strip().split("/")[-1]
	#trim and select contigs
	oF=open("Align/%s" %(name),"w")
	for x in SeqIO.parse(i,"fasta"):
		if len(str(x.seq))>=contig_length+2*borders:
			x.seq=x.seq[borders:-borders]
			SeqIO.write(x,oF,"fasta")


	oF.close()



### INCOMPLETED HERE ???

#for i in glob.glob("Align"):
#os.system("cp %s/*.fna Align" %(ABS_query_folder))

###

#check_res_vir(threads)

Mummer_snp_call(threads,ref)

print("Performing phylogenetic reconstruction...\n")

cmd="raxmlHPC-PTHREADS -f a -x 12345 -p 12345 -# 100 -m ASC_GTRGAMMA -s SNP_alignment.core.fasta -n coreSNPs_Phylo.nwk -T %i --no-bfgs --asc-corr=lewis" %(threads)
os.system(cmd)


os.system("Rscript ../Snpbreaker.R SNP_alignment.core.fasta %s" %(args.snp_threshold))
os.system("Rscript ../annotated_tree.R RAxML_bipartitionsBranchLabels.coreSNPs_Phylo.nwk")



if args.meta is not None:
	
	os.system("Rscript ../contact_network.R ../%s" %(args.meta))





os.chdir("..")


time=float(time.time() - start_time)
hours=int(time/3600)
minutes=((time/3600)-hours)*60

#outF="Results_%s" + datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
#os.makedirs(outF)




"""
os.makedirs("phylogeny outbreak_detection")
os.system("mv RAxML_* phylogeny/")
os.system("mv annotated_tree.svg phylogeny/"

os.system("mv *csv outbreak_detection/")
os.system("mv *svg outbreak_detection/")
"""


print ("\n\n\n\n\n\n\n\n\n")
print ("####################\n")
print ("Analysis completed in %i hours and %f minutes" %(hours,minutes))
print ("####################\n")

