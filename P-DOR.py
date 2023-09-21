#!/usr/bin/env python3.8

## LIBRARIES
import argparse
from Bio import SeqIO
import sys
import time
import datetime
import glob
import os
from PDOR_lib import core_snps_2_fasta
from PDOR_lib import core_snps_list_path

## FUNCTIONS
class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def parse_args():

    parser=MyParser(usage='%(prog)s -q QUERY_FOLDER -sd SD_SKETCH -ref REFERENCE_GENOME -snp_thr SNP_THRESHOLD [options]',
        prog='P-DOR.py',
        argument_default=argparse.SUPPRESS,
        epilog="According to the legend, P-dor is the Son of K-mer, but it also likes SNPs")
    requiredNamed = parser.add_argument_group('Input data (all required, except when in TEST mode)')
    requiredNamed.add_argument('-q', help='query folder containing genomes in .fna format',metavar="<dirname>",dest="query_folder",type=str,required=False,default=None)
    requiredNamed.add_argument("-sd", help="Source Dataset (SD) sketch file",metavar="<dirname>",dest="db_sketch", required=False,default=None)
    requiredNamed.add_argument("-ref", help="reference genome", dest="ref", metavar="<filename>",required=False,default=None)
    requiredNamed.add_argument("-snp_thr", dest="snp_threshold", help="Threshold number of SNPs to define an epidemic cluster",required=False,metavar="<int>",default=None)

    optional = parser.add_argument_group('Additional arguments')
    optional.add_argument("-amrf", help="if selected, P-DOR uses amrfinder-plus to search for antimicrobial resistance and virulence genes in the entire Analysis Dataset [Default: %(default)s]",action="store_true",required=False)
    optional.add_argument("-meta", help="metadata file; see example file for formatting [Default: %(default)s]",default=None,required=False)
    optional.add_argument("-sd_folder", help="folder containing the genomes from which the Source Dataset (SD) sketch was created [Default: P-DOR downloads the background genomes from BV-BRC]",required=False,default="",dest="bkg_folder",nargs="?")
    optional.add_argument("-borders", type=int, help="length of the regions at the contig extremities from which SNPs are not called [Default: %(default)s]",metavar="<int>",default=20)
    optional.add_argument("-min_contig_length", type=int, help="minimum contig length to be retained for SNPs analysis [Default: %(default)s]",metavar="<int>",default=500)
    optional.add_argument("-snp_spacing", type=int, help="Number of bases surrounding the mutated position to call a SNP [Default: %(default)s]",metavar="<int>",default=10)
    optional.add_argument('-species', help="full species name of the genomes in the query folder [Default: %(default)s]\nNeeds to be written within single quotes, e.g.: -species 'Klebsiella pneumoniae'.\nThis option is used for a quicker and more precise gene search with AMRFinderPlus", default="")
    optional.add_argument('-n', type=int,help='Maximum closest genomes from Source Dataset (SD) [Default: %(default)s]',metavar="<int>",dest="near",default=20)
    optional.add_argument('-t', type=int,help='number of threads [Default: %(default)s]',metavar="<int>",dest="threads",default=2)
    parser.add_argument('-TEST', action= 'store_true',help='TEST MODE. Warning: this option overrides all other arguments', required=False)
    parser.add_argument('-v', '--version', action='version', version='P-DOR v1.1')

    parser.set_defaults(amrf=False)
    parser.set_defaults(TEST=False)
    args = parser.parse_args()
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()
    elif (not args.TEST and (args.query_folder is None or args.db_sketch is None or args.ref is None or args.snp_threshold is None)):
        parser.error("\n\nERROR!! If not in TEST mode, the following arguments are required: -q, -sd, -ref, -snp_thr\n\n")
    return args

def logfile():
	with open("PDOR.log","w") as p:
		for arg, value in sorted(vars(args).items()):
			p.write("%s\t%s\n" %(arg,value))

def parse_species(sp):
    orglist=["Acinetobacter_baumannii", "Burkholderia_cepacia", "Burkholderia_pseudomallei", \
    "Campylobacter", "Clostridioides_difficile", "Enterococcus_faecalis", "Enterococcus_faecium", \
    "Escherichia", "Klebsiella_oxytoca", "Klebsiella_pneumoniae", "Neisseria_gonorrhoeae", \
    "Neisseria_meningitidis", "Pseudomonas_aeruginosa", "Salmonella", "Staphylococcus_aureus", \
    "Staphylococcus_pseudintermedius", "Streptococcus_agalactiae", "Streptococcus_pneumoniae", \
    "Streptococcus_pyogenes", "Vibrio_cholerae"]
    organism=""
    if sp.strip()!="":
        test="_".join(sp.strip().split())
        if test in orglist:
            organism=test
        else:
            test=sp.strip().split()[0]
            if test in orglist:
                organism=test
    return organism

def check_res_vir(threads,AD_folder,i_ref):
	print ("Checking for resistance and virulence genes...")
	path_res=path_dir+"/"+Results_folder_name+"/"+AD_folder
	os.chdir(path_res)
	os.system("cp ../%s ." %(i_ref))
	os.system("amrfinder --update")
	organism=parse_species(args.species)
	if organism=="":
		os.system("ls *fna | parallel -j %i 'amrfinder -n {} --plus --threads 1 -o {}.txt --name {}'" %(threads))
	else:
		os.system("ls *fna | parallel -j %i 'amrfinder -n {} --plus --threads 1 -o {}.txt -O %s --name {}'" %(threads,organism))
	os.system("cat *fna.txt >summary_resistance_virulence.txt")
	os.system("rm *fna.txt")
	os.system("rm %s/%s" %(path_res,i_ref.split("/")[-1]))

	os.system("mv summary_resistance_virulence.txt %s/%s" %(path_dir,Results_folder_name))
	os.chdir("..")

def Mummer_snp_call(threads,i_ref,Align_folder):
	print("Aligning genomes with Mummer4...\n")
	os.chdir(path_dir+"/"+Results_folder_name+"/"+Align_folder)
	os.system('ls *fna | parallel -j %i "nucmer --maxgap=500 -p {} %s {}; delta-filter -1 {}.delta > {}_filtered.delta; show-snps -H -C -I -r -T {}_filtered.delta | cut -f1,2,3 >{}.snp"' %(threads,i_ref))
	os.system('ls *_filtered.delta | parallel -j %i "dnadiff -d {} -p {}_info >/dev/null 2>&1"' %(threads))
	reportFiles=glob.glob("*.report")
	rFdict={}
	for rF in reportFiles:
		rFcontent=open(rF,"r").read().split("\n")
		rFdict[rF.strip().replace(".fna_filtered.delta_info.report","")]=float([t.strip().split()[-1].split("(")[1].split("%")[0] for t in rFcontent if "AlignedBase" in t][0])
	with open("Coverage_report.txt","w") as rFout:
		sleep=False
		for t in rFdict.keys():
			if rFdict[t]<70:
				sleep=True
				print("\n\nWARNING! Genome %s aligns poorly to the REFERENCE GENOME (%f %% horizontal coverage). Please, check!\n\n" %(t,rFdict[t]))
			rFout.write("%s\t%f\n" %(t,rFdict[t]))
		if sleep==True:
			time.sleep(10)
	#sys.exit()
	os.chdir(path_dir+"/"+Results_folder_name)
	os.system("ls %s/*.snp >SNP_genomes.list" %(Align_folder))
	core_snps_list_path("SNP_genomes.list", snp_spacing, "SNP_positions")
	core_snps_2_fasta("SNP_genomes.list", "SNP_positions.core.list",i_ref,"SNP_alignment")
	if os.stat("SNP_alignment.core.fasta").st_size > 10:
		print ("Core-SNPs aligmment detected!\n")
	else:
                sys.exit("\nERROR: Core-SNPs alignment not present, check your input files!\n")


### MAIN
print('\n\n      (/T\)  _____\n      (-,-) ()____)\n      \(o)/   -|3       BEHOLD! This is P-DOR!\n     /="""===//||\n,OOO//|   |    ""\nO:O:O LLLLL\n\OOO/ || ||\n     C_) (_D\n\n')
args = parse_args()
time_now=datetime.datetime.now()
datestamp=time_now.strftime("%Y-%m-%d_%H-%M-%S")
start_time = time.time()
prog_dir = os.path.abspath( os.path.dirname( sys.argv[0]))
Results_folder_name="Results_%s" %(datestamp)

##CHECK IF TEST MODE
if args.TEST==True:
    Results_folder_name="Test_%s" %(datestamp)
    args.amrf=False
    args.bkg_folder='%s/data/test_DB/' %(prog_dir)
    args.borders=20
    args.db_sketch='%s/data/sketches.msh' %(prog_dir)
    args.meta='%s/data/sample_metadata_table.txt' %(prog_dir)
    args.min_contig_length=500
    args.near=5
    args.query_folder='%s/data/test_query/' %(prog_dir)
    args.ref='%s/data/NJST258.fna' %(prog_dir)
    args.snp_spacing=10
    args.snp_threshold='20'
    args.species='Klebsiella pneumoniae'
    args.threads=2
    print("\n\nTesting AMRFinderPlus function on a single genome")
    os.system("amrfinder --update")
    organism=parse_species(args.species)
    os.system("amrfinder -n %s --plus --threads %i -o test.amrftest -O %s --name REF" %(args.ref,args.threads,organism))
    os.system("head -n2 test.amrftest")
    os.system("rm test.amrftest")

os.mkdir(Results_folder_name)
logfile()
path_dir = os.path.abspath( os.path.dirname( Results_folder_name))


os.system("mv PDOR.log %s" %(Results_folder_name))

query_folder=args.query_folder
db_sketch=os.path.abspath(args.db_sketch)
ref=args.ref
threads=args.threads
nearest=args.near
borders=args.borders
contig_length=args.min_contig_length
snp_spacing=args.snp_spacing
bkg_folder=args.bkg_folder
if args.meta is not None:
    ABSmeta=os.path.abspath(args.meta)

ref=os.path.abspath(ref)
ABS_query_folder=os.path.abspath(query_folder)

exts=['*.fasta', '*.fna', '*.fa']

files = [f for ext in exts for f in glob.glob(os.path.join(ABS_query_folder, ext))]

genomes_query = [i.strip().split("/")[-1] for i in files]
#genomes_query=[i.replace(".fasta",".fna").replace(".fa",".fna") for i in genomes_query]

print ("\nChecking input formats...\n")
### CHECK N 1
# REF has 1 contig
num_cont_c1=int(os.popen('grep -c ">" %s' %(ref)).read().strip())
if num_cont_c1>1:
	sys.exit('\nERROR: the reference file needs to contain only one sequence!')
### CHECK N 2
# REF is in fasta format
elif num_cont_c1==0:
	sys.exit('\nERROR: the reference file needs to be a FASTA file!')
### CHECK N 3
if db_sketch.split(".")[-1]!="msh":
	sys.exit("\n\nERROR: please check you sketch file!!!")
### CHECK N 4
for gen in genomes_query:
	if gen.endswith(".fasta") or gen.endswith(".fna") or gen.endswith(".fa"):
		n_genomes=len(genomes_query)
		print ("Detected %i genomes in the query folder\n" %n_genomes)
		break
	else:
		sys.exit("\n\nERROR: your query genome folder does not contain fasta files!!!")

### format REF
alt_ref=SeqIO.read(ref,"fasta")
alt_ref.id="REF_%s" %(str(alt_ref.id).split(".")[0])
SeqIO.write(alt_ref,"%s/%s.fna" %(Results_folder_name,alt_ref.id),"fasta")

for gen_names in genomes_query:
	ext=gen_names.strip().split(".")[-1]
	fname=".".join(gen_names.strip().split(".")[:-1])
	if (ext == "fa" or ext == "fasta") :
		os.rename ("%s/%s" %(ABS_query_folder,gen_names), "%s/%s.fna" %(ABS_query_folder,fname))

NEAREST=[]
for i in genomes_query:
	cmd="mash dist %s/%s %s -p %i 2>/dev/null | sort -gr -k5 | head -n %i | cut -f2"  %(ABS_query_folder,i,db_sketch,threads,nearest)
	NEAREST=NEAREST+(os.popen(cmd).read().strip().strip('\n').split('\n'))

NEAREST=list(set(NEAREST))
NEAREST=[i.split("/")[-1] for i in NEAREST]
print ("Sketches completed!\n")

os.chdir(Results_folder_name)
with open("query_genome_list.txt","w") as q:
	for i in genomes_query:
		q.write("%s\n" %(i))
with open("backgroud_genome_list.txt","w") as q:
	for i in NEAREST:
		q.write("%s\n" %(i))
os.mkdir("Background")

if args.bkg_folder:
	print ("Retrieving nearest genomes from the %s folder\n" %args.bkg_folder)
	os.chdir(path_dir)
	background_folder=os.path.abspath(bkg_folder)
	os.chdir(Results_folder_name)
	for N in NEAREST:
		cmd="cp %s/%s Background/" %(background_folder,N)
		os.system(cmd)
else:
	os.chdir("Background")
	print ("Retrieving nearest genomes from the BV-BRC Database...\n")
	for N in NEAREST:
		N=N.replace(".fna","").replace(".fasta","").replace(".fa","")
		cmd = "wget -qN ftp://ftp.bvbrc.org/genomes/%s/%s.fna" %(N,N)
		os.system(cmd)
	os.chdir(path_dir+"/"+Results_folder_name)

os.mkdir("Align")
os.mkdir("Analysis_Dataset")

for i in glob.glob("Background/*"):
	name=i.strip().split("/")[-1]
	oF=open("Align/BD_%s" %(name),"w")
	origsize=0
	trimmedsize=0
	for x in SeqIO.parse(i,"fasta"):
		origsize+=len(str(x.seq))
		if len(str(x.seq))>=contig_length+2*borders:
			x.seq=x.seq[borders:-borders]
			trimmedsize+=len(str(x.seq))
			SeqIO.write(x,oF,"fasta")
	oF.close()
	if trimmedsize<float(origsize)*0.90:
		print("\nDatabase genome %s was too short after contig polishing and was excluded from the Analysis Dataset (AD)" %(name))
		os.system("rm Align/BD_%s" %(name))
	else:
		os.system("cp %s Analysis_Dataset/BD_%s" %(i,name))

for i in glob.glob("%s/*" %(ABS_query_folder)):
	name=i.strip().split("/")[-1]
	oF=open("Align/%s" %(name),"w")
	origsize=0
	trimmedsize=0
	for x in SeqIO.parse(i,"fasta"):
		origsize+=len(str(x.seq))
		if len(str(x.seq))>=contig_length+2*borders:
			x.seq=x.seq[borders:-borders]
			trimmedsize+=len(str(x.seq))
			SeqIO.write(x,oF,"fasta")
	oF.close()
	if trimmedsize<float(origsize)*0.90:
		sys.exit("\n\nERROR: your query genome %s was too short after contig polishing, please check your genome quality! EXECUTION ABORTED" %(name))
	else:
		os.system("cp %s Analysis_Dataset/%s" %(i,name))

aln_folder_size=os.popen("ls Align/ | wc -l").read().strip()
os.system("rm -rf Background")

if args.amrf==True:
    check_res_vir(threads,"Analysis_Dataset","%s.fna" %(alt_ref.id))

alt_ref_abs=os.path.abspath("%s.fna" %(alt_ref.id))
Mummer_snp_call(threads,alt_ref_abs,"Align")
os.mkdir("SNPs")
os.system("mv Align/*.snp SNPs/")
os.system("mv Align/Coverage_report.txt .")
snp_num=os.popen("ls SNPs/*.snp | wc -l").read().strip()
if int(snp_num)!=int(aln_folder_size):
    sys.exit("\n\nERROR! Something went wrong! Not all .snp files were produced")
os.system("rm -rf Align")

print("Performing phylogenetic reconstruction...\n")
cmd="iqtree -s SNP_alignment.core.fasta -pre SNP_alignment -m MFP+GTR+ASC -bb 1000 -nt %i" %(threads)
os.system(cmd)

os.system("$CONDA_PREFIX/bin/Rscript %s/Snpbreaker.R SNP_alignment.core.fasta %s" %(prog_dir,args.snp_threshold))


if args.amrf==True:
    emptyAMRF=os.popen("cat summary_resistance_virulence.txt | sort | uniq | wc -l").read().strip()
    if int(emptyAMRF)<=1:
        os.system("rm summary_resistance_virulence.txt")

os.system("$CONDA_PREFIX/bin/Rscript %s/annotated_tree.R SNP_alignment.treefile" %(prog_dir))

if args.meta is not None:
	os.system("$CONDA_PREFIX/bin/Rscript %s/contact_network.R %s" %(prog_dir,ABSmeta))

###ORGANIZE OUTPUT
os.mkdir("Other_output_files")
os.system("mv SNPs Other_output_files/SNPs")
os.system("mv Analysis_Dataset Other_output_files/Analysis_Dataset/")
os.system("mv SNP_alignment* Other_output_files/")

os.system("mv backgroud_genome_list.txt Other_output_files/")
os.system("mv clusters_manual_threshold_%s_snps.csv Other_output_files/" %(args.snp_threshold))
os.system("mv query_genome_list.txt Other_output_files/")
os.system("mv snp_distance_matrix.tsv Other_output_files/")
os.system("mv SNP_genomes.list Other_output_files/")
os.system("mv SNP_positions.core.list Other_output_files/")
if args.amrf==True:
    os.system("mv summary_resistance_virulence.txt Other_output_files/")

os.chdir(path_dir)
os.system("rm %s" %(alt_ref_abs))

time=float(time.time() - start_time)
hours=int(time/3600)
minutes=((time/3600)-hours)*60
print ("\n\n\n\n")
print ("####################\n")
print ("Analysis completed in %i hours and %f minutes" %(hours,minutes))
print ("Results are stored in folder %s/%s" %(path_dir,Results_folder_name))
print ("####################\n")
