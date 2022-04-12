## impostare cosi: cartella genomes (sono i query) || path assoluta della cartella con il db || cartella Purple e una cartella ref in cui si inserisce la reference
### a parte la cartella database il resto deve essere tutto all'interno dello stesso folder

### esempio ls

### genomes ref Purple
### mettere script nella cartella genomes e lanciarlo da li con ### CMD= python miniskynet.py . /mnt/pinkydisk/gherard/Skynet_rib/prova/db/ 20 3 20 20


import sys
import time
import datetime
import os
import glob
from Bio import SeqIO
import shutil
import argparse
import subprocess

from Skynet_SNPs_library_4_0 import core_snps_2_fasta
from Skynet_SNPs_library_4_0 import core_snps_two_files
from Skynet_SNPs_library_4_0 import core_snps_list_path
from Skynet_SNPs_library_4_0 import snp_calling_mauve
from Skynet_SNPs_library_4_0 import card_resistance_genes


def rename_input_files(query_folder):
	Qglobber="%s/*" %(query_folder) #cambiare la desinenza in modo da aggiungere poi mashed
	genomes_query = glob.glob(Qglobber)
	for q in genomes_query:
		if q.split(".")[-1] in ("fasta","fa"):
			new=".".join(q.split(".")[:-1])+".fna"
			cmd="mv %s %s" %(q,new)
			os.system(cmd)


timefunct=datetime.datetime.now()
Results_folder_name="Results_%s-%s-%s_%s-%s" %(str(timefunct.year),str(timefunct.month),str(timefunct.day),str(timefunct.hour),str(timefunct.minute)) ####FIX zeros in p-dor
os.mkdir(Results_folder_name)

start_time = time.time()

query_folder=sys.argv[1].strip() #path con i genomi query
db_folder=sys.argv[2].strip() #path database genomi target
align_ref=sys.argv[3].strip()
threads=int(sys.argv[4].strip()) #MASH threads
Nnearest=int(sys.argv[5].strip()) #numero massimo di genomi vicini da identificare per ogni genoma query


rename_input_files(query_folder)

Qglobber="%s/*.fna" %(query_folder)
genomes_query = glob.glob(Qglobber)

out=[]
query_list=[]
for i in genomes_query:
	i=i.strip().split("/")[-1]
	cmd="mash dist %s/%s %s/*.msh -p %i | sort -gr -k5 | head -n %i | cut -f2" %(query_folder,i,db_folder,threads,Nnearest)
	output = subprocess.check_output(cmd,shell=True)
	out+=str(output).strip("\n'").strip("b'").replace("\\n","\n").replace("\\t","\t").strip("\n").split("\n")
out=list(set(out))
print(out)


os.system("rm -rf tmp")
os.mkdir("tmp")

for i in out:
	cmd="rsync %s/%s tmp" %(db_folder,i)
	os.system(cmd) 

os.system("mmv 'tmp/*' 'tmp/DB_#1'")

for i in genomes_query:
	cmd="rsync %s tmp" %(i)
	os.system(cmd) 


#######################################################################################

#LANCIA MUMMER,SKynet,RAxML

os.chdir("tmp/")
cmd='ls | parallel -j %i "nucmer -p {} ../%s {}; show-snps -H -C -I -r -T {}.delta | cut -f1,2,3 > {}.snp "' %(threads,align_ref)
os.system(cmd)
os.chdir("..")
os.system("ls tmp/*.snp > SNP_genomes.list")


############SKYNET LAUNCH


core_snps_list_path("SNP_genomes.list", 2, "SNP_positions")

core_snps_2_fasta("SNP_genomes.list", "SNP_positions.core.list",align_ref,"SNP_alignment")



#####RAxML LAUNCH

os.system("rm *.SNP_Phylo")

cmd="raxmlHPC-PTHREADS -f a -x 12345 -p 12345 -# 100 -m ASC_GTRGAMMA -s SNP_alignment.core.fasta -n SNP_Phylo -T %i --no-bfgs --asc-corr=lewis" %(threads)
os.system(cmd)



############ R Plot

cmd="Rscript Snp_breaker.R"
os.system(cmd)

########## Folders

FF="%s/Phylo" %(Results_folder_name)
os.mkdir(FF)
FF="%s/SNPs_Plots" %(Results_folder_name)
os.mkdir(FF)

cmd="mv RAxML* %s/Phylo" %(Results_folder_name)
os.system(cmd)

cmd="mv SNP_alignment.core.fasta %s/Phylo" %(Results_folder_name)
os.system(cmd)

cmd="mv *.eps %s/SNPs_Plots" %(Results_folder_name)
os.system(cmd)

cmd="mv *snp_distance_tab.xlsx %s/SNPs_Plots" %(Results_folder_name)
os.system(cmd)

cmd="mv *.png %s/SNPs_Plots" %(Results_folder_name)
os.system(cmd)

cmd="mv *.pdf %s/SNPs_Plots" %(Results_folder_name)
os.system(cmd)

cmd="cp %s/Phylo/RAxML_bipartitionsBranchLabels.* %s/" %(Results_folder_name,Results_folder_name)
os.system(cmd)

cmd="cp %s/SNPs_Plots/overlap_plot.png %s/" %(Results_folder_name,Results_folder_name)
os.system(cmd)
