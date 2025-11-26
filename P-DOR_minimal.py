
## LIBRARIES
import argparse
from Bio import SeqIO
import time
import datetime
import glob
import sys
import re, sys, io, os

parser = argparse.ArgumentParser()

## Mandatory

parser.add_argument("-f", "--genomes_folder", help="Genomes folder (fastA)", required=True)
parser.add_argument("-r", "--reference", help="Reference genome", required=True)
parser.add_argument("-o", "--output", help="Output name", required=True)

## Non-mandatory

# Analyses

parser.add_argument("-ws", "--window_size", help="Half conserved window size", required=False, default = 5)
parser.add_argument("-pr", "--print_reference", help="Print reference in the output file", required=False, default = "y")
parser.add_argument("-pos", "--snp_pos_list", help="List of the SNPs positions", required=False, default = "auto")
parser.add_argument("-mask", "--mask_ref_pos", help="Mask reference genome positions to SNPs calling", required=False, default = "none")
parser.add_argument("-thr", "--threads", help="Number of threads", required=False, default = 5)

args = parser.parse_args()

genomes_folder = args.genomes_folder
reference = args.reference
output = args.output

half_window_size = args.window_size
print_ref = args.print_reference
snp_pos_list = args.snp_pos_list
to_mask_pos = args.mask_ref_pos
threads = int(args.threads)

outname_snp_fasta="%s.assembly_core_snp.fasta" %(output.strip())

#######
# DEF
#######

def core_snps_list_path(list_path, window, output):

	# caricare un input (lista file .snp)

	list_p=open(list_path)

	# crea il dizionario vuoto e un diz della reference

	diz_pos={}
	diz_pos_ref={}

	# scorre riga per riga
	for path_line in list_p:
		path = path_line.split('\t')[1]
		p=open(path.strip())
	
		# Reading text
	
		for x in p:
   
		# rimuovere le posizioni con il .
			if not re.search(r'\.',x):
				pos=int(x.split('\t')[0])
				base_ref=x.split('\t')[1]
				base_org=x.split('\t')[2]
				letters = set('atcgATCG')   
			
	# se il diz non presenta una posizione la aggiunge ed assegna valore 0

				try:
					diz_pos[pos]
				except:
					diz_pos.update({int(pos):'0'})
   
	# altrimenti verifica che base_ref e base_org siano != ATCG in quella posizione, se verificato mette in quella posizione valore 1  
   
			if not ((letters & set(base_org)) and (letters & set(base_ref))):
				diz_pos.update({int(pos):'1'})

	# Close opened file
		p.close()

	keys_sort=sorted(diz_pos.keys())
	core=[]

	if ((keys_sort[1]-keys_sort[0]) > int(window)) and (diz_pos[keys_sort[0]] == "0"):
			core.append(keys_sort[0])

	for i in range(1,len(keys_sort)-2):
		if (keys_sort[i]-keys_sort[i-1]) > int(window) and (keys_sort[i+1]-keys_sort[i]) > int(window) and (diz_pos[keys_sort[i]] == "0"):
			core.append(keys_sort[i])

	if (keys_sort[-1]-keys_sort[-2]) > int(window) and (diz_pos[keys_sort[-1]] == "0"):
			core.append(keys_sort[-1])

	# stampare in un file la lista core

	outname_snp_list="%s.assembly_core_snp.pos_list" %(output.strip())

	with open(outname_snp_list, "w") as file:
		for line in core:
			file.write(f'{line}\n')


def core_snps_2_fasta(snp_paths_list, core_snp_pos, reference_fasta, output, print_ref):

# Carica posizioni core

	core_file=open(core_snp_pos,'r')
	core=[]
	for i in core_file:
		i2=i.replace("\n","")
		core.append(i2)


# Carica le basi della reference

	count=0
	for seq_record in SeqIO.parse(reference_fasta,"fasta"):
		list_pos_ref=(list(seq_record.seq.lower()))
		reference=seq_record.id
		count=count+1   

		if count>1: 
			print ('Error')
			break	

# Prerara files output

	outname_snp_fasta="%s.assembly_core_snp.fasta" %(output.strip())
	outfasta=open(outname_snp_fasta,"w")

# Print fasta reference

	if (print_ref == "y"):
		s=''
		for i in core:
			list_pos_ref[int(i)-1].strip()
			s+=list_pos_ref[int(i)-1]

		outfasta.write(">%s\n" %(reference.strip()))
		outfasta.write("%s\n" %(s.strip()))

# Per ogni file, carica basi organismo e stampa il fasta

	list_path2=open(snp_paths_list,'r')

	for path2_line in list_path2:

		basi={}
		path2 = path2_line.split('\t')[1]
		organism_name = path2_line.split('\t')[0]
		p2=open(path2.strip())

		for x2 in p2: 
			if not re.search(r'\.',x2):
				pos2=int(x2.split('\t')[0])
				base_org2=x2.split('\t')[2]
				basi.update({int(pos2):base_org2.strip()})
        
		# Stampa fasta
		s=''
		for i in core:
			try:
				basi[int(i)]
				s+=basi[int(i)]
			except:
				list_pos_ref[int(i)-1].strip(),
				s+=list_pos_ref[int(i)-1]

		outfasta.write(">%s\n" %(organism_name.strip()))
		outfasta.write("%s\n" %(s.strip()))

		p2.close()

	outfasta.close()

#############
# Run NUCmer
#############

os.system('ls %s/* | parallel -j %i "nucmer --maxgap=500 -p {} %s {}; delta-filter -1 {}.delta > {}_filtered.delta; show-snps -H -C -I -r -T {}_filtered.delta | cut -f1,2,3 > {}.snp"' %(genomes_folder, threads, reference))

# Create the snp list file

os.system('ls %s/*.snp > %s.snp.paths' %(genomes_folder, output))
os.system('cp %s.snp.paths %s.snp.names' %(output, output))

names = open(f'{output}.snp.names','r')
out = open(f'{output}.snp.names2','w')

for i in names:
	i2 = os.path.basename(i)
	i3 = re.sub(".snp$", "", i2)
	out.write(i3)

out.close()

os.system('paste %s.snp.names2 %s.snp.paths > %s.snp.tab' %(output, output, output))

snp_tab_list = f'{output}.snp.tab' # tab of snp paths

if (snp_pos_list == "auto"):

	core_snps_list_path(snp_tab_list, half_window_size, output)
	# Inputs: [lista path file .snps] [lista coresnps] [path fasta reference] [output]

	snp_pos_list="%s.assembly_core_snp.pos_list" %(output.strip())

#### MASCHERO QUELLE NON BUONE

if (to_mask_pos != "none"):
	
	called_pos_file = open(snp_pos_list,'r')
	called_pos = called_pos_file.readlines()	

	positions_to_mask_file = open(to_mask_pos,'r')
	mask_pos = positions_to_mask_file.readlines()

	selected_pos = list(set(called_pos) - set(mask_pos))
	selected_pos.sort(key=int)
	
	selected_pos_file="%s.assembly_core_snp.selected.pos_list" %(output.strip())

	with open(selected_pos_file, 'w') as sorted_list:
    		for item in selected_pos:
        		sorted_list.write("%s" % item)

	snp_pos_list = selected_pos_file


#### PRODUCO IL FASTA SNP

core_snps_2_fasta(snp_tab_list, snp_pos_list, reference, output, print_ref)


### ORGANISE THE OUTPUTs

os.system('mkdir %s' %(output))
os.system('mkdir %s/tmp' %(output))

os.system('mv %s %s/' %(outname_snp_fasta, output))
os.system('mv %s %s/tmp/' %(snp_pos_list, output))
os.system('mv %s %s/tmp/' %(snp_tab_list, output))

os.system('mv %s/*.delta %s/tmp/' %(genomes_folder, output))
os.system('mv %s/*.snp %s/tmp/' %(genomes_folder, output))
os.system('mv %s.snp.names %s/tmp/' %(output, output))
os.system('mv %s.snp.names2 %s/tmp/' %(output, output))
os.system('mv %s.snp.paths %s/tmp/' %(output, output))


