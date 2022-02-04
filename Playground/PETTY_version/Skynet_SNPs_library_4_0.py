#!/usr/bin/python
# -*- coding: utf-8 -*- 

########## CREARE UN FILE CORE_SNP DA DUE FILES DI INPUT

import re, sys,os

# crea tre diz (diz_ref, diz_1 e diz_2). diz_ref sono le pos_ref (sono le keys dei dizionari) e le basi (sono i value assegnati alle keys). gli altri due diz contengono le pos dei due files di input in cui non sono presenti "." e le rispettive basi.  

def core_snps_two_files(file1, file2, window):
	count=0
	diz_ref={}
	files=[file1,file2]
	for arg in files:
		diz_base={}
		f=open(arg,"r")
		lines=f.readlines()
		for x in lines:
			if not re.search(r'\.',x):
				diz_base.update({int(x.split(' ')[0]):x.split(' ')[2]}) 
				diz_ref.update({int(x.split(' ')[0]):x.split(' ')[1]})   
				f.close()
		if count==0:
			diz_1={}
			diz_1=diz_base 
			count=count+1
		else:
			diz_2={}
			diz_2=diz_base

	# crea due liste, keys_1 e keys_2 e la lista intersect.

	keys_1=diz_1.keys() 
	keys_2=diz_2.keys()
	intersec= set(keys_1).intersection(keys_2)

	# crea una lista snv contenente le chiavi di intersect che presentano basi diverse nei dizionari diz_1 e diz_2 nella posizione esaminata.

	snv=[]
	for i in intersec:
		if diz_1[i] != diz_2[i]:
			snv.append(i)

	# crea lista only_k1. Aggiunge a diz_2 queste keys assegnandovi come value la bas_ref in quella posizione.

	only_k1= list(set(keys_1)-set(keys_2))  
	for i in only_k1:   
		diz_2.update({i:diz_ref[i]})

	# crea lista only_k2. Aggiunge a diz_1 queste keys assegnandovi come value la bas_ref in quella posizione.

	only_k2= list(set(keys_2)-set(keys_1))  
	for i in only_k2:
		diz_1.update({i:diz_ref[i]})

	# crea una lista all_pos contenenti tutte le posizioni.

	all_pos=[]  
	all_pos= snv+only_k1+only_k2

	# crea lista no_dup_pos contenente all_pos senza pos doppie e sortata.

	no_dup_pos=sorted(list(set(all_pos)))   

	# crea una lista core contenenti i coreSNP. Questa lista contiene le posizioni di no_dup_pos che distano almeno il valore della finestra (window). Viene verificato che ogni base sia 'atcgATCG'.

	letters = set('atcgATCG')
	core=[]

	if (no_dup_pos[1]-no_dup_pos[0])<int(window) and letters & set(diz_2[no_dup_pos[0]]) and letters & set(diz_1[no_dup_pos[0]]) and diz_1[no_dup_pos[0]] != diz_2[no_dup_pos[0]] :
			core.append(no_dup_pos[0])

	for i in range(1,len(no_dup_pos)-1):
		if (no_dup_pos[i]-no_dup_pos[i-1])<int(window) and (no_dup_pos[i+1]-no_dup_pos[i])<int(window) and letters & set(diz_1[no_dup_pos[i]]) and letters & set(diz_2[no_dup_pos[i]]) and diz_1[no_dup_pos[i]] != diz_2[no_dup_pos[i]]:
			core.append(no_dup_pos[i])

	if (no_dup_pos[-1]-no_dup_pos[-2])<int(window) and letters & set(diz_1[no_dup_pos[-1]]) and letters & set(diz_2[no_dup_pos[-1]]) and diz_1[no_dup_pos[-1]] != diz_1[no_dup_pos[-1]]:
			core.append(no_dup_pos[i])

	# stampa a schermo una tabella formata dalle colonne: primo file di input, secondo file di input, finestra e lunghezza della lista core; con uno dei due comandi sotto riportati.

	#outname="%s.%s.core.list" %(file1.replase(".XMFA.Mauvealn.fasta.snp", ""), file2.replase(".XMFA.Mauvealn.fasta.snp", ""))
	

	file1_nopath=os.path.basename(file1)
	file2_nopath=os.path.basename(file2)

	file1_rename= file1_nopath.replace(".XMFA.Mauvealn.fasta.snp", "")
	file2_rename= file2_nopath.replace(".XMFA.Mauvealn.fasta.snp", "")
    
	outname= "%s.%s.core.list" %(file1_rename, file2_rename)

	output=open(outname,"w")
	output.write("%s \t %s \t %d \t %d\n" % (file1_rename, file2_rename, (int(window)*2)+1, len(core)))


########## CREARE UN FILE CORE_SNP DA UNA LISTA DI FILE .SNP

import re, sys, io
from Bio import SeqIO

def core_snps_list_path(list_path, window, output):

	# caricare un input (lista file .snp)

	list_p=open(list_path,'r')

	# crea il dizionario vuoto e un diz della reference

	diz_pos={}
	diz_pos_ref={}

	# scorre riga per riga

	for path in list_p:
		p=open(path.strip(),'r')
	
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
					diz_pos_ref.update({int(pos):base_ref})
   
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

	outname_snp_list="%s.core.list" %(output.strip())
	#print(core)
	with open(outname_snp_list, "w") as file:
		file.write("\n".join(map(lambda x: str(x), core)))




########## CREARE IL FASTA DALLA LISTA DEI CORE_SNP

def core_snps_2_fasta(list_path,core_path,reference_fasta,output):

# Carica posizioni core

	core_file=open(core_path,'r')
	core=[int(i.strip()) for i in core_file.read().strip().split("\n")]

# Carica le basi della reference

	count=0
	for seq_record in SeqIO.parse(reference_fasta,"fasta"):
		list_pos_ref=(list(seq_record.seq.lower()))
		reference=seq_record.id
		count=count+1   

		assert count==1, "ERROR! Alignment reference contains more than one sequence"

# Prerara files output

	outname_snp_fasta="%s.core.fasta" %(output.strip())
	outfasta=open(outname_snp_fasta,"w")


# Print fasta reference

	s=''
	for i in core:
		list_pos_ref[int(i)-1].strip()
		s+=list_pos_ref[int(i)-1]

	outfasta.write(">REF_%s\n" %(reference.strip().replace(".snp","").replace(".fasta","").replace(".fna","")))
	outfasta.write("%s\n" %(s.strip().upper()))

# Per ogni file, carica basi organismo e stampa il fasta

	list_path2=open(list_path,'r')

	for path2 in list_path2:

		basi={}
		p2=open(path2.strip(),'r')

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

		org_name_fasta=path2.split('/')[-1]
		org_name_fasta2=org_name_fasta.replace(".XMFA.Mauvealn.fasta.snp", "")

		outfasta.write(">%s\n" %(org_name_fasta2.strip().replace(".snp","").replace(".fasta","").replace(".fna","")))
		outfasta.write("%s\n" %(s.strip().upper()))

		p2.close()

	outfasta.close()


########## CREA IL FILE CON I .SNP USANDO MAUVE

import os, sys

def snp_calling_mauve(reference, genomes):

	os.system("rm " + genomes + ".sslist")
	os.system("rm " + genomes + ".xmfa")
	os.system("rm " + genomes + ".xmfa.backbone")
	os.system("rm " + genomes + ".xmfa.bbcols")
	os.system("rm " + reference + ".sslist")

	os.system("./mauve_snapshot_2015-02-13/linux-x64/progressiveMauve --output=" + genomes + ".xmfa " + reference + " " + genomes)
	os.system("java -cp mauve_snapshot_2015-02-13/Mauve.jar org.gel.mauve.analysis.SnpExporter -f " + genomes + ".xmfa -o " + genomes + ".mauve_snps")

	file_mauve_snp=open(genomes+".mauve_snps", "r")
	snp_output=open(genomes+".snps","w")

	os.system("rm " + genomes + ".sslist")
	os.system("rm " + genomes + ".xmfa")
	os.system("rm " + genomes + ".xmfa.backbone")
	os.system("rm " + genomes + ".xmfa.bbcols")

	os.system("rm " + reference + ".sslist")

	for i in file_mauve_snp.readlines()[1:]:
		snp=i.split('\t')[0]
		snp_ref=list(snp)[0].lower()
		snp_org=list(snp)[1].lower()
		pos=i.split('\t')[2]

		snp_output.write("%s %s %s\n" % (pos, snp_ref, snp_org))

	file_mauve_snp.close()
	snp_output.close()

	os.system("rm " + genomes + ".mauve_snps")


########## CARD

import sys, os

def card_resistance_genes(genome):

	cmd='python rgi.py -i {} -n 1 -o {}'.format(genome, genome+'.card')
	os.system(cmd)

	cmd='python convertJsonToTSV.py -i {} -o {}'.format(genome+'.card.json',genome+'.genetab')
	os.system(cmd)

	cmd="python clean.py"
	os.system(cmd)

	cmd='rm {}.card.json'.format(genome)
	os.system(cmd)



