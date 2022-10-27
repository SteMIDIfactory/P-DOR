#!/usr/bin/python
# -*- coding: utf-8 -*-

import re, sys, io, os
from Bio import SeqIO

def core_snps_list_path(list_path, window, output):

	list_p=open(list_path,'r')
	diz_pos={}
	diz_pos_ref={}
	for path in list_p:
		p=open(path.strip(),'r')
		for x in p:
			if not re.search(r'\.',x):
				pos=int(x.split('\t')[0])
				base_ref=x.split('\t')[1]
				base_org=x.split('\t')[2]
				letters = set('atcgATCG')

				try:
					diz_pos[pos]
				except:
					diz_pos.update({int(pos):'0'})
					diz_pos_ref.update({int(pos):base_ref})

			if not ((letters & set(base_org)) and (letters & set(base_ref))):
				diz_pos.update({int(pos):'1'})
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


	outname_snp_list="%s.core.list" %(output.strip())
	with open(outname_snp_list, "w") as file:
		file.write("\n".join(map(lambda x: str(x), core)))


def core_snps_2_fasta(list_path,core_path,reference_fasta,output):

	core_file=open(core_path,'r')
	core=[]
	for i in core_file:
		i2=i.replace("\n","")
		core.append(i2)

	count=0
	for seq_record in SeqIO.parse(reference_fasta,"fasta"):
		list_pos_ref=(list(seq_record.seq.lower()))
		reference=seq_record.id
		count=count+1

		if count>1:
			print('Error')
			break


	outname_snp_fasta="%s.core.fasta" %(output.strip())
	outfasta=open(outname_snp_fasta,"w")

	s=''
	for i in core:
		list_pos_ref[int(i)-1].strip()
		s+=list_pos_ref[int(i)-1]

	outfasta.write(">%s\n" %(reference.strip().replace(".snp","").replace(".fasta","").replace(".fna","")))
	outfasta.write("%s\n" %(s.strip().upper()))

	list_path2=open(list_path,'r')

	for path2 in list_path2:

		basi={}
		p2=open(path2.strip(),'r')

		for x2 in p2:
			if not re.search(r'\.',x2):
				pos2=int(x2.split('\t')[0])
				base_org2=x2.split('\t')[2]
				basi.update({int(pos2):base_org2.strip()})

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
