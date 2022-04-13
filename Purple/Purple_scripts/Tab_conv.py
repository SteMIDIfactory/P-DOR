import sys
from Bio import SeqIO

inF=open(sys.argv[1].strip(),"r")
nnnum=int(sys.argv[2].strip())



tot=0
for x in SeqIO.parse(inF,"fasta"):
	for i in range(1,len(str(x.seq))+1):
		tot+=1
		string="%s\t%i\t%i" %(str(x.id),i,tot)
		print(string)
	tot+=nnnum


inF.close()
