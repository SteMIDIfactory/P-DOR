import sys

inF=open(sys.argv[1].strip(),"r")

oname="%s_sorted" %(sys.argv[1].strip())
outsnp=open(oname,"w")

oldnum=0
c=0

out=""
for line in inF.readlines():
	if c==0:
		
		outsnp.write(line)
		
	else:
		num=float(line.strip().split()[0])
		if num-oldnum<2:
			out+=line

		else:
			out+="\n"			
			out+=line

		oldnum=num
	c+=1

inF.close()


List=out.split("\n\n")

c=0

for lll in List:
	if len(lll.strip().split("\n"))>1:
		c+=1
		outsnp.write("%s Stretch_%i_start\n" %(lll.strip().split("\n")[0].strip(),c))
		pattern=lll.strip().split("\n")[0].strip().split()
		pattern.pop(0)
		pattern="".join(pattern)
		ppp=""

		for h in pattern:
			if h!="-":
				ppp+="#"
			else:
				ppp+="-"
		akg=0
		for i in lll.strip().split("\n"):
			pattern=i.strip().split()
			pattern.pop(0)
			pattern="".join(pattern)

			ppp2=""
			for h in pattern:
				if h!="-":
					ppp2+="#"
				else:
					ppp2+="-"

			if ppp2!=ppp:
				akg+=1
				outsnp.write("%s Stretch_%i_mid_%i\n" %(i.strip(),c,akg))
				ppp=ppp2
		
		outsnp.write("%s Stretch_%i_end\n" %(lll.strip().split("\n")[-1].strip(),c))
	else:
		for i in lll.strip().split("\n"):
			c+=1
			outsnp.write("%s SNV_%i\n" %(i,c))
outsnp.close()


