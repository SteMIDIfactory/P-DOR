
$m=shift;
$tree=shift;
$name=shift;

open(TTT,"<$tree");
@alb=<TTT>;
$albero=$alb[0];
chomp $albero;

open(MMM,"<$m");
@mat=<MMM>;

chomp $mat[0];
@t = split("\t",$mat[0]);

print "#NEXUS\n\nBEGIN TAXA;\n	TITLE Taxa;\n	DIMENSIONS NTAX=$#t;	TAXLABELS\n		",join(" ",@t[1..$#t]),";\n	BLOCKID WM1148988015a107;\n\nEND;\n\n\n";
print "BEGIN CHARACTERS;\n	TITLE  'Matrix in file \"matrix\"';\n	DIMENSIONS  NCHAR=",$#mat,";\n	FORMAT DATATYPE = STANDARD GAP = k MISSING = ? SYMBOLS = \"a t g c -\"\;\n	CHARSTATELABELS\n		"; 

for $i(1..$#mat)
{
	@s=split("\t",$mat[$i]);
	print $i," ",$s[0];

	if ($i < $#mat){print ", ";}
	else {print ";\n";}
}


print "	MATRIX\n";

for $i(0..$#mat)
{

	chomp $mat[$i];
	@s=split("\t",$mat[$i]);


	for $i2(0..$#s)
	{
		$tot[$i][$i2]=$s[$i2];
	}

}


for $o(1..$#s)
{

	print "	";
	
	for $o2(0..$#mat)
	{
		print lc($tot[$o2][$o])," ";
	}

	print "\n";

}


print "\n;\n		BLOCKID WM1148988018a0;\n\n\nEND;\nBEGIN TREES;\n	Title 'Trees from \"01-DataMatrix.nex\"';\n	LINK Taxa = Taxa;\n	TRANSLATE\n";


for $to(1..$#t)
{
	print "		$to $t[$to]";

	$albero =~ s/,$t[$to]\:/,$to\:/;
	$albero =~ s/\($t[$to]\:/\($to\:/;


	if ($to < $#t){print ",\n";}
	else {print ";\n";}

}


print "	TREE $name = $albero\n";











