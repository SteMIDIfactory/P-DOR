
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


for $i(1..$#mat)
{
	@s=split("\t",$mat[$i]);
}


for $i(0..$#mat)
{

	chomp $mat[$i];
	@s=split("\t",$mat[$i]);


	for $i2(0..$#s)
	{
		$tot[$i][$i2]=$s[$i2];
	}

}


for $to(1..$#t)
{
	$albero =~ s/,$t[$to]\:/,$to\:/;
	$albero =~ s/\($t[$to]\:/\($to\:/;
}


print "	TREE $name = $albero\n";











