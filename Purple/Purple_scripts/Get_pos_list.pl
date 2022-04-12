

$all=shift;
open(AAA,"<$all");
@a=<AAA>;

open(OUT,">$all.poslist");

for $i(1..$#a)
{
	@s=split(" ",$a[$i]);
	print OUT $i," ",$s[0],"\n";
}
