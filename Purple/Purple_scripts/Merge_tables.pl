

$tab1=shift;
$tab2=shift;

open(UNO,"<$tab1"); # patterns
open(DUE,"<$tab2"); # numbers


$c=0;
while(<UNO>)
{
	chomp $_;

	if ($c==0)
	{
		@org=split("\t",$_);
		$c++;
	}

	else
	{
		@s=split("\t",$_);

		$hash_row{$s[0]}="";

		for $j(1..$#s)
		{
			$hash_tab{$s[0]}{$org[$j]}=$s[$j];
		}
	}

}



$c2=0;
while(<DUE>)
{
	chomp $_;

	if ($c2==0)
	{
		@org2=split("\t",$_);
		$c2++;
	}

	else
	{
		@s2=split("\t",$_);

		$hash_row{$s2[0]}="";

		for $j(1..$#s2)
		{
			if ($s2[$j] == 1){$hash_tab{$s2[0]}{$org2[$j]}="c";}
			elsif ($s2[$j] == 0){$hash_tab{$s2[0]}{$org2[$j]}="a";}
			else {$hash_tab{$s2[0]}{$org2[$j]}="-";}
		}
	}

}



print "Features\t";

for $o(1..$#org)
{
	chomp $org[$o];
	print $org[$o],"\t";
}

print "\n";


$count=1;
foreach $k(sort {$a <=> $b} keys %hash_row)
{
	chomp $k;
	print $count."_".$k,"\t";
	$count++;

	foreach $o2(1..$#org)
	{
		print $hash_tab{$k}{$org[$o2]},"\t";		
	}

	print "\n";
}


















