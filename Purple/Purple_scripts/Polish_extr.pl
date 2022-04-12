

$all=shift;
$snp=shift;
$num=shift;

$nn= "N" x $num;

open(AAA,"<$all");
open(TTT,">$all.tmp");

$cc2=0;
$cc=0;
while(<AAA>)
{

if($cc==0)
{$cc++}

else
{

	$cc2++;
	@s_pos=split(" ",$_);
	$hash_pos{$s_pos[0]}=$cc2;
	
	if($_ =~ /N/){print TTT "N";}
	elsif ($_ =~ /-/){print TTT "0";}
	else {print TTT "1";}

}

}

close(TTT);

open(BBB,"<$all.tmp") or die;
@zu=<BBB>;
$b=$zu[0];

$b=$nn.$zu[0];

while ($b=~ /(0+$nn)/)
{
	$b=~ /(0+$nn)/;
	$c=$1;
	$d=$c;
	$c=~s/0/2/g;

	$b =~ s/$d/$c/;
}


$br=reverse $b;

while ($br=~ /(0+$nn)/)
{
	$br=~ /(0+$nn)/;
	$cr=$1;
	$dr=$cr;
	$cr=~s/0/2/g;

	$br =~ s/$dr/$cr/;
}


$b2=reverse $br;
$b2=~ s/$nn//;
@b3=split("",$b2);

open(SSS,"<$snp");
@snp_array=<SSS>;

open(OUT,">$snp.selected");
print OUT $snp_array[0];

for $sa(1..$#snp_array)
{
	@snp_line=split(" ",$snp_array[$sa]);

	if (($b3[$hash_pos{$snp_line[0]}-1] eq "0") or ($b3[$hash_pos{$snp_line[0]}-1] eq "1"))
	{print OUT $snp_array[$sa];}

}

`rm $all.tmp`;






