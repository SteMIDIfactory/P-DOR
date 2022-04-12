
$tab=shift;
$org=shift;

open(TTT,"<$tab");
open(OUT,">$tab.formatted");


$c=0;
while(<TTT>)
{
	chomp $_;
	@s=split("\t",$_);

if ($c==0)
{
	
	print OUT join("\t",@s[0..6]),"\t",join("\t",@s[7+$org..$org+14]),"\t","TAG\tGene_start\tGene_end\tGene_Annotation","\t",join("\t",@s[7..$org+6]),"\t",join("\t",@s[$org+15..$#s-9]),"\n";

$c++;

}

else
{	

	if ($s[$org+7] eq "INTERGENIC")
	{print OUT join("\t",@s[0..6]),"\t",$s[$org+7],"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t",join("\t",@s[7..$org+6]),"\t",@s[$org+14..$#s],"\n";}

	else
	{
	print OUT join("\t",@s[0..6]),"\t",join("\t",@s[7+$org..$org+14]),"\t",$s[$#s-6],"\t",$s[$#s-5],"\t",$s[$#s-4],"\t",$s[$#s],"\t",join("\t",@s[7..$org+6]),"\t",join("\t",@s[$org+15..$#s-9]),"\n";
	}

}

	
}

