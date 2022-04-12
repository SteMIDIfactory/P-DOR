

$fasta_snp=shift;
$annot_tab=shift;

@s=split('\.',$fasta_snp);
$out=join('.',@s[0..$#s-1]);

`sed '/^#/ d' $out\_sta.gr | sed '/^&/ d' - > $out\_sta.gr\_mod`;

open(HOM,"<$out\_sta.gr\_mod");
open(ANN,"<$annot_tab");


while(<HOM>)
{
	chomp $_;
	@arr_hom=split(" ",$_);
	$hom_values[$arr_hom[0]]=$arr_hom[1];
}


$c=0;
$count=1;
while(<ANN>)
{
	chomp $_;

	if ($c==0)
	{
		@st=split("\t",$_);
		print $st[0],"\t",$st[1],"\tHomoplasy_value\t",join("\t",@st[2..$#st]),"\n";
		$c++;
	}

	else
	{
		@st2=split("\t",$_);
		
		if ($st2[1] eq "CORE")
		{
			print $st2[0],"\t",$st2[1],"\t",$hom_values[$count],"\t",join("\t",@st2[2..$#st2]),"\n";
			$count++;			
		}

		else
		{
			print $st2[0],"\t",$st2[1],"\tNA\t",join("\t",@st2[2..$#st2]),"\n";			
		}


	}

}
