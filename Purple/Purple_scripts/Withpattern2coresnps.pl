

$with_pattern=shift;
$n_org=shift;

open(WWW,"<$with_pattern");

$c=0;

while(<WWW>)
{
	chomp $_;	
	@s=split("\t",$_);

	if ($c==0)
	{	
		for $i(6..6+$n_org-1)
		{
			$hash_snps{$i}="";
			$hash_name{$i}=$s[$i];
		}
		$c++,
	}


	else
	{
		$str_snps=join("",@s[6..6+$n_org-1]);
		#$sf=split("_",$s[0]);

		if (($str_snps !~ /-/) and ($str_snps !~ /[qweryuiopsdfhjklzxvbnm]/i) and ($s[0] =~ /SNV_(.+?)_1/) and ($s[$#s-13] eq "GENIC"))
		{
			for $i(6..6+$n_org-1)
			{
				$hash_snps{$i}=$hash_snps{$i}.$s[$i];
			}
		}

		elsif (($str_snps !~ /-/) and ($str_snps !~ /[qweryuiopsdfhjklzxvbnm]/i) and ($s[$#s] eq "INTERGENIC") and ($s[0] =~ /SNV_/))
		{
			for $i(6..6+$n_org-1)
			{
				$hash_snps{$i}=$hash_snps{$i}.$s[$i];
			}
		}
		
	}

}


for $i(6..6+$n_org-1)
{

	print ">",$hash_name{$i},"\n",$hash_snps{$i},"\n";



}


