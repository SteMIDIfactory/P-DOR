

$with_pattern=shift;
$n_org=shift;

open(WWW,"<$with_pattern");

$c=0;

while(<WWW>)
{
	chomp $_;	
	@s=split(" ",$_);

	if ($c==0)
	{	
		for $i(1..$n_org)
		{
			$hash_snps{$i}="";
			$hash_name{$i}=$s[$i];
		}
		$c++,
	}


	else
	{
		$str_snps=join("",@s[1..$n_org]);
		
		if ($str_snps !~ /-/ and $str_snps !~ /[qweryuiopsdfhjklzxvbnm]/i and $s[$#s] =~ /SNV_/)
		{

			for $i(1..$n_org)
			{
				$hash_snps{$i}=$hash_snps{$i}.$s[$i];
			}

		}
	}

}


for $i(1..$n_org)
{

	print ">",$hash_name{$i},"\n",$hash_snps{$i},"\n";



}




