

# lancia lo scrit con l'annotation tab finale per ottenere i fasta degli SNPs sinonimi e non


$with_pattern=shift;
$n_org=shift;

open(WWW,"<$with_pattern");

open(SIN,">$with_pattern.Syn.fasta");
open(NOT,">$with_pattern.NotSyn.fasta");

$c=0;
while(<WWW>)
{
	chomp $_;	
	@s=split("\t",$_);

	if ($c==0)
	{	
		for $i(6..6+$n_org-1)
		{
			$hash_snps_syn{$i}="";
			$hash_snps_notsyn{$i}="";
			$hash_name{$i}=$s[$i];
		}
		$c++,
	}


	else
	{
		$str_snps=join("",@s[6..6+$n_org-1]);
		
		if ($str_snps !~ /-/ and $str_snps !~ /[qweryuiopsdfhjklzxvbnm]/i and $s[0] =~ /SNV_/)
		{


			if ($s[$n_org+10] eq "Syn")
			{

				for $i(6..6+$n_org-1)
				{
					$hash_snps_syn{$i}=$hash_snps_syn{$i}.$s[$i];
				}

			}


			elsif ($s[$n_org+10] eq "NotSyn")
			{

				for $i(6..6+$n_org-1)
				{
					$hash_snps_notsyn{$i}=$hash_snps_notsyn{$i}.$s[$i];
				}

			}


			else
			{}



		}
	}

}


for $i(6..6+$n_org-1)
{

	print SIN ">",$hash_name{$i},"\n",$hash_snps_syn{$i},"\n";
	print NOT ">",$hash_name{$i},"\n",$hash_snps_notsyn{$i},"\n";

}


