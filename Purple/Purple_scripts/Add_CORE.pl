
$tab=shift;   # annotation tab
$org=shift;   # number of organisms


	use Bio::Perl;
	use Bio::SeqIO;
	use Bio::Seq;


open(OUT,">$tab.withCore.tab");

open(TTT,"<$tab");
$c=0;
while(<TTT>)
{
	chomp $_;
	@s=split("\t",$_);

		# Table header

	if ($c==0)
	{
		$c++;
		print OUT $s[0],"\tCORE\t",join("\t",@s[1..$#s]),"\n";
	}	


	else
	{

		if ($s[0] =~ /SNV/)
		{	
	          $sb=join("",@s[6..6+$org-1]);
	          if ($sb =~ /[qweryuiopsdfhjklzxvbnm-]/i){print OUT $s[0],"\tNotCORE\t",join("\t",@s[1..$#s]),"\n";}
	          else                                    {print OUT $s[0],"\tCORE\t",join("\t",@s[1..$#s]),"\n";}
		}

		else
		{ print OUT $s[0],"\tNA\t",join("\t",@s[1..$#s]),"\n";}


	}

} # while
