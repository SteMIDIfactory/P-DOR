

$in=shift;


open(INN,"<$in");
open(OUT,">$in.sel");

$c=0;
while(<INN>)
{

	if ($c==0)
	{print OUT $_;$c++;}

	else
	{
		chomp $_;
		@s=split("\t",$_);

		$not=0;
		$one=0;		
	
		for $i(1..$#s)
		{
			if ($s[$i] == 1){$one++;}
		}


		if (($one < $#s))
		{print OUT $_,"\n";}
	}
}

close(OUT);
close(INN);

