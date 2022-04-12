
# This script add the frame to each SNP


$withpattern=shift;  #  .withpattern table
$n_org=shift;     # number of organisms

open(WIT,"<$withpattern");
open(OUT,">$withpattern.withframes");


$c=0;
while(<WIT>)
{

	chomp $_;
	@s=split("\t",$_);

	if ($c==0)
	{
		print OUT $_,"\tLOC\tPOS\tFRAME\tSTRAND\tCONS\tGFF\n";
		$c++;
	}


	else
	{

	if ($s[$n_org+6] eq "GENIC")
	{		


		$pos="bho";
		$frame="error";
		
		if ($s[$n_org+14] eq "+")
		{

			$pos=$s[3]-$s[$n_org+11]+1;

			$frame_pos = sprintf("%.1f", ($pos/3)-int(($pos/3)));
			chomp $frame_pos;

			if ($frame_pos == 0.0){$frame=3;}
			elsif ($frame_pos == 0.7){$frame=2;}
			elsif ($frame_pos == 0.3){$frame=1;}

		}

		else
		{
			$pos=abs($s[3]-$s[$n_org+12])+1;

			$frame_pos = sprintf("%.1f", ($pos/3)-int(($pos/3)));
			chomp $frame_pos;

			if ($frame_pos == 0.0){$frame=3;}
			elsif ($frame_pos == 0.7){$frame=2;}
			elsif ($frame_pos == 0.3){$frame=1;}

		}

		$strand=$s[$n_org+14];
		$_=~ s/\tGENIC\t/\tGENIC\t$pos\t$frame\t$strand\t/;

		print OUT $_,"\n";



	} # genic


	else
	{print OUT $_,"\n";}




	} # $c




}
