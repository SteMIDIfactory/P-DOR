
$folder=shift;


%hash_pos=();
%hash_cont=();

open(POS,"<$folder/Results/$folder.mergedsnps.tab.allbases.poslist");
while(<POS>)
{
	chomp $_;
	@spos=split(" ",$_);
	$hash_pos{$spos[1]}=$spos[0];

}


open(CON,"<$folder/Results/$folder.mergedsnps.tab.allbases.contigposlist");
while(<CON>)
{
	chomp $_;
	@scon=split("\t",$_);
	$hash_con{$scon[2]}=$scon[0]."\t".$scon[1];

}


open(OUU,">$folder/Results/$folder.mergedsnps.tab.allbases.position");



foreach $kk(sort {$a <=> $b} keys %hash_pos)
{

		if ($hash_con{int($kk)} eq ""){print OUU "N\tN","\t",$kk,"\t",$hash_pos{$kk},"\n";}
		else{print OUU $hash_con{int($kk)},"\t",$kk,"\t",$hash_pos{$kk},"\n";}

}

