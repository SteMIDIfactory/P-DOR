
$folder=shift;
$refe=shift;

@ls=`ls $folder/*`;


foreach $file(@ls)
{
	chomp $file;

	@fff=split("/",$file);

	$fff[$#fff] =~ s/\.XMFA\.Mauvealn\.fasta\.sis\.Mauvealn\.fasta\.snp//g;
	$fff[$#fff] =~ s/\.XMFA\.Mauvealn\.fasta\.snp//g;

	push(@files,$fff[$#fff]);

	open(SNP,"<$file");

	while(<SNP>)
	{
		chomp $_;
		@s=split(" ",$_);

		$hash{$fff[$#fff]}{$s[0]}=$s[2];
		$hash{"Ref"}{$s[0]}=$s[1];	

		$hash_pos{$s[0]}="";

	}

}


@k=sort {$a <=> $b} keys %hash_pos;

print "Pos Ref_$refe ";

foreach $file(@files)
{chomp $file; 
print $file," ";
}

print "\n";


foreach $kk(@k)
{

	print $kk," ",$hash{"Ref"}{$kk}," ";

	foreach $file(@files)
	{
		chomp $file;

		if ($hash{$file}{$kk} eq "")
		{print $hash{"Ref"}{$kk}," ";}
		else
		{print $hash{$file}{$kk}," ";}

	}
	
	print "\n";

}
