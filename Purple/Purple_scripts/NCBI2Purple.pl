
# This script converts the .fna and .gff file available on NCBI database in Purple inputs

	use Bio::Perl;
	use Bio::SeqIO;
	use Bio::Seq;


$pa1=shift;
$sh1=shift;
$pa2=shift;
$sh2=shift;

$hash{$pa1}=$sh1;
$hash{$pa2}=$sh2;

my $ref = $hash{'-r'};
my $gff = $hash{'-g'};


if (($ref eq "") or ($gff eq ""))
{
	print "\nUsage: perl NCBI2Purple.pl -r [reference fasta] -g [reference gff]\n\n";exit;
}



	$in  = Bio::SeqIO->new(-file => $ref , '-format' => 'Fasta');

	open(FAS,">$ref.2Purple.fasta");
	open(GFA,">$gff.2Purple.gff");
	
	while($seq=$in->next_seq)
	{
		$name=$seq->display_id;
		@s=split('\|',$name);
		chomp $s[3];
		print FAS ">",$s[3],"\n",$seq->seq,"\n";
	}

open(GGG,"<$gff");

while(<GGG>)
{

	if ($_ !~ /^#/)
	{	
		chomp $_;
		@ss=split("\t",$_);

		if ($ss[2] ne "region"){ $hash_line{$ss[3]."_".$ss[4]}=$_; }
	}
}


	foreach $line(keys %hash_line)
	{
		print GFA $hash_line{$line},"\n";
	}







