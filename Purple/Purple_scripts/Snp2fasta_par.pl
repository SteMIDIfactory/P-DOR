

	use Bio::Perl;
	use Bio::SeqIO;
	use Bio::Seq;

$snp=shift; # snps table
$ref=shift; # concatenate

my $seqobj  = Bio::SeqIO->new( -format => 'fasta',-file => $ref) or die;
$seq=$seqobj->next_seq;
$conc=$seq->seq;
$conc =~ s/-//g;

@pos=split("",$conc);


open(SNP,"<$snp");

while(<SNP>)
{
  chomp $_;
  @s=split(" ",$_);
  $num=$#s;
  $hash{int($s[0])}="OK";

  $hash_out{$s[0]}= $_;

}

for $pp(0..$#pos)
{
   $pp2=$pp+1;

  if ($hash{$pp2} eq "")
  {
     $basi=$pos[$pp]." ";
     $basi2=$basi x $num;

   $hash_out{$pp2}=$pp2." ".$basi2;

  }

}


open(OUT,">$snp.allbases");

foreach $kk(sort {$a <=> $b} keys %hash_out)
{
  print OUT $hash_out{$kk},"\n";
}


close(OUT);


$line=`head -n 1 $snp.allbases`;
@s_line=split(" ",$line);

$out_file=$snp;
$out_file =~ s/\.mergedsnps\.tab/\.multifaconc\.fasta/;

open(OUT,">$out_file");

for $i(1..$#s_line)
{
	$i2=$i+1;

	$seq=`cut -d ' ' -f$i2 $snp.allbases`;

	@s_seq=split("\n",$seq);

	print OUT ">",$s_seq[0],"\n",join("",@s_seq[1..$#s_seq]),"\n";
}



