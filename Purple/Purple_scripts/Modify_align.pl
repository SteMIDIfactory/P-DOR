

$fasta=shift;

	use Bio::Perl;
	use Bio::SeqIO;
	use Bio::Seq;

my $seqobj  = Bio::SeqIO->new( -format => 'fasta',-file => $fasta) or die;


$c=0;
while($seq=$seqobj->next_seq)
{
	$array_seq[$c]=$seq->seq;
	$array_name[$c]=$seq->display_id;	
	$c++;
}


$s=$array_seq[0];


	$error=0;

	while ($s =~ /([N]+[-]+[N]+)/)
	{
		$s =~ /([N]+[-]+[N]+)/;
		$s2=$1;
		$s3=$s2;
		$s2=~s/\-/2/g;
		$s =~ s/$s3/$s2/;

		$error++;
	}


if ($error > 0)
{

open(OUT,">$fasta.sis.Mauvealn.fasta");

@s_s=split("",$s);
@s_array_seq_ref=split("",$array_seq[0]);
@s_array_seq_aln=split("",$array_seq[1]);


# print reference

print OUT ">",$array_name[0],"\n";

for $i(0..$#s_s)
{

	if ($s_s[$i] ne "2")
	{print OUT $s_array_seq_ref[$i];}
}

print OUT "\n";


# print aln

print OUT ">",$array_name[1],"\n";

for $i(0..$#s_s)
{
	if ($s_s[$i] ne "2")
	{print OUT $s_array_seq_aln[$i];}
}

print OUT "\n";



`mv $fasta $fasta.alnproblem`;

}





