


$fasta=shift;

open(OUT,">$fasta.snp");


	use Bio::Perl;
	use Bio::SeqIO;
	use Bio::Seq;

my $seqobj  = Bio::SeqIO->new( -format => 'fasta',-file => $fasta) or die;


$ref=$seqobj->next_seq;
$seq_ref=$ref->seq;

$org=$seqobj->next_seq;
$seq_org=$org->seq;

@s_ref=split("",$seq_ref);
@s_org=split("",$seq_org);


$pos=0;
$under=0;

for $b(0..$#s_ref)
{
	
	if ($s_ref[$b] eq "-")
	{$under=$under+0.00000001;

		if ($s_ref[$b] ne $s_org[$b])
		{print OUT $pos+$under," ",$s_ref[$b]," ",$s_org[$b],"\n";}

	}

	else
	{$pos++;

	$under=0;

		if ($s_ref[$b] ne $s_org[$b])
		{print OUT $pos," ",$s_ref[$b]," ",$s_org[$b],"\n";}

	}

}









