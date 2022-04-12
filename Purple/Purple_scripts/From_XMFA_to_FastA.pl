
	use Bio::Perl;
	use Bio::SeqIO;
	use Bio::Seq;


$xmfa=shift;
open(XXX,"<$xmfa");

open(TMP,">$xmfa.tmp");

$org=0;
while(<XXX>)
{
	if ($_ !~ /#/)
	{print TMP $_;}	
	else
	{
		if (($_ =~ /File\t/) and ($_ =~ /Sequence/))
		{@line_s=split("\t",$_);
		 $org++;
		 $org_name[$org]=$line_s[1];
		}

	}

}



close(TMP);

open(INN,"<$xmfa.tmp");
$tmp_str=join("",<INN>);
$tmp_str=~ s/=\n/=/g;
close(INN);

@tmp_s=split("=",$tmp_str);
mkdir("$xmfa\_Mauve_block_fasta");

foreach $block(@tmp_s)
{

  if ($block =~ /> 1:/)
    {
	@blo = $block =~ /> 1:(.+?) /;
	@block_s=split("-",$blo[0]);
	$ff{$block_s[0]}="";
	
	open(BLK,">$xmfa\_Mauve_block_fasta/$block_s[0].fasta");
	print BLK $block;
	close(BLK);
    }

}

mkdir("$xmfa\_tmp_conc");
for $u(1..$org)
{open(CON,">>$xmfa\_tmp_conc/$u.fasta");print CON ">$org_name[$u]\n";close(CON);}

$nn="N" x 0;

$co=0;

foreach $fas(sort {$a <=> $b} keys %ff)
{	

	`sed -i \"s/\\n>/>/\" $xmfa\_Mauve_block_fasta/$fas.fasta`;


	my $seqobj  = Bio::SeqIO->new( -format => 'fasta',-file => "$xmfa\_Mauve_block_fasta/$fas.fasta") or die;
	%hash=();

	while($seq = $seqobj -> next_seq)
	{
		@id_s=split("\:",$seq->display_id);
		$hash{$id_s[0]}=$seq->seq;

		if ($id_s[0] == 1){$l=$seq->length;}
	}

	
	$j="?" x $l;

if ($co == 0)
{

	for $i(1..$org)
	{

		if ($hash{$i} eq "")
		{open (CON,">>$xmfa\_tmp_conc/$i.fasta");print CON $j;close(CON);}

		else
		{open (CON,">>$xmfa\_tmp_conc/$i.fasta");print CON $hash{$i};close(CON);}

	}

}

else
{

	for $i(1..$org)
	{

		if ($hash{$i} eq "")
		{open (CON,">>$xmfa\_tmp_conc/$i.fasta");print CON $nn.$j;close(CON);}

		else
		{open (CON,">>$xmfa\_tmp_conc/$i.fasta");print CON $nn.$hash{$i};close(CON);}

	}




}

$co++;

}

	for $y(1..$org)
	{open (CON,">>$xmfa\_tmp_conc/$y.fasta");print CON "\n";close(CON);}



`cat $xmfa\_tmp_conc/1.fasta $xmfa\_tmp_conc/3.fasta > $xmfa.Mauvealn.fasta`;

`rm -r $xmfa\_tmp_conc`;
`rm -r $xmfa\_Mauve_block_fasta`;
`rm $xmfa.tmp`;





