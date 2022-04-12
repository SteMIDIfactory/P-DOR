use Bio::SeqIO;

$inp=shift;

    $in  = Bio::SeqIO->new(-file => $inp , '-format' => 'Fasta');



$n_seq=0;
$str="";

while (my $seq=$in->next_seq)
{
	$n_seq++;
	$len=$seq->length;
	
	$str=$str.$seq->display_id."          ".$seq->seq."\n";	
}

print " ",$n_seq," ",$len,"\n",$str;








