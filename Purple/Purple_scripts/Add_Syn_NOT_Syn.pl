

$fasta=shift; # multifaconc
$tab=shift;   # annotation tab
$org=shift;   # number of organisms


	use Bio::Perl;
	use Bio::SeqIO;
	use Bio::Seq;

my $seqobj  = Bio::SeqIO->new( -format => 'fasta',-file => $fasta) or die;
while($seq=$seqobj->next_seq)
{push(@seq_array,$seq->seq);}

open(OUT,">$tab.with_syn.tab");

open(TTT,"<$tab");
$c=0;
while(<TTT>)
{
	chomp $_;
	@s=split("\t",$_);

		# Table header

	if ($c==0)
	{
		$c++;
		print OUT join("\t",@s[0..$org+9]),"\tSyn_NOTsyn\tSNPs_count\tShannon_index\t",join("\t",@s[$org+10..$#s]),"\n";
	}	


	else
	{

		# Intergenic positions
		if ($s[$org+6] eq "INTERGENIC")
		{print OUT join("\t",@s[0..$org+9]),"\t",join("\t",@s[$org+10..$#s]),"\n";}

		# conserved genes have only synonimous mutations
	#	elsif ($s[$org+10] eq "CONSERVED")
	#	{print OUT join("\t",@s[0..$org+9]),"\tSyn\t",join("\t",@s[$org+10..$#s]),"\n";}
		
		# not conserved genes
		else
		{
			$pos=$s[5]-1; # substr start from zero
			$frame=$s[$org+8];
			$strand=$s[$org+9];

			$gap=0;
			%hash_aa=();
			$shannon=0;

				if ($frame == 2)
				{

						
						foreach $d(@seq_array)
						{
							$aa="";
							$sb=substr($d,$pos-1,3);
							if ($sb =~ /[qweryuiopsdfhjklzxvbnm-]/i){$gap=1;}
							else{

							$codon = Bio::Seq->new(-seq => $sb,-alphabet => 'dna' );

							if ($strand eq "+")
							{$aa=$codon->translate->seq;}
							else
							{$aa=$codon->revcom->translate->seq;}

							$hash_aa{$aa}=$hash_aa{$aa}+1;

						            }
						}

							@k=reverse sort { $hash_aa{$a} <=> $hash_aa{$b} } keys(%hash_aa);
							
							if ($gap == 0)
							{
								if ($#k == 0)
								{print OUT join("\t",@s[0..$org+9]),"\tSyn\t",$k[0],"(",$hash_aa{$k[0]},")\tNA\t",join("\t",@s[$org+10..$#s]),"\n";}

								else
								{print OUT join("\t",@s[0..$org+9]),"\tNotSyn\t";

								 for $ii(0..$#k-1)
								 {print OUT $k[$ii],"(",$hash_aa{$k[$ii]},"),";
								  $shannon=$shannon+(($hash_aa{$k[$ii]}/$org)*(log($hash_aa{$k[$ii]}/$org)));
								  }

								  print OUT $k[$#k],"(",$hash_aa{$k[$#k]},")\t";

								  $shannon=$shannon+(($hash_aa{$k[$#k]}/$org)*(log($hash_aa{$k[$#k]}/$org)));

								  print OUT -$shannon,"\t",join("\t",@s[$org+10..$#s]),"\n";}


							}


							else
							{
								print OUT join("\t",@s[0..$org+9]),"\tNA\tNA\tNA\t",join("\t",@s[$org+10..$#s]),"\n";	
							}


				}# frame2			

######

				elsif ($frame == 1)
				{

						
						foreach $d(@seq_array)
						{
							$aa="";

							if ($strand eq "+"){$sb=substr($d,$pos,3);}
							else{$sb=substr($d,$pos-2,3);}

							if ($sb =~ /[qweryuiopsdfhjklzxvbnm-]/i){$gap=1;}
							else{

							$codon = Bio::Seq->new(-seq => $sb,-alphabet => 'dna' );

							if ($strand eq "+")
							{$aa=$codon->translate->seq;}
							else
							{$aa=$codon->revcom->translate->seq;}

							$hash_aa{$aa}=$hash_aa{$aa}+1;

						            }
						}

							@k=reverse sort { $hash_aa{$a} <=> $hash_aa{$b} } keys(%hash_aa);

							if ($gap == 0)
							{
								if ($#k == 0)
								{print OUT join("\t",@s[0..$org+9]),"\tSyn\t",$k[0],"(",$hash_aa{$k[0]},")\tNA\t",join("\t",@s[$org+10..$#s]),"\n";}

								else
								{print OUT join("\t",@s[0..$org+9]),"\tNotSyn\t";

								 for $ii(0..$#k-1)
								 {print OUT $k[$ii],"(",$hash_aa{$k[$ii]},"),";
								  $shannon=$shannon+(($hash_aa{$k[$ii]}/$org)*(log($hash_aa{$k[$ii]}/$org)));
								  }

								  print OUT $k[$#k],"(",$hash_aa{$k[$#k]},")\t";

								  $shannon=$shannon+(($hash_aa{$k[$#k]}/$org)*(log($hash_aa{$k[$#k]}/$org)));

								  print OUT -$shannon,"\t",join("\t",@s[$org+10..$#s]),"\n";}

							}


							else
							{
								print OUT join("\t",@s[0..$org+9]),"\tNA\tNA\tNA\t",join("\t",@s[$org+10..$#s]),"\n";	
							}


				}# frame 1			


######


				elsif ($frame == 3)
				{

						
						foreach $d(@seq_array)
						{
							$aa="";

							if ($strand eq "+"){$sb=substr($d,$pos-2,3);}
							else{$sb=substr($d,$pos,3);}

							if ($sb =~ /[qweryuiopsdfhjklzxvbnm-]/i){$gap=1;}
							else{

							$codon = Bio::Seq->new(-seq => $sb,-alphabet => 'dna' );

							if ($strand eq "+")
							{$aa=$codon->translate->seq;}
							else
							{$aa=$codon->revcom->translate->seq;}

							$hash_aa{$aa}=$hash_aa{$aa}+1;

						            }
						}

							@k=reverse sort { $hash_aa{$a} <=> $hash_aa{$b} } keys(%hash_aa);

							if ($gap == 0)
							{
								if ($#k == 0)
								{print OUT join("\t",@s[0..$org+9]),"\tSyn\t",$k[0],"(",$hash_aa{$k[0]},")\tNA\t",join("\t",@s[$org+10..$#s]),"\n";}

								else
								{print OUT join("\t",@s[0..$org+9]),"\tNotSyn\t";

								 for $ii(0..$#k-1)
								 {print OUT $k[$ii],"(",$hash_aa{$k[$ii]},"),";
								  $shannon=$shannon+(($hash_aa{$k[$ii]}/$org)*(log($hash_aa{$k[$ii]}/$org)));
								  }

								  print OUT $k[$#k],"(",$hash_aa{$k[$#k]},")\t";

								  $shannon=$shannon+(($hash_aa{$k[$#k]}/$org)*(log($hash_aa{$k[$#k]}/$org)));

								  print OUT -$shannon,"\t",join("\t",@s[$org+10..$#s]),"\n";}

							}


							else
							{
								print OUT join("\t",@s[0..$org+9]),"\tNA\tNA\tNA\t",join("\t",@s[$org+10..$#s]),"\n";	
							}


				}# frame 3


######


		} # else





	} # else










} # while



