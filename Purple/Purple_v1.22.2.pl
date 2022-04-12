
	use Bio::Perl;
	use Bio::SeqIO;
	use Bio::Seq;

$pa1=shift;
$sh1=shift;
$pa2=shift;
$sh2=shift;
$pa3=shift;
$sh3=shift;
$pa4=shift;
$sh4=shift;
$pa5=shift;
$sh5=shift;
$pa6=shift;
$sh6=shift;
$pa7=shift;
$sh7=shift;
$pa8=shift;
$sh8=shift;
$pa9=shift;
$sh9=shift;
$pa10=shift;
$sh10=shift;
$pa11=shift;
$sh11=shift;
$pa12=shift;
$sh12=shift;
$pa13=shift;
$sh13=shift;

$hash{$pa1}=$sh1;
$hash{$pa2}=$sh2;
$hash{$pa3}=$sh3;
$hash{$pa4}=$sh4;
$hash{$pa5}=$sh5;
$hash{$pa6}=$sh6;
$hash{$pa7}=$sh7;
$hash{$pa8}=$sh8;
$hash{$pa9}=$sh9;
$hash{$pa10}=$sh10;
$hash{$pa11}=$sh11;
$hash{$pa12}=$sh12;
$hash{$pa13}=$sh13;

my $folder = $hash{'-f'};
my $ref = $hash{'-r'};
my $gff = $hash{'-g'};
my $tree = $hash{'-t'};
my $make_fasta = $hash{'-fa'};
my $num_nn = $hash{'-n'};
my $cpu = $hash{'-cpu'};
my $steps = $hash{'-s'};
my $clus = $hash{'-c'};
my $indel = $hash{'-indel'};
my $homo = $hash{'-h'};
my $ncbi = $hash{'-ncbi'};

# Default

if ($tree eq ""){$tree = "n";}
if ($num_nn eq ""){$num_nn=1000;}
if($cpu eq ""){$cpu=1;}
if ($make_fasta eq ""){$make_fasta="n";}
if($steps eq ""){$steps="all";} # Steps: can be only "all", "align" or "call"
if($clus eq ""){$clus="n";}
if($indel eq ""){$indel = "n";}
if ($homo eq ""){$homo = "y";}
if($ncbi eq ""){$ncbi = "n";}

if ($ref eq "" or $folder eq "" or $gff eq ""){

print "
###############################################################################
			Welcome! Enjoy your SNP analysis!
###############################################################################

This is the quickest command line:
perl Purple.pl -r [reference genome] -g [annotation file] -f [genomes folder] 

Inputs:
-r	Reference genome in FastA format			(required)
-g	Reference annotation in GFF format			(required)
-f	Folder including the genomes to align			(required)


Other Parameters:
-t	Perform phylogenetic analyses (y/n) 			def: n
-fa	Extract gene alignments (time consuming) (y/n)  	def: n
-cpu 	Number of threats 					def: 1
-indel	Include indels in the multiple genomes alignment (y/n)	def: n
-s	Steps to perfom (all/align/call) 			def: all
-h 	Perform analysis of homoplasy (y/n)			def: y
-ncbi format the reference .fna and .gff files retrieved from NCBI  (y/n)   def: n 

Expert:
-n 	Number of N between contigs		 		def: 1000

";exit;}


# Inputs check

if (-e $ref){}else{print "\nI can't find the $ref file, please check the file name (or position) and try again.\n
If you need help information, please launch the script without any parameters.\n\n";exit;}
if (-d $folder){}else{print "\nI can't find the $folder folder, please check the folder name (or position) and try again.\n
If you need help information, please launch the script without any parameters.\n\n";exit;}
if (-d "$folder\_OUT"){print "\nI have found an existing folder with name \"$folder\_OUT\" (the name that I will use for the output folder),\nplease rename or delete it and try again.\n
If you need help information, please launch the script without any parameters.\n\n";exit;}
if (-e $gff){}else{print "\nI can't find the $gff file, please check the file name (or position) and try again.\n
If you need help information, please launch the script without any parameters.\n\n";exit;
}

# Parameters check

if (($tree ne "y") and ($tree ne "n")){print "\nThe \"$tree\" is not an allowed value for the -t flag. The only allowed values are y or n.\nPlease check your input line and try again.\n
If you need help information, please launch the script without any parameters.\n\n";exit;}

if (($make_fasta ne "y") and ($make_fasta ne "n")){print "\nThe \"$make_fasta\" is not an allowed value for the -fa flag. The only allowed values are y or n.\nPlease check your input line and try again.\n
If you need help information, please launch the script without any parameters.\n\n";exit;}

if (($clus ne "y") and ($clus ne "n")){print "\nThe \"$clus\" is not an allowed value for the -c flag. The only allowed values are y or n.\nPlease check your input line and try again.\n
If you need help information, please launch the script without any parameters.\n\n";exit;}

if (($indel ne "y") and ($indel ne "n")){print "\nThe \"$indel\" is not an allowed value for the -indel flag. The only allowed values are y or n.\nPlease check your input line and try again.\n
If you need help information, please launch the script without any parameters.\n\n";exit;}

if (($homo ne "y") and ($homo ne "n")){print "\nThe \"$homo\" is not an allowed value for the -h flag. The only allowed values are y or n.\nPlease check your input line and try again.\n
If you need help information, please launch the script without any parameters.\n\n";exit;}

if (($num_nn ne $num_nn+0) and (int($num_nn) eq $num_nn) and ($num_nn > 0)){print "\nThe \"$num_nn\" is not an allowed value for the -n flag. The only allowed values are not zero natural numbers.\nPlease check your input line and try again.\n
If you need help information, please launch the script without any parameters.\n\n";exit;}

if (($num_nn ne $num_nn+0) and (int($num_nn) eq $num_nn) and ($num_nn > 0)){print "\nThe \"$cpu\" is not an allowed value for the -cpu flag. The only allowed values are not zero natural numbers.\nPlease check your input line and try again.\n
If you need help information, please launch the script without any parameters.\n\n";exit;}

if (($steps ne "call") and ($steps ne "align") and ($steps ne "all")){print "\nThe \"$steps\" is not an allowed value for the -s flag. The only allowed values are all, align and call.\nPlease check your input line and try again.\n
If you need help information, please launch the script without any parameters.\n\n";exit;}

################### NCBI

if ($ncbi eq "y")
{
`perl Purple_scripts/NCBI2Purple.pl -g $gff -r $ref`;
$gff=$gff.".2Purple.gff";
$ref=$ref.".2Purple.fasta";
}

else
{}

####################

if ($steps eq "call")
{
@ls=`ls $folder/snp_tab/*`;
$num_org=$#ls+2;
}

if ($steps eq "call")
{

# Concatenate reference contigs

my $seqobj  = Bio::SeqIO->new( -format => 'fasta',-file => $ref) or die;

$conc="";
$nn="N" x $num_nn;

while($seq=$seqobj->next_seq)
{
	$conc=$conc.lc($seq->seq).$nn;
	$name=$seq->display_id;
}

open(OUT,">$ref.conc.fasta");
print OUT ">",$name,"\n",$conc;
close(OUT);

}

# how many organims?

if (($steps eq "all") or ($steps eq "align") or ($steps eq "call"))
{

@ls=`ls $folder/*`;
$num_org=$#ls+2;

##############################################     CHECK and FORMAT INPUTS

print "\n\n0. CHECK and FORMAT INPUTS ... ";

# Read contig names from .GFF

open(GFF,"<$gff");

while(<GFF>)
{
	if ($_ !~ /^#/)
	{
		@s=split("\t",$_);
		chomp $s[0];
		$hash_gff_contigs{$s[0]}="gff";		
	}	
}

close(GFF);

# Read contig names from FastA

my $seqobj_ref  = Bio::SeqIO->new( -format => 'fasta',-file => $ref) or die;

while($seq_ref=$seqobj_ref->next_seq)
{
	$hash_fasta_contigs{$seq_ref->display_id}="fasta";
}

# lists comparison

@k_gff=keys %hash_gff_contigs;

$error=0;
for $i(0..$#k_gff)
{
	if ($hash_fasta_contigs{$k_gff[$i]} eq "")
	{$error=1;}
}


if ($error == 1)
{
	print "ERROR: Contig names in GFF and in FastA file don't correspond!\n\n";
	exit;
}


}



if (($steps eq "all") or ($steps eq "align"))
{


foreach $file(@ls)
{
	chomp $file;

	`cp $file $file.tomin`;

	my $seqobj_org  = Bio::SeqIO->new( -format => 'fasta',-file => "$file.tomin") or die;

	open(OUT,">$file");
	
	while($seq_org=$seqobj_org->next_seq)
	{
		print OUT ">",$seq_org->display_id,"\n",lc($seq_org->seq),"\n";
	}

	close(OUT);

	`rm $file.tomin`;

}

# Concatenate reference contigs

my $seqobj  = Bio::SeqIO->new( -format => 'fasta',-file => $ref) or die;

$conc="";
$nn="N" x $num_nn;

while($seq=$seqobj->next_seq)
{

	$conc=$conc.lc($seq->seq).$nn;
	$name=$seq->display_id;
}

open(OUT,">$ref.conc.fasta");
print OUT ">",$name,"\n",$conc;
close(OUT);

	`cp $ref.conc.fasta $ref.conc.REF.fasta`;
	`sed -i "s/>/>REF_/g" $ref.conc.REF.fasta`;

$dat=`date`;
chomp $dat;
print "Done! ($dat)\n";


##############################################     ALIGN GENOMES WITH MAUVE

print "1. ALIGN GENOMES WITH MAUVE ... "; 

	`ls $folder/* | parallel -j $cpu "progressiveMauve --output={}.XMFA $ref.conc.fasta $ref.conc.REF.fasta {}"`;

	`ls $folder/*.XMFA | parallel -j $cpu "perl Purple_scripts/From_XMFA_to_FastA.pl {}"`;

	`ls $folder/*.Mauvealn.fasta | parallel -j $cpu "perl Purple_scripts/Modify_align.pl {}"`;

$dat=`date`;
chomp $dat;
print "Done! ($dat)\n";


##############################################     SNVs EXTRACTION

print "2. SNVs CALLING ... ";

  	`ls $folder/*.Mauvealn.fasta | parallel -j $cpu "perl Purple_scripts/Fasta2svps.pl {}"`;

$dat=`date`;
chomp $dat;
print "Done! ($dat)\n";

    if ($steps ne "align")
    {`mkdir $folder/snp_tab`;`mv $folder/*.snp $folder/snp_tab`;}

    else
    {
        `mkdir $folder/tmp`;
        `mv $folder/*.sslist $folder/tmp`;
        `mv $folder/*.XMFA $folder/tmp`;
        `mv $folder/*.backbone $folder/tmp`;
        `mv $folder/*.bbcols $folder/tmp`;
        `mv $folder/*.XMFA.Mauvealn.fasta $folder/tmp`;

        `mkdir $folder/snp_tab`;
        `mv $folder/*.snp $folder/snp_tab`;
    
        `mv $ref.conc.fasta $folder/tmp/`;
        `mv $ref.conc.fasta.sslist $folder/tmp/`;
        `mv $ref.conc.REF.fasta $folder/tmp/`;
        `mv $ref.conc.REF.fasta.sslist $folder/tmp/`;

        `mkdir $folder\_OUT`;
        `mv $folder/snp_tab $folder\_OUT/`;
        `mv $folder/tmp $folder\_OUT/`;

        exit;
    }


} # if steps 

else{}

##############################################     SNVs CALLING and ANALYSES

if (($steps eq "all") or ($steps eq "call"))
{

if ($indel eq "n")
{
	`mv $folder/snp_tab/ $folder/snp_tab_with_indel/`;
	`mkdir $folder/snp_tab`;

	`ls $folder/snp_tab_with_indel/* | parallel -j $cpu "perl Purple_scripts/Remove_indels.pl {} $folder/snp_tab"`;
}

print "3. MERGE SNVs ... ";

  	`perl Purple_scripts/Merge_snps.pl $folder/snp_tab $ref > $folder.mergedsnps.tab`;
  	`perl Purple_scripts/Snp2fasta_par.pl $folder.mergedsnps.tab $ref.conc.fasta`;

  	`mkdir $folder/reference`;
  	`cp $ref*  $folder/reference/`;
  	`cp $folder/reference/$ref .`; # modifica
  	`mkdir $folder/out_MAUVE`;

  	
if ($steps eq "all")
{
  	`mv $folder/*.sslist $folder/out_MAUVE`;
  	`mv $folder/*.XMFA* $folder/out_MAUVE`;
}
 
  	`mkdir $folder/Results`;
  	`mv *.mergedsnps.tab* $folder/Results`;
  	`mv *.multifaconc.fasta $folder/Results`;

  	`perl Purple_scripts/Polish_extr.pl $folder/Results/$folder.mergedsnps.tab.allbases $folder/Results/$folder.mergedsnps.tab $num_nn`;

	`perl Purple_scripts/Get_pos_list.pl $folder/Results/$folder.mergedsnps.tab.allbases`;

	`python Purple_scripts/Tab_conv.py $folder/reference/$ref $num_nn > $folder/Results/$folder.mergedsnps.tab.allbases.contigposlist`;

	`perl Purple_scripts/Merge_pos_tab.pl $folder`;

	`rm $folder/Results/$folder.mergedsnps.tab.allbases.contigposlist`;
	`rm $folder/Results/$folder.mergedsnps.tab.allbases.poslist`;

	`python Purple_scripts/Sort_snps.py $folder/Results/$folder.mergedsnps.tab.selected`;

	`cp $gff $folder/reference/`;

$dat=`date`;
print "Done! ($dat)\n";

print "4. SNVs ANNOTATION ... ";

	`perl Purple_scripts/Split_table.pl $folder/Results/$folder.mergedsnps.tab.selected_sorted $cpu`;

	`ls $folder/Results/$folder.mergedsnps.tab.selected_sorted.split_* | parallel -j $cpu "perl Purple_scripts/Annotation_SNPS.pl $folder/reference/$gff $folder {} $make_fasta $cpu"`;

	`perl Purple_scripts/Merge_splitted_results.pl $folder`;

	`mkdir $folder/tmp`;
	`mkdir $folder/tmp/tables`;

	`perl Purple_scripts/Select_genes.pl  $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.gene_matrix`;

	
	`perl Purple_scripts/Annotated2pattern.pl $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated $num_org`;

	`perl Purple_scripts/Merge_tables.pl $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.pattern_list $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.gene_matrix.sel > $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.gene_matrix_pattern_list`;

	`mv $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.gene_matrix $folder/tmp/tables/`;

	`mkdir $folder/out_MAUVE/reference`;

	`mv $folder/Results/$folder.mergedsnps.tab $folder/tmp/tables/`;
	`mv $folder/Results/$folder.mergedsnps.tab.allbases $folder/tmp/tables/`;
	`mv $folder/Results/$folder.mergedsnps.tab.allbases.position $folder/tmp/tables/`;
	`mv $folder/Results/$folder.mergedsnps.tab.selected $folder/tmp/tables/`;
	`mv $folder/Results/$folder.mergedsnps.tab.selected_sorted $folder/tmp/tables/`;

	`mv $folder/out_MAUVE/ $folder/tmp/`;
	`mv $folder/snp_tab/ $folder/tmp/`;

	`perl Purple_scripts/Add_codon_frame2.pl $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns $num_org`;
	`perl Purple_scripts/Add_Syn_NOT_Syn.pl $folder/Results/$folder.multifaconc.fasta $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes $num_org`;

	`perl Purple_scripts/Add_CORE.pl $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.with_syn.tab $num_org`; # komma
	`perl Purple_scripts/Format_table2.2.pl $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.with_syn.tab.withCore.tab $num_org`;

$dat=`date`;
chomp $dat;
print "Done! ($dat)\n";

print "5. GET SNPs MULTIFASTA ...";

# intergenic SNPs

$int_num=$num_org+7;

	`head -n 1 $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes > $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.intergenic`;

	`awk '\$$int_num == "INTERGENIC"' $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes >> $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.intergenic`;

# SNPS in position 1,2,3

$pos_num=$num_org+9;

	`head -n 1 $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes > $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.pos1`;

	`awk '\$$pos_num == 1' $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes >> $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.pos1`;


	`head -n 1 $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes > $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.pos2`;

	`awk '\$$pos_num == 2' $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes >> $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.pos2`;


	`head -n 1 $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes > $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.pos3`;

	`awk '\$$pos_num == 3' $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes >> $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.pos3`;


# neutral SNPs

	`head -n 1 $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes > $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.neutrali`;

	`awk '\$$int_num == "INTERGENIC"' $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes >> $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.neutrali`;

	`awk '\$$pos_num == 3' $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes >> $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.neutrali`;


# Make fastas

	`perl Purple_scripts/Withpattern2coresnps.pl $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.neutrali $num_org > $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.neutrali.fasta`;

	`perl Purple_scripts/Withpattern2coresnps.pl $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.intergenic $num_org > $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.intergenic.fasta`;

	`perl Purple_scripts/Withpattern2coresnps.pl $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.pos1 $num_org > $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.pos1.fasta`;

	`perl Purple_scripts/Withpattern2coresnps.pl $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.pos2 $num_org > $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.pos2.fasta`;

	`perl Purple_scripts/Withpattern2coresnps.pl $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.pos3 $num_org > $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.pos3.fasta`;

	`perl Purple_scripts/Withpattern2coresnps.pl $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes $num_org > $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.core.fasta`;

	`perl Purple_scripts/WithSyn2FastaSnp.pl $folder/Results/*.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.with_syn.tab $num_org`;


$dat=`date`;
chomp $dat;
print "Done! ($dat)\n";


if ($clus eq "y")
{

	open(RRR,">$folder/Results/$folder.core.R");
	print RRR "library(seqinr)\n";
	print RRR "a=read.alignment(\"$folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.core.fasta\",format=\"fasta\")\n";
	print RRR "d=dist.alignment(a,matrix = \"identity\")\n";
	print RRR "pv=hclust(as.dist(d))\n";
	print RRR "pdf(\"$folder/Results/$folder.core_SNPs.hclust.pdf\")\n";
	print RRR "plot(pv)\n";
	print RRR "dev.off()\n";

	close(RRR);

	system("Rscript $folder/Results/$folder.core.R"); 

	`rm $folder/Results/$folder.core.R`;	
}

else
{}

if ($tree eq "y")
{

print "6. PHYLOGENETIC ANALYSES ...";

# fasttree

`fasttree -nt -gtr -verbose 0 $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.intergenic.fasta > $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.intergenic.tree`;
`fasttree -nt -gtr -verbose 0 $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.neutrali.fasta > $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.neutrali.tree`;
`fasttree -nt -gtr -verbose 0 $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.pos1.fasta > $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.pos1.tree`;
`fasttree -nt -gtr -verbose 0 $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.pos2.fasta > $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.pos2.tree`;
`fasttree -nt -gtr -verbose 0 $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.pos3.fasta > $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.pos3.tree`;
`fasttree -nt -gtr -verbose 0 $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.core.fasta > $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.core.tree`;
`fasttree -nt -gtr -verbose 0 $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.with_syn.tab.Syn.fasta > $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.with_syn.tab.Syn.tree`;
`fasttree -nt -gtr -verbose 0 $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.with_syn.tab.NotSyn.fasta > $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.with_syn.tab.NotSyn.tree`;


	`perl Purple_scripts/Patternlist2nex.pl $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.gene_matrix_pattern_list $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.core.tree CORE> $folder/Results/$folder.nex`;

	`perl Purple_scripts/Patternlist2nex_append.pl $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.gene_matrix_pattern_list $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.intergenic.tree INTERGENIC>> $folder/Results/$folder.nex`;

	`perl Purple_scripts/Patternlist2nex_append.pl $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.gene_matrix_pattern_list $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.neutrali.tree NEUTRAL>> $folder/Results/$folder.nex`;

	`perl Purple_scripts/Patternlist2nex_append.pl $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.gene_matrix_pattern_list $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.pos1.tree POS1>> $folder/Results/$folder.nex`;

	`perl Purple_scripts/Patternlist2nex_append.pl $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.gene_matrix_pattern_list $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.pos2.tree POS2>> $folder/Results/$folder.nex`;

	`perl Purple_scripts/Patternlist2nex_append.pl $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.gene_matrix_pattern_list $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.pos3.tree POS3>> $folder/Results/$folder.nex`;

	`perl Purple_scripts/Patternlist2nex_append.pl $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.gene_matrix_pattern_list $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.with_syn.tab.Syn.tree Syn>> $folder/Results/$folder.nex`;

	`perl Purple_scripts/Patternlist2nex_append.pl $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.gene_matrix_pattern_list $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.with_syn.tab.NotSyn.tree NotSyn>> $folder/Results/$folder.nex`;

$dat=`date`;
chomp $dat;
print "Done! ($dat)\n";

}

`mv $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated $folder/tmp/tables/`;
`mv $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns $folder/tmp/tables/`;
`mv $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.intergenic $folder/tmp/tables/`;
`mv $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.neutrali $folder/tmp/tables/`;
`mv $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.pos1 $folder/tmp/tables/`;
`mv $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.pos2 $folder/tmp/tables/`;
`mv $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.pos3 $folder/tmp/tables/`;

`mv $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.with_syn.tab $folder/tmp/tables/`;
`mv $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.with_syn.tab.withCore.tab $folder/tmp/tables/`;

# rename outputs

`mv $folder/Results/*.mergedsnps.tab.selected_sorted.annotated.gene_matrix_pattern_list $folder/Results/$folder.Allpatterns`;
`mv $folder/Results/*.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes $folder/tmp/tables`;
`mv $folder/Results/*.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.with_syn.tab.withCore.tab.formatted $folder/Results/$folder.SNPs.annotation_tab`;

`mkdir $folder/Results/Fasta_SNPs`;

`mv $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.core.fasta $folder/Results/Fasta_SNPs/$folder.core_SNPs.fasta`;
`mv $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.intergenic.fasta $folder/Results/Fasta_SNPs/$folder.intergenic_SNPs.fasta`;
`mv $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.neutrali.fasta $folder/Results/Fasta_SNPs/$folder.neutral_SNPs.fasta`;
`mv $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.pos1.fasta $folder/Results/Fasta_SNPs/$folder.pos1_SNPs.fasta`;
`mv $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.pos2.fasta $folder/Results/Fasta_SNPs/$folder.pos2_SNPs.fasta`;
`mv $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.pos3.fasta $folder/Results/Fasta_SNPs/$folder.pos3_SNPs.fasta`;
`mv $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.with_syn.tab.Syn.fasta $folder/Results/Fasta_SNPs/$folder.Syn_SNPs.fasta`;
`mv $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.with_syn.tab.NotSyn.fasta $folder/Results/Fasta_SNPs/$folder.NotSyn_SNPs.fasta`;
`mv $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.gene_matrix.sel $folder/Results/$folder.Gene_patterns`;
`mv $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.pattern_list $folder/Results/$folder.SNV_patterns`;


if($tree eq "y")
{

`mkdir $folder/Results/trees`;
`mv $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.core.tree $folder/Results/trees/$folder.core_SNPs.tree`;
`mv $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.intergenic.tree $folder/Results/trees/$folder.intergenic_SNPs.tree`;
`mv $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.neutrali.tree $folder/Results/trees/$folder.neutral_SNPs.tree`;
`mv $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.pos1.tree $folder/Results/trees/$folder.pos1_SNPs.tree`;
`mv $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.pos2.tree $folder/Results/trees/$folder.pos2_SNPs.tree`;
`mv $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.pos3.tree $folder/Results/trees/$folder.pos3_SNPs.tree`;

`mv $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.with_syn.tab.Syn.tree $folder/Results/trees/$folder.Syn_SNPs.tree`;
`mv $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.with_patterns.withframes.with_syn.tab.NotSyn.tree $folder/Results/trees/$folder.NotSyn_SNPs.tree`;
}


} # if steps

else
{}

# MOVE FOLDERS

`mkdir $folder\_OUT`;
`mv $folder/Results/ $folder\_OUT/`;
`mv $folder/tmp/ $folder\_OUT/`;
`mv $folder/reference/ $folder\_OUT/tmp/`;

# ANALYSIS of HOMOPLASY 

if ($homo eq "y")
{
	
`cp $folder\_OUT/Results/Fasta_SNPs/$folder.core_SNPs.fasta \.`;
`mv $folder\_OUT/Results/$folder.SNPs.annotation_tab \.`;
`mv $folder.SNPs.annotation_tab $folder.SNPs.annotation_tab_nohomoplasy`;
`cp Purple_scripts/Add_homoplasy.pl \.`;
`noisy $folder.core_SNPs.fasta`;
`perl Add_homoplasy.pl $folder.core_SNPs.fasta $folder.SNPs.annotation_tab_nohomoplasy > $folder.SNPs.annotation_tab`;
`rm Add_homoplasy.pl`;
`rm $folder.core_SNPs_sta.gr*`;
#`rm $folder.core.SNPs_sta.gr_mod`;
`rm $folder.core_SNPs_typ.eps`;
`rm $folder.core_SNPs_idx.txt`;

`mv $folder.core_SNPs_out.fas $folder.core_SNPs.NoHomoplasy.fasta`;
`sed -i "s/ //g" $folder.core_SNPs.NoHomoplasy.fasta`;

if ($tree eq "y")
{
`fasttree -nt -gtr -verbose 0 $folder.core_SNPs.NoHomoplasy.fasta > $folder.core_SNPs.NoHomoplasy.tree`;
`perl Purple_scripts/Patternlist2nex_append.pl $folder\_OUT/Results/$folder.Allpatterns $folder.core_SNPs.NoHomoplasy.tree NoHomoplasy >> $folder\_OUT/Results/$folder.nex`;
`mv $folder.core_SNPs.NoHomoplasy.tree $folder\_OUT/Results/trees/`;
}
else
{}

`rm $folder.core_SNPs.fasta`;
`mv $folder.core_SNPs.NoHomoplasy.fasta $folder\_OUT/Results/Fasta_SNPs/`;
`mv $folder.SNPs.annotation_tab $folder\_OUT/Results/`;
`mv $folder.SNPs.annotation_tab_nohomoplasy $folder\_OUT/tmp/tables/`;

}


if ($steps eq "call")
{`mv $folder\_OUT/tmp/snp_tab $folder/`;}

if ($indel eq "n")
{
  `rm -r $folder\_OUT/tmp/snp_tab`;
  `rm -r $folder/snp_tab`;
  `mv $folder/snp_tab_with_indel $folder/snp_tab`
}


# multifasta conc

`mv $folder\_OUT/Results/$folder.multifaconc.fasta $folder\_OUT/Results/$folder.multifaconc.fasta.tmp`;
open(FFF,"<$folder\_OUT/Results/$folder.multifaconc.fasta.tmp");
open(OOO,">$folder\_OUT/Results/$folder.multifaconc.fasta");

$c_nn=0;
while(<FFF>)
{
	if ($c_nn == 0)
	{
		print OOO $_;
		$c_nn = 1;
	}

	else
	{
		$subsub=substr($_,0,length($_)-$num_nn-1);
		print OOO $subsub,"\n";
		$c_nn = 0;	
	}

}

`rm $folder\_OUT/Results/$folder.multifaconc.fasta.tmp`;

# Convert Annotation tab to VCF 

`perl Purple_scripts/Annotation_tab2VCF.pl $folder\_OUT/Results/$folder.SNPs.annotation_tab $homo`;

`rm $ref.conc.fasta`;



