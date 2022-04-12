

$folder=shift;

# merge gene_matrix

@ls_gene=`ls $folder/Results/*.split_*.annotated.gene_matrix`;

chomp $ls_gene[0];
`head -n 1 $ls_gene[0] > $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.gene_matrix`;
`sed -i 1d $folder/Results/*.split_*.annotated.gene_matrix`;
`cat $folder/Results/*.split_*.annotated.gene_matrix >> $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated.gene_matrix`;

foreach $g(@ls_gene)
{
	chomp $g;
	`rm $g`;
}

# merge annotated

@ls_ann=`ls $folder/Results/*.split_*.annotated`;

chomp $ls_ann[0];
`head -n 1 $ls_ann[0] > $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated`;
`sed -i 1d $folder/Results/*.split_*.annotated`;
`cat $folder/Results/*.split_*.annotated >> $folder/Results/$folder.mergedsnps.tab.selected_sorted.annotated`;

foreach $a(@ls_ann)
{
	chomp $a;
	`rm $a`;
}

# remove splitted files

`rm $folder/Results/*.selected_sorted.split_*`;
