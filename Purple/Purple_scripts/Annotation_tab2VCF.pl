
$annotation_tab=shift;
$homo=shift;


#`awk '\$2 == \"CORE\"' Salm.SNPs.annotation_tab > $annotation_tab.CORE`;

open(INN,"<$annotation_tab");
open(VCF,">$annotation_tab.vcf");

# print header
print VCF "##fileformat=VCFv4.0\n";
print VCF "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	";

$c=0;
while(<INN>)
{
    chomp $_;
	@s=split("\t",$_);

    if ($homo eq "y"){$ref_col=20;}
    else{$ref_col=19;}

    if ($c > 0)
     {
	 	 
	  %hash_SNP=();
	  $hash_SNP{$s[$ref_col]}=0;
	  $base=0;
	  $base_str=$s[$ref_col];
	  	  
	  $vcf_str="";
	  
	  $vcf_str=$vcf_str.$hash_SNP{$s[$ref_col]}."\t";
	  
	  for $i($ref_col+1..$#s)
	  {

	    chomp $s[$i];
	  
	    if (! exists $hash_SNP{$s[$i]}){$base++;$hash_SNP{$s[$i]}=$base;$base_str=$base_str.",".$s[$i];}
	    else{}
	    
	    $vcf_str=$vcf_str.$hash_SNP{$s[$i]}."\t";
	  
	  }

	       # line, position is "padded"

	  chomp $base_str;     
	  chomp $vcf_str;
	  
	  if ($s[1] eq "NA"){$s[1]="Stretch";}
	  else{}



        if ($homo eq "y")
        {
          print VCF $s[4],"\t",$s[7],"\t",$s[0],"\t",$s[20],"\t",$base_str,"\t",".","\t",$s[1],"\t",$s[10],",",$s[11],",",$s[12],",",$s[13],",",$s[19],"\t","GT","\t",$vcf_str,"\n";
        }

        else
        {
          print VCF $s[3],"\t",$s[6],"\t",$s[0],"\t",$s[19],"\t",$base_str,"\t",".","\t",$s[1],"\t",$s[9],",",$s[10],",",$s[11],",",$s[12],",",$s[18],"\t","GT","\t",$vcf_str,"\n";
        }

     }
    
    else
     {$c++;
       # print header
       print VCF join("\t",@s[$ref_col..$#s]),"\n";   
     }

}






