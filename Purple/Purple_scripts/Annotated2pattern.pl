

$sort=shift; # file _sorted.annotated
$org=shift; # number of organisms

open(SOR,"<$sort");


$pattern=1;
%hash_pat=();

open(TAB,">$sort.pattern_list");
open(ANN,">$sort.with_patterns");

$c=0;
while(<SOR>)
{
			chomp $_;
			@s=split("\t",$_);

	if ($c == 0){
			print ANN join("\t",$s[0]),"\tPattern\t",join("\t",@s[1..$#s]),"\n";$c++;

			print TAB "Pattern_name\t",join("\t",@s[5..5+$org-1]),"\n";
		    }

	else
		    {


		     $snp = join("\t",@s[5..5+$org-1]),"\n";
		     $snp =~ s/[qweryuiopsdfhjklzxvbnm]/\?/gi;

		     	if ($hash_pat{$snp} eq "")
		     	{
				$hash_pat{$snp}="Pattern_$pattern";
				print TAB $hash_pat{$snp},"\t",$snp,"\n";
				$pattern++;
			}

		    print ANN join("\t",$s[0]),"\t",$hash_pat{$snp},"\t",join("\t",@s[1..$#s]),"\n";

		    }

}
