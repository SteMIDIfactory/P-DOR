

$snp=shift;
$out_folder=shift;

open(SNN,"<$snp");
chomp $snp;
@snp_name=split("/",$snp);

open(OUT,">$out_folder/$snp_name[$#snp_name]");

while(<SNN>)
{
	if ($_ !~ /[\.N]/)
	{print OUT $_;}
}
