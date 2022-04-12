

$in=shift;
$n=shift;

open(TAB,"<$in");
@arr=<TAB>;

$head=$arr[0];
@arr2=@arr[1..$#arr];

$lines=int(($#arr2+1)/$n);

$c=0;

for $i(1..$n-1)
{

	open(OUT,">$in.split_$i");

	print OUT $head;
	print OUT join("",@arr2[$c..($c+$lines)]);
	close(OUT);
	

	$c=$c+$lines+1;


}


	$i=$n;

	open(OUT,">$in.split_$i");

	print OUT $head;
	print OUT join("",@arr2[$c..$#arr2]);
	close(OUT);

