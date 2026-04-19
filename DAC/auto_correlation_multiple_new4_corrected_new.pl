
#! perl -w

open (INPUT1, "<$ARGV[0]");
open (INPUT2, "<$ARGV[1]");
open (OUTPUT1, ">$ARGV[2]");

$max = 1000;

while (<INPUT1>) {
	$line = $_;
	chomp($line);

	@data = split(/\t/, $line);

	push(@all_data1, $data[1]);
	push(@all_count1, $data[2]);
}

while (<INPUT2>) {
	$line = $_;
	chomp($line);

	@data = split(/\t/, $line);

	push(@all_data2, $data[1]);
	push(@all_count2, $data[2]);
}

#print "here\n";

%occurrence = ();

for ($i = 0; $i <= $max; $i++) {
	$occurrence{$i} = 0;
}

#for ($i = 1; $i <= $max; $i++) {
#	$occurrence_sum = 0;
#	$length_sum = 0;
#	$Pn = 0;

	for ($j = 0; $j < @all_data1; $j++) {

		for ($m = $j; $m < @all_data2; $m++) {
			$dist = $all_data2[$m] - $all_data1[$j];
			$abs_dist = abs($dist);

			if ($abs_dist == 0) {
				$occurrence{$abs_dist} = $occurrence{$abs_dist} + ($all_count1[$j] * $all_count2[$m]);
			}
			elsif ($abs_dist > 0 && $abs_dist <= $max) {
				$occurrence{$abs_dist} = $occurrence{$abs_dist} + $all_count1[$j] * $all_count2[$m];
			}
			else {
				last;
			}
		}
	}

foreach $key (sort numeric keys %occurrence) {
	print "$key\t$occurrence{$key}\n";
	print OUTPUT1 "$key\t$occurrence{$key}\n";
}

sub numeric { $a <=> $b }