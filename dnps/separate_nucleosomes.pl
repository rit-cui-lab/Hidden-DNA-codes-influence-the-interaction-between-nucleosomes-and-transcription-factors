#! perl -w

open (INPUT, "<$ARGV[0]");
open (OUTPUT, ">$ARGV[1]");

@minor_position_147 = (	[qw(5 6 7)],
			[qw(15 16 17)],
			[qw(26 27 28)],
			[qw(37 38 39)],
			[qw(47 48 49)],
			[qw(57 58 59)],
			[qw(88 89 90)],
			[qw(98 99 100)],
			[qw(108 109 110)],
			[qw(119 120 121)],
			[qw(130 131 132)],
			[qw(140 141 142)],
		);

@major_position_147 = (	[qw(10 11 12)],
			[qw(21 22)],
			[qw(32 33)],
			[qw(42 43 44)],
			[qw(52 53 54)],
			[qw(62 63 64)],
			[qw(83 84 85)],
			[qw(93 94 95)],
			[qw(103 104 105)],
			[qw(114 115)],
			[qw(125 126)],
			[qw(135 136 137)],
		);


#if ($ARGV[0] =~ /^>(\w+)/) {
#	$chr = $1;
#	$chr = "chr" . $chr;

#	$output1 = $chr . "_type1_fasta.txt";
#	open(OUTPUT1, ">$output1");

#	$output2 = $chr . "_type2_fasta.txt";
#	open(OUTPUT2, ">$output2");

#	$output3 = $chr . "_type3_fasta.txt";
#	open(OUTPUT3, ">$output3");

#	$output4 = $chr . "_type4_fasta.txt";
#	open(OUTPUT4, ">$output4");
#}

$type1_count = 0;
$type2_count = 0;
$type3_count = 0;
$type4_count = 0;

$coefficient = 1.0909;

while(<INPUT>) {
	$name = $_;
	chomp($name);
	$seq = <INPUT>;
	chomp($seq);
	$seq = uc($seq);

	$WW_minor = 0;
	$WW_major = 0;
	$SS_minor = 0;
	$SS_major = 0;

	for ($i=0; $i<= $#minor_position_147; $i++) {
		$pos = $minor_position_147[$i][0];
		$length = @{$minor_position_147[$i]};
		$frag = substr($seq, $pos - 1, $length + 1);

		if ($frag =~ /AAAA|AAAT|AATA|AATT|ATAA|ATAT|ATTA|ATTT|TAAA|TAAT|TATA|TATT|TTAA|TTAT|TTTA|TTTT/) {
			$WW_minor = $WW_minor + 3;
		}
		elsif ($frag =~ /AAA|AAT|ATA|ATT|TAA|TAT|TTA|TTT/) {
			$WW_minor = $WW_minor + 2;
		}
		elsif ($frag =~ /AA|AT|TA|TT/) {
			$WW_minor = $WW_minor + 1;
		}

		if ($frag =~ /GGGG|GGGC|GGCG|GGCC|GCGG|GCGC|GCCG|GCCC|CGGG|CGGC|CGCG|CGCC|CCGG|CCGC|CCCG|CCCC/) {
			$SS_minor = $SS_minor + 3;
		}
		elsif ($frag =~ /GGG|GGC|GCG|GCC|CGG|CGC|CCG|CCC/) {
			$SS_minor = $SS_minor + 2;
		}
		elsif ($frag =~ /GG|GC|CC|CG/){
			$SS_minor = $SS_minor + 1;
		}
	}

	for ($j=0; $j<= $#major_position_147; $j++) {
		$pos = $major_position_147[$j][0];
		$length = @{$major_position_147[$j]};
		$frag = substr($seq, $pos - 1, $length + 1);

		if (length($frag) == 4) {
			if ($frag =~ /AAAA|AAAT|AATA|AATT|ATAA|ATAT|ATTA|ATTT|TAAA|TAAT|TATA|TATT|TTAA|TTAT|TTTA|TTTT/) {
				$WW_major = $WW_major + 3;
			}
			elsif ($frag =~ /AAA|AAT|ATA|ATT|TAA|TAT|TTA|TTT/) {
				$WW_major = $WW_major + 2;
			}
			elsif ($frag =~ /AA|AT|TA|TT/) {
				$WW_major = $WW_major + 1;
			}

			if ($frag =~ /GGGG|GGGC|GGCG|GGCC|GCGG|GCGC|GCCG|GCCC|CGGG|CGGC|CGCG|CGCC|CCGG|CCGC|CCCG|CCCC/) {
				$SS_major = $SS_major + 3;
			}
			elsif ($frag =~ /GGG|GGC|GCG|GCC|CGG|CGC|CCG|CCC/) {
				$SS_major = $SS_major + 2;
			}
			elsif ($frag =~ /GG|GC|CC|CG/){
				$SS_major = $SS_major + 1;
			}
		}
		elsif (length($frag) == 3) {
			if ($frag =~ /AAA|AAT|ATA|ATT|TAA|TAT|TTA|TTT/) {
				$WW_major = $WW_major + 2;
			}
			elsif ($frag =~ /AA|AT|TA|TT/) {
				$WW_major = $WW_major + 1;
			}

			if ($frag =~ /GGG|GGC|GCG|GCC|CGG|CGC|CCG|CCC/) {
				$SS_major = $SS_major + 2;
			}
			elsif ($frag =~ /GG|GC|CC|CG/){
				$SS_major = $SS_major + 1;
			}
		}
	}

	if ($WW_minor >= ($WW_major * $coefficient) && $SS_minor <= ($SS_major * $coefficient)) {
#		print OUTPUT1 "$name\n";
#		print OUTPUT1 "$seq\n";
		$type1_count++;
	}
	elsif ($WW_minor >= ($WW_major * $coefficient) && $SS_minor > ($SS_major * $coefficient)) {
#		print OUTPUT2 "$name\n";
#		print OUTPUT2 "$seq\n";
		$type2_count++;
	}
	elsif ($WW_minor < ($WW_major * $coefficient) && $SS_minor <= ($SS_major * $coefficient))	{
#		print OUTPUT3 "$name\n";
#		print OUTPUT3 "$seq\n";
		$type3_count++;
	}
	elsif ($WW_minor < ($WW_major * $coefficient) && $SS_minor > ($SS_major * $coefficient)) {
#		print OUTPUT4 "$name\n";
#		print OUTPUT4 "$seq\n";
		$type4_count++;
	}	
}
print OUTPUT "$type1_count\t$type2_count\t$type3_count\t$type4_count\n";

close(INPUT);
close(OUTPUT);
#close(OUTPUT1);
#close(OUTPUT2);
#close(OUTPUT3);
#close(OUTPUT4);