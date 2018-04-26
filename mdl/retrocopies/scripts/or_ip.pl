#!/usr/bin/perl
#####################################
# Program: or_ip.pl  -  Date: Thu Jul 31 15:20:45 BRT 2014
# Autor: Fabio C. P. Navarro - Ludwig
# Goal:
#
# Input:
#
# Output:
#
#####################################


use Getopt::Long;
use Statistics::Basic qw(:all nofill);
use strict;
use warnings;

my ( $desc, $file ) = ( "",  "" );
GetOptions( "d|desc=s"      => \$desc,
			"f|file=s"      => \$file  )
    || die "Error while parsing comand line arguments";

#chr2   198368040  198368072  32       -  chrX   964764     965317     553   521      12     HSPE1_protein_coding_6787           IN  IN  IN     IN
my @tokens = split ( /[ \t\cI]+/,$desc );
my $chr_par = $tokens[0];
$chr_par =~ s/chr//g;
my $start_par = $tokens[1];
my $end_par = $tokens[2];

my $chr_ip = $tokens[5];
$chr_ip =~ s/chr//g;
my $start_ip = $tokens[6];
my $end_ip = $tokens[7];

if ( $chr_ip eq $chr_par ) {
	$chr_par = "=";
}

#HG00107.ERR229778.37654657 163 15 40854235 60 101M = 40854530 395 CAATATGCTCCCATTCTCAACAATCAATCTATTTATGTAAGTTTTTCAAACTCCAGCATCAGAAATCCACATGTAGCCCTGGGGCCATCACATTTTCTTGC
#HG00107.ERR229778.99998091 177 15 40854240 37 101M 7 26246001 0 TGCTCCCATTCTCAACAATCAATCTATTTATGTAAGTTTTTCAAACTCCAGCATCAGAAATCCACATGTAGCCCTGGGGCCATCACATTTTCTTGCAGGAT
#HG00107.ERR229778.74449444 99 15 40854241 60 101M = 40854494 352 GCTCCCATTCTCAACAATCAATCTATTTATGTAAGTTTTTCAAACTCCAGCATCAGAAATCCACATGTAGCCCTGGGGCCATCACATTTTCTTGCAGGATG
my $min_up = 3000000000;
my $max_down = -1;
my $strand = "";
my $color = "";
my $SUP = 5;
my $support = 0;
my $support_exonic = 0;
my $output = "";
my $up_count = 0;
my $down_count = 0;

my (@cov, @up_cov, @down_cov) = ((), (), ());

my %exons = ();
open(IN0,"<$file");
while (<IN0>) {
	chomp $_;
	my @tokens = split ( /[ \t\cI]+/,$_ );
	$exons{$tokens[3]} = \@tokens;
}

my $min_par = 3000000000;
my $max_par = -1;
foreach my $exon ( keys %exons ) {
	$min_par = $exons{$exon}->[1] if ( $exons{$exon}->[1] < $min_par );
	$max_par = $exons{$exon}->[2] if ( $exons{$exon}->[2] > $max_par );
}

while ( <> ) {
	chomp $_;
	my @tokens = split ( /[ \t\cI]+/,$_ );

	if ( ($tokens[7] eq $chr_par || $tokens[7] eq "chr".$chr_par) && $tokens[8] >= $min_par-1000 && $tokens[8] <= $max_par+1000 ) {
			$support++;

			my @cigar = split (//,$tokens[6]);
			my $sum = 0;
			my $num = "";
			foreach my $chr ( @cigar ) {
				if ( $chr =~ /[0-9]/ ) {
					$num .= $chr;
				}
				else {
					if ( $chr eq "M" || $chr eq "D" || $chr eq "N" || $chr eq "X" ) {
						$sum += int($num);
					}
					$num = "";
				}
			}

			my $exonic = 0;
			foreach my $exon ( keys %exons ) {
				if ( $tokens[8] >= $exons{$exon}->[1]-15 && $tokens[8] <= $exons{$exon}->[2]+15 ) {
					$exonic = 1;
				}
			}
			if ( $exonic ) {
				$support_exonic++;
			}

			if ( $tokens[2] & 16 ) { #DOWN_STREAM
				$min_up = $tokens[4] if ( $tokens[4] < $min_up );
				$down_count++;
				$strand = "-";
				if ( $exonic ) {
					$color = "250,101,54";
				}
				else {
					$color= "255,192,173";
				}

				for ( my $i = $tokens[4]; $i <= $tokens[4]+$sum; $i++ ) {
					push(@down_cov,$i);
				}
			}
			else { #UP_STREAM
				$max_down = $tokens[4]+$sum if ( $tokens[4]+$sum > $max_down );
				$up_count++;
				$strand = "+";
				if ( $exonic ) {
					$color = "46,211,164";
				} else {
					$color = "162,239,215";
				}	

				for ( my $i = $tokens[4]; $i <= $tokens[4]+$sum; $i++ ) {
					push(@up_cov,$i);
				}
			}

			#track name="retroCNV" description="" itemRgb="On"
			#chr15 40854290 40854357 HG00097SRR741384.28540818 1000 - 40854290 40854357 255,0,0
			#chr7 26245985 26245995 HG00097SRR741384.28540818 1000 - 26245985 26245995 255,0,0
			$output .= "chr$tokens[3] $tokens[4] ".($tokens[4]+$sum)." $tokens[1] 1000 $strand $tokens[4] ".($tokens[4]+$sum)." $color\n";
			if ( $chr_par eq "=" ) {
				$output .= "chr$tokens[3] $tokens[8] ".( $tokens[8]+10 )." $tokens[1] 1000 $strand $tokens[8] ".( $tokens[8]+10 )." $color\n";
			} else {
				$output .= "chr$tokens[7] $tokens[8] ".( $tokens[8]+10 )." $tokens[1] 1000 $strand $tokens[8] ".( $tokens[8]+10 )." $color\n";
			}
	}
}

@up_cov = sort {$a <=> $b} @up_cov;
pop(@up_cov);
shift(@up_cov);

@down_cov = sort {$a <=> $b} @down_cov;
pop(@down_cov);
shift(@down_cov);

my ($up_median,$up_sd,$down_median,$down_sd) = (median(@up_cov),stddev(@up_cov),median(@down_cov),stddev(@down_cov));
my ($up_start,$up_end,$down_start,$down_end) = (0,0,0,0);
my ( $overlap_start, $overlap_end ) = ( 0, 0 );


if ( $support >= $SUP && $support_exonic >= $SUP && $min_up != 3000000000 && $max_down != -1 ) {
	($up_start,$up_end,$down_start,$down_end) = ($up_median-2*$up_sd, $up_median+2*$up_sd, $down_median-2*$down_sd, $down_median+2*$down_sd);

	if ( $up_start > $down_start ) {
		$overlap_start = $up_start;
	} else {
		$overlap_start = $down_start;	
	}
	
	if ( $up_end < $down_end ) {
		$overlap_end = $up_end;
	} else {
		$overlap_end = $down_end;	
	}
}


#open(OUT0,">cov.temp");
#print OUT0 join("\n", @cov);
#close(OUT0);

print "track name=\"retroCNV\" description=\"\" itemRgb=\"On\"\n";
print "$output";

if ( $support >= $SUP ) {
	if ( $support_exonic >= $SUP && $support_exonic/$support >= 0.85 ) {
		if ( ( ($up_end-$up_start) == 0 || ($down_end-$down_start) == 0 || $overlap_end-$overlap_start < 0 ) ||  ( ($overlap_end-$overlap_start)/($up_end-$up_start) <= 0.7 && ($overlap_end-$overlap_start)/($down_end-$down_start) <= 0.7) ) {
			my $size = $max_down-$min_up;
			if ( $min_up < $max_down ) {
				print STDERR "$desc IN_$size\_chr$chr_ip:$min_up-$max_down\_$up_count\_$down_count\n";
			}
			else {
				print STDERR "$desc IN_$size\_chr$chr_ip:$max_down-$min_up\_$up_count\_$down_count\n";
			}
		} 
		else {
			print STDERR "$desc OUT_OVERLAP\n";
		}
	}
	else {
		print STDERR "$desc OUT_EXON\n";
	}
}
else {
	print STDERR "$desc OUT_SUP\n";
}

sub usage {
    die "Usage: $0 [options]
Available Options :
   -f  | --file         : Filename in.
   -h  | --help --usage : Print this message.
";

}
