#!/usr/bin/perl
#####################################
# Program: remove_rep.pl  -  Date: Mon Sep  9 14:06:16 BRT 2013
# Autor: Fabio C. P. Navarro - Ludwig
# Goal:
#
# Input:
#
# Output:
#
#####################################


use Getopt::Long;
use strict;
use warnings;
use File::Basename;

my ( $file, $file2, $usage ) = ( "", "", 0 );
GetOptions( "f|file=s"      => \$file,
           "f2|file2=s"      => \$file2,
           "h|help|usage"  => \$usage )
    || die "Error while parsing comand line arguments";

usage() if( $usage );

if( $file eq "" ) {
    usage();
}

#chr1    16777160        16777470        SINE
#chr1    25165800        25166089        SINE

my $WINDOW = 10000;
my $RAND = int(rand(20000));

my %reps = ();
my $count = 0;
open ( IN0, "<$file" ) or die "Can't open '$file': $!";
while ( <IN0> ) {
	chomp $_;
	my @tokens = split ( /[ \t\cI]+/,$_ );
	for (my $i=int($tokens[1]/$WINDOW); $i<=int($tokens[2]/$WINDOW); $i++ ) {
		push(@{$reps{$tokens[0]}->{$i}},\@tokens);
	}
}

#		my @files = glob( './*.max_ins_size' );

#foreach my $file2 ( @files ) {
	my $output = "";
	my $output_final = "";
	open (IN1, "<$file2");
	#chr7 26242115 26252719 (10604) - chr2 33141304 33141592 (288) 10316 11
	#chr7 26242594 26251970 (9376) - chr12 2895872 2896604 (732) 8644 17
	while ( <IN1> ) {
		$count++;
		chomp $_;
		my $event = "$_";
		my @tokens = split ( /[ \t\cI]+/,$_ );
		my $print=1;
		#FOREACH REPETITIVE ELEMENT
		my $window = int($tokens[6]/$WINDOW);
		my $window_ip = int($tokens[1]/$WINDOW);
		my $rep_inside = "";
		my @reps_inside = "";
		my $status = "";

		foreach my $rep ( @{$reps{$tokens[0]}->{$window_ip-1}},@{$reps{$tokens[0]}->{$window_ip}},@{$reps{$tokens[0]}->{$window_ip+1}} ) {
			if ( ( $tokens[1] >= $rep->[1] && $tokens[1] <= $rep->[2] ) || ( $tokens[2] >= $rep->[1] && $tokens[2] <= $rep->[2] ) || 
                 ( $tokens[1] <= $rep->[1] && $tokens[2] >= $rep->[2] ) ) {

				#MY PARENTAL LOCI IS INSIDE THIS REP
				my $min = 0;
				if ( $tokens[1] > $rep->[1] ) {
					$min = $tokens[1];
				}
				else {
					$min = $rep->[1];
				}
				my $max = 0;
				if ( $tokens[2] < $rep->[2] ) { 
					$max = $tokens[2];
				}
				else {
					$max = $rep->[2];
				}

				if ( ($max-$min)/(($tokens[2]-$tokens[1])) >= 0.4 ) {
					$print = 0;
					$status = "PAR";
					last;
				}
			}
		}

		foreach my $rep ( @{$reps{$tokens[5]}->{$window-1}},@{$reps{$tokens[5]}->{$window}},@{$reps{$tokens[5]}->{$window+1}} ) {
			if ( ( $tokens[6] >= $rep->[1] && $tokens[6] <= $rep->[2] ) || ( $tokens[7] >= $rep->[1] && $tokens[7] <= $rep->[2] ) || 
                 ( $tokens[6] <= $rep->[1] && $tokens[7] >= $rep->[2] ) ) {

				if ( $rep->[3] eq "Low_complexity" || $rep->[3] eq "Simple_repeat" || $rep->[3] eq "Satellite" ) {
					push (@reps_inside,$rep);
				}

				#IM INSIDE THIS REP
				my $min = 0;
				if ( $tokens[6] > $rep->[1] ) {
					$min = $tokens[6];
				}
				else {
					$min = $rep->[1];
				}
				my $max = 0;
				if ( $tokens[7] < $rep->[2] ) { 
					$max = $tokens[7];
				}
				else {
					$max = $rep->[2];
				}

				if ( ($max-$min)/($tokens[7]-$tokens[6]) >= 0.9 ) {
					$print = 0;
					$status = "IP";
					last;
				}
			}
		}

			
		if ( $print ) {
			#There isnt a rep element that overlaps 90% of the insertion point. The sum of rep elements overlaps > 90% of the insertoin point?
			#$output .= "$_ OUT_REP\n";
			
			shift(@reps_inside);
			if ( scalar(@reps_inside) >= 2 ) {
				@reps_inside = sort {$a->[1] <=> $b->[1]} @reps_inside;
				my $fn = basename($file2);
				open (OUT1,">/dev/shm/$fn.$RAND.temp_rep_el");
				foreach my $rep ( @reps_inside ) {
					print OUT1 "$rep->[0]\t$rep->[1]\t$rep->[2]\t$rep->[3]\n";
				}
				close(OUT1);
				`mergeBed -i /dev/shm/$fn.$RAND.temp_rep_el > /dev/shm/$fn.$RAND.temp_rep_el.merged`;
				`rm /dev/shm/$fn.$RAND.temp_rep_el`;
				open (IN2,"< /dev/shm/$fn.$RAND.temp_rep_el.merged");
				my @buffer = ();
				while( <IN2> ) {
					my @tokens = split ( /[ \t\cI]+/,$_ );
					push(@buffer,\@tokens);
				}
				close(IN2);
				`rm /dev/shm/$fn.$RAND.temp_rep_el.merged`;
	
				my $sum = 0;
	
				foreach my $rep ( @buffer ) {
					if ( ( $tokens[6] >= $rep->[1] && $tokens[6] <= $rep->[2] ) || ( $tokens[7] >= $rep->[1] && $tokens[7] <= $rep->[2] ) || 
	                 ($tokens[6] <= $rep->[1] && $tokens[7] >= $rep->[2] ) ) {
	
						#IM INSIDE THIS REP
						my $min = 0;
						if ( $tokens[6] > $rep->[1] ) {
							$min = $tokens[6];
						}
						else {
							$min = $rep->[1];
						}
						my $max = 0;
						if ( $tokens[7] < $rep->[2] ) { 
							$max = $tokens[7];
						}
						else {
							$max = $rep->[2];
						}

						$sum += $max-$min;
					}
				}

				if ( ($sum)/($tokens[7]-$tokens[6]) >= 0.85 ) {
					$print = 0;
				}
				if ( $print ) {
					$output .= "$event IN\n";
					$output_final .= "$event IN\n";
				}
				else {
					$output .= "$event OUT_REP_MULT\n";
				}
			}
			else {
					$output .= "$event IN\n";
					$output_final .= "$event IN\n";
			}
		}
		else {
			$output .= "$event OUT_REP_ONE_$status\n";
		}
		if ( $count %10000 == 0 ) {
			print "$count\n";
		}
	}
	close(IN1);
	
	open (OUT1,">$file2.norep_debug");
	print OUT1 "$output";
	close(OUT1);
	open (OUT1,">$file2.norep");
	print OUT1 "$output_final";
	close(OUT1);
#}

#open ( IN1, "<$file2" );
#while ( <IN1> ) {
#	chomp $_;
#	my @tokens = split ( /[ \t\cI]+/,$_ );
#}



sub usage {
    die "Usage: $0 [options]
Available Options :
   -f  | --file         : Filename in.
   -h  | --help --usage : Print this message.
";

}
