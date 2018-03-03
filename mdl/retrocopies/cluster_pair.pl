#!/usr/bin/perl
#####################################
# Program: clean.pl  -  Date: Thu Sep 15 22:53:28 BRT 2011
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

my ( $window, $support, $usage ) = ( 4000, 10, 0 );
GetOptions( "s|support=i"      => \$support,
            "w|window=i"      => \$window,
            "h|help|usage"  => \$usage )
    || die  "Error while parsing comand line arguments";

usage() if( $usage );

#7  26251857  11  27828252
#7  26250582  15  30723341

my $line = 0;
my @buffer = ();
my @final = ();
while ( <> ) {
	chomp $_;
	my @tokens = split ( /[ \t\cI]+/,$_ );
	my @tmp = split ( /[\_\-]/,$tokens[0] );

	my @del = ();
	my $inserted = 0;
	for ( my $i = 0; $i < scalar(@buffer); $i++ ) {
		if ( $buffer[$i] && $tokens[3] <= $buffer[$i]->[2]+$window && $tokens[2] eq $buffer[$i]->[7] ) {
				if ( $tokens[3] > $buffer[$i]->[3] ) {
					$buffer[$i]->[3] = $tokens[3];
				}
				elsif ( $tokens[3] < $buffer[$i]->[2] ) {
					$buffer[$i]->[2] = $tokens[3];
				}

				if ( $tokens[1] > $buffer[$i]->[1] ) {
					$buffer[$i]->[1] = $tokens[1];
				}
				elsif ( $tokens[1] < $buffer[$i]->[0] ) {
					$buffer[$i]->[0] = $tokens[1];
				}
				$buffer[$i]->[4]++;
				$inserted = 1;
		}
		else {
			push ( @del, $i );
		}
	}
	foreach my $i ( @del ) {
		if ( $buffer[$i] && $buffer[$i]->[4] >= $support && ($buffer[$i]->[3]-$buffer[$i]->[2]) >= 60 ) {
			my @tmp = @{$buffer[$i]};
			push ( @final, \@tmp );
		}
		delete( $buffer[$i] );
	}
	if ( ! $inserted ) {
		my @tmp = ($tokens[1],$tokens[1],$tokens[3],$tokens[3],1,"","$tokens[0]","$tokens[2]");
		push ( @buffer,\@tmp );
	}
	$line++;
}
foreach ( my $i = 0; $i < scalar(@buffer); $i++ ) {
	if ( $buffer[$i] && $buffer[$i]->[4] >= $support && ($buffer[$i]->[3]-$buffer[$i]->[2]) >= 60 ) {
		my @tmp = @{$buffer[$i]};
		push ( @final, \@tmp );
	}
	delete( $buffer[$i] );
}

if ( 1 ) {
my $size = scalar(@final);
for ( my $i=0; $i < $size-1; $i++ ) {
	for ( my $j = $i+1; $j < $size; $j++ ) {
		if ( $final[$i] && $final[$j] && $final[$i]->[6] eq $final[$j]->[6] && $final[$i]->[7] eq $final[$j]->[7]) {
			if ( (	( $final[$i]->[0] >= $final[$j]->[0] && $final[$i]->[0] <= $final[$j]->[1] ) || 
				( $final[$j]->[0] >= $final[$i]->[0] && $final[$j]->[0] <= $final[$i]->[1] ) || 
				( $final[$i]->[1] >= $final[$j]->[0] && $final[$i]->[1] <= $final[$j]->[1] ) ||
				( $final[$j]->[1] >= $final[$i]->[0] && $final[$j]->[1] <= $final[$i]->[1] ) ) &&
			
			     (	( $final[$i]->[2] >= $final[$j]->[2] && $final[$i]->[2] <= $final[$j]->[3] ) || 
				( $final[$j]->[2] >= $final[$i]->[2] && $final[$j]->[2] <= $final[$i]->[3] ) || 
				( $final[$i]->[3] >= $final[$j]->[2] && $final[$i]->[3] <= $final[$j]->[3] ) ||
				( $final[$j]->[3] >= $final[$i]->[2] && $final[$j]->[3] <= $final[$i]->[3] ) ) ) {
				$final[$i]->[0] = $final[$j]->[0] if ( $final[$j]->[0] < $final[$i]->[0] );
				$final[$i]->[1] = $final[$j]->[1] if ( $final[$j]->[1] > $final[$i]->[1] );
				$final[$i]->[2] = $final[$j]->[2] if ( $final[$j]->[2] < $final[$i]->[2] );
				$final[$i]->[3] = $final[$j]->[3] if ( $final[$j]->[3] > $final[$i]->[3] );

				$final[$i]->[4] += $final[$j]->[4];
				$final[$i]->[5] .= $final[$j]->[5];
				delete ($final[$j]);
				$j--;
				$size--;
			}
		}
	}
}
}

my $DEBUG=0;

for ( my $i=0; $i < scalar(@final); $i++ ) {
	if ( $final[$i] ) {
		print "$final[$i]->[6] $final[$i]->[0] $final[$i]->[1] (".($final[$i]->[1]-$final[$i]->[0]).") - $final[$i]->[7] $final[$i]->[2] $final[$i]->[3] (".($final[$i]->[3]-$final[$i]->[2]).") ".(abs(($final[$i]->[3]-$final[$i]->[2]) -($final[$i]->[1]-$final[$i]->[0])))." $final[$i]->[4]\n";
		if ( $DEBUG ) {
			print "$final[$i]->[5]\n";
		}
	}
}

sub usage {
    die "Usage: $0 [options]
Available Options :
   -s  | --support         : Min cluster support.
   -w  | --window         : Distance to close clusters (Max distance between reads from a cluster).
   -h  | --help --usage : Print this message.
";

}
