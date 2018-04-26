#!/usr/bin/perl
#####################################
# Program: similarity_filter.pl  -  Date: Sat Dec 17 11:55:02 BRST 2011
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

my ( $file, $genome, $pd, $usage ) = ( "", "", "/tmp", 0 );
GetOptions( "f|file=s"      => \$file,
           "g|genome=s"      => \$genome,
		   "p|processing_dir=s" => \$pd,
           "h|help|usage"  => \$usage )
    || die "Error while parsing comand line arguments";

usage() if( $usage );

if( $file eq "" || $genome eq "") {
    usage();
}
#if( $file2 eq "" ) {
#    usage();
#}


#chr22 39355953 39358192 (2239) - chr22 39479455 39483340 (3885) 1646 6

my %clusters = ();
my $id = "";
my $RAND = int(rand(20000));

open ( IN0, "<$file" );
while ( <IN0> ) {
	chomp $_;
	my @tokens = split ( /[ \t\cI]+/,$_ );
	
	$id = "$tokens[0]_$tokens[1]_$tokens[2]_$tokens[5]_$tokens[6]_$tokens[7]";
	$clusters{$id}->[0] .= "$_ ";

	#chr first
	$clusters{$id}->[1] = $tokens[0];
	#ADD A CONDITION TO VERIFY THE REFERENCE GENOME FORMAT
	#$clusters{$id}->[1] =~ s/chr//g;
	#min first
	$clusters{$id}->[2] = $tokens[1];
	#max first
	$clusters{$id}->[3] = $tokens[2];

	#second chr
	$clusters{$id}->[4] = $tokens[5];
	#ADD A CONDITION TO VERIFY THE REFERENCE GENOME FORMAT
	#$clusters{$id}->[4] =~ s/chr//g;
	#min second
	$clusters{$id}->[5] = $tokens[6];
	#max second
	$clusters{$id}->[6] = $tokens[7];
}


my $border = 500;
foreach my $cluster ( keys %clusters ) {


	open (OUT1,">$pd/a.$RAND.bed");
	print OUT1 "$clusters{$cluster}->[1]\t".($clusters{$cluster}->[2]-$border)."\t".($clusters{$cluster}->[3]+$border)."\t$clusters{$cluster}->[1]_$clusters{$cluster}->[2]_$clusters{$cluster}->[3]\n";
	close(OUT1);

	open (OUT2,">$pd/b.$RAND.bed");
	print OUT2 "$clusters{$cluster}->[4]\t".($clusters{$cluster}->[5]-$border)."\t".($clusters{$cluster}->[6]+$border)."\t$clusters{$cluster}->[4]_$clusters{$cluster}->[5]_$clusters{$cluster}->[6]\n";
	close(OUT2);

	`bedtools getfasta -fi $genome -bed $pd/a.$RAND.bed -fo $pd/a.$RAND.fa`;
	`bedtools getfasta -fi $genome -bed $pd/b.$RAND.bed -fo $pd/b.$RAND.fa`;

	`blat -mask=lower -tileSize=12 -minIdentity=75 -minScore=30 -out=pslx $genome $pd/a.$RAND.fa $pd/ab.$RAND.pslx`; 
	`blat -mask=lower -tileSize=12 -minIdentity=75 -minScore=30 -out=pslx $genome $pd/b.$RAND.fa $pd/ba.$RAND.pslx`; 
	my $remove = 0;
	my @a_alignments = [];
	my @a_coord = [];
	my @b_alignments = [];
	my @b_coord = [];

	open(IN1,"< $pd/ab.$RAND.pslx");
	while (<IN1>) {
		chomp $_;
		my @tokens = split ( /[ \t\cI]+/,$_ );

		if ( $_ =~ /^[0-9]/ ) {
			my $len = $tokens[12];
			my $qgap = $tokens[5];
			my $matches = $tokens[0]; 
			@a_coord = split( /[:-]/,$tokens[9] );
			if ( $qgap/$len < 2/3 && $matches/$len > 0.4 ) {
				push(@a_alignments,[$tokens[13],$tokens[18],$tokens[20]]);
			}
		}
	}
	close(IN1);

	open(IN1,"< $pd/ba.$RAND.pslx");
	while (<IN1>) {
		chomp $_;
		my @tokens = split ( /[ \t\cI]+/,$_ );

		if ( $_ =~ /^[0-9]/ ) {
			my $len = $tokens[12];
			my $qgap = $tokens[7];
			my $matches = $tokens[0];
			@b_coord = split( /[:-]/,$tokens[9] );
			if ( $qgap/$len < 2/3 && $matches/$len > 0.4 ) {
				push(@b_alignments,[$tokens[13],$tokens[18],$tokens[20]]);
			}
		}
	}
	close(IN1);

	foreach ( my $i = 1; $i < scalar(@a_alignments); $i++ ) {
		if ( $remove == 1 ) {
			last;
		}
		if ( $a_alignments[$i]->[0] eq $b_coord[0] ) {
			my @t_len = split( /[,]/,$a_alignments[$i]->[1] );
			my @t_start = split( /[,]/,$a_alignments[$i]->[2] );
			
			for ( my $j = 0; $j < scalar(@t_start); $j++ ) {
				if ( ($t_start[$j] && $t_len[$j]) && ( $b_coord[1] >= $t_start[$j] && $b_coord[1] <= ($t_start[$j]+$t_len[$j]) ||
				     $b_coord[2] >= $t_start[$j] && $b_coord[2] <= ($t_start[$j]+$t_len[$j]) ||
				     $b_coord[1] <= $t_start[$j] && $b_coord[2] >= ($t_start[$j]+$t_len[$j]) ) ) {
#Debug
#					print "$a_alignments[$i]->[0] = $b_coord[0];  $b_coord[1] : $b_coord[2] =  $t_start[$j] ".($t_start[$j]+$t_len[$j])."\n";
					$remove = 1;
					last;
				}
			}
		}
	}

	if ( $remove == 0 ) {
		foreach ( my $i = 1; $i < scalar(@b_alignments); $i++ ) {
			if ( $b_alignments[$i]->[0] eq $a_coord[0] ) {
				my @t_len = split( /[,]/,$b_alignments[$i]->[1] );
				my @t_start = split( /[,]/,$b_alignments[$i]->[2] );
				
				for ( my $j = 0; $j < scalar(@t_start); $j++ ) {
					if ( ($t_start[$j] && $t_len[$j]) && ($a_coord[1] >= $t_start[$j] && $a_coord[1] <= ($t_start[$j]+$t_len[$j]) ||
					     $a_coord[2] >= $t_start[$j] && $a_coord[2] <= ($t_start[$j]+$t_len[$j]) ||
					     $a_coord[1] <= $t_start[$j] && $a_coord[2] >= ($t_start[$j]+$t_len[$j]) ) ) {
#Debug
#						print "$b_alignments[$i]->[0] = $a_coord[0];  $a_coord[1] : $a_coord[2] =  $t_start[$j] ".($t_start[$j]+$t_len[$j])."\n";	
						$remove = 1
					}
				}
			}
		}
	}	
	

	if ( ! $remove ) { 
		print "$clusters{$cluster}->[0]IN\n"
		#print "$clusters{$cluster}->[1] $clusters{$cluster}->[2] $clusters{$cluster}->[3] $clusters{$cluster}->[4] $clusters{$cluster}->[5] $clusters{$cluster}->[6]\n";
	}
	else {
		print "$clusters{$cluster}->[0]OUT_SIMILARITY\n"
	}
	`rm $pd/a.$RAND.bed`;
	`rm $pd/b.$RAND.bed`;
	`rm $pd/a.$RAND.fa`;
	`rm $pd/b.$RAND.fa`;
	`rm $pd/ab.$RAND.pslx`;
	`rm $pd/ba.$RAND.pslx`;


}
#foreach cluster
#  fetch seq
#  blat against each other
#  mask

#  see if they align on each other
#  if true flag it
#
#}

#print all unflagged
#open ( IN1, "<$file2" );
#while ( <IN1> ) {
#	chomp $_;
#	my @tokens = split ( /[ \t\cI]+/,$_ );
#}



sub usage {
    die "Usage: $0 [options]
Available Options :
   -f  | --file         : Filename in.
   -g  | --genome       : Reference genome
   -h  | --help --usage : Print this message.
";

}
