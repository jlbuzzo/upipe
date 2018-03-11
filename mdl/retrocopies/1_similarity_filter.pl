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

my ( $file, $genome, $usage ) = ( "", "", 0 );
GetOptions( "f|file=s"      => \$file,
           "g|genome=s"      => \$genome,
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
	$clusters{$id}->[1] =~ s/chr//g;
	#min first
	$clusters{$id}->[2] = $tokens[1];
	#max first
	$clusters{$id}->[3] = $tokens[2];

	#second chr
	$clusters{$id}->[4] = $tokens[5];
	#ADD A CONDITION TO VERIFY THE REFERENCE GENOME FORMAT
	$clusters{$id}->[4] =~ s/chr//g;
	#min second
	$clusters{$id}->[5] = $tokens[6];
	#max second
	$clusters{$id}->[6] = $tokens[7];
}


my $border = 500;
foreach my $cluster ( keys %clusters ) {

	#Legacy
    #my $seq = BioInfoGenome::extract2bit($genome,"$clusters{$cluster}->[1]:".($clusters{$cluster}->[2]-$border)."-".($clusters{$cluster}->[3]+$border));
    #my $seq2 = BioInfoGenome::extract2bit($genome,"$clusters{$cluster}->[4]:".($clusters{$cluster}->[5]-$border)."-".($clusters{$cluster}->[6]+$border));

	open (OUT1,">/dev/shm/a.$RAND.bed");
	print OUT1 "$clusters{$cluster}->[1]\t".($clusters{$cluster}->[2]-$border)."\t".($clusters{$cluster}->[3]+$border)."\t$clusters{$cluster}->[1]_$clusters{$cluster}->[2]_$clusters{$cluster}->[3]\n";
	close(OUT1);

	open (OUT2,">/dev/shm/b.$RAND.bed");
	print OUT2 "$clusters{$cluster}->[4]\t".($clusters{$cluster}->[5]-$border)."\t".($clusters{$cluster}->[6]+$border)."\t$clusters{$cluster}->[4]_$clusters{$cluster}->[5]_$clusters{$cluster}->[6]\n";
	close(OUT2);

<<<<<<< HEAD
	`bedtools getfasta -fi $genome -bed /dev/shm/a.bed -fo /dev/shm/a.fa`;
	`bedtools getfasta -fi $genome -bed /dev/shm/b.bed -fo /dev/shm/b.fa`;
=======
	`bedtools getfasta -fi $genome -bed /dev/shm/a.$RAND.bed -fo /dev/shm/a.$RAND.fa`;
	`bedtools getfasta -fi $genome -bed /dev/shm/b.$RAND.bed -fo /dev/shm/b.$RAND.fa`;
>>>>>>> 08e82b936125bc1b3fd8cb75a695b4e02f6fae7c

	`blat -mask=lower -tileSize=12 -minIdentity=75 -minScore=30 -out=pslx /dev/shm/a.$RAND.fa /dev/shm/b.$RAND.fa /dev/shm/ab.$RAND.pslx`; 
	`blat -mask=lower -tileSize=12 -minIdentity=75 -minScore=30 -out=pslx /dev/shm/b.$RAND.fa /dev/shm/a.$RAND.fa /dev/shm/ba.$RAND.pslx`; 
	my $remove = 0;

	open(IN1,"< /dev/shm/ab.$RAND.pslx");
	while (<IN1>) {
		chomp $_;
		my @tokens = split ( /[ \t\cI]+/,$_ );
		if ( $_ =~ /^[0-9]/ ) {
			if ( $tokens[0]/$tokens[12] > 0.8 ) {
				$remove = 1;
			}
		}
	}
	close(IN1);

	open(IN1,"< /dev/shm/ba.$RAND.pslx");
	while (<IN1>) {
		chomp $_;
		my @tokens = split ( /[ \t\cI]+/,$_ );
		if ( $_ =~ /^[0-9]/ ) {
			if ( $tokens[0]/$tokens[12] > 0.8 ) {
				$remove = 1;
			}
		}
	}
	close(IN1);

	if ( ! $remove ) { 
		print "$clusters{$cluster}->[0]IN\n"
		#print "$clusters{$cluster}->[1] $clusters{$cluster}->[2] $clusters{$cluster}->[3] $clusters{$cluster}->[4] $clusters{$cluster}->[5] $clusters{$cluster}->[6]\n";
	}
	else {
		print "$clusters{$cluster}->[0]OUT_SIMILARITY\n"
	}
	`rm /dev/shm/a.$RAND.bed`;
	`rm /dev/shm/b.$RAND.bed`;
	`rm /dev/shm/a.$RAND.fa`;
	`rm /dev/shm/b.$RAND.fa`;
	`rm /dev/shm/ab.$RAND.pslx`;
	`rm /dev/shm/ba.$RAND.pslx`;


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
