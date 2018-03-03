#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: gtf2bed.pl
#
#        USAGE: ./gtf2bed.pl [gene|exon] <path to gencode>
#
#  DESCRIPTION: grep [gene|exon] from gencode with the following columns:
#  			chr
#  			start
#  			end
#  			gene_name
#  			gene_type
#  			level
#  			gene_id
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Thiago Miller (tmiller), tmiller@mochsl.org.br
# ORGANIZATION: Group of Bioinformatics
#      VERSION: 1.0
#      CREATED: 03/20/2017 11:23:28 PM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use autodie;

die usage() unless scalar @ARGV == 2;
my $type = shift;
die usage() unless $type eq 'gene' or $type eq 'exon';
my $gencode = shift;
die usage() . ".gtf required\n" unless $gencode =~ /\.gtf$/;

my $gene_id_rg = qr/gene_id\s+\"(.*?)\"\;/;
my $transcript_id_rg = qr/transcript_id\s+\"(.*?)\"\;/;
my $gene_name_rg = qr/gene_name\s+\"(.*?)\"\;/;
my $gene_type_rg = qr/gene_type\s+\"(.*?)\"\;/;
my $level_rg = qr/level\s+([0-9]*?)\;/;

use constant {
	CHR  => 0,
	TYPE => 2,
	FROM => 3,
	TO   => 4,
};

open my $fh, "<" => $gencode;

while (my $line = <$fh>) {
	chomp $line;
	next if $line =~ /^#/;

	my @fields = split /\t/ => $line;
	next if $fields[TYPE] ne $type;

	my $gene_type = $line =~ $gene_type_rg ? $1 : undef;
	next if $gene_type ne 'protein_coding';
	my $gene_name = $line =~ $gene_name_rg ? $1 : undef;
	my $level = $line =~ $level_rg ? $1 : undef;
	if ($type eq 'gene') {
	my $gene_id = $line =~ $gene_id_rg ? $1 : undef;
	print "$fields[CHR]\t$fields[FROM]\t$fields[TO]\t$gene_name\@\@$gene_type\@\@$level\@\@$gene_id\n";
	 }
	 else {
	my $transcript_id = $line =~ $transcript_id_rg ? $1 : undef;
	print "$fields[CHR]\t$fields[FROM]\t$fields[TO]\t$gene_name\_$gene_type\_$level\_$transcript_id\n";
	 }
}

close $fh;

sub usage {
	return "perl gtf2bed.pl [gene|exon] <path to gencode>\n";
}
