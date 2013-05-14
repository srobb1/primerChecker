#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

my $file = shift;
if (!defined $file){
  print "Usage:

	./fasta_splitter.pl FASTA_FILE

FASTA file will be split into many smaller files with one sequence per file


"
}

my $seqIO_obj = Bio::SeqIO->new (-file=> $file, -format=>'fasta');
my $count = 0;
while (my $seqObj = $seqIO_obj->next_seq){
  my ($base) = $file =~ /(.+)\.fa|fasta%/;
  my $outFile = $base.part_$count.fasta;
  my $out_seqIO_obj = Bio::SeqIO->new (-file=> ">$outFile", -format=>'fasta');
  $out_seqIO_obj->write_seq($seqObj);
  $count++;
} 
