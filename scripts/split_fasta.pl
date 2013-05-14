#!/usr/bin/perl -w
use strict;

my $file = shift;
if (!defined $file){
  die "Usage:

	./split_fasta.pl FASTA_FILE

FASTA file will be split into many smaller files with one sequence per file


"
}

open FASTA, $file or die "Can't open $file\n";
my ($base) = $file =~ /(.+)\.fa|fasta%/;
my $count = 0;

my $id;
my $seq; 
my $first_line = <FASTA>;
chomp $first_line;
if ($first_line =~ /^>(\S+)/){
  $id = $1;
  while (my $seq_line = <FASTA>){
    chomp $seq_line;
    if ($seq_line =~ /^>(\S+)/){
      printSeq(\$seq);
      $seq='';
      ($id) = $seq_line =~ /^>(\S+)/;
      $count++;
    }else{
      $seq.=$seq_line;
    }
  }
}
## for last seq
printSeq(\$seq);

sub printSeq {
  my $ref2seq = shift; 
  my $outFASTA = "$base.part_$count.fasta";
  open OUT, ">$outFASTA" or die "Can't open $outFASTA\n";
  ${$ref2seq} =~ s/(.{80})/$1\n/g;
  print OUT ">$id\n${$ref2seq}\n";
  close OUT;
}
