#!/usr/bin/perl
use warnings;
use strict;
use Bio::DB::SeqFeature::Store;
use Data::Dumper;

my $GFF = shift;
die "Please provide a GFF3 file with mRNA and exon features" if !defined $GFF;

#############################################################################
## make a Bio::DB::SeqFeature::Store object (contains info about  organism in
## SQLite database and all the methods needed to interact with the database
#############################################################################

#Open the sequence database
open GFF, $GFF or die "Can't open $GFF\n";
open OUTEXONS, ">transcript_exons_info.txt" or die "Can't open transcript_exons_info.txt";
my %genes;
while (my $line = <GFF>){
  last if $line =~ /FASTA/;
  next if $line =~ /^#/;
  chomp $line;
  my ($ref,$source,$type,$start,$end,$five,$strand,$seven,$nine) = split /\t/ , $line;
  #print "($ref,$source,$type,$start,$end,$five,$strand,$seven,$eight,$nine)\n";
  if ($type eq 'mRNA'){
    my ($id) = $nine =~ /ID=(.+?);/;
    my ($name) = $nine =~ /Name=(.+?);/;
    $genes{$ref}{$id}{name}=$name;
    $genes{$ref}{$id}{start}=$start;
    $genes{$ref}{$id}{end}=$end;
    #$genes{$ref}{$id}{strand}=$strand;
  }elsif($type eq 'exon'){
    my ($parent_id) = $nine =~ /Parent=(.+?);/;
    push @{$genes{$ref}{$parent_id}{exons}} , [$start,$end];
  }
}
foreach my $ref (sort keys %genes){
  foreach my $gene (sort keys %{$genes{$ref}}){
    my $name =  $genes{$ref}{$gene}{name};
    my $start =  $genes{$ref}{$gene}{start};
    my $end =  $genes{$ref}{$gene}{end};
    my $exons =  $genes{$ref}{$gene}{exons};
    my @exons;
    foreach my $exon (sort {${$a}[0]<=>${$b}[0]} @{$exons}){
      push @exons, "$$exon[0],$$exon[1]";
    }
    print OUTEXONS "$name\t$ref\t$start,$end\t",join(';',@exons),"\n";
  }
} 
