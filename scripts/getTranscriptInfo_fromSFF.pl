#!/usr/bin/perl
use warnings;
use strict;
use Bio::DB::SeqFeature::Store;
use Data::Dumper;

my $sqlite = shift;
die "Please provide a SQLite datafile of a seqfeature db" if !defined $sqlite;

#############################################################################
## make a Bio::DB::SeqFeature::Store object (contains info about  organism in
## SQLite database and all the methods needed to interact with the database
#############################################################################

#Open the sequence database
my $db_obj = Bio::DB::SeqFeature::Store->new(
  -adaptor => 'DBI::SQLite',
  -dsn     => $sqlite
);

## get_feature_by_type can get any type from your original GFF, anything in Col3
my @features_type_gene = sort {$a->ref cmp $b->ref } $db_obj->get_features_by_type('mRNA');
my %exons;
open OUTEXONS, ">transcript_exons_info.txt" or die "Can't open transcript_exons_info.txt";
foreach my $feature ( @features_type_gene ) {
  my $f_name  = $feature->name;
  my $f_start = $feature->start;
  my $f_end   = $feature->end;
  my $ref     = $feature->ref;
  my $strand = $feature->strand;
  $strand = $strand > 0 ? '+' : '-';
  my @features_exons = $db_obj->features(
            -type => 'exon',
            -seq_id => $ref,
            -start  => $f_start,
            -end    => $f_end
        );
  my @exons;
  my @coords;
  foreach my $f (sort {$a->start <=> $b->start} @features_exons) {
    my %attr = $f->attributes;
    my $parent_id = ${$attr{parent_id}}[0];
    next unless $parent_id eq $f_name; 
    my $e_start = $f->start;
    my $e_end   = $f->end;
    push @exons, "$e_start,$e_end";
    push @coords, $e_start,$e_end;
  }
  my @sorted_coords =  sort {$a <=> $b} @coords;
  my $start = shift @coords;
  my $end = pop @coords;
  print OUTEXONS "$f_name\t$ref\t$start,$end\t".join(';',@exons)."\n";
}
