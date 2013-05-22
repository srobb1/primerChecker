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
my @features_type_exon = sort {$a->ref cmp $b->ref } $db_obj->get_features_by_type('exon');
my %exons;
my $last_ref;
my $ref_count = 0;
my $exon_count = 0;
open OUTEXONS, ">exon_lengths.txt" or die "Can't open exon_lengths.txt";
foreach my $feature ( @features_type_exon ) {
  my $f_name  = $feature->name;
  my $f_start = $feature->start;
  my $f_end   = $feature->end;
  my $ref     = $feature->ref;
  if (!defined $last_ref){## first time through
    $last_ref = $ref;
  }elsif ( $last_ref ne $ref ) {
    print_seq( \%exons, $ref_count );
    %exons = {};
    $last_ref = $ref ;
    $ref_count++;
  }
  my $strand = $feature->strand;
  $strand = $strand > 0 ? '+' : '-';
  my %attr = $feature->attributes;
  my ( $g_name, $exon_id ) = $f_name =~ /(.+)[_.](E.+)$/i;
  my $f_seq = $db_obj->fetch_sequence(
    -seq_id => $ref,
    -start  => $f_start,
    -end    => $f_end
  );
  $exons{$g_name}{$f_start}{$f_name}{len}    = length $f_seq;
  $exons{$g_name}{$f_start}{$f_name}{header} = "$f_name($strand)|$ref:$f_start..$f_end";
  $exons{$g_name}{$f_start}{$f_name}{seq}    = $f_seq;
}
## for last feature
print_seq( \%exons, $ref_count );

sub print_seq {
  my $ref2exons = shift;
  my $ref_count = shift;
  open OUT, ">exons_part$ref_count.fasta"
    or die "Can't open exons_part$ref_count.fasta\n";
  foreach my $g_name ( sort keys %{$ref2exons} ) {
    my $count = 0;
    foreach my $exon_start (
      sort { $a <=> $b } keys %{${$ref2exons}{ $g_name }})
    {
      foreach my $exon_name ( sort keys %{${$ref2exons}{$g_name}{$exon_start}} ){
        $count++;
        my $header = ${$ref2exons}{$g_name}{$exon_start}{$exon_name}{header};
        my $seq = ${$ref2exons}{$g_name}{$exon_start}{$exon_name}{seq};
        my $len = ${$ref2exons}{$g_name}{$exon_start}{$exon_name}{len};
        print OUT ">$g_name:$count|$header|$exon_count\n$seq\n";
        print OUTEXONS "$exon_count\t$g_name:$count\t$len\n";
        $exon_count++;
      }
    }
  }
  close OUT;
}
