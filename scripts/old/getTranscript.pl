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
my $gene_count = 0;
open OUTEXONS, ">transcript_exons.txt" or die "Can't open transcript_exons.txt";
foreach my $feature ( @features_type_exon ) {
  my $f_name  = $feature->name;
  my $f_start = $feature->start;
  my $f_end   = $feature->end;
  my $ref     = $feature->ref;
  if (!defined $last_ref){## first time through
    $last_ref = $ref;
  }elsif ( $last_ref ne $ref ) {
    print_seq( \%exons, $ref_count );
    %exons = ();
    $last_ref = $ref ;
    $ref_count++;
  }
  my $strand = $feature->strand;
  $strand = $strand > 0 ? '+' : '-';
  my %attr = $feature->attributes;
  my $g_name = ${$attr{parent_id}}[0];
  if (ref $g_name =~ /HASH|ARRAY/){
    my ( $g_name ) = $f_name =~ /(.+)[_.]E/i;
  }
  my $f_seq = $db_obj->fetch_sequence(
    -seq_id => $ref,
    -start  => $f_start,
    -end    => $f_end
  );
  $exons{$g_name}{$f_start}{$f_name}{seq}    = $f_seq;
  $exons{$g_name}{$f_start}{$f_name}{end}    = $f_end;
  $exons{$g_name}{$f_start}{$f_name}{ref}    = $ref;
}
## for last feature
print_seq( \%exons, $ref_count );

sub print_seq {
  my $ref2exons = shift;
  my $ref_count = shift;
  my $ref_name;
  open OUT, ">transcripts_part$ref_count.fasta"
    or die "Can't open transcripts_part$ref_count.fasta\n";
  foreach my $g_name ( sort keys %{$ref2exons} ) {
    my $seq;
    my @exons;
    my @coords;
    foreach my $exon_start (
      sort { $a <=> $b } keys %{${$ref2exons}{ $g_name }})
    {
      foreach my $exon_name ( sort keys %{${$ref2exons}{$g_name}{$exon_start}} ){
        $seq .= ${$ref2exons}{$g_name}{$exon_start}{$exon_name}{seq};
        my $exon_end = ${$ref2exons}{$g_name}{$exon_start}{$exon_name}{end};
        $ref_name = ${$ref2exons}{$g_name}{$exon_start}{$exon_name}{ref};
        push @exons, "$ref_name:$exon_start..$exon_end"; 
        push @coords, $exon_start,$exon_end; 
      }
    }
    my @sorted_coords =  sort {$a <=> $b} @coords;
    my $start = shift @coords;
    my $end = pop @coords;
    print OUT ">$g_name|$ref_name:$start..$end|$gene_count\n$seq\n";
    print OUTEXONS "$g_name\t".join(' ',@exons)."\n";
    $gene_count++;
  }
  close OUT;
}
