#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use File::Spec;
use Tie::File;

use Getopt::Std;
use CGI':standard';
use IO::String;


if (!param){

print header;
print
    start_html('Run Primer Checker'),
    h1('Run PrimerChecker'),

    start_multipart_form,

#textbox for seq
    "Input one or more Primer sets in the following format:
NAME<tab>forwardPrimer<tab>ReversePrimer 
NAME<tab>forwardPrimer<tab>ReversePrimer", 

        textarea(-name=>'primerSets',-rows=>15,-cols=>90),

    br,
        #end textbox
    br,
    "<u>Organism</u>: ( maize or rice only)", "&nbsp;",

         radio_group(-name=>'organism', -value=>{maize=>'Zea mays(corn)',rice=>'Oryza sativa ssp nipponbare(rice)'},-linebreak=>'true',-default=>'maize'),

    br,
    br,
        "<u>Primers for Genomic or cDNA</u> ",
        "&nbsp;",
        "&nbsp;",
         radio_group(-name=>'type', -value=>{genomic=>'Genomic',cDNA=>'cDNA'},-linebreak=>'true',-default=>'genomic'),
        "&nbsp;",
    submit(-name=>'primer3', -value=>'Check Primers'),
    end_form,
    end_html;

}
if ( param ){
    print header;
    print start_html('Your Primer Results');

    my $primerSets = param('primerSets');
    my $organism = param('organism');
    my $type = param('type'); 
    my @primerSets = split /\n/ , $primerSets;
my $db_dir = "/rhome/project/robb/maize_mites/FALL2013/primers/dbs/$organism";    ## genome file
my %exons;
##transcript_exons_info.txt
open EXONS, "$db_dir/transcript_exons_info.txt"
  or die "Can't open exon info file: transcript_exons_info.txt\n";
while ( my $line = <EXONS> ) {
  chomp $line;
  my ( $t_name, $ref, $t_coord, $e_coords ) = split /\t/, $line;
  my @e_coords = split ';', $e_coords;
  my ($first_digit) = $t_coord =~ /^(\d)/;
  my ( $t_s, $t_e ) = split ',', $t_coord;
  if ( $ref =~ /^\d/ ) {
    $ref = "chr$ref";
  }
  $exons{$ref}{$first_digit}{$t_s}{$t_name}{end} = $t_e;
  foreach my $e_coord (@e_coords) {
    my ( $e_s, $e_e ) = split ',', $e_coord;
    push @{ $exons{$ref}{$first_digit}{$t_s}{$t_name}{exons} }, [ $e_s, $e_e ];
  }
}
 my $time = time();
open GBROWSE, ">$time.forBrowser.txt"
  or die "Can't Open $time.forBrowser\n";
print GBROWSE "[primers]
glyph = segments
feature = primers
\n\n";
open FASTA,   ">$time.fa" or die "Can't Open $time.fa\n";

my %primers;
my ( $id, $p1, $p2 );
print "inputed primer sets
==========================================\n";
  foreach my $line( @primerSets ) {
    next if $line =~ /^#/;
    next if $line =~ /^\s*$/;
    ( $id, $p1, $p2 ) = split /\t/, $line;
    $id =~ s/\s+//g;
    $p1 =~ s/\s+//g;
    $p2 =~ s/\s+//g;
    print "$id,$p1,$p2\n";

    $primers{$id}{p1}{seq} = $p1;
    $primers{$id}{p2}{seq} = $p2;
    $primers{$id}{p1}{good} = 0;
    $primers{$id}{p2}{good} = 0 ;
    print FASTA ">$id", "_p1\n$p1\n";
    print FASTA ">$id", "_p2\n$p2\n";
  }

#print "\nGood Blat hits: No mismatches: Full Primer Length
#==========================================\n";
#   print
#     "qName\tqLen\ttName\ttStart\ttEnd\tstrand\tmatches\tmismatches\n";
############ RUN BLAT ################
my @db_files = <$db_dir/*fasta>;

#print "db:$db_dir:@db_files\n";
my $count = 0;
`rm $time.blatout` if -e "$time.blatout";
foreach my $db (@db_files) {

  #   print "blat against $db\n";
`blat -noHead -tileSize=7 -minScore=10 $db $time.fa $time.$count.blatout`;
  `cat $time.$count.blatout >> $time.blatout`;
  `rm -f $time.*.blatout`;
  $count++;
}

#my @file_path = split '/', $blat_file;
#my $file_name = pop @file_path;
##blat parser
open INBLAT, "$time.blatout",
  or die "Please provide a blat output file\n";

while ( my $line = <INBLAT> ) {
  my @line        = split /\t/, $line;
  my $matches     = $line[0];
  my $mismatches  = $line[1];
  my $qBaseInsert = $line[5];
  my $tBaseInsert = $line[7];
  my $strand      = $line[8];
  my $qName       = $line[9];
  my $qLen        = $line[10];
  my $tName       = $line[13];
  my $tLen        = $line[14];
  my $tStart      = $line[15] + 1;
  my $tEnd        = $line[16];
  my $aln_bp      = $matches + $qBaseInsert + $mismatches;

  ## throw out if alignment is too small
  next if ( $matches + $mismatches ) < $qLen;
  ## throw out if there are too many MM
  next if $mismatches > 0;
  my $id = $qName;
  $id =~ s/_(p\d+)//;
  my ($pair) = $qName =~ /_(p\d)$/;
  $primers{$id}{qLen}{$pair} = $qLen;
  ## foreach id: foreach target: 1. are there 2 pairs? 2. is there 1 pair foreach strand
  $primers{$id}{hit}{$tName}{$pair}{$tStart}{strand}     = $strand;
  $primers{$id}{hit}{$tName}{$pair}{$tStart}{tEnd}       = $tEnd;
  $primers{$id}{hit}{$tName}{$pair}{$tStart}{matches}    = $matches;
  $primers{$id}{hit}{$tName}{$pair}{$tStart}{mismatches} = $mismatches;
}

print "\nGood Primer Pairs with product Size
==========================================\n";
my %results;

print "id\tproduct_size\ttName\tp1_range\tp2_range\n";
foreach my $id ( sort keys %primers ) {
  foreach my $tName ( sort keys %{ $primers{$id}{hit} } ) {
    my $pairs_per_target = scalar keys %{ $primers{$id}{hit}{$tName} };
    next if $pairs_per_target < 2;
    my %hits;
    foreach my $pair ( sort keys %{ $primers{$id}{hit}{$tName} } ) {
      foreach my $tStart (
        sort {
          $primers{$id}{hit}{$tName}{$pair}{$a} <=> $primers{$id}{hit}{$tName}
            {$pair}{$b}
        } keys %{ $primers{$id}{hit}{$tName}{$pair} }
        )
      {
        my $strand = $primers{$id}{hit}{$tName}{$pair}{$tStart}{strand};
        my $tEnd   = $primers{$id}{hit}{$tName}{$pair}{$tStart}{tEnd};
        push @{ $hits{$strand}{$pair} }, [ $tStart, $tEnd ];
      }
    }
    ## next if we dont have hits on both strands
    next if keys %hits < 2;
    my %pairs;    ## pairs
    foreach my $strand (%hits) {
      foreach my $p ( keys %{ $hits{$strand} } ) {
        $pairs{$p}++;
      }
    }
    my ( $p1, $p2 ) = sort keys %pairs;
    next if keys %pairs != 2;
    ## next if we dont have one pair on - and one on +
    next if ( keys %{ $hits{'-'} } < 1 or keys %{ $hits{'+'} } < 1 );
    foreach my $strand (%hits) {
      next if !exists $hits{$strand}{$p1};
      ## sort by smallest to biggest start
      foreach
        my $p1_range ( sort { ${$a}[0] <=> ${$b}[0] } @{ $hits{$strand}{$p1} } )
      {
        my $p1_strand = $strand;
        my $p2_strand;
        my ( $p1_s, $p1_e ) = sort { $a <=> $b } @{$p1_range};
        my $other_strand = $strand eq '+' ? '-' : '+';
        next if !exists $hits{$other_strand}{$p2};
        foreach my $p2_range ( sort { ${$a}[0] <=> ${$b}[0] }
          @{ $hits{$other_strand}{$p2} } )
        {
          my $p2_strand = $other_strand;
          my ( $p2_s, $p2_e ) = sort { $a <=> $b } @{$p2_range};
          my @sorted_starts = sort { $a <=> $b } ( $p1_s, $p2_s );
          my @sorted_ends   = sort { $b <=> $a } ( $p1_e, $p2_e );
          my $hit_start     = $sorted_starts[0];
          my $hit_end       = $sorted_ends[0];
          my $product_size = $hit_end - $hit_start + 1;
          my $p1_range_str = join( '..', @$p1_range );
          my $p2_range_str = join( '..', @$p2_range );

          $primers{$id}{p1}{good} = 1;
          $primers{$id}{p2}{good} = 1 ;
          if ( $p1_s > $p2_s ) {
            $results{$id}{$product_size}{gbrowse} =
              "primers\t$id($product_size)\t$tName:$p2_range_str,$p1_range_str";
          }
          else {
            $results{$id}{$product_size}{gbrowse} =
              "primers\t$id($product_size)\t$tName:$p1_range_str,$p2_range_str";
          }
          $results{$id}{$product_size}{info} =
"$id\t$product_size\tp1|$p1_strand;p2|$p2_strand\t$tName:$p1_range_str\t$tName:$p2_range_str";
          if ( $type eq 'cDNA' and $product_size < 10000 ) {
            my ($hit_first_digit) = $hit_start =~ /^(\d)/;
            my $last_start        = 0;
            my $done              = 0;
            foreach my $gene_start (
              sort { $a <=> $b }
              keys %{ $exons{$tName}{$hit_first_digit} }
              )
            {
              if ( $hit_start <= $gene_start ) {
                foreach my $gene_name (
                  sort keys %{ $exons{$tName}{$hit_first_digit}{$last_start} } )
                {
                  my $gene_end =
                    $exons{$tName}{$hit_first_digit}{$last_start}{$gene_name}
                    {end};
                  if ( $hit_end <= $gene_end ) {
                    $done = 1;
                    my @exon_ranges =
                      @{ $exons{$tName}{$hit_first_digit}{$last_start}
                        {$gene_name}{exons} };
                    my %cDNA;
                    my $exon_count = 0;
                    ## pre-ordered smallest to biggest in ref
                    foreach my $exon_range (@exon_ranges) {
                      my $p1_len = length $primers{$id}{p1}{seq};
                      my $p2_len = length $primers{$id}{p2}{seq};
                      my $p1_info =
                        getExonInfo( $p1_len, $exon_range, $p1_range );
                      my $p2_info =
                        getExonInfo( $p2_len, $exon_range, $p2_range );
                      my $exon_size =
                        range_get_end($exon_range) -
                        range_get_start($exon_range) + 1;

                      $cDNA{$exon_count}{size}  = $exon_size;
                      $cDNA{$exon_count}{range} = $exon_range;
                      $cDNA{$exon_count}{p1}    = $p1_info;
                      $cDNA{$exon_count}{p2}    = $p2_info;
                      $exon_count++;
                    }
                    my $exon_p1_product_size = 0;
                    my $exon_p2_product_size = 0;
                    my $p1_exon;
                    my $p2_exon;
                    my $exon_score;
                    foreach my $exon ( sort { $a <=> $b } keys %cDNA ) {
                      next
                        if !exists $cDNA{$exon}{p1}
                          and !exists $cDNA{$exon}{p2};
                      my $exon_size  = $cDNA{$exon}{size};
                      my $exon_range = $cDNA{$exon}{range};
                      if ( $cDNA{$exon}{p1}) {
                        $p1_exon = $exon;
                        $exon_score+=$cDNA{$exon}{p1};
                        if ( $p1_s < $p2_s ) {
                          $exon_p1_product_size =
                            range_get_end($exon_range) - $p1_s + 1;
                        }
                        else {
                          $exon_p1_product_size =
                            $p1_e - range_get_start($exon_range) + 1;
                        }
                      }
                      elsif ( $cDNA{$exon}{p2}) {
                        $exon_score+=$cDNA{$exon}{p2};
                        $p2_exon = $exon;
                        if ( $p1_s < $p2_s ) {
                          $exon_p2_product_size =
                            $p2_e - range_get_start($exon_range) + 1;
                        }
                        else {
                          $exon_p2_product_size =
                            range_get_end($exon_range) - $p2_s + 1;
                        }
                      }
                    }
                    my $inner_exon_size = 0;
                    my @sorted_exons    = ( $p1_exon, $p2_exon );
                    my $smallest        = 'p1';
                    next if !defined $p1_exon or !defined $p2_exon;
                    if ( $p1_exon > $p2_exon ) {
                      @sorted_exons = ( $p2_exon, $p1_exon );
                      $smallest = 'p2';
                    }
                    if ( ( $sorted_exons[1] - $sorted_exons[0] ) > 1 ) {
                      for (
                        my $i = $sorted_exons[0] + 1 ;
                        $i < $sorted_exons[1] ;
                        $i++
                        )
                      {
                        $inner_exon_size += $cDNA{$i}{size};
                      }
                    }

                    my $cDNA_product_size =
                      $exon_p1_product_size +
                      $inner_exon_size +
                      $exon_p2_product_size;
                    $results{$id}{$product_size}{info} =
"$id\t$product_size|$cDNA_product_size\t$gene_name;$exon_score;p1|$p1_strand|$cDNA{$p1_exon}{p1};p2|$p2_strand|$cDNA{$p2_exon}{p2}\t$tName:$p1_range_str\t$tName:$p2_range_str";
                    if ( $p1_s > $p2_s ) {
                      $results{$id}{$product_size}{gbrowse} =
"primers\t$id($product_size|$cDNA_product_size)\t$tName:$p2_range_str,$p1_range_str";
                    }
                    else {
                      $results{$id}{$product_size}{gbrowse} =
"primers\t$id($product_size|$cDNA_product_size)\t$tName:$p1_range_str,$p2_range_str";
                    }
                  }
                }
              }
              $last_start = $gene_start;
              last if $done;
            }
          }
        }
      }
    }
  }
}
foreach my $id ( sort keys %results ) {
  foreach my $product_size ( sort { $a <=> $b } keys %{ $results{$id} } ) {
    print $results{$id}{$product_size}{info}, "\n";
    print GBROWSE $results{$id}{$product_size}{gbrowse}, "\n";
  }
}
foreach my $id (sort keys %primers){
  my $good = $primers{$id}{p1}{good};
  if (!$good){
    print "$id\tnoHits\n";
  }
}


}
#####SUBROUTINES########
sub getExonInfo {
  my $p_len      = shift;
  my $exon_range = shift;
  my $p_range    = shift;

  my $overlap = range_overlap( $exon_range, $p_range );
  if ( $overlap == $p_len ) {
    ## primer completely within 1 exon
    return 1;
  }
  elsif ( $overlap > 0 ) {
    ## primer overlap junction, in more than 1 exon
    return 0.5;
  }
  else {
    ## not found in this exon
    return 0;
  }
}

sub range_create {
## range needs to be s<=e
## range is in 1 base notation
## takes 2 numbers and returns an anonymous array
  my @range = ( $_[0], $_[1] );
  if ( $range[0] !~ /^\d+$/ or $range[1] !~ /^\d+$/ ) {
    die "Ranges be numbers but \n";
  }
  elsif ( $range[0] == 0 or $range[1] == 0 ) {
    die "Ranges are in 1 base notation not 0, numbers must be > 0\n";
  }
  @range = sort { $a <=> $b } @range;
  return \@range;
}

sub range_check {
  my $range = shift;
  if ( ref $range !~ /ARRAY/ ) {
    die "range functions need to be given an array ref\n";
  }
  if ( scalar @$range != 2 ) {
    die "range functions need to have 2 values\n";
  }
  return $range;
}

sub range_get_start {
  my $r = shift;
  return $$r[0];
}

sub range_get_end {
  my $r = shift;
  return $$r[1];
}

sub range_overlap {
  ## returns the size of the overlap
  my ( $r1, $r2 ) = @_;
  $r1 = range_check($r1);
  $r2 = range_check($r2);

  my ( $s,  $e )  = @$r1;
  my ( $s2, $e2 ) = @$r2;
  ## r     |-----------|
  ## r2 |-->
  ## r2    |
  if ( $s2 <= $s and $e2 >= $s ) {
    if ( $e2 <= $e ) {
      return ( $e2 - $s + 1 );
    }
    elsif ( $e2 > $e ) {
      return ( $e - $s + 1 );
    }
    else {
      die "range: error1\n";
    }
  }
  ## r  |---------------|
  ## r2 |--->
  ## r2    |--->
  elsif ( $s2 >= $s and $s2 <= $e and $e2 >= $s ) {
    if ( $e2 <= $e ) {
      return ( $e2 - $s2 + 1 );
    }
    elsif ( $e2 > $e ) {
      return ( $e - $s2 + 1 );
    }
    else {
      die "range: error2\n";
    }
  }
  else {
    return 0;
  }
}