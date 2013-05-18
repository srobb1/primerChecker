#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use File::Spec;
use Getopt::Long;
use Tie::File;

if ( !defined @ARGV ) {
  &getHelp();
}
my ( $organism, $primersFile, $type, $format );
GetOptions(
  'o|organism:s'    => \$organism,
  'f|primersFile:s' => \$primersFile,
  't|type:s'        => \$type,
  'i|format:i'      => \$format,
  'h|help'          => \&getHelp,
);

sub getHelp {
  print ' 
usage:
./primerChecker.pl [-o organism][-f primerFile][-t type][-i inputFormat][-h] 

options:
-o | --organism		string		rice,maize
-f | --primersFile      file		file with primer id and sequences 
-t | --type		int		1=Genomic; 2=cDNA
-i | --format		int		1=l line with both primers; 2=2 lines, 1 for each primers
-h | --help				this help message


## format=1 for 1 line for each primer set: CommonName\tseq1\tseq2\n
## foramt=2 for 2 lines for each primer set: Primer1Name\tseq1\nPrimer2Name\tseq1\n

';
  exit 1;
}
$type = $type eq 1 ? 'genomic' : 'cDNA'; 
my $db_dir = "dbs/$organism/$type";    ## genome file
tie my @exon_lengths, 'Tie::File', "dbs/$organism/cDNA/exon_lengths.txt"
  or die "Can't open exon length file: exon_lengths.txt\n";

open PRIMERS, $primersFile       or die "Can't Open $primersFile\n";
open FASTA,   ">$primersFile.fa" or die "Can't Open $primersFile.fa\n";

my %primers;
my ( $id, $p1, $p2 );
print "inputed primer sets
==========================================\n";
if ( $format == 1 ) {
  while ( my $line = <PRIMERS> ) {
    chomp $line;
    next if $line =~ /^#/;
    ( $id, $p1, $p2 ) = split /\t/, $line;
    print "$id,$p1,$p2\n";
    $id =~ s/\s+//g;
    $p1 =~ s/\s+//g;
    $p2 =~ s/\s+//g;

    $primers{$id}{p1}{seq} = $p1;
    $primers{$id}{p2}{seq} = $p2;
    print FASTA ">$id", "_p1\n$p1\n";
    print FASTA ">$id", "_p2\n$p2\n";
  }
}
elsif ( $format == 2 ) {
  while ( my $line = <PRIMERS> ) {
    chomp $line;
    my ( $id1, $id2 );
    ( $id1, $p1 ) = split /\t/, $line;
    $line = <PRIMERS>;
    chomp $line;
    ( $id2, $p2 ) = split /\t/, $line;
    my $copy_1 = $id1;
    my $copy_2 = $id2;

    $copy_1 =~ s/_?(F|R)$//i;
    $copy_2 =~ s/_?(F|R)$//i;
    if ( $copy_1 eq $copy_2 ) {
      $id = $copy_1;
    }
    else {
      $id = join( '-', $id1, $id2 );
    }
    print "$id,$p1,$p2\n";
    $id =~ s/\s+//g;
    $p1 =~ s/\s+//g;
    $p2 =~ s/\s+//g;

    $primers{$id}{p1}{seq} = $p1;
    $primers{$id}{p2}{seq} = $p2;
    print FASTA ">$id", "_p1\n$p1\n";
    print FASTA ">$id", "_p2\n$p2\n";
  }
}

#print "\nGood Blat hits: No mismatches: Full Primer Length
#==========================================\n";
#   print
#     "qName\tqLen\ttName\ttStart\ttEnd\tstrand\tmatches\tmismatches\n";
############ RUN BLAT ################
my @db_files = <$db_dir/*fasta>;
#print "db:$db_dir:@db_files\n";
my $count    = 0;
`rm $primersFile.blatout` if -e "$primersFile.blatout";
foreach my $db (@db_files) {

  #   print "blat against $db\n";
`blat -noHead -tileSize=7 -minScore=10 $db $primersFile.fa $primersFile.$count.blatout`;
`cat $primersFile.$count.blatout >> $primersFile.blatout`;
`rm -f $primersFile.*.blatout`;
  $count++;
}

#my @file_path = split '/', $blat_file;
#my $file_name = pop @file_path;
##blat parser
open INBLAT, "$primersFile.blatout",
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
  if ( $type eq 'genomic' ) {
    ## foreach id: foreach target: 1. are there 2 pairs? 2. is there 1 pair foreach strand
    $primers{$id}{hit}{$tName}{$pair}{$tStart}{strand}     = $strand;
    $primers{$id}{hit}{$tName}{$pair}{$tStart}{tEnd}       = $tEnd;
    $primers{$id}{hit}{$tName}{$pair}{$tStart}{matches}    = $matches;
    $primers{$id}{hit}{$tName}{$pair}{$tStart}{mismatches} = $mismatches;
  }
  else {    ##cDNA
    ## AC148152.3_FG001:3|AC148152.3_FG001.exon4(-)|2:231360108..231360555|86071
    my ( $exon_id, $exon_name, $loc, $exon_count ) = split /\|/, $tName;
    my ( $gene, $gene_exon_count ) = split /:/, $exon_id;
    $primers{$id}{hit}{$gene}{$tName}{$pair}{$tStart}{strand}     = $strand;
    $primers{$id}{hit}{$gene}{$tName}{$pair}{$tStart}{tEnd}       = $tEnd;
    $primers{$id}{hit}{$gene}{$tName}{$pair}{$tStart}{matches}    = $matches;
    $primers{$id}{hit}{$gene}{$tName}{$pair}{$tStart}{mismatches} = $mismatches;
    $primers{$id}{hit}{$gene}{$tName}{$pair}{$tStart}{tLen}       = $tLen;
  }

#   print
#     "$qName\t$qLen\t$tName\t$tStart\t$tEnd\t$strand\t$matches\t$mismatches\n";

}

#print Dumper \%primers;
print "\nGood Primer Pairs with product Size
==========================================\n";
my %results;

if ( $type eq 'genomic' ) {
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
          my $strand  = $primers{$id}{hit}{$tName}{$pair}{$tStart}{strand};
          my $tEnd    = $primers{$id}{hit}{$tName}{$pair}{$tStart}{tEnd};
          my $matches = $primers{$id}{hit}{$tName}{$pair}{$tStart}{matches};
          my $qLen    = $primers{$id}{qLen}{$pair};
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
        foreach my $p1_range ( sort { ${$a}[0] <=> ${$b}[0] }
          @{ $hits{$strand}{$p1} } )
        {
          my ( $p1_s, $p1_e ) = sort { $a <=> $b } @{$p1_range};
          my $other_strand = $strand eq '+' ? '-' : '+';
          next if !exists $hits{$other_strand}{$p2};
          foreach my $p2_range ( sort { ${$a}[0] <=> ${$b}[0] }
            @{ $hits{$other_strand}{$p2} } )
          {
            my ( $p2_s, $p2_e ) = sort { $a <=> $b } @{$p2_range};
            my @sorted_starts = sort { $a <=> $b } ( $p1_s, $p2_s );
            my @sorted_ends   = sort { $b <=> $a } ( $p1_e, $p2_e );
            my $product_size = $sorted_ends[0] - $sorted_starts[0] + 1;

            #print "$id\t$product_size\t$tName\t@$p1_range\t@$p2_range\n";
            $results{$id}{$product_size} =
              "$id\t$product_size\t$tName\t@$p1_range\t@$p2_range";
          }
        }
      }
    }
  }

}
else {
print "id\tproduct_size\texon_sum\tgene\tP1(exon_range;genomic_range)\tP2(exon_range;genomic_range)\n";

  foreach my $id ( sort keys %primers ) {
    foreach my $gene ( sort keys %{ $primers{$id}{hit} } ) {
      my $pairs_per_gene = scalar keys %{ $primers{$id}{hit} };
      my %hits;
      foreach my $exon ( sort keys %{ $primers{$id}{hit}{$gene} } ) {
        my $pairs_per_exon = scalar keys %{ $primers{$id}{hit}{$gene} };
        next if $pairs_per_gene < 2;
        #next if $pairs_per_exon > 1;
        foreach my $pair ( sort keys %{ $primers{$id}{hit}{$gene}{$exon} } ) {
          foreach my $tStart (
            sort {
              $primers{$id}{hit}{$gene}{$exon}{$pair}{$a} <=> $primers{$id}{hit}
                {$gene}{$exon}{$pair}{$b}
            } keys %{ $primers{$id}{hit}{$gene}{$exon}{$pair} }
            )
          {
            my $strand =
              $primers{$id}{hit}{$gene}{$exon}{$pair}{$tStart}{strand};
            my $tEnd = $primers{$id}{hit}{$gene}{$exon}{$pair}{$tStart}{tEnd};
            my $matches =
              $primers{$id}{hit}{$gene}{$exon}{$pair}{$tStart}{matches};
            my $qLen = $primers{$id}{qLen}{$pair};
            my $tLen = $primers{$id}{hit}{$gene}{$exon}{$pair}{$tStart}{tLen};
            push @{ $hits{$strand}{$pair} }, [ $tStart, $tEnd, $tLen, $exon ];
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
          foreach my $p1_range ( sort { ${$a}[0] <=> ${$b}[0] }
            @{ $hits{$strand}{$p1} } )
          {
            my $exon_p1_long = pop @{$p1_range};
            my $exon_p1_len  = pop @{$p1_range};
            my ( $p1_s, $p1_e ) = sort { $a <=> $b } @{$p1_range};
            my $other_strand = $strand eq '+' ? '-' : '+';
            next if !exists $hits{$other_strand}{$p2};
            foreach my $p2_range ( sort { ${$a}[0] <=> ${$b}[0] }
              @{ $hits{$other_strand}{$p2} } )
            {
              my $exon_p2_long = pop @{$p2_range};
              my $exon_p2_len  = pop @{$p2_range};
              my ( $p2_s, $p2_e ) = sort { $a <=> $b } @{$p2_range};
              my @sorted_starts = sort { $a <=> $b } ( $p1_s, $p2_s );
              ## exons are numbered left to right in genome ## 5'-e1-e2-e3-3' 3'-e1-e2-e3-5'
              my $exon_p1_product_size;
              my $exon_p2_product_size;
              ## GRMZM2G152059:1|GRMZM2G152059_E01(-)|7:302416..304114|12345
              my ( $exon_id_p1, $exon_name_p1, $loc_p1, $exon_count_p1 ) =
                split /\|/, $exon_p1_long;
              my ( $gene_p1, $exon_p1 ) = split /:/, $exon_id_p1;
              my ( $exon_id_p2, $exon_name_p2, $loc_p2, $exon_count_p2 ) =
                split /\|/, $exon_p2_long;
              my ( $gene_p2, $exon_p2 ) = split /:/, $exon_id_p2;
              my $inner_exon_size = 0;
              if ( $exon_p1 < $exon_p2 ) {    ## plus strand
                $exon_p1_product_size = $exon_p1_len - $p1_s + 1;
                $exon_p2_product_size = $p2_e;
                if ( ( $exon_p2 - $exon_p1 ) > 1 ) {
                  for ( my $i = $exon_count_p1 + 1 ;
                    $i < $exon_count_p2 ; $i++ )
                  {
                    my ( $exon_num, $exon_name, $exon_len ) = split /\t/,
                      $exon_lengths[$i];
                    $inner_exon_size += $exon_len;
                  }
                }
              }
              else {
                $exon_p1_product_size = $p1_e;
                $exon_p2_product_size = $exon_p2_len - $p2_s + 1;
              }
              my $product_size = $exon_p1_product_size + $exon_p2_product_size + $inner_exon_size;
              my $exon_sum = $exon_p1_product_size . '+' . $inner_exon_size . '+' . $exon_p2_product_size; 

              $results{$id}{$product_size} =
"$id\t$product_size\t$exon_sum\t$gene\tP1(exon_$exon_p1:".join('..',@$p1_range).";$loc_p1)\tP2(exon_$exon_p2:".join('..',@$p2_range).";$loc_p2)\n";
            }
          }
        }
      }
    }
  }
}

foreach my $id ( sort keys %results ) {
  foreach my $product_size ( sort { $a <=> $b } keys %{ $results{$id} } ) {
    print $results{$id}{$product_size}, "\n";
  }
}

