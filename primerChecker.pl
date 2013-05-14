#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use File::Spec;

my $db_dir      = shift;    ## genome file
my $primersFile = shift;
my $format      = shift;
## format=1 for 1 line for each primer set: CommonName\tseq1\tseq2\n
## foramt=2 for 2 lines for each primer set: Primer1Name\tseq1\nPrimer2Name\tseq1\n

open PRIMERS, $primersFile       or die "Can't Open $primersFile\n";
open FASTA,   ">$primersFile.fa" or die "Can't Open $primersFile.fa\n";

my %primers;
my ( $id, $p1, $p2 );
print "inputed primer sets
==========================================\n";
if ( $format == 1 ) {
  while ( my $line = <PRIMERS> ) {
    chomp $line;
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


print "\nGood Blat hits: No mismatches: Full Primer Length
==========================================\n";
   print
     "qName\tqLen\ttName\ttStart\ttEnd\tstrand\tmatches\tmismatches\n";
############ RUN BLAT ################
my @db_files = <$db_dir/*fasta>;
my $count = 0;
`rm $primersFile.blatout` if -e "$primersFile.blatout";
foreach my $db (@db_files){
#   print "blat against $db\n";
  `blat -noHead -tileSize=7 -minScore=10 $db $primersFile.fa $primersFile.$count.blatout` ;
   `cat $primersFile.$count.blatout >> $primersFile.blatout`;
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
  $id =~ s/_(p\d)$//;
  my ($pair) = $qName =~ /_(p\d)$/;
  ## foreach id: foreach target: 1. are there 2 pairs? 2. is there 1 pair foreach strand
  $primers{$id}{hit}{$tName}{$pair}{$tStart}{strand}     = $strand;
  $primers{$id}{hit}{$tName}{$pair}{$tStart}{tEnd}       = $tEnd;
  $primers{$id}{hit}{$tName}{$pair}{$tStart}{matches}    = $matches;
  $primers{$id}{hit}{$tName}{$pair}{$tStart}{mismatches} = $mismatches;
  $primers{$id}{qLen}{$pair}                             = $qLen;

   print
     "$qName\t$qLen\t$tName\t$tStart\t$tEnd\t$strand\t$matches\t$mismatches\n";

}
#print Dumper \%primers;
print "\nPrimer Pairs with product Size
==========================================\n";
print "id\tproduct_size\ttName\tp1_range\tp2_range\n";
my %results;
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
    my %pairs; ## pairs
    foreach my $strand ( %hits ) {
      foreach my $p ( keys %{ $hits{$strand} }){
        $pairs{$p}++;
      }
    }
    my ($p1, $p2) = sort keys %pairs;
    next if keys %pairs != 2;;
    ## next if we dont have one pair on - and one on +
    next if ( keys %{ $hits{'-'} } < 1 or keys %{ $hits{'+'} } < 1 );
    foreach my $strand ( %hits ) {
      next if !exists $hits{$strand}{$p1} ;
      ## sort by smallest to biggest start
      foreach
        my $p1_range ( sort { ${$a}[0] <=> ${$b}[0] } @{ $hits{$strand}{$p1} } )
      {
        my ( $p1_s, $p1_e ) = sort { $a <=> $b } @{$p1_range};
        my $other_strand = $strand eq '+' ? '-' : '+';
        next if !exists $hits{$other_strand}{$p2} ;
        foreach my $p2_range ( sort { ${$a}[0] <=> ${$b}[0] }
          @{ $hits{$other_strand}{$p2} } )
        {
          my ( $p2_s, $p2_e ) = sort { $a <=> $b } @{$p2_range};
          my @sorted_starts = sort { $a <=> $b } ( $p1_s, $p2_s );
          my @sorted_ends   = sort { $b <=> $a } ( $p1_e, $p2_e );
          my $product_size = $sorted_ends[0] - $sorted_starts[0] + 1;
          #print "$id\t$product_size\t$tName\t@$p1_range\t@$p2_range\n";
          $results{$id}{$product_size} = "$id\t$product_size\t$tName\t@$p1_range\t@$p2_range";
        }
      }
    }
  }
}
foreach my $id (sort keys %results){
  foreach my $product_size (sort { $a <=> $b } keys %{$results{$id}}){
    print $results{$id}{$product_size},"\n";
  }
}
