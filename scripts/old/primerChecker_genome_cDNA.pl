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
tie my @exons, 'Tie::File', "dbs/$organism/cDNA/transcript_exons.txt"
  or die "Can't open exon length file: transcript_exons.txt\n";

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
    ## AC148152.3_FG001|2:231360108..231360555|12345
    my ( $transcript_id, $loc , $tie_array_id) = split /\|/, $tName;
    $primers{$id}{hit}{$tName}{$pair}{$tStart}{strand}     = $strand;
    $primers{$id}{hit}{$tName}{$pair}{$tStart}{tEnd}       = $tEnd;
    $primers{$id}{hit}{$tName}{$pair}{$tStart}{matches}    = $matches;
    $primers{$id}{hit}{$tName}{$pair}{$tStart}{mismatches} = $mismatches;
    $primers{$id}{hit}{$tName}{$pair}{$tStart}{tLen}       = $tLen;
    $primers{$id}{hit}{$tName}{$pair}{$tStart}{loc}       = $loc;
    $primers{$id}{hit}{$tName}{$pair}{$tStart}{tie_array_id}       = $tie_array_id;
  }

#   print
#     "$qName\t$qLen\t$tName\t$tStart\t$tEnd\t$strand\t$matches\t$mismatches\n";

}

#print Dumper \%primers;
print "\nGood Primer Pairs with product Size
==========================================\n";
my %results;

#if ( $type eq 'genomic' ) {
print "id\tproduct_size\ttName\tp1_range\tp2_range\texon_info\n";

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
          my $tie_array_id = $primers{$id}{hit}{$tName}{$pair}{$tStart}{tie_array_id};
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
            my $exon_info_sum=0;
            my @exon_info;
            if ($type eq 'cDNA'){
              my ($t_id,$loc,$tie_array_id) = split /\|/ ,$tName ;
              $tName = $t_id;
              my ($genomic_ref, $genomic_s, $genomic_e) = $loc =~ /(.+):(\d+)\.\.(\d+)/;
              
              my $p1_genomic_s = $genomic_s + range_get_start($p1_range) + 1 ; 
              my $p1_genomic_e = $genomic_s + range_get_end($p1_range) + 1 ; 
              my $p1_len    = $primers{$id}{qLen}{$p1};

              my $p2_genomic_s = $genomic_s +  range_get_start($p2_range) + 1 ; 
              my $p2_genomic_e = $genomic_s + range_get_end($p2_range) + 1 ; 
              my $p2_len    = $primers{$id}{qLen}{$p2};

              $p1_range = range_create($p1_genomic_s,$p1_genomic_e);
              $p2_range = range_create($p2_genomic_s,$p2_genomic_e);
              my $exon_ranges =  $exons[ $tie_array_id ];
              my @exon_ranges =  split /\s+/, $exon_ranges;
              my $gene_name = shift @exon_ranges; 
              my $exon_count = 0;
              foreach my $exon_range (@exon_ranges){ 
                $exon_count++;
                my ($e_s,$e_e) = $exon_range =~ /.+:(\d+)\.\.(\d+)/;
                my $e_range = range_create ($e_s,$e_e);
                my $exon_size = $e_e - $e_s + 1;
                ## something about different exons or spanning exon junctions
                
                my $p1_exon_info = getExonInfo( $p1_len, $e_range , $p1_range );
                my $p2_exon_info    = getExonInfo( $p2_len, $e_range , $p2_range ) ;
                $exon_info_sum += ($p1_exon_info + $p2_exon_info);
                push @exon_info, "ex:$exon_count($exon_size bp)|p1:$p1_exon_info|p2:$p2_exon_info";
              }
              $p1_range =  "$genomic_ref:".join('..',@$p1_range);
              $p2_range =  "$genomic_ref:".join('..',@$p2_range);
            }else{
              $p1_range =  "$tName:".join('..',@$p1_range);
              $p2_range =  "$tName:".join('..',@$p2_range);
            }
            #print "$id\t$product_size\t$tName\t@$p1_range\t@$p2_range\n";
            my $exon_info = join(";",@exon_info);
            $results{$id}{$product_size} =
              "$id\t$product_size\t$tName\t$p1_range\t$p2_range\texonSpanScore=$exon_info_sum|$exon_info\n";
          }
        }
      }
    }
  }
#}

foreach my $id ( sort keys %results ) {
  foreach my $product_size ( sort { $a <=> $b } keys %{ $results{$id} } ) {
    print $results{$id}{$product_size}, "\n";
  }
}

#####SUBROUTINES########
sub getExonInfo {
  my $p_len = shift;
  my $exon_range = shift;
  my $p_range = shift;
  
  my $overlap = range_overlap ($exon_range,$p_range);
  if ($overlap == $p_len){
  ## primer completely within 1 exon
    return 1;
  }elsif ($overlap > 0){
  ## primer overlap junction, in more than 1 exon
    return 0.5;
  }else { 
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
