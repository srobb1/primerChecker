#!/usr/bin/perl -w
use strict;
use Getopt::Long;
if ( !defined @ARGV ) {
  &getHelp();
}
my $file ;
my $organism ;
GetOptions(
  'f|fasta:s'         => \$file,
  'o|organism:s'           => \$organism,
  'h|help' => \&getHelp,
);

if (!defined $file){
  getHelp();
}
sub getHelp {
  die "Usage:

	./split_fasta.pl -f FASTA_FILE -o Organism

1. FASTA file will be split into many smaller files with one sequence per file
2. A directory named dbs/{ORGANISM} will be created in the working path.
3. All smaller FASTA files will be placed in this new directory
4. primerChecker.pl will look for dbs/{ORGANISM} at the same level as itself

required options:
-f | --fasta       file     fasta containing genomic sequences  
-o | --organism    str      name of organism, ex: maize
-h | --help 		    this message
 
";
  exit 1;
}

my $dir = "dbs/$organism";
system ("mkdir -p $dir");
open FASTA, $file or die "Can't open $file\n";
my @path = split /\|/ , $file;
my $fasta = pop @path;
my ($base) = $fasta =~ /(.+)\.fa|fasta%/;
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
  open OUT, ">$dir/$outFASTA" or die "Can't open $outFASTA\n";
  ${$ref2seq} =~ s/(.{80})/$1\n/g;
  print OUT ">$id\n${$ref2seq}\n";
  close OUT;
}
