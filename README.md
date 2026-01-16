#PrimerChecker
check primer mapping to genome or cDNA

1. need to run getTranscriptInfo_fromSFF.pl or getTranscriptInfo_fromGFF.pl to generate a file of transcript exon info. This is necessary for cDNA primer checking.
2. need to make a dir in your cgi directory called 'dbs'. In 'dbs' you need to have a directory for each organim that you have as an option in the checkbox. Would be nice to auto generate the check box based on the direcotries in dbs
3. make sure to modify the organims checkbox to only include organisms you have in dbs/
4. include a genome fasta in dbs/{organism}
5. need to have blat installed and have its complete path in the code.
6. split the fasta into smaller fastas or the website will timeout while  blat is running. or fix the timeout parameter

#my $db_dir = "dbs/$organism";
#open EXONS, "$db_dir/transcript_exons_info.txt" 
#  or die "Can't open exon info file: transcript_exons_info.txt\n" if $type eq 'cDNA';
#my @db_files = <$db_dir/*fasta>;
