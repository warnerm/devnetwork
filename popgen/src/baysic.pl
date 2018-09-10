#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use List::Util qw(max sum);
use File::Temp; 
use File::Basename;
use File::Which; 

my $usage = "example usage:\n\t " . 'baysic.pl --statsOutFile combined.stats --pvalCutoff 0.8 --vcf file1.vcf --vcf file2.vcf --countsOutFile combined.cts --vcfOutFile combined.vcf' . "\n";
my @files;
my %sites;
my $ctr = 0;
my $countsOutFile;
my $vcfOutFile;
my $statsOutFile,
my $pvalCutoff = 0.8;

# by default look in directory where $0 is
my $basedir = dirname($0);
my $lcaScript = join("/", $basedir, "lca.R");

my $opt = GetOptions (
    "vcf=s" => \@files, 
    "countsOutFile=s" => \$countsOutFile,
    "vcfOutFile=s" => \$vcfOutFile,
    "statsOutFile=s" => \$statsOutFile,
    "pvalCutoff=f" => \$pvalCutoff,
    );

die "--lcaScript $lcaScript isn't executable!\n$usage" unless -x $lcaScript;
die "no --vcf arguments given!\n$usage" if @files < 1;
die "need more than 1 --vcf argument!\n$usage" if @files < 2;
die "no --countsOutFile arguments given!\n$usage" if ! defined $countsOutFile;
die "no --vcfOutFile arguments given!\n$usage" if ! defined $vcfOutFile;
die "no --statsOutFile arguments given!\n$usage" if ! defined $statsOutFile;

my $nfiles = scalar(@files);

# make sure we can find tabix
my $tabix_exec;
unless ( $tabix_exec = which("tabix")){
    die "can't find an executable tabix, is it installed? (See https://bitbucket.org/jtr4v/baysic for more installation advice.)";
}

# make sure we can find bgzip
my $bgzip_exec;
unless ( $bgzip_exec = which("bgzip")){
    die "can't find an executable bgzip, is it installed? (See https://bitbucket.org/jtr4v/baysic for more installation advice.)";
}

my $vcf_merge_exec;
# make sure we can find vcf-merge
unless ( $vcf_merge_exec = which("vcf-merge") ){
    die "can't find an executable vcf-merge, is vcftools installed? (See https://bitbucket.org/jtr4v/baysic for more installation advice.)";
}

open my $COUNTS, ">$countsOutFile" || die "can't open counts outfile $countsOutFile: $!\n";
open my $VCF, ">$vcfOutFile" || die "can't open vcf outfile $vcfOutFile: $!\n";
open my $VCFPOS, ">$vcfOutFile.pos" || die "can't open vcf outfile $vcfOutFile: $!\n";
open my $STATS, ">$statsOutFile" || die "can't open statsOutFile $statsOutFile: $!\n";

# write position file header
print $VCFPOS "chromosome\tposition\n";

for my $file (@files) {
  open(FILE, $file) or die "problem accessing file $file: $!";
  while (my $line = <FILE>) {
    if ( $line =~ m/^#/) {
	next;
    }
    chomp $line;
    my ($chr, $pos) = split(" ", $line);
    $sites{$chr}->{$pos} ||= [ (0) x $nfiles ];
    $sites{$chr}->{$pos}->[$ctr]++;
  }
  $ctr++;
}

#######################
# write out counts file 
#######################
writeCountsFile();

##########################################
# run lca.R to get posterior probabilities
##########################################
system("$lcaScript --countsFile $countsOutFile --statsOutFile $statsOutFile > /dev/null 2>&1");

# suck in statistics
my $stats = "";
open my $INSTATS, "$statsOutFile" || die "can't read in stats file $statsOutFile!\n";
while ( my $line = <$INSTATS> ){
    $stats .= $line;
}
warn "WARNING: stats file seems to be empty, so it's likely that statistics were not generated - are R and the required R packages installed correctly?\nSee here for installation hints: https://bitbucket.org/jtr4v/baysic\n" if length($stats) == 0;
close $INSTATS;

# stats summary comes out like this:
# postprobs[1] 9.999619e-01 5.134753e-06 4.895795e-08   4.895795e-08
# postprobs[2] 1.356514e-06 9.548155e-07 9.103808e-09   9.325555e-09
# postprobs[3] 1.344586e-06 9.527360e-07 9.083981e-09   8.905481e-09
# postprobs[4] 1.829956e-12 2.037073e-12 1.942272e-14   0.000000e+00
# .. [m]

# $postProbArray is the parsed stats summary, arranged like this:
# $postProbArray->[3]->[1]
# where 3 is the index of the SNP caller
# and item number 1 is the mean, 2 is the SD, 3 is the Naive SE, 4 is the Time-series SE
my $postProbArray = retrievePostProbs( $stats );

# make a grid of 0/1 values to match up with probabilities for the sake of my own sanity
# (may or may not use this)
my $ContingencyArrayRef =  makeContingencyArray();

# make a hash to turn contingency string into an index I can use to look up post. probabilities
# using postProbArray
# eg. '010' -> 5
#  --> $postProbArray->[5]->[1] == 0.02 
my $contigency2Index = makeContigency2IndexHash();

# go through file again and write out those SNPs that pass posterior probability threshold
# to get 0/1 status of each site for all callers: 
#       join("", @{$sites{1}->{1046829}})

# for each input file, write out new VCF with only lines that pass the pvalue threshold
my @filteredVCFFiles; # keep track of where we write out the passing files
my %seenPos;
for my $file (@files) {
  open(FILE, $file) or die "problem accessing file $file: $!";

  my $filename_wo_path = fileparse( $file );
  my $outTempFile = $filename_wo_path . getTempFile();
  push @filteredVCFFiles, $outTempFile;
  open(my $FILTERED, ">$outTempFile") or die "problem open tempfile $outTempFile while writing out filtered VCF: $!";
  while (my $line = <FILE>) {
      chomp $line;
      if ($line =~ m/^#/){
	  print $FILTERED $line . "\n";
	  next;
      }
      my ($chr, $pos) = split(/\t/, $line);
      my $contingency_string = getContingencyString($chr, $pos);
      warn "Didn't get $nfiles items in contingency string!" if (length $contingency_string != $nfiles);
      if ( getPostProb( $contingency_string ) > $pvalCutoff ){
	  print $FILTERED $line . "\n";
	  print $VCFPOS join("\t", $chr, $pos) . "\n" unless exists $seenPos{$chr}->{$pos};
	  $seenPos{$chr}->{$pos}++;
      }
  }
}

# bgzip files and tabix them
# warn "bgzipping and tabix'ing files...";
foreach my $file ( @filteredVCFFiles ){
    system("$bgzip_exec $file");
    system("$tabix_exec -f -p vcf ${file}.gz");
}
# warn "done.";

# warn "merging vcf files";
my $vcf_file_arg = join(" ", map {$_ . ".gz"} @filteredVCFFiles);
my $cmd = "$vcf_merge_exec $vcf_file_arg 1> $vcfOutFile 2> /dev/null";
my $mergeReturn = system("$vcf_merge_exec $vcf_file_arg 1> $vcfOutFile");
if ( $mergeReturn != 0){
    warn "The VCF merge command:\n\t$cmd\n seems to have failed ($mergeReturn)! The vcf is likely to be empty (the position file should still be useful though).\n";
}

# remove tempfiles in @filteredVCFFiles
foreach my $file ( @filteredVCFFiles ){
    unlink $file;
    unlink $file . ".gz";
    unlink $file . ".gz.tbi";
}

sub getPostProb {
    my $contingency_string = shift;
    my $indexInPostProbArray;
    $indexInPostProbArray = $contigency2Index->{$contingency_string};
    if ( exists $contigency2Index->{$contingency_string} && exists $postProbArray->[$indexInPostProbArray]->[1] ){
	return $postProbArray->[$indexInPostProbArray]->[1];
    }
    else {
	warn "couldn't find postProb for $contingency_string!\n";
	return undef;
    }
}

sub writeCountsFile { # and also decide which SNVs are in each VCF file
    my @counts = (0) x (2 ** $ctr);

    my $sum = 0;
    while (my ($chr, $chrsites) = each %sites) {
	$sum += max(grep {$_ =~ m/\d+/} keys %$chrsites); # for some reason "position" and "chr" were showing up here and screwing up the sum
	while (my ($pos, $ctrs) = each %$chrsites) {
	    my @status = map { $_ > 0 ? 1 : 0 } @$ctrs;
#	    print $VCF join("\t", $chr, $pos, ".", "N", "N", ".", ".", join(";", @status, sum(@status)) ) . "\n";
	    my $idx = 0;
	    for (my $i = 0 ; $i < $nfiles ; $i++) {
		if ($ctrs->[$i] > 0) {
		    $idx += 2 ** $i;
		}
	    }
	    $counts[$idx]++;
	}
    }
    
    print $COUNTS "Estimated sum: $sum\n";
    
    $counts[0] = $sum - sum(@counts[1..$#counts]);
    
    for (my $idx = 0 ; $idx < @counts ; $idx++) {
	my @labels = (0) x $nfiles;
	for (my $i = 0 ; $i < @labels ; $i++) {
	    if ($idx & (1 << $i)) {
		$labels[$i] = 1;
	    }
	}
	print $COUNTS join("\t", @labels, $counts[$idx]), "\n";
    }
    
}

sub getContingencyString {
    my $chr = shift;
    my $pos = shift;
    warn "Can't find any counts for $chr:$pos! (This shouldn't happen!)\n" unless (exists $sites{$chr}->{$pos});
    # the following makes a string of 0s and 1s of length $nfiles, such that any integer > 1 becomes a 1
    return join("", map {$_ > 0 ? 1 : 0} @{$sites{$chr}->{$pos}});
}

sub getTempFile {
    my $tempfile_object = new File::Temp( 
	TEMPLATE => '_XXXXXX',
	SUFFIX => '.temp'
	);
    return $tempfile_object->filename;
}

sub makeContingencyArray { # just a convenience method so I can check match up 0/1 labels with posterior probabilities from R
    my @zero_one_labels;
    my $num_permutations = 2**$nfiles;
    for (my $idx = 0 ; $idx < $num_permutations; $idx++) {
	my @labels = (0) x $nfiles;
	for (my $i = 0 ; $i < @labels ; $i++) {
	    if ($idx & (1 << $i)) {
		$labels[$i] = 1;
	    }
	}
	push @zero_one_labels, join("", @labels);
    }
    return @zero_one_labels;
}

sub makeContigency2IndexHash { # just a convenience method so I can check match up 0/1 labels with posterior probabilities from R
    my %contingency2index;
    my $num_permutations = 2**$nfiles;
    my $count = 0;
    for (my $idx = 0 ; $idx < $num_permutations; $idx++) {
	my @labels = (0) x $nfiles;
	for (my $i = 0 ; $i < @labels ; $i++) {
	    if ($idx & (1 << $i)) {
		$labels[$i] = 1;
	    }
	}
	$contingency2index{join("", @labels)} = $count;
	$count++;
    }
    return \%contingency2index;
}

sub retrievePostProbs {
    my $stats = shift;
    my $postProbs;
    foreach my $line ( split /\n/, $stats ){
	next if $line !~ /^postprobs\[\d/;
	my @fields = split /\s+/, $line;
	push @{$postProbs}, \@fields;
    }
    return $postProbs;
}
