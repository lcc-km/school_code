#!/usr/bin/perl
use strict;
use warnings;

###############################################################
## usage :  perl fasta_hash1_in_keep.pl <file1: fasta > <file2: geneid & protein id> >out.file
##
## eg:   perl fasta_in_keep.pl prtein.fa keep_input.txt >result.fa
##
###############################################################

my $file1=$ARGV[0];   # input fasta file name
my %hash1=();
my ($name,$line);


open(FASTA,"$file1") or  die "cannot open input file: $!";

while (<FASTA>) {
    chomp;
    my $line=$_;
    if ( $line=~/^>./) {            ## match line begin with ">"
        $hash1{$line}="";        
        $name=$line;
    
    } else{
        if ($hash1{$name} eq "") {
        $hash1{$name}=$line
        }  else {
        $hash1{$name}.=$line;
        
        }
    
    }
        
}
close FASTA;
#print "$_\n$hash1{$_}\n" foreach (keys %hash1);  # print content of hash1
#print keys %hash1;


my $file2=$ARGV[1];   # input longest protein ID name
open(KEEP, $file2) or  die "cannot open input file: $!";

while (<KEEP>) {
    chomp;
	my @itm=split /\t/,$_;
    my $pname=$itm[0];
	my $gname=$itm[1];
    foreach  (keys %hash1) {
        if ($_=~/.$pname./) {            # mach protein ID name with fasta name line
			my $fname= substr ($_, 1, length($_));			# delele > in fasta seq name
            print  ">$gname $fname\n$hash1{$_}\n";    # add geneid before fasta
        }
    }
}