#!/usr/bin/perl
#use strict;
use warnings;

print "Hello, World...\n";

my $file1=$ARGV[0];
my %hash1=();
#my ($name,$line);
#my $name;

open(FASTA,"$file1") or  die "cannot open input file: $!";

while (<FASTA>) {
	chomp;
	my $line=$_;
	if ( $line=~/^>./) {    #match line begin with ">"

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
print "$_\n$hash1{$_}\n" foreach (keys %hash1);