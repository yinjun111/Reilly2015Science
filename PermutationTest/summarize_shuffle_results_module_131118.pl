#!/usr/bin/perl -w
use strict;

my ($infile,$shufflefiles,$resultfile)=@ARGV;

#original file
my %ori2num;
my %ori2all;
my @title;
open(IN,$infile) || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;	
	if ($_=~/^Feature/) {
		@title=@array;
		next;
	}
	$ori2num{$array[0]}=$array[1];
	$ori2all{$array[0]}=$array[2];#for binom test
}
close IN;

my @files=glob($shufflefiles);
my %shuffle2num;
foreach my $file (@files) {
	open(IN,$file) || die $!;
	while(<IN>) {
		tr/\r\n//d;
		next if $_=~/^Feature/;
		my @array=split/\t/;
		push @{$shuffle2num{$array[0]}},$array[1]; #Edited 03.26
	}
	close IN;
}

open(OUT,">$resultfile") || die $!;
print OUT "Name\t$title[2]\t$title[1]\tAverage of random shuffling\tNumber of higher or equal count in ",scalar(@files)," times shuffling\tChance of higher or equal count\n";
foreach my $motif (sort keys %ori2num) {
	my $rank=rank_num($ori2num{$motif}, @{$shuffle2num{$motif}}); #Edited 03.26
	my $avg=average( @{$shuffle2num{$motif}});  #Edited 03.26
	print OUT $motif,"\t",$ori2all{$motif},"\t",$ori2num{$motif},"\t",$avg,"\t",$rank,"\t",$rank/@files,"\n"; #10.03
}
close OUT;

sub rank_num {
	my ($num,@nums)=@_;
	my $rank=0;
	foreach my $n (sort {$b<=>$a} @nums) {
		if($n>=$num) { #Edited 03.26
			$rank++;
		}
		else {
			return $rank;
		}
	}
	return $rank;
}

sub average {
	my @nums=@_;
	#if(@nums==0) {
	#	return 0;
	#}
	#else {
		my $sum=0;
		foreach my $num (@nums) {
			$sum+=$num;
		}
		return sprintf("%0.0f",$sum/@nums);
	#}
}

