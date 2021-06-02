#!/usr/bin/perl -w
use strict;
use List::Util qw(shuffle);

my($bedfile,$annofile,$outfile) = @ARGV;

my %bed2coord;
open(IN,$bedfile) || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	$bed2coord{$array[3]}=[@array[0..2]];
}
close IN;

#get all the hs_select (3-way ortho) beds
my %hssels;
my %bed2cate;
open(IN,$annofile) || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	foreach my $cate (split(/;/,$array[1])) {
		if($cate=~/select/) {
			$hssels{$cate}{$array[0]}++;
			$bed2cate{$array[0]}=$cate;
		}
	}
}
close IN;

#shuffle every hs_select list
#reassign the link
foreach my $cate (sort keys %hssels) {
	my @beds=sort keys %{$hssels{$cate}};
	my @shuffled_beds=shuffle @beds;
	for(my $num=0;$num<@beds;$num++) {
		$hssels{$cate}{$beds[$num]}=$shuffled_beds[$num];
	}
}

#output
open(OUT,">$outfile") || die $!;
foreach my $bed (sort keys %bed2coord) {
	if(defined $bed2cate{$bed}) {
		print OUT join("\t",@{$bed2coord{$bed}}),"\t",$hssels{$bed2cate{$bed}}{$bed},"\n";
	}
	else {
		print OUT join("\t",@{$bed2coord{$bed}}),"\t",$bed,"\n";
	}
}
close OUT;


