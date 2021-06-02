#!/usr/bin/perl -w
use strict;

my ($sumfile,$peakannofile,$moduleannofile,$outfile,$cate1,$cate2,$cate1col,$cate2col)=@ARGV;
#file human, out count, e44, hs_gain, all, enhancer
#$cate1col=$cate2col=1;

my %seleles1;
my %seleles2;
open(IN,$peakannofile) || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	
	my %anno1=map {$_,1} split(";",$array[$cate1col]);
	my %anno2=map {$_,1} split(";",$array[$cate2col]);
	if($cate1 eq "all") {
		$seleles1{$array[0]}++;
	}
	else {
		if(defined $anno1{$cate1} ) {
				$seleles1{$array[0]}++; 
		}
	}
	if($cate2 eq "all") {
		$seleles2{$array[0]}++;
	}
	else {
		if(defined $anno2{$cate2}) {
				$seleles2{$array[0]}++;
		}       
	}
}
close IN;

open(IN,$moduleannofile) || die $!;
my %gene2module;
while(<IN>) {
	tr/\r\n//d;
	next unless $_=~/^ENS/;
	my @array=split/\t/;
	$gene2module{$array[0]}=$array[1];
}
close IN;	

#matrix of bed 2 module
open(IN,$sumfile) || die $!;
my %modules;
my %module2bed1;
my %module2bed2;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	my ($id,$name)=split(":",$array[0]);
	if($gene2module{$id}) {
		$modules{$gene2module{$id}}++;
	}
	if(defined $seleles1{$array[1]}) {
		if(defined $gene2module{$id}) {
			$module2bed1{$gene2module{$id}}{$array[1]}++;  #counting how many hs_gains are regulating the module
		}
	}
	if(defined $seleles2{$array[1]}) {
		if(defined $gene2module{$id}) {	
			$module2bed2{$gene2module{$id}}{$array[1]}++;
		}
	}
}
close IN;

#outfile
open(OUT,">$outfile") || die $!;
#$cate1,$cate2
print OUT "Feature\tNo. of $cate1 with feature\tAll $cate1\tNo. of $cate2 with feature\tAll $cate2\n";

foreach my $module (sort keys %modules) {
	print OUT $module,"\t";
	if(defined $module2bed1{$module}) {
		print OUT scalar(keys %{$module2bed1{$module}}),"\t";
	}
	else {
		print OUT "0\t";
	}
	print OUT scalar(keys %seleles1),"\t";
	if(defined $module2bed2{$module}) {
		print OUT scalar(keys %{$module2bed2{$module}}),"\t";
	}
	else {
		print OUT "0\t";
	}
	print OUT scalar(keys %seleles2),"\n";	
}
close OUT;


    

