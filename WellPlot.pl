#!/usr/bin/perl

use strict;
use warnings;

chdir "Software";
my $WellName;
my $Mode;
if($#ARGV !=1){
  print " perl WellPlot.pl WellName AnalyteID\n";
  exit;
}else{
  $WellName = $ARGV[0];
  $Mode = $ARGV[1];
}

mkdir "../Plots";
`Rscript WellPlot.R $WellName $Mode`;
if($Mode <11 || $Mode >500){
  `tar zcvf ../plots.tar ../Plots`;
}

