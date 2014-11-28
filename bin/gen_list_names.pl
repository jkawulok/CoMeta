#!/usr/bin/env perl

#	This file is a part of CoMeta software.
#	The homepage of the CoMeta project is http://sun.aei.polsl.pl/cometa

use strict;
use Benchmark;
use warnings;
use File::Spec;
use vars qw($SEP);
use Getopt::Long;

my $pathScores="";
my $nameFile=""; 
my $pathFile=""; 

GetOptions(
			'ps:s'      => \$pathScores,
			'nf:s'		=> \$nameFile,
			'pf:s'		=> \$pathFile,
       );

my $allPath=$pathFile."/".$nameFile;
open(QUERY, "<$allPath") or die "Can't open $allPath:$!, \n";    

my $old_name="";
my $allPath_OUT;
my $start=1;

while (<QUERY>) 
{
        chomp $_;

        if (m/^\s+$/) 
		{		# Avoid blank lines
            next;
        }
		
		my @Tab = split( /_NT_/, $_);
		
		if ($old_name eq $Tab[0])
		{
			print FILE_OUT $_ ,"\n"; # name
		}
		else
		{
			if ($start==0)
			{
				close(FILE_OUT)
			}
			
			$old_name=$Tab[0];
			$allPath_OUT=$pathScores."/NW_".$old_name.".txt";
			open(FILE_OUT, ">$allPath_OUT") or die "Can't open $allPath:$!, \n";  
			print FILE_OUT $_ ,"\n"; # name			
			# print STDERR "$old_name \n";
			$start=0;
		}

		
}
 
 close(QUERY);