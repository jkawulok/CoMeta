# A program for dividing reference sequences into the group rank

#	This file is a part of CoMeta software.
#	The homepage of the CoMeta project is http://sun.aei.polsl.pl/cometa


#!/usr/bin/env perl
use lib qw( /mnt/auto/people/plgjkawulok/perl_modules/lib/perl5/site_perl/5.8.8/);
use strict;
use Benchmark;
use Bio::LITE::Taxonomy::NCBI;
use warnings;
use File::Spec;
use vars qw($SEP);
use Getopt::Long;
# use Bio::DB::Taxonomy;
use Benchmark;

my %dict_phylum; 
my %namefile;
my %num_seq;
my $num_seqALL=0;
my $file_in; 
my $prefix_out='';
my $suffix_out='';
my $path_in='';
my $path_out='';
my $path_NCBItax='';
my $queryid;
my $indTAX = 0;
my $indTAXfile = 0;
my $indTAXfilePREV = 0;
my $nameTAX = '';
my $nameTAX1 = '';
my $nameTAX2 = '';
my $nameTAX_old = '';
my $first = 1;


my $name;

#my $dbINT = new Bio::DB::Taxonomy(-source => 'entrez');
#my %tabidx;

my ($t0read, $t1read, $tdread);
my ($t0read2, $t1read2, $tdread2);
my ($t0ALL, $t1ALL, $tdALL);
my ($level_first, $level_second);


GetOptions(
		'Pin:s'			=> \$path_in,
		'Pout:s'          		=> \$path_out,
		'Pncbi:s'			=> \$path_NCBItax,
		'LF:s'          		=> \$level_first,
		'LS:s'          		=> \$level_second,
		'i|in:s'          		=> \$file_in,
		'op|outprefix:s'        	=> \$prefix_out,
		'os|outsuffix:s'        	=> \$suffix_out,
		# 'h|help'          		=> sub { system('perldoc', $0); exit }
       );

	   
#---------------------------------
print "\n";
$t0ALL= new Benchmark;	


print "Input file: $file_in \n";



# Open files with sequences
my $path_file_in = $path_in . "/" . $file_in;
my ($pin, $fin) = f_path_file($path_file_in);
my $fnamesh = f_ending($fin);
open(QUERY, "<$path_file_in") or die "Can't open $path_file_in $! \n";



$t0read= new Benchmark;	
my $dbLITE  = Bio::LITE::Taxonomy::NCBI->new (
						db=>"NCBI",
						names=> $path_NCBItax . "/names.dmp",
						nodes=> $path_NCBItax . "/nodes.dmp"
                                              );
											    $t1read = new Benchmark;
$tdread = timediff($t1read, $t0read);
print "Time read nodes and names was taken  ", timestr($tdread, 'all'), " seconds \n\n";






while (<QUERY>) 
{
	if (m/^\s+$/) 
	{		# Avoid blank lines
		next;
	}

	s/\r//g;			# removes takM
	if (/>/)  #
	{
		chomp($queryid = $_); 
		
		# $indGI = f_take_GI($queryid); #number GI
		$indTAX = f_take_TAX($queryid); #number TAX
		
		# ($nameTAX, $nameTAXorg) = f_take_name_tax($indTAX);
		
		$nameTAX1 = $dbLITE->get_term_at_level($indTAX, $level_first);
		$nameTAX2 = $dbLITE->get_term_at_level($indTAX, $level_second);
		
		$nameTAX1 =~ s{\/}{=}g;
		$nameTAX2 =~ s{\/}{=}g;
		
		if ($nameTAX2 eq 'undef')
		{
		
			$indTAXfile = f_take_TAXnamefile($path_file_in);
			if (!$indTAXfile)
			{	
				$indTAXfile = $dbLITE->get_taxid_from_name($nameTAX1);	
			}	
		
		}
		else
		{
			$indTAXfile = $dbLITE->get_taxid_from_name($nameTAX2);
			
			if ($nameTAX1 eq 'undef')
			{
				$indTAXfilePREV = f_take_TAXnamefile($path_file_in);
				# print "indTAXfilePREV $indTAXfilePREV \n";				
			}		
		}


		
		if (!$nameTAX1)
		{
			$nameTAX="notknow_".$fnamesh;
		}
		else
		{
			if ( ($nameTAX1 eq 'undef') && ($nameTAX2 ne 'undef') )
			{
				$nameTAX=$nameTAX1.'_idPREV'.$indTAXfilePREV.'_'.$nameTAX2.'_idTAX'.$indTAXfile;
				#print "nameTAX: $nameTAX \n";
			}
			else
			{		
				$nameTAX=$nameTAX1.'_'.$nameTAX2.'_idTAX'.$indTAXfile;
			}
		}
		# print "nameTAX: $nameTAX \n";
		
		unless (exists($namefile{$nameTAX}))
		{
			start_file($nameTAX);
		}
		
		$num_seq{$nameTAX}++;
		$num_seqALL++;
		
		# if ($num_seqALL%100000 == 0){
		#	print "Nr seq: ", $num_seqALL, "\n";}
		
		#print "namefile\{nameTAX\}: $namefile{$nameTAX} \n";
		if ($first == 0)
		{					
			if ($nameTAX_old ne $nameTAX)
			{
				close (FILE_OUT);
				open (FILE_OUT, ">>$namefile{$nameTAX}") or die "Can't write to $nameTAX $! \n";			
			}			
		}
		else
		{		
			open (FILE_OUT, ">>$namefile{$nameTAX}") or die "Can't write to $nameTAX $! \n";# open output file 		
			$first = 0;		
		}
		
		
		
		#print "queryid: $queryid \n";		
		#if ($indGI){print "indGI: $indGI\n";}		
		#if ($indTAX){print "indTAX: $indTAX\n";}	
		
		print FILE_OUT $queryid ,"\n"; # name
		$nameTAX_old = $nameTAX;	

	} 
	else #save sequence
	{
		print FILE_OUT $_;  
	}	
}

close (FILE_OUT);
close (QUERY);


$t1ALL = new Benchmark;
$tdALL = timediff($t1ALL, $t0ALL);

Write_file_sum($prefix_out, $tdALL, $suffix_out);
print "nALL time taken  ", timestr($tdALL, 'all'), " seconds \n";
print "-------------------------------------------------------- \n\n";


###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################

#--------------------------------------------

sub Write_file_sum
{
	my $path_file_out = $path_out . "/" . $_[0] . "_summary" . $_[2] . ".txt";
	open (SUM, ">>$path_file_out") or die "Can't write to: $path_file_out $!\n";
	foreach (sort mysort keys %num_seq )
	{
		print "Key $_ num = $num_seq{$_} \n";
		print SUM "Key $_ num = $num_seq{$_} \n";
	}
	
	print "\nALL seq: $num_seqALL\n";
	print SUM "\nALL seq: $num_seqALL\n";
	
	print SUM "\nTIME: ", timestr($_[1]), "\n\n";
	close (SUM);	
}

sub start_file
{
	
	$namefile{$_[0]} = ($path_out . "/" . $prefix_out . $_[0] . $suffix_out . ".fa");	
	$num_seq{$_[0]} = 0;
	open (FILE_OUT_start, ">>$namefile{$_[0]}") or die "Can't write to $namefile{$_[0]} $! \n";
	close(FILE_OUT_start);	
	# print "\nFILE: $namefile{$_[0]} \n";
}

#GI = f_take_GI($queryid); #number GI
sub	f_take_GI #number GI
{
	$_[0] =~ m/gi\|(.*)\|(.*)/;
	my @Tab = split( /\|/, $1);
	return $Tab[0];
}

#TAX =f_take_GI($queryid); #number GI
sub	f_take_TAX #number GI
{
	$_[0] =~ m/<idTAX\|(.*)\|/;
	return $1;
}

sub	f_take_TAXnamefile #number GI from file name
{
	$_[0] =~ m/idTAX(.*).fa/;
	return $1;
}


sub mysort 
{
	return $a cmp $b ;
}

sub	f_path_file #without path
{
	$_[0] =~ m/(.*)\/(.*)/;

	return $1, $2;
}

sub	f_ending#
{
	$_[0] =~ m/(.*).fa/;
	
	return $1;
}


