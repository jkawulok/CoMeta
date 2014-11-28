/*
	This file is a part of CoMeta software.
	The homepage of the CoMeta project is http://sun.aei.polsl.pl/cometa
 
	Authors: Jolanta Kawulok and Sebastian Deorowicz
	
	Version: 0.3
	Date   : 2014-11-26	
	
	The source codes are based on KMC 0.1 codes written by Sebastian Deorowicz and 
	Agnieszka Debudaj-Grabysz and published:
    Deorowicz S, Debudaj-Grabysz A, Grabowski S (2013) Disk-based k-mer 
	counting on a PC. BMC Bioinformatics 14.
*/

#ifndef _FASTQ_READER_H
#define _FASTQ_READER_H

#include "defs.h"
#include <stdio.h>
#include <iostream>
#include <vector>
#include "queues.h"
using namespace std;

class CFastqReader {
	CMemoryMonitor *mm;

	list<string> input_file_name;
	bool is_fasta;
	bool is_read; //
	int kmer;

	FILE *in;
	uint64 part_size;
	
	uchar *part;
	char *part_seq;
	uint64 part_filled;
	
	bool SkipNextEOL(uchar *part, int64 &pos, int64 max_pos, char *&part_seq, int kmer); //

public:
	CFastqReader(CMemoryMonitor *_mm, bool _is_fasta, bool _is_read, int _kmer);
	~CFastqReader();

	static uint64 OVERHEAD_SIZE;

	bool SetNames(list<string> _input_file_name);
	bool SetPartSize(uint64 _part_size);
	bool OpenFiles();

	bool GetPart(uchar *&_part, uint64 &_size, char *&_part_seq, int &kmer); //
};

#endif
