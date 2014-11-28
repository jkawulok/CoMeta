/*
	This file is a part of CoMeta software.
	The homepage of the CoMeta project is http://sun.aei.polsl.pl/cometa
 
	Authors: Jolanta Kawulok and Sebastian Deorowicz
	
	Version: 0.3
	Date   : 2014-11-26

*/

#ifndef _KMER_READER_H
#define _KMER_READER_H

#include "defs.h"
#include <stdio.h>
#include <iostream>
#include <vector>
#include "queues.h"
using namespace std;



class CKmerReader {

	CMemoryMonitor *mm;
	string input_file_name;
	int kmer;
		
	FILE *in;
	uint64 part_size;
	uint64* start_prefix;
	
	uchar *part;
	uint64 number_suff;
	uchar number_suff_char[8];

	list<uint64> *list_num_suff;
	
	uchar kmer_len;
	uchar pref_old;
	
	int num_byte;

public:
	CKmerReader(CMemoryMonitor *_mm,  uchar &_kmer_len, uint64* _start_prefix, uint64 _part_size);
	~CKmerReader();

	static uint64 OVERHEAD_SIZE;

	bool SetNames(string _input_file_name);
	bool SetPartSize(uint64 _part_size);
	bool OpenFiles(uchar &_kmer_len);

	bool GetPart(uchar *&_part, uint64 &_size, list<uint64> *&_list_num_suff, int &_size_list);

};


#endif













