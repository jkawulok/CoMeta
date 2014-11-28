/*
	This file is a part of CoMeta software.
	The homepage of the CoMeta project is http://sun.aei.polsl.pl/cometa
 
	Authors: Jolanta Kawulok and Sebastian Deorowicz
	
	Version: 0.3
	Date   : 2014-11-26
	
*/
#ifndef _TR_WRAPERS_H
#define _TR_WRAPERS_H

#include <stdio.h>
#include <iostream>
#include <queue>
//#include <tuple>
#include <list>
#include <map>
#include "defs.h"
#include <string>
#include "kmer_reader.h"
#include "queues.h"
#include "splitter_kmer.h"

using namespace boost; //////////////////////////
using namespace std;

//--------------------- Load k-mer database ------------------------------------
class CWKmersReader {
	CMemoryMonitor *mm;

	CKmerReader *fkm;
	string file_name;
	uint64 part_size;
	CPartKMERQueue *part_kmer;
	uchar kmer_len;
	uint64* start_prefix;

public:
	CWKmersReader(CMemoryMonitor *_mm, string _file_name, uint64 _part_size,  CPartKMERQueue *_part_kmer, uchar &_kmer_len, uint64* start_prefix);
	~CWKmersReader();

	uchar get_kmer_len(){
		return kmer_len;
	};

	void operator()();
};



//-------------------- Separation loaded k-mers----------------------------
class CWKmersSplitter_new{
	CMemoryMonitor *mm;

	CSplitter_KMER_new *spl;
	CBinKmers_all_new *bin_kmers; 

	CPartKMERQueue *pkmer;
	//kosze

public:
	CWKmersSplitter_new(CMemoryMonitor *_mm, CPartKMERQueue *_pkmer, CBinKmers_all_new *_bin_kmers, uchar _kmer_len, uint32 _prefix_len);
	~CWKmersSplitter_new();

	void operator()();

};



#endif
