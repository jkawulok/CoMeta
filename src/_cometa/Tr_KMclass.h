/*
	This file is a part of CoMeta software.
	The homepage of the CoMeta project is http://sun.aei.polsl.pl/cometa
 
	Authors: Jolanta Kawulok and Sebastian Deorowicz
	
	Version: 0.3
	Date   : 2014-11-26

*/

#ifndef _TR_KM_CLASS_H
#define _TR_KM_CLASS_H

#include "defs.h"
#include <iostream>
#include <string>
#include <vector>
#include "queues.h"
#include "wrapers.h"
//#include "timer.h"
#include <boost/thread.hpp>

using namespace std;

class TrKMclass {

	bool initialized;

	string data_file_name;
	string output_file_name;


	bool is_fasta;

	uchar kmer_len;			// kmer length
	uint32 prefix_len;			// prefix length; 4
	uint64 num_kmer; //number of unique kmer

	int n_bins;				// number of bins; 
	int bin_part_size;		// size of a bin part; fixed:
	int fastq_buffer_size;	// size of FASTQ file buffer; fixed: 
		uint32  n_cores;
	uint32 kmers_buffer_size;	// size of kmers file buffer; fixed: 

	int n_splitters;		// number of splitters; default: 1

	
	uint64 max_mem_size;		// maximum amount of memory to be used in GBs;	default: 30GB
	
	uint64 max_mem_fastq;		// maximum amount of memory for FASTQ parts queue
	uint64 max_mem_reads;		// maximum amount of memory for reads parts queue
	uint64 max_mem_seqKMER;		
	uint64 max_mem_binsKMER;



	// Memory monitor
	CMemoryMonitor *mm;


	// Queues
	//CPartQueue *pq;
	CPartKMERQueue *pkmer;



	CBinKmers_all_new *bin_kmers; // all sorted kmers  from file new --------------------------------------------

	FILE *outDatabase;

	// Boost thread groups
	boost::thread_group gr1, gr2, gr3, gr4;		// thread groups for 1st stage
	boost::thread_group  gr5, gr6;		// thread groups for 2nd stage

	uint64 tmp_size;

	// Class wrapers
	//CWFastqReader *w_fastq;
	CWKmersReader *w_kmers;

	vector<CWKmersSplitter_new*> w_splittersKMER; //insert kmers into bins new -------------------------------
	
	void AdjustToHardware(uint32 cores);
	bool AdjustMemoryLimits();		

public:
	TrKMclass();
	~TrKMclass();

	void SetFileNames(string _data_file_name, string _output_file_name);
	void SetParams(int _n_threads, uint32 _max_mem_size_in_gb);

	bool Process();
	bool New_kmer_len(CWKmersReader *w_kmers);

	//void GetStats(double &_time1, double &_time2, uint64 &_num_alig_reads, uint64 &_num_REValig_reads, uint64 &_num_NOalig_reads,  uint64 &_n_reads, uchar &_kmer_len, double &_matchcutoff, uint32 &_n_cores, uint32 &_n_splitters, uint32 &_n_classification, uint32 &_max_mem_size_in_gb);
	void GetStats();

	void SaveOutFile();

};

#endif