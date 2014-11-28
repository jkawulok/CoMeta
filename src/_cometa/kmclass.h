/*
	This file is a part of CoMeta software.
	The homepage of the CoMeta project is http://sun.aei.polsl.pl/cometa
 
	Authors: Jolanta Kawulok and Sebastian Deorowicz
	
	Version: 0.3
	Date   : 2014-11-26	
*/



#ifndef _KMclass_H
#define _KMclass_H

#include "defs.h"
#include <iostream>
#include <string>
#include <vector>
#include "queues.h"
#include "wrapers.h"
#include "timer.h"
#include <boost/thread.hpp>

using namespace std;

class CKMclass {

	bool initialized;

	list<string> read_file_name;
	string data_file_name;
	string output_file_nameMatch;
	string output_file_nameMisMatch;

	bool is_fasta;
	bool b_sortKmer;

	uchar kmer_len;			// kmer length
	uchar step_k;
	uint32 prefix_len;			// prefix length; 4
	uint64 num_kmer; //number of unique kmer

	int n_disks;
	int n_bins;				// number of bins; 
	int bin_part_size;		// size of a bin part; fixed:
	int fastq_buffer_size;	// size of FASTQ file buffer; fixed: 
	uint32  n_cores;
	uint32 kmers_buffer_size;	// size of kmers file buffer; fixed: 
	double matchcutoff;

	int n_splitters;		// number of splitters; default: 1
	int n_classification;	// number of classification
	
	
	uint64 max_mem_size;		// maximum amount of memory to be used in GBs;	default: 30GB
	
	uint64 max_mem_fastq;		// maximum amount of memory for FASTQ parts queue
	uint64 max_mem_reads;		// maximum amount of memory for reads parts queue
	uint64 max_mem_seqKMER;		
	uint64 max_mem_binsKMER;


	CStopWatch w1, w2;

	// Memory monitor
	CMemoryMonitor *mm;


	// Queues
	CPartQueue *pq;
	CPartKMERQueue *pkmer;
	CReadsQueue *pread;


	//CBinKmers_all *bin_kmers; // all sorted kmers  from file old -------------------------------------------
	CBinKmers_all_new *bin_kmers; // all sorted kmers  from file new --------------------------------------------

	FILE *outMatch;
	FILE *outMisMatch;


	// Boost thread groups
	boost::thread_group gr1, gr2, gr3, gr4;		// thread groups for 1st stage
	boost::thread_group  gr5, gr6;		// thread groups for 2nd stage

	uint64 n_reads, tmp_size;

	uint64 num_alig_reads, num_REValig_reads, num_NOalig_reads;

	// Class wrapers
	CWFastqReader *w_fastq;
	CWKmersReader *w_kmers;

	vector<CWKmersSplitter_new*> w_splittersKMER; //insert kmers into bins new -------------------------------
	
	CWSplitter2class *w_splitters2class; // split reads
	
	vector<CWClassifier_new*> w_classifier;   //new -----------------------------------------------
	
	void AdjustToHardware(uint32 cores);
	bool AdjustMemoryLimits();		

public:
	CKMclass();
	~CKMclass();

	void SetFileNames(list<string> _read_file_name, string _data_file_name, string _output_file_nameM, string _output_file_nameMM, double _matchcutoff, bool _b_sortKmer);
	void SetParams(int _step_k, int _n_splitters,  int _n_classification, uint32 _max_mem_size_in_gb);
	void SetParams(int _step_k, int _n_threads, uint32 _max_mem_size_in_gb);

	bool Process();
	bool New_kmer_len(CWKmersReader *w_kmers);

	void GetStats(double &_time1, double &_time2, uint64 &_num_alig_reads, uint64 &_num_REValig_reads, uint64 &_num_NOalig_reads,  uint64 &_n_reads, uchar &_kmer_len, double &_matchcutoff, uint32 &_n_cores, uint32 &_n_splitters, uint32 &_n_classification, uint32 &_max_mem_size_in_gb);

	void OpenOutFile();

};

#endif
