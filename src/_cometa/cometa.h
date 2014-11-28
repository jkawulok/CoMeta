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

#ifndef _COMETA_H
#define _COMETA_H

#include "defs.h"
#include <iostream>
#include <string>
#include <vector>
#include "queues.h"
#include "wrapers.h"
#include "timer.h"

using namespace std;

class CCOMETA {
	bool initialized;

	list<string> input_file_name;
	string output_file_name;
	vector<string> working_directories;
	bool is_fasta;
	bool is_read;

	int kmer_len;			// kmer length
	int prefix_len;			// prefix length; fixed: 4 or 5

	int n_disks;
	int n_bins;				// number of bins; fixed: 448
	int bin_part_size;		// size of a bin part; fixed: 2^16
	int fastq_buffer_size;	// size of FASTQ file buffer; fixed: 2^25

	int n_splitters;		// number of splitters; default: 1
	int n_sorters;			// number of sorters; default: 1
	int n_omp_threads;		// number of OMP threads per each sorter

	uint64 max_mem_size;		// maximum amount of memory to be used in GBs;	default: 30GB
	uint64 max_mem_fastq;		// maximum amount of memory for FASTQ parts queue
	uint64 max_mem_bin_part;	// maximum amount of memory for parts of bins queue
	uint64 max_mem_storer;		// maximum amount of memory for internal buffers of KmerStorer
	uint64 max_mem_stage2;		// maximum amount of memory in stage 2

	CStopWatch w1, w2;

	// Memory monitor
	CMemoryMonitor *mm;
	CDiskMonitor *dm;
	CBinOrdering *bo;

	// Queues
	CPartQueue *pq;
	CBinPartQueue *bpq;
	CBinDesc *bd;
	CBinQueue *bq;
	CKmerQueue *kq;

	// Boost thread groups
	boost::thread_group gr1, gr2, gr3;		// thread groups for 1st stage
	boost::thread_group gr4, gr5, gr6;		// thread groups for 2nd stage

	uint64 n_unique, n_singletons, n_total, n_reads, tmp_size;
	uint32 n_cores;

	// Class wrapers
	CWFastqReader *w_fastq;
	
	vector<CWSplitter*> w_splitters;

	CWKmerBinStorer *w_storer;

	vector<CWKmerBinReader*> w_readers;
	vector<CWKmerBinSorter*> w_sorters;
	CWKmerBinCompleter *w_completer;

	void AdjustToHardware(uint32 cores);
	bool AdjustMemoryLimits();

public:
	CCOMETA();
	~CCOMETA();

	void SetFileNamesList(string file_name_FILEsequence, list<string> &file_name_sequence, string path_seq);
	void SetFileNames(list<string> _input_file_name, string _output_file_name, vector<string> _working_directories, bool _is_fasta, bool _is_read);
	void SetFileNames(list<string> _input_file_name, string _output_file_name, string _working_directory, bool _is_fasta, bool _is_read);
	void SetParams(int _kmer_len, int _n_splitters, int _n_sorters, int _n_omp_threads, uint32 _max_mem_size_in_gb);
	void SetParams(int _kmer_len, int _n_threads, uint32 _max_mem_size_in_gb);

	bool Process();

	void GetStats(double &time1, double &time2, uint64 &_n_unique, uint64 &_n_singletons, uint64 &_n_total, uint64 &_n_reads, uint64 &_tmp_size, uchar &_kmer_len, uint32 &_n_cores, uint32 &_n_splitters, uint32 &_n_sorters, uint32 &_n_omp_th, uint32 &_max_mem_size_in_gb);
};

#endif
