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

#ifndef _WRAPERS_H
#define _WRAPERS_H

#include <stdio.h>
#include <iostream>
#include <queue>
//#include <tuple>
#include <list>
#include <map>
#include "defs.h"
#include <string>
#include "fastq_reader.h"
#include "kmer_reader.h"
#include "queues.h"
#include "splitter.h"
#include "splitter2class.h"
#include "splitter_kmer.h"
#include "classifier.h"
#include "Tr_wrapers.h"

#include "boost/tuple/tuple.hpp" ////////////////
using namespace boost; //////////////////////////


using namespace std;

//----------------------------------------------------------------------------------
class CWFastqReader {
	CMemoryMonitor *mm;

	CFastqReader *fqr;
	list<string> file_name;
	uint64 part_size;
	CPartQueue *part_queue;

public:
	CWFastqReader(CMemoryMonitor *_mm, list<string> _file_name, uint64 _part_size, CPartQueue *_part_queue, bool _is_fasta, bool _is_read, int _kmer);
	CWFastqReader(CMemoryMonitor *_mm, string _file_name, uint64 _part_size, CPartQueue *_part_queue, bool _is_fasta, bool _is_read, int _kmer);
	~CWFastqReader();

	void operator()();
};


//----------------------------------------------------------------------------------
class CWSplitter {
	CMemoryMonitor *mm;

	CSplitter *spl;

	CPartQueue *pq;
	CBinPartQueue *bpq;
	CBinDesc *bd;

public:
	CWSplitter(CMemoryMonitor *_mm, CPartQueue *_pq, CBinPartQueue *_bpq, CBinDesc *_bd, uint32 _buffer_size, int _kmer_len, int _prefix_len, int _n_bins, bool _is_fasta, bool _is_read);
	~CWSplitter();

	void operator()();
	void GetTotal(uint64 &_n_reads);
};


//-------------------- Separation loaded reads to the classification----------------------------
class CWSplitter2class {
	CMemoryMonitor *mm;

	CSplitter2class *spl;

	CPartQueue *pq;
	CReadsQueue *pread;


public:
	CWSplitter2class(CMemoryMonitor *_mm, CPartQueue *_pq, CReadsQueue *_pread,  bool _is_fasta);
	~CWSplitter2class();

	void operator()();
	void GetTotal(uint64 &_n_reads);
};

//----------------------------------------------------------------------------------
class CWKmerBinStorer {
	CMemoryMonitor *mm;

	CKmerBinStorer *kbs;

	int n_bins;
	CBinPartQueue *q;
	CBinDesc *bd;
	vector<string> working_directories;
	uint32 prefix_len;
	uint32 kmer_len;

public:
	CWKmerBinStorer(CMemoryMonitor *_mm, int _n_bins, CBinPartQueue *_q, CBinDesc *_bd, vector<string> _working_directories, 
		uint32 _prefix_len, uint32 _kmer_len, uint64 _max_mem_buffer);
	~CWKmerBinStorer();

	void operator()();
};

//----------------------------------------------------------------------------------
class CWKmerBinReader {
	CMemoryMonitor *mm;
	CDiskMonitor *dm;
	CBinOrdering *bo;

	CKmerBinReader *kbr;

	CBinDesc *bd;
	CBinQueue *bq;

public:
	CWKmerBinReader(CMemoryMonitor *_mm, CDiskMonitor *_dm, CBinOrdering *_bo, CBinDesc *_bd, CBinQueue *_bq);
	~CWKmerBinReader();

	void operator()();
};

//----------------------------------------------------------------------------------
class CWKmerBinSorter {
	CMemoryMonitor *mm;

	CKmerBinSorter *kbs;

	int32 n_bins;
	CBinDesc *bd;
	CBinQueue *bq;
	CKmerQueue *kq;
	uint32 kmer_len;

public:
	CWKmerBinSorter(CMemoryMonitor *_mm, int32 _n_bins, CBinDesc *_bd, CBinQueue *_bq, CKmerQueue *_kq, int _n_omp_threads, uint32 _kmer_len);
	~CWKmerBinSorter();

	void operator()();
};

//----------------------------------------------------------------------------------
class CWKmerBinCompleter {
	CMemoryMonitor *mm;

	CKmerBinCompleter *kbc;

	string file_name;
	CKmerQueue *kq;
	CBinDesc *bd;
	uint32 kmer_len;

public:
	CWKmerBinCompleter(CMemoryMonitor *_mm, string _file_name, CBinDesc *_bd, CKmerQueue *_kq, uint32 _kmer_len);
	~CWKmerBinCompleter();

	void operator()();

	void GetTotal(uint64 &_n_unique, uint64 &_n_singletons, uint64 &_n_total);
};



//--------------------- Classification (new bin) -------------------------------
class CWClassifier_new{
	
	CMemoryMonitor *mm;

	CClassifier_new *clr;

	CBinKmers_all_new *bin_kmers;
	CReadsQueue *pread;
	double matchcutoff;



public:
	CWClassifier_new(CMemoryMonitor *_mm, CReadsQueue *_pread, CBinKmers_all_new *_bin_kmers, uchar _kmer_len, uchar _step_k, double _matchcutoff, FILE* _outMatch, FILE* _outMisMatch);
	
	CWClassifier_new(CWClassifier_new& obj)
	{
		printf("Operation forbidden!\nCWClassifier cannot be copied...");
		getchar();
	}
	~CWClassifier_new();

	void operator()();


	//CWClassifier_new& operator = (CWClassifier_new& obj)
	//{
	//	printf("Operation forbidden!\nCWClassifier cannot be assigned...");
	//	getchar();
	//}


	void GetTotal(uint64 &_num_alig_reads, uint64 &_num_REValig_reads, uint64 &_num_NOalig_reads);

};


#endif
