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

#ifndef _KMER_BIN_H
#define _KMER_BIN_H

#include "defs.h"
#include "queues.h"
#include <string>
#include <algorithm>
//#include <array>
#include <vector>
#include <stdio.h>

using namespace std;

class CKmerBinCollector {
	CMemoryMonitor *mm;

	uint32 prefix;
	uint32 buffer_size; 
	uint32 prefix_len; 
	uint32 kmer_len; // 

	uint32 suffix_len;
	uint32 count_len; 
	uint64 kmer_mask; 
	uint64 kmer_pck_mask;

	uint32 buffer_counter, total_counter;
	uint64 *buffer, *raw_buffer;
	uchar *file_buffer;
	uint32 *pref_pack_cnts;
	uint32 bucket_shift;  // bucket_shift=p2*2 

	void Start();
	void Release();

public:
	CKmerBinCollector(CMemoryMonitor *_mm, uint32 _prefix, uint32 _buffer_size, uint32 _prefix_len, uint32 _kmer_len);
	~CKmerBinCollector();

	inline bool PutKmer(uint64 x);
	void GetCompressionParams(uint32 &_count_len)
	{
		_count_len  = count_len;
	}
	bool GetBinPart(uchar *&_bin_part, uint64 &_bin_part_size, uint64 &_n_rec);
};

inline bool CKmerBinCollector::PutKmer(uint64 x)
{
	buffer[buffer_counter++] = x & kmer_mask; // x - kmer //inserting edgings (without p1 - without number bin)

	return buffer_counter == buffer_size;
}


//**************************************************************************************************
class CKmerBinStorer {
	CMemoryMonitor *mm;

	vector<string> working_directories;
	int n_bins;
	CBinPartQueue *q;
	CBinDesc *bd;
	uint64 buffer_size_bytes;
	uint64 max_mem_buffer;
	uchar *write_buffer;
	uint32 write_buffer_size;

	FILE** files;
	string* files_name;  //
	vector<vector<pair<uchar *, uint64> > > buffer;

	void Release();
	string GetName(int n, uint32 &c_disk);
	void CheckBuffer();
	void WriteBin(uint32 n);

public:
	CKmerBinStorer(CMemoryMonitor *_mm, int _n_bins, CBinPartQueue *_q, CBinDesc *_bd, vector<string> _working_directories, uint32 _prefix_len, uint32 _kmer_len, uint64 _max_mem_buffer);
	~CKmerBinStorer();

	bool OpenFiles();
	bool CloseFiles();
	void ProcessQueue();

};


//**************************************************************************************************
class CKmerBinReader {
	CMemoryMonitor *mm;
	CDiskMonitor *dm;
	CBinOrdering *bo;

	CBinDesc *bd;
	CBinQueue *bq;

public:
	CKmerBinReader(CMemoryMonitor *_mm, CDiskMonitor *_dm, CBinOrdering *_bo, CBinDesc *_bd, CBinQueue *_bq);
	~CKmerBinReader();

	void ProcessBins();
};

//************************************************************************************************************
class CKmerBinSorter {
	CMemoryMonitor *mm; 
	CBinDesc *bd;
	CBinQueue *bq;
	CKmerQueue *kq;

	int32 n_bins;
	int32 bin_id;
	uchar *data;
	uint64 size;
	uint64 n_rec;
	string desc;
	uint32 buffer_size;
	uint32 prefix_len;
	uint32 kmer_len;
	uint32 count_len;
	int n_omp_threads;

	uint64 n_unique, n_singletons, n_total;

	uint64 *raw_buffer, *buffer;

	void Expand();
	void Sort();
	void Compact();

public:
	CKmerBinSorter(CMemoryMonitor *_mm, int32 _n_bins, CBinDesc *_bd, CBinQueue *_bq, CKmerQueue *_kq, int _n_omp_threads, uint32 _kmer_len);
	~CKmerBinSorter();

	void ProcessBins();
};


//**************************************************************************************************
class CKmerBinCompleter {
	CMemoryMonitor *mm;
	string file_name;
	CKmerQueue *kq;
	CBinDesc *bd;
	uint32 kmer_len;

	uint64 n_unique, n_singletons, n_total;

public:
	CKmerBinCompleter(CMemoryMonitor *_mm, string _file_name, CBinDesc *_bd, CKmerQueue *_kq, uint32 _kmer_len);
	~CKmerBinCompleter();

	void ProcessBins();
	void GetTotal(uint64 &_n_unique, uint64 &_n_singletons, uint64 &_n_total);
};

#endif
