/*
	This file is a part of CoMeta software.
	The homepage of the CoMeta project is http://sun.aei.polsl.pl/cometa
 
	Authors: Jolanta Kawulok and Sebastian Deorowicz
	
	Version: 0.3
	Date   : 2014-11-26	
*/


#ifndef _CLASSIFIER_H
#define _CLASSIFIER_H

#include "defs.h"
#include "queues.h"
#include <stdio.h>
#include <iostream>
#include <vector>
#include "splitter_kmer.h"
using namespace std;

//--------------------------------------------
//--------------------------------------------
//--------------------------------------------
class CClassifier_new{

	CMemoryMonitor *mm;

	CBinKmers_all_new *bin_kmers;
	CReadsQueue *pread;
	double matchcutoff;	
	double MISmatchcutoff;

	uchar kmer_len;
	uchar step_k;
	
	char* name;
	char* seq;
	uint64 size;

	uint32 match_nucl;
	uint32 trackmatch;
	uint32 trackposition;
	uint32 min_error;
	uint64 num_alig_reads;
	uint64 num_REValig_reads;
	uint64 num_NOalig_reads;

	uint64 kmer_mask;
	uint32 second_mask;
	uint64 ending_mask;

	FILE *outMatch;
	FILE *outMisMatch;


public:
	CClassifier_new(CMemoryMonitor *_mm, CBinKmers_all_new *_bin_kmers, CReadsQueue *_pread,  uchar _kmer_len, uchar _step_k, double _matchcutoff, FILE* _outMatch, FILE* _outMisMatch);
	~CClassifier_new();

	bool ProcessClass(char* _name, char* _seq, uint64 _size);
	void Check_seqence();
	void Reverse_Comp();
	uint64 Check_kmer_BIN(uchar prefix, uchar second, uint64 ending);
	uint64 Check_kmer_INTER(uchar prefix, uchar second, uint64 ending);
	void Write_score(int s);

	void GetTotal(uint64 &_num_alig_reads, uint64 &_num_REValig_reads, uint64 &_num_NOalig_reads);
};


#endif