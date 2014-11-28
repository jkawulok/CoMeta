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
#include "stdafx.h"
#include <algorithm>
#include "defs.h"
#include "splitter.h"

uint32 CSplitter::MAX_LINE_SIZE = 1 << 14;

//----------------------------------------------------------------------------------
CSplitter::CSplitter(CMemoryMonitor *_mm, CBinPartQueue *_q, CBinDesc *_bd, uint32 _buffer_size, uint32 _kmer_len, uint32 _prefix_len, uint32 _n_bins, bool _is_fasta, bool _is_read)
{
	mm = _mm;
	is_fasta = _is_fasta;
	is_read = _is_read;

	q  = _q;
	bd = _bd;
	kmer_len   = _kmer_len;
	prefix_len = _prefix_len;
	n_bins     = _n_bins;

	part = NULL;

	for(int i = 0; i < 256; ++i)
		codes[i] = -1;
	codes['A'] = codes['a'] = 0;
	codes['C'] = codes['c'] = 1;
	codes['G'] = codes['g'] = 2;
	codes['T'] = codes['t'] = 3;

	uint32 buffer_size = _buffer_size;
	uint32 count_len;

	bins = new CKmerBinCollector*[n_bins];
	for(uint32 i = 0; i < n_bins; ++i)  // change bins
	{
		//if(i < n_bins / 4)
		//{
			bins[i] = new CKmerBinCollector(mm, i, buffer_size, prefix_len+1, kmer_len);
			bins[i]->GetCompressionParams(count_len);
			bd->insert(i, "", 0, 0, 0, buffer_size, prefix_len+1, kmer_len, count_len);
		//}
		//else if(i % 4 == 0)
		//{
		//	bins[i] = new CKmerBinCollector(mm, i, buffer_size, prefix_len, kmer_len);
		//	bins[i]->GetCompressionParams(count_len);
		//	bd->insert(i, "", 0, 0, 0, buffer_size, prefix_len, kmer_len, count_len);
		//}
		//else
		//	bins[i] = NULL;
	}

	count_all_A = 0;	
	count_all_C = 0;
	count_all_A_end_G = 0;
	count_all_T = 0;

	n_reads = 0;
}

//----------------------------------------------------------------------------------
CSplitter::~CSplitter()
{
	if(bins)
		for(uint32 i = 0; i < n_bins; ++i)
			if(bins[i])
				delete bins[i];
	delete[] bins;
}

//----------------------------------------------------------------------------------
void CSplitter::Complete()  // 
{
	uint64 i;
	uint32 bin_shift = (2*kmer_len - 2*prefix_len);

	count_all_A = min<uint64>(256ULL, count_all_A);
	count_all_C = min<uint64>(256ULL, count_all_C);
	
	count_all_T = min<uint64>(256ULL, count_all_T);
	count_all_A_end_G = min<uint64>(256ULL, count_all_A_end_G);

	for(i = 0; i < count_all_A; ++i) // inserting only AAAA (max 256)
		if(bins[0]->PutKmer(0))
			PushBinToQueue(0);

	for(i = 0; i < count_all_T; ++i) // inserting only  TTTTT (max 256) 		
		if(bins[(all_T >> (bin_shift-2))]->PutKmer(all_T)) //		
			PushBinToQueue((int32)(all_T >> (bin_shift-2)));
	
	for(i = 0; i < count_all_C; ++i) // inserting only  CCCC (max 256)// 		
		if(bins[(all_C >> (bin_shift-2))]->PutKmer(all_C)) 
			PushBinToQueue((int32)(all_C >> (bin_shift-2)));

	for(i = 0; i < count_all_A_end_G; ++i) // inserting only  AAAA...AAG (max 256)
		if(bins[0]->PutKmer(firstA_endG))
			PushBinToQueue(0);

	if(bins)
		for(uint32 i = 0; i < n_bins; ++i) 
		{
			
			if(bins[i]) // bin exsits
				PushBinToQueue(i);
		}


}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
bool CSplitter::ProcessReads(uchar *_part, uint64 _part_size)
{
	part	  = _part;
	part_size = _part_size;
	part_pos  = 0;

	char *seq;
	
	if (is_read)
		seq = new char[MAX_LINE_SIZE];
	else
		seq = new char[part_size+10];
		
	uint32 seq_size;
	int omit_next_n_kmers = 0;
	kmer_t kmer_str, kmer_rev, kmer_can;
	uint32 i;
	kmer_t kmer_mask = 0;
	uint32 bin_shift = (2*kmer_len - 2*prefix_len);
	
	all_C = 0;
	for(i = 0; i < 2*kmer_len; ++i)
	{
		kmer_mask = (kmer_mask << 1) | 1;
		if(i % 2 == 0)
			all_C |= ((uint64) 1) << i;
	}
	first_A = (1ULL << (2*kmer_len-2)) - 1;
	firstA_endG = 0ULL + 2;
	all_T = (1ULL << (2*kmer_len)) - 1;
	
	uint32 kmer_len_shift = (kmer_len-1)*2;

	//------
	while(GetSeq(seq, seq_size))
	{
		if (seq_size < kmer_len)
			continue;

		// Init k-mer
		kmer_str = 0;
		kmer_rev = 0;
		for(i = 0; i < kmer_len-1; ++i)
		{
			if(seq[i] < 0)
			{
				seq[i] = 0;
				omit_next_n_kmers = i+1;
			}
			kmer_str = (kmer_str << 2) + seq[i];
			//kmer_rev = kmer_rev + (((kmer_t) 3 - seq[i]) << (i+i));
		}
		//kmer_rev <<= 2;
		for(; i < seq_size; ++i)
		{
			if(seq[i] < 0)
			{
				seq[i] = 0;
				omit_next_n_kmers = kmer_len;
			}
			kmer_str = ((kmer_str << 2) + seq[i]) & kmer_mask;
			//kmer_rev = (kmer_rev >> 2) + (((kmer_t) 3 - seq[i]) << kmer_len_shift);
			
			if(omit_next_n_kmers > 0)
			{
				omit_next_n_kmers--;
				continue;
			}

			//if(kmer_str < kmer_rev)
			kmer_can = kmer_str;
			//else
				//kmer_can = kmer_rev;

			if(kmer_can == 0ULL) 
			{
				count_all_A++;
				continue;
			}


			if(kmer_can == all_T)
			{
				count_all_T++;
				continue;
			}


			if(kmer_can == firstA_endG)
			{
				count_all_A_end_G++;
				continue;
			}
			
			if(kmer_can == all_C)
			{
				count_all_C++;
				continue;
			}

			uint32 bin_no = (uint32) (kmer_can >> (bin_shift-2));
			//if(bin_no >= (1u << (prefix_len*2)))  // change bins
			//	bin_no &= ~3U;
			if(bins[bin_no]->PutKmer(kmer_can))
				PushBinToQueue(bin_no);
		}
	}

	putchar('.');
	fflush(stdout);

	delete[] seq;

	return true;
}

//----------------------------------------------------------------------------------
void CSplitter::GetTotal(uint64 &_n_reads)
{
	_n_reads = n_reads;
}

// ***** EOF
