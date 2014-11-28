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

#ifndef _SPLITTER_H
#define _SPLITTER_H

#include "defs.h"
#include "kmer_bin.h"
#include "queues.h"
#include <stdio.h>
#include <iostream>
#include <vector>
#include "queues.h"
using namespace std;

class CSplitter {
	CMemoryMonitor *mm;

	uchar *part;
	uint64 part_size, part_pos;
	CKmerBinCollector **bins;
	CBinPartQueue *q;
	CBinDesc *bd;
	char codes[256];
	bool is_fasta;
	bool is_read;

	uint32 kmer_len;
	uint32 prefix_len;
	uint32 n_bins;

	uint64 all_C, first_A, firstA_endG, all_T;
	uint64 count_all_A;
	uint64 count_all_T;
	uint64 count_all_C;
	uint64 count_all_A_end_G;

	uint64 n_reads;

	inline bool GetSeq(char *&seq, uint32 &seq_size);
	inline void PushBinToQueue(int bin_no);

public:
	static uint32 MAX_LINE_SIZE;

	CSplitter(CMemoryMonitor *_mm, CBinPartQueue *_q, CBinDesc *_bd, uint32 _buffer_size, uint32 _kmer_len, uint32 _prefix_len, uint32 _n_bins, bool _is_fasta, bool _is_read);
	~CSplitter();

	bool ProcessReads(uchar *_part, uint64 _part_size);
	void Complete();

	void GetTotal(uint64 &_n_reads);
};


//----------------------------------------------------------------------------------
inline bool CSplitter::GetSeq(char *&seq, uint32 &seq_size)
{
	uchar c;
	uint32 pos = 0;

	if(is_fasta)
	{
		// Title
		if(part_pos >= part_size)
			return false;
		c = part[part_pos++];
		//if(c != '>') // 
			//return false;
		
		if(c == '>') // 
		{
			n_reads++;
			for(; part_pos < part_size;)
			{
				c = part[part_pos++];
				if(c < 32)					// newliners
					break;
			}
		}
		else
			part_pos--;
		
		if(part_pos >= part_size)
			return false;

		c = part[part_pos++];
		if(c >= 32)
			part_pos--;
		else if(part_pos >= part_size)
			return false;

		// Sequence
		for(; part_pos < part_size;)
		{
			c = part[part_pos++];
		
			if (c=='N' || c=='n')
				int ui = 0;
			
			if(c < 32)					// newliners
				continue;					

			if(c == '>')			
			{
				part_pos--;
				c = part[part_pos-1];
				break;
			}

			
			
			if (c=='N' || c=='n')
			{
				if( part[part_pos] != 'N' && part[part_pos] != 'n' && part[part_pos] != '\n' )
					break;
				else
					continue;
			}


			seq[pos++] = codes[c];
		}
		seq_size = pos;

		if(part_pos >= part_size)
			return true;

		if(part[part_pos++] >= 32)
			part_pos--;
		else if(part_pos >= part_size)
			return true;
	}
	

	return (c == '\n' || c == '\r' || c == 'N'  || c == 'n'); 
}

//----------------------------------------------------------------------------------
inline void CSplitter::PushBinToQueue(int bin_no)
{
	uchar *bin_part;
	uint64 bin_part_size;
	uint64 n_rec;

	if(bins[bin_no]->GetBinPart(bin_part, bin_part_size, n_rec))
	{
		q->push(bin_no, bin_part, bin_part_size);
		bd->insert(bin_no, "", 0, bin_part_size, n_rec); // n_rec - number in bin
	}
}

#endif
