/*
	This file is a part of CoMeta software.
	The homepage of the CoMeta project is http://sun.aei.polsl.pl/cometa
 
	Authors: Jolanta Kawulok and Sebastian Deorowicz
	
	Version: 0.3
	Date   : 2014-11-26

*/

#ifndef _SPLITTER2class_H
#define _SPLITTER2class_H

#include "defs.h"
#include "kmer_bin.h"
#include "queues.h"
#include <stdio.h>
#include <iostream>
#include <vector>
#include "queues.h"
using namespace std;

class CSplitter2class {
	CMemoryMonitor *mm;

	uchar *part;
	uint64 part_size, part_pos;
	CReadsQueue *q;
	
	//CBinPartQueue *q;
	//CBinDesc *bd;
	char codes[256];
	bool is_fasta;
	bool is_read;

	//uint32 kmer_len;
	//uint32 prefix_len;
	//uint32 n_bins;

	//uint64 all_C, first_A;
	//uint64 count_all_A;
	//uint64 count_all_C;

	uint64 n_reads;

	inline bool GetSeq(char *&title, char *&seq, uint64 &seq_size, uint32 &title_size);
	//inline void PushBinToQueue(int bin_no);

public:
	static uint32 MAX_LINE_SIZE;

	//CSplitter2class(CMemoryMonitor *_mm, CBinPartQueue *_q,  uint32 _buffer_size, uint32 _kmer_len, uint32 _prefix_len, uint32 _n_bins, bool _is_fasta);
	CSplitter2class(CMemoryMonitor *_mm, CReadsQueue *_q,  bool _is_fasta);
	~CSplitter2class();

	bool ProcessReads(uchar *_part, uint64 _part_size);
	void Complete();

	void GetTotal(uint64 &_n_reads);
};


//----------------------------------------------------------------------------------
inline bool CSplitter2class::GetSeq(char *&title, char *&seq, uint64 &seq_size, uint32 &title_size)
{
	uchar c;
	uint32 pos = 0;
	uint32 posT = 0;

	if(is_fasta)
	{
		// Title
		if(part_pos >= part_size)
			return false;
		c = part[part_pos++];
		if(c != '>') 
			return false;
		
		if(c == '>') 
		{
			n_reads++;
			for(; part_pos < part_size;)
			{
				c = part[part_pos++];
				if(c < 32)					// newliners
				{
					title[posT] = 0;
					break;
				}
				title[posT++] = c;
			}
		}
		else
			part_pos--;
		
		title_size = posT;



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
		
			if(c < 32)					// newliners
				continue;					

			if(c == '>')			
			{
				part_pos--;
				c = part[part_pos-1];
				break;
			}

			seq[pos++] = codes[c];
		}
		seq_size = pos;

		if (part_pos>1015870)
			int yyy=3;

		if(part_pos >= part_size)
			return true;

		if(part[part_pos++] >= 32)
			part_pos--;
		else if(part_pos >= part_size)
			return true;
	}
	
	if (n_reads > 2200 && n_reads%1 == 0)
		int il =0;

	return (c == '\n' || c == '\r'); 
}

//----------------------------------------------------------------------------------



#endif