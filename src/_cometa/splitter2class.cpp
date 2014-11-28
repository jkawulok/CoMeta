/*
	This file is a part of CoMeta software.
	The homepage of the CoMeta project is http://sun.aei.polsl.pl/cometa
 
	Authors: Jolanta Kawulok and Sebastian Deorowicz
	
	Version: 0.3
	Date   : 2014-11-26
	
*/

#include "stdafx.h"

#include <algorithm>
#include "defs.h"
#include "splitter2class.h"
 
uint32 CSplitter2class::MAX_LINE_SIZE = (1 << 14);  

//----------------------------------------------------------------------------------
CSplitter2class::CSplitter2class(CMemoryMonitor *_mm, CReadsQueue *_q,  bool _is_fasta)
{
	mm = _mm;
	is_fasta = _is_fasta;

	q  = _q;

	part = NULL;

	for(int i = 0; i < 256; ++i)
		codes[i] = -1;
	codes['A'] = codes['a'] = 0;
	codes['C'] = codes['c'] = 1;
	codes['G'] = codes['g'] = 2;
	codes['T'] = codes['t'] = 3;

	uint32 count_len;
	n_reads = 0;
}

//----------------------------------------------------------------------------------
CSplitter2class::~CSplitter2class()
{

}


//----------------------------------------------------------------------------------
bool CSplitter2class::ProcessReads(uchar *_part, uint64 _part_size)
{
	part	  = _part;
	part_size = _part_size;
	part_pos  = 0;
	
	char *seq = new char[MAX_LINE_SIZE];
	char *title = new char[MAX_LINE_SIZE/3];
	
	

	uint64 seq_size = 1;
	uint32 title_size;
	uint32 i;

	while(GetSeq(title, seq, seq_size, title_size))
	{
		//cout << title << "\n";
		q->push(title, seq, seq_size, title_size);
		title = new char[MAX_LINE_SIZE/3];
		seq = new char[MAX_LINE_SIZE];
	}
	
	//cout << part_size << "\n";
	//cout << seq_size << "\n";
	putchar('-');	
	fflush(stdout);

	delete[] seq;
	delete[] title;

	return true;
}

//----------------------------------------------------------------------------------
void CSplitter2class::GetTotal(uint64 &_n_reads)
{
	_n_reads = n_reads;
}

// ***** EOF
