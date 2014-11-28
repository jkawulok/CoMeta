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
#include "Tr_wrapers.h"

#ifdef _DEBUG
   #ifndef DBG_NEW
      #define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
      #define new DBG_NEW
   #endif
#endif  // _DEBUG


//************************************************************************************************************
// CWKmersReader
//************************************************************************************************************
CWKmersReader::CWKmersReader(CMemoryMonitor *_mm, string _file_name, uint64 _part_size,  CPartKMERQueue *_part_kmer, uchar &_kmer_len, uint64* _start_prefix)
{
	mm = _mm;

	file_name = _file_name;
	part_size = _part_size;
	part_kmer = _part_kmer;
	start_prefix = _start_prefix;
	

	fkm = new CKmerReader(mm,  _kmer_len, _start_prefix, _part_size);
	fkm->SetNames(file_name);
	fkm->SetPartSize(part_size);
}

//----------------------------------------------------------------------------------
CWKmersReader::~CWKmersReader()
{
}

//----------------------------------------------------------------------------------
void CWKmersReader::operator()()
{
	uchar *part;
	uint64 part_filled = 0;
	uchar _kmer_len;
	
	list<uint64>* list_num_suff;
	int size_list = 0;

	fkm->OpenFiles(_kmer_len);
	kmer_len = _kmer_len;
	long sumALL = 0; //// added
	long sum_part_filled = 0; //// added

	// Reading kmer parts
	while(fkm->GetPart(part, part_filled, list_num_suff, size_list)) 
	{
		if (size_list == 0)
			break;
		 
		for (list<uint64>::iterator si = (*list_num_suff).begin(); si != (*list_num_suff).end(); si++) /// ------------------ added
			sumALL += *si;

		sum_part_filled += part_filled; // ------------------ added

		part_kmer->push(part, part_filled, list_num_suff, size_list);
	}
	delete fkm;
	part_kmer->mark_completed();
	cout << "\nEnd read kmer file\n";
	cout << "Check:   No. of kmer: " << setw(12) << sumALL << "\n"; // ------------------ added
	cout << "Check:   sum_part_filled: " << setw(12) << sum_part_filled+9 << "\n"; // ------------------ added
	if (sumALL == 0)
	{
		cout << "\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ERROR"; 
		cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
	}
}


//************************************************************************************************************
//CWKmersSplitter
//************************************************************************************************************
CWKmersSplitter_new::CWKmersSplitter_new(CMemoryMonitor *_mm, CPartKMERQueue *_pkmer, CBinKmers_all_new *_bin_kmers, uchar _kmer_len, uint32 _prefix_len)
{
	mm = _mm;
	pkmer  = _pkmer;
	bin_kmers = _bin_kmers;

	spl = new CSplitter_KMER_new(mm, _pkmer, _bin_kmers,  _kmer_len);
}

//----------------------------------------------------------------------------------
CWKmersSplitter_new::~CWKmersSplitter_new()
{
	if (spl)
		delete spl;
	spl = NULL;
}

//----------------------------------------------------------------------------------
void CWKmersSplitter_new::operator()()
{
	// Splitting parts
	while(!pkmer->completed())
	{
		uchar *part;
		uint64 size;
		list<uint64>* list_num_suff;
		int size_list;


		if(pkmer->pop(part, size, list_num_suff, size_list))
		{
			spl->ProcessReads(part, size, list_num_suff);
			delete[] part;
		}
	}
	bin_kmers->mark_completed();
	cout << "End split kmer file\n";

}






// ***** EOF
