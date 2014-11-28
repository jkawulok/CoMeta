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
#include "wrapers.h"
#include "Tr_wrapers.h"

#ifdef _DEBUG
   #ifndef DBG_NEW
      #define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
      #define new DBG_NEW
   #endif
#endif  // _DEBUG

//************************************************************************************************************
// CWFastqReader
//************************************************************************************************************
CWFastqReader::CWFastqReader(CMemoryMonitor *_mm, list<string> _file_name, uint64 _part_size, CPartQueue *_part_queue, bool _is_fasta, bool _is_read, int _kmer) 
{
	mm = _mm;

	file_name = _file_name;
	part_size = _part_size;
	part_queue = _part_queue;

	fqr = new CFastqReader(mm, _is_fasta, _is_read, _kmer);
	fqr->SetNames(file_name);
	fqr->SetPartSize(part_size);
}

//----------------------------------------------------------------------------------
CWFastqReader::~CWFastqReader()
{
}

//----------------------------------------------------------------------------------
void CWFastqReader::operator()()
{
	uchar *part;
	char *part_seq = NULL;
	int kmer;

	uint64 part_filled;
	int nrF = 0;
	cout << "Open files\n"; 
	
	while (	fqr->OpenFiles())
	{
		printf("\nOpen file number %d\n", nrF++);
		// Reading Fastq parts
		while(fqr->GetPart(part, part_filled, part_seq, kmer) )
		{
			part_queue->push(part, part_filled);
		
			if (part_seq != NULL)
				if (strlen(part_seq) > kmer) 
				{
					part_queue->push( (uchar*)part_seq, strlen(part_seq));
					part_seq = NULL;
				}
		
		}
	}
	delete fqr;
	part_queue->mark_completed();
	cout << "\n End of the fasta file reading";	fflush(stdout);
}

//************************************************************************************************************
// CWSplitter
//************************************************************************************************************
CWSplitter::CWSplitter(CMemoryMonitor *_mm, CPartQueue *_pq, CBinPartQueue *_bpq, CBinDesc *_bd, uint32 _buffer_size, int _kmer_len, int _prefix_len, int _n_bins, bool _is_fasta, bool _is_read)
{
	mm = _mm;

	pq  = _pq;
	bpq = _bpq;
	bd  = _bd;
	spl = new CSplitter(mm, bpq, bd, _buffer_size, _kmer_len, _prefix_len, _n_bins, _is_fasta, _is_read);
}

//----------------------------------------------------------------------------------
CWSplitter::~CWSplitter()
{
	if (spl)
		delete spl;
}

//----------------------------------------------------------------------------------
void CWSplitter::operator()()
{
	// Splitting parts
	while(!pq->completed())
	{
		uchar *part;
		uint64 size;
		if(pq->pop(part, size))
		{
			spl->ProcessReads(part, size);
			delete[] part;
		}
	}
	spl->Complete(); // addef  AAA, CCC, //build bins (bin) in memory
	bpq->mark_completed();
}

//----------------------------------------------------------------------------------
void CWSplitter::GetTotal(uint64 &_n_reads)
{
	spl->GetTotal(_n_reads);
}


//************************************************************************************************************
// CWSplitter2class
//************************************************************************************************************
//CWSplitter2class::CWSplitter2class(CMemoryMonitor *_mm, CPartQueue *_pq, CBinPartQueue *_bpq, uint32 _buffer_size, int _kmer_len, int _prefix_len, int _n_bins, bool _is_fasta)
CWSplitter2class::CWSplitter2class(CMemoryMonitor *_mm, CPartQueue *_pq, CReadsQueue *_pread,  bool _is_fasta)
{
	mm = _mm;

	pq  = _pq;
	pread= _pread;

	spl = new CSplitter2class(mm, pread,  _is_fasta);

}

//----------------------------------------------------------------------------------
CWSplitter2class::~CWSplitter2class()
{
	delete spl;
}

//----------------------------------------------------------------------------------
void CWSplitter2class::operator()()
{
	// Splitting parts to reads
	while(!pq->completed())
	{
		uchar *part;
		uint64 size;
		if(pq->pop(part, size))
		{
			spl->ProcessReads(part, size);
			delete[] part;
		}
	}
	pread->mark_completed();
	cout << "\n End of the 'reads' spliting";	fflush(stdout);
}

//----------------------------------------------------------------------------------
void CWSplitter2class::GetTotal(uint64 &_n_reads)
{
	spl->GetTotal(_n_reads);
}



//************************************************************************************************************
// CWKmerBinStorer
//************************************************************************************************************
CWKmerBinStorer::CWKmerBinStorer(CMemoryMonitor *_mm, int _n_bins, CBinPartQueue *_q, CBinDesc *_bd, vector<string> _working_directories, 
	uint32 _prefix_len, uint32 _kmer_len, uint64 _max_mem_buffer)
{
	mm = _mm;

	n_bins = _n_bins;
	q = _q;
	bd = _bd;
	working_directories = _working_directories;
	prefix_len = _prefix_len;
	kmer_len = _kmer_len;

	kbs = new CKmerBinStorer(mm, n_bins, q, bd, working_directories, prefix_len, kmer_len, _max_mem_buffer);
	kbs->OpenFiles(); // tworzenie koszy (bin) na frag k-mer o danych prefixach dl p1 - moj kom
}

//----------------------------------------------------------------------------------
CWKmerBinStorer::~CWKmerBinStorer()
{
	kbs->CloseFiles();
		delete kbs;
}

//----------------------------------------------------------------------------------
void CWKmerBinStorer::operator()()
{
	kbs->ProcessQueue();
}


//************************************************************************************************************
// CWKmerBinReader
//************************************************************************************************************
CWKmerBinReader::CWKmerBinReader(CMemoryMonitor *_mm, CDiskMonitor *_dm, CBinOrdering *_bo, CBinDesc *_bd, CBinQueue *_bq)
{
	mm = _mm;
	dm = _dm;
	bo = _bo;

	bd = _bd;
	bq = _bq;

	kbr = new CKmerBinReader(mm, dm, bo, bd, bq);
}

//----------------------------------------------------------------------------------
CWKmerBinReader::~CWKmerBinReader()
{
	delete kbr;
}

//----------------------------------------------------------------------------------
void CWKmerBinReader::operator()()
{
	kbr->ProcessBins();
}


//************************************************************************************************************
// CWKmerBinSorter
//************************************************************************************************************
CWKmerBinSorter::CWKmerBinSorter(CMemoryMonitor *_mm, int32 _n_bins, CBinDesc *_bd, CBinQueue *_bq, CKmerQueue *_kq, int _n_omp_threads, uint32 _kmer_len)
{
	mm = _mm;

	n_bins = _n_bins;
	bd = _bd;
	bq = _bq;
	kq = _kq;
	kmer_len = _kmer_len;

	kbs = new CKmerBinSorter(mm, n_bins, bd, bq, kq, _n_omp_threads, kmer_len);
}

//----------------------------------------------------------------------------------
CWKmerBinSorter::~CWKmerBinSorter()
{
	delete kbs;
}

//----------------------------------------------------------------------------------
void CWKmerBinSorter::operator()()
{
	kbs->ProcessBins();
}

//************************************************************************************************************
// CWKmerBinCompleter
//************************************************************************************************************
CWKmerBinCompleter::CWKmerBinCompleter(CMemoryMonitor *_mm, string _file_name, CBinDesc *_bd, CKmerQueue *_kq, uint32 _kmer_len)
{
	mm = _mm;

	file_name = _file_name;
	kq		  = _kq;
	bd		  = _bd;
	kmer_len  = _kmer_len;

	kbc = new CKmerBinCompleter(mm, file_name, bd, kq, kmer_len);
}

//----------------------------------------------------------------------------------
CWKmerBinCompleter::~CWKmerBinCompleter()
{
	if (kbc)
		delete kbc;
}

//----------------------------------------------------------------------------------
void CWKmerBinCompleter::operator()()
{
	kbc->ProcessBins();
}

//----------------------------------------------------------------------------------
void CWKmerBinCompleter::GetTotal(uint64 &_n_unique, uint64 &_n_singletons, uint64 &_n_total)
{
	if(kbc)
		kbc->GetTotal(_n_unique, _n_singletons, _n_total);
}

//************************************************************************************************************
//CWClassifier_new
//************************************************************************************************************
CWClassifier_new::CWClassifier_new(CMemoryMonitor *_mm, CReadsQueue *_pread, CBinKmers_all_new *_bin_kmers, uchar _kmer_len, uchar _step_k, double _matchcutoff, FILE* _outMatch, FILE* _outMisMatch)
{
	mm = _mm;
	pread  = _pread;
	bin_kmers = _bin_kmers;
	matchcutoff = _matchcutoff;
	
	clr = new CClassifier_new(mm, bin_kmers, pread, _kmer_len, _step_k, matchcutoff, _outMatch, _outMisMatch);

}

CWClassifier_new::~CWClassifier_new()
{
	if(clr)
		delete clr;
	clr = NULL;
}

void CWClassifier_new::operator()()
{
	while(!pread->completed())
	{
		char *name; 
		char *seq;
		uint64 size; 
		uint32 size_name;

		if(pread->pop(name, seq, size, size_name) )
		{
			clr->ProcessClass(name, seq, size);
			
			delete[] name;
			delete[] seq;
		}
	}
	cout << "\n End of the classifing";	fflush(stdout);
}

void CWClassifier_new::GetTotal(uint64 &_num_alig_reads, uint64 &_num_REValig_reads, uint64 &_num_NOalig_reads)
{
	if(clr)
		clr->GetTotal(_num_alig_reads, _num_REValig_reads,_num_NOalig_reads);
}



// ***** EOF
