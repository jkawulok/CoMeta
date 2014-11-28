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
#include "cometa.h"
#include "asmlib.h"
#include <boost/filesystem.hpp>


//----------------------------------------------------------------------------------
CCOMETA::CCOMETA()
{
#if !defined(_OPENMP)
	BOOST_STATIC_ASSERT_MSG(false, "You need to use OpenMP");
#endif

	initialized = false;
	kmer_len = 0;
	n_splitters = 1;
	n_sorters = 1;
	n_omp_threads = 1;
}

//----------------------------------------------------------------------------------
CCOMETA::~CCOMETA()
{

}

void CCOMETA::SetFileNamesList(string file_name_FILEsequence, list<string>& file_name_sequence, string path_seq)
{
	FILE *in;
	if((in = fopen(file_name_FILEsequence.c_str(), "rt")) == NULL)
	{
		perror("fopen");
		cout << file_name_FILEsequence.c_str() << " does not exist \n";		
	}

	char* pszCursor = new char[10000];
	char* pszCursorADD = new char[10000];
	

	while(fscanf(in,  "%[^\n]", pszCursor) != EOF) //whole line
	{
		file_name_sequence.push_back( path_seq+string(pszCursor));
		fscanf(in,  "%[\n]", pszCursorADD);
	}

	delete[] pszCursor;
	delete[] pszCursorADD;


}


//----------------------------------------------------------------------------------
void CCOMETA::SetFileNames(list<string> _input_file_name, string _output_file_name, vector<string> _working_directories, bool _is_fasta, bool _is_read)
{
	input_file_name     = _input_file_name;
	output_file_name    = _output_file_name;
	working_directories = _working_directories;
	n_disks				= (uint32) working_directories.size();
	is_fasta			= _is_fasta;
	is_read				= _is_read;

	initialized = true;

	SetMemcpyCacheLimit(1);
}



//----------------------------------------------------------------------------------
void CCOMETA::SetFileNames(list<string> _input_file_name, string _output_file_name, string _working_directory, bool _is_fasta, bool _is_read)
{
	vector<string> vs;
	vs.push_back(_working_directory);
	SetFileNames(_input_file_name, _output_file_name, vs, _is_fasta, _is_read);
}

//----------------------------------------------------------------------------------
void CCOMETA::SetParams(int _kmer_len, int _n_splitters, int _n_sorters, int _n_omp_threads, uint32 _max_mem_size_in_gb)
{
	kmer_len = NORM(_kmer_len, 1, 32);

	prefix_len    = 3;  //
	n_bins        = (1 << (2 * prefix_len)) * 4;
	bin_part_size = 1 << 16;
	fastq_buffer_size = 1 << 22;

	n_splitters   = NORM(_n_splitters, 1, 32);
	n_sorters     = NORM(_n_sorters, 1, 32);
	n_omp_threads = NORM(_n_omp_threads, 1, 32);
	max_mem_size  = NORM(((uint64) _max_mem_size_in_gb) << 30, 2ull << 30, 1024ull << 30);
	//max_mem_size *= (1.4/2.0); /////////////////////////////////////////////////////////////////////////////

	initialized = true; 
}

//----------------------------------------------------------------------------------
void CCOMETA::SetParams(int _kmer_len, int _n_threads, uint32 _max_mem_size_in_gb)
{
	kmer_len = NORM(_kmer_len, 1, 32);

	prefix_len    = 3; // 
	n_bins        = (1 << (2 * prefix_len)) * 4;
	bin_part_size = 1 << 16;
	fastq_buffer_size = 1 << 22;

	AdjustToHardware(_n_threads);

	max_mem_size  = NORM(((uint64) _max_mem_size_in_gb) << 30, 2ull << 30, 1024ull << 30);

	initialized = true; 
}

//----------------------------------------------------------------------------------
void CCOMETA::AdjustToHardware(uint32 cores)
{
	if(!cores)
		cores = boost::thread::hardware_concurrency();

	if(cores < 6)
	{
		n_splitters   = NORM(cores-2, 1, 8);
		n_sorters     = 1;
		n_omp_threads = NORM(cores-2, 1, 8);
	}
	if(cores <= 8)
	{
		n_splitters   = NORM(cores-2, 1, 8);
		n_sorters     = 2;
		n_omp_threads = NORM(cores/2, 1, 8);
	}
	else
	{
		n_splitters   = NORM(cores/2, 8, 16);
		n_sorters     = cores/8;
		n_omp_threads = cores/n_sorters;
	}	
	n_cores = cores;
}

//----------------------------------------------------------------------------------
bool CCOMETA::AdjustMemoryLimits()
{
	uint64 m_splitters = n_bins * 7 / 16 * bin_part_size * sizeof(uint64) * n_splitters;
	uint64 m_rest = max_mem_size - m_splitters;
	
	if(max_mem_size >= 16ull << 30)
		max_mem_storer = 2 * m_rest / 3;
	else
		max_mem_storer    = m_rest / 2;
	m_rest			 -= max_mem_storer;
	max_mem_fastq     = m_rest * 1 / 4;
	m_rest			 -= max_mem_fastq;
	max_mem_bin_part  = m_rest;
	max_mem_stage2	  = max_mem_size;

#ifdef DEBUG_MODE	
	cout << "Splitters: " << m_splitters << "\n";
	cout << "Storer   : " << max_mem_storer << "\n";
	cout << "FASTQ    : " << max_mem_fastq << "\n";
	cout << "Bin part : " << max_mem_bin_part << "\n";
	cout << "Stage 2  : " << max_mem_stage2 << "\n";
#endif

	if(max_mem_storer < (1 << 28))
		return false;
	if(max_mem_fastq < 3 * (uint32)fastq_buffer_size)
		return false;
	if(max_mem_bin_part < bin_part_size * sizeof(uint64) * 2)
		return false;

	return true;
}

//----------------------------------------------------------------------------------
bool CCOMETA::Process()
{
	if(!initialized)
		return false;

	if(!AdjustMemoryLimits())
		return false;
		
	w1.startTimer();

	// Create monitors
	mm = new CMemoryMonitor(max_mem_stage2);
	dm = new CDiskMonitor(n_disks);
	bo = new CBinOrdering(n_bins);

	// Create queues
	pq  = new CPartQueue(max_mem_fastq);
	bpq = new CBinPartQueue(n_splitters, max_mem_bin_part);
	bd  = new CBinDesc;
	bq  = new CBinQueue(n_disks);
	kq  = new CKmerQueue(n_bins, n_sorters);

	// ***** Stage 1 *****

	w_splitters.resize(n_splitters);
	for(int i = 0; i < n_splitters; ++i)
	{
		w_splitters[i] = new CWSplitter(mm, pq, bpq, bd, bin_part_size, kmer_len, prefix_len, n_bins, is_fasta, is_read);
		gr2.create_thread(boost::ref(*w_splitters[i]));
	}

	w_storer = new CWKmerBinStorer(mm, n_bins, bpq, bd, working_directories, prefix_len, kmer_len, max_mem_storer); // tworzenie koszy (bin) 
	gr3.create_thread(boost::ref(*w_storer));

	w_fastq = new CWFastqReader(mm, input_file_name, fastq_buffer_size, pq, is_fasta, is_read, kmer_len);
	gr1.create_thread(boost::ref(*w_fastq));

	gr1.join_all(); // read sequence	
	gr2.join_all(); // split sequence
	gr3.join_all(); // insert into bins

	n_reads = 0;

	delete w_fastq;
	for(int i = 0; i < n_splitters; ++i)
	{
		uint64 _n_reads;
		w_splitters[i]->GetTotal(_n_reads);
		n_reads += _n_reads;
		delete w_splitters[i];
	}

	delete w_storer;
	w1.stopTimer();
	w2.startTimer();	
	

	// ***** End of Stage 1 *****
	
	// ***** Stage 2 *****
	w_readers.resize(n_disks);
	for(int i = 0; i < n_disks; ++i)
	{
		w_readers[i] = new CWKmerBinReader(mm, dm, bo, bd, bq);
		gr4.create_thread(boost::ref(*w_readers[i]));
	}

	//---------------
	w_sorters.resize(n_sorters);
	for(int i = 0; i < n_sorters; ++i)
	{
		w_sorters[i] = new CWKmerBinSorter(mm, n_bins, bd, bq, kq, n_omp_threads, kmer_len);
		gr5.create_thread(boost::ref(*w_sorters[i]));
	}

	//----------------
	w_completer = new CWKmerBinCompleter(mm, output_file_name, bd, kq, kmer_len);
	gr6.create_thread(boost::ref(*w_completer));


	//--------------------------------------------------------------------------------
	gr4.join_all(); // read bins
	gr5.join_all();  // sort
	gr6.join_all();

	// ***** End of Stage 2 *****

	w_completer->GetTotal(n_unique, n_singletons, n_total);

	for(int i = 0; i < n_disks; ++i)
		delete w_readers[i];
	for(int i = 0; i < n_sorters; ++i)
		delete w_sorters[i];
	delete w_completer;

	delete mm;
	delete dm;
	delete bo;

	// ***** Removing temporary files *****
	int32 bin_id;
	string name;
	uint32 c_disk;
	uint64 size;
	uint64 n_rec;

	bd->reset_reading();
	while((bin_id = bd->get_next_bin()) >= 0)
	{
		bd->read(bin_id, name, c_disk, size, n_rec);
		boost::filesystem::remove(boost::filesystem::path(name));
		tmp_size += size;
	}

	delete bd;
	delete bq;
	delete pq;
	delete bpq;
	delete kq;

	w2.stopTimer();

	return true;
}

//----------------------------------------------------------------------------------
void CCOMETA::GetStats(double &time1, double &time2, uint64 &_n_unique, uint64 &_n_singletons, uint64 &_n_total, uint64 &_n_reads, uint64 &_tmp_size, uchar &_kmer_len, uint32 &_n_cores, uint32 &_n_splitters, uint32 &_n_sorters, uint32 &_n_omp_th, uint32 &_max_mem_size_in_gb)
{
	time1 = w1.getElapsedTime();
	time2 = w2.getElapsedTime();

	_n_unique	  = n_unique;
	_n_singletons = n_singletons;
	_n_total      = n_total;
	_n_reads      = n_reads;
	_tmp_size     = tmp_size;
	_kmer_len     = kmer_len;
	_n_cores     = n_cores;
	_n_splitters = n_splitters;
	_n_sorters = n_sorters;
	_n_omp_th = n_omp_threads; //  number of threads per single sorter\n";
	_max_mem_size_in_gb = max_mem_size >> 30; 

}


// ***** EOF
