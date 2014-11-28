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
#include "Tr_KMclass.h"
#include "asmlib.h"
#include <boost/filesystem.hpp>

#ifdef _DEBUG
   #ifndef DBG_NEW
      #define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
      #define new DBG_NEW
   #endif
#endif  // _DEBUG



//----------------------------------------------------------------------------------
TrKMclass::TrKMclass()
{
#if !defined(_OPENMP)
	BOOST_STATIC_ASSERT_MSG(false, "You need to use OpenMP");
#endif

	initialized = false;
	kmer_len = 0;
	n_splitters = 1;

}

//----------------------------------------------------------------------------------
TrKMclass::~TrKMclass()
{

}


//----------------------------------------------------------------------------------
void TrKMclass::SetFileNames(string _data_file_name, string _output_file_name)
{
	data_file_name		= _data_file_name;
	output_file_name    = _output_file_name;

	initialized = true;


	SetMemcpyCacheLimit(1);
}

//----------------------------------------------------------------------------------
void TrKMclass::SetParams(int _n_threads, uint32 _max_mem_size_in_gb)
{
	kmer_len = 0;

	prefix_len    = 4;
	n_bins        = (1 << (2 * prefix_len));
	bin_part_size = 1 << 18; //??
	fastq_buffer_size = 1 << 20; //?

	if (_max_mem_size_in_gb < 40)
		kmers_buffer_size = 1 << 29; //??
	else
		kmers_buffer_size = 1 << 31; //??

	AdjustToHardware(_n_threads);

	max_mem_size  = NORM(((uint64) _max_mem_size_in_gb) << 30, 2ull << 30, 1024ull << 30);

	initialized = true; 
}

//----------------------------------------------------------------------------------
void TrKMclass::AdjustToHardware(uint32 cores)
{
	if(!cores)
		cores = boost::thread::hardware_concurrency();
	cout << "n_threads " << cores << "\n";

	n_cores = cores;
}

//----------------------------------------------------------------------------------
bool TrKMclass::AdjustMemoryLimits()
{
	FILE* pFile = NULL;
	if((pFile = fopen(data_file_name.c_str(), "rb")) == NULL)
	{
		perror("fopen");
		cout << data_file_name.c_str() << "\n";
		return false;
	}
	fread(&kmer_len, 1, 1, pFile);
	
	
	#ifdef WIN32
		_fseeki64(pFile, 0L, SEEK_END); //Windows
		uint64 len_all = _ftelli64(pFile) - 8;// windows
		_fseeki64(pFile, len_all, SEEK_SET); // windows		
	#else
		fseeko(pFile, 0L, SEEK_END);   // linux
		uint64 len_all = ftello(pFile) - 8; // linux
		fseeko(pFile, len_all, SEEK_SET); // linux	 
	#endif
	
	//max_mem_binsKMER = ftell(pFile) + 256*256*8 + 100; // pamiec dla struktury z wszystkimi kmerami wczytanymi
	 
	num_kmer=0;
	fread(&num_kmer, 1, sizeof(uint64), pFile); 
	max_mem_binsKMER = 256 * (1 + 8 + 8*256 + 8) + 8 * num_kmer;
	cout << "   No. of kmer: " << setw(12) << num_kmer << "\n";

	fclose(pFile);
	
	if (max_mem_size <= max_mem_binsKMER)
	{
		cout << "max_mem_size <= max_mem_binsKMER\n";
		return false;
	}
	uint64 m_rest = max_mem_size - max_mem_binsKMER;
	max_mem_seqKMER   = m_rest; //wczytywane sekwencje z kmerami
	

	//part 2
	m_rest = max_mem_size - max_mem_binsKMER;


#ifdef DEBUG_MODE	
	cout << "Splitters: " << m_splitters << "\n";
	//cout << "Storer   : " << max_mem_storer << "\n";
	cout << "FASTQ    : " << max_mem_fastq << "\n";
	cout << "READS    : " << max_mem_reads << "\n";
	//cout << "Bin part : " << max_mem_bin_part << "\n";
	cout << "Stage 2  : " << max_mem_stage2 << "\n";
#endif


	//if(max_mem_fastq < 2 * (uint32)fastq_buffer_size)
	//{
	//	cout << "max_mem_fastq < 2 * (uint32)fastq_buffer_size \n";
	//	return false;
	//}

	
	return true;
}

bool TrKMclass::New_kmer_len(CWKmersReader *w_kmers)
{
	kmer_len = w_kmers->get_kmer_len();
	return true;
}


void TrKMclass::SaveOutFile()
{

	outDatabase	= fopen(output_file_name.c_str(), "wb");


	if(!outDatabase)
	{
		cout << "Error: Cannot create " << output_file_name << "\n";
		return;
	}

	fwrite(&kmer_len, 1, 1, outDatabase);	
	bin_kmers->save_outKMERS(outDatabase);
}


//void TrKMclass::GetStats(double &time1, double &time2, uint64 &_num_alig_reads, uint64 &_num_REValig_reads, uint64 &_num_NOalig_reads, uint64 &_n_reads, uchar &_kmer_len, double &_matchcutoff, uint32 &_n_cores, uint32 &_n_splitters, uint32 &_n_classification, uint32 &_max_mem_size_in_gb)
void TrKMclass::GetStats()
{
	fclose(outDatabase);


}


//----------------------------------------------------------------------------------
bool TrKMclass::Process()
{
	if(!initialized)
		return false;

	if(!AdjustMemoryLimits())
		return false;
		
	//w1.startTimer();

	//// Create monitors
	mm = new CMemoryMonitor(max_mem_size);
	//bo = new CBinOrdering(n_bins);

	//// Create queues
	int n_splitters_read = 1;
	pkmer  = new CPartKMERQueue(max_mem_seqKMER); // wczytany plik z kmerami	


		

	// ---------------------------------- kmers -------------

	//new -----------------
	bin_kmers = new CBinKmers_all_new(n_splitters, kmer_len, num_kmer, max_mem_binsKMER);
	w_splittersKMER.resize(n_splitters); 
	for(int i = 0; i < n_splitters; ++i)
	{
		w_splittersKMER[i] = new CWKmersSplitter_new(mm, pkmer, bin_kmers, kmer_len, prefix_len); //wstawianie kmerow do koszy
		gr2.create_thread(boost::ref(*w_splittersKMER[i]));
	}
	

	w_kmers = new CWKmersReader(mm, data_file_name, kmers_buffer_size, pkmer, kmer_len, bin_kmers->fun_start_prefix());
	gr1.create_thread(boost::ref(*w_kmers));
	
	gr1.join_all(); //  read kmers
	gr2.join_all();	//	split kmers
	
	delete w_kmers;
	for(int i = 0; i < n_splitters; ++i)
	{
		delete w_splittersKMER[i];
	}

	//w1.stopTimer();
	cout << "\n";

	//w2.startTimer();
	SaveOutFile();


	
	// ----------------------------------------------------
	//delete pread;
	delete pkmer;
	delete bin_kmers;

	//delete bpq;
	//delete bd;
	//delete bq;	
	//delete kq;

	//w2.stopTimer();

	return true;
}