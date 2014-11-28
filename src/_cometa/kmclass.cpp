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
#include "kmclass.h"
#include "asmlib.h"
#include <boost/filesystem.hpp>

#ifdef _DEBUG
   #ifndef DBG_NEW
      #define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
      #define new DBG_NEW
   #endif
#endif  // _DEBUG



//----------------------------------------------------------------------------------
CKMclass::CKMclass()
{
#if !defined(_OPENMP)
	BOOST_STATIC_ASSERT_MSG(false, "You need to use OpenMP");
#endif

	initialized = false;
	kmer_len = 0;
	step_k = 0;
	n_splitters = 1;
	n_classification = 1;

}

//----------------------------------------------------------------------------------
CKMclass::~CKMclass()
{

}


//----------------------------------------------------------------------------------
void CKMclass::SetFileNames(list<string> _read_file_name, string _data_file_name, string _output_file_nameM, string _output_file_nameMM, double _matchcutoff, bool _p_sortKmer)
{
	read_file_name      = _read_file_name;
	data_file_name		= _data_file_name;
	output_file_nameMatch    = _output_file_nameM;
	output_file_nameMisMatch    = _output_file_nameMM;
	b_sortKmer = _p_sortKmer;

	matchcutoff			= _matchcutoff;

	initialized = true;
	is_fasta = true;

	SetMemcpyCacheLimit(1);
}

//----------------------------------------------------------------------------------
void CKMclass::SetParams(int _p_stepk, int _n_splitters,  int _n_classification, uint32 _max_mem_size_in_gb)
{
	kmer_len = 0;
	step_k = _p_stepk;
	n_cores = 0;

	prefix_len    = 4;
	n_bins        = (1 << (2 * prefix_len)) * 4;
	bin_part_size = 1 << 15; // ??
	fastq_buffer_size = 1 << 18; //??
	
	if (_max_mem_size_in_gb < 40)
		kmers_buffer_size = 1 << 29; //??
	else
		kmers_buffer_size = 1 << 31; //??

	n_splitters		 = NORM(_n_splitters, 1, 32);
	n_classification = NORM(_n_classification, 1, 32);

	max_mem_size  = NORM(((uint64) _max_mem_size_in_gb) << 30, 2ull << 30, 1024ull << 30);

	initialized = true; 
}
//----------------------------------------------------------------------------------
void CKMclass::SetParams(int _p_stepk, int _n_threads, uint32 _max_mem_size_in_gb)
{
	kmer_len = 0;
	step_k = _p_stepk;

	prefix_len    = 4;
	n_bins        = (1 << (2 * prefix_len));
	bin_part_size = 1 << 15; //??
	fastq_buffer_size = 1 << 18; //?

	if (_max_mem_size_in_gb < 30)
		kmers_buffer_size = 1 << 29; //??
	else
		kmers_buffer_size = 1 << 31; //??

	AdjustToHardware(_n_threads);

	max_mem_size  = NORM(((uint64) _max_mem_size_in_gb) << 30, 2ull << 30, 1024ull << 30);

	initialized = true; 
}

//----------------------------------------------------------------------------------
void CKMclass::AdjustToHardware(uint32 cores)
{
	if(!cores)
		cores = boost::thread::hardware_concurrency();

	//if(cores <= 4)
	//{
	//	n_splitters			= NORM(cores/2, 1, 8);
	//	n_classification	= NORM(cores/2, 1, 8);
	//}
	//else
	//{
		n_splitters			= NORM(cores-1, 1, 16);

		if (cores>1)
			n_classification	= NORM(cores-2, 1, 16);
	//}
	
	n_cores = cores;
}

//----------------------------------------------------------------------------------
bool CKMclass::AdjustMemoryLimits()
{
	FILE* pFile = NULL;
	if((pFile = fopen(data_file_name.c_str(), "rb")) == NULL)
	{
		perror("fopen");
		cout << data_file_name.c_str() << " does not exist \n";
		return false;
	}
	fread(&kmer_len, 1, 1, pFile);

	if (step_k == 0)
		step_k = kmer_len;

	num_kmer=0;
	if (b_sortKmer)
	{
		fread(&num_kmer, 1, sizeof(uint64), pFile); 
	}
	else
	{	
		#ifdef WIN32
			_fseeki64(pFile, 0L, SEEK_END); //Windows
			uint64 len_all = _ftelli64(pFile) - 8;// windows
			_fseeki64(pFile, len_all, SEEK_SET); // windows		
		#else
			fseeko(pFile, 0L, SEEK_END);   // linux
			uint64 len_all = ftello(pFile) - 8; // linux
			fseeko(pFile, len_all, SEEK_SET); // linux	 
		#endif
	 
		
		fread(&num_kmer, 1, sizeof(uint64), pFile); 
	}


	max_mem_binsKMER = 256 * (1 + 8 + 8*256 + 8) + 8 * num_kmer;
	cout << "   No. of kmer: " << setw(12) << num_kmer << "\n";

	fclose(pFile);
	
	if (max_mem_size <= max_mem_binsKMER)
		return false;

	uint64 m_rest = max_mem_size - max_mem_binsKMER;
	max_mem_seqKMER   = m_rest; //loaded sequnces with k-mers
	

	//part 2
	m_rest = max_mem_size - max_mem_binsKMER;

	if (m_rest > 5 * fastq_buffer_size)
	{
		max_mem_fastq = 3 * fastq_buffer_size;
		max_mem_reads = 2 * fastq_buffer_size;

	}
	else
	{
		max_mem_fastq   = m_rest * 2 / 3;  //loaded sequnces with reads 
		max_mem_reads	= m_rest * 1 / 3;  //memory for 'reads'
	}
		

#ifdef DEBUG_MODE	
	cout << "Splitters: " << m_splitters << "\n";
	//cout << "Storer   : " << max_mem_storer << "\n";
	cout << "FASTQ    : " << max_mem_fastq << "\n";
	cout << "READS    : " << max_mem_reads << "\n";
	//cout << "Bin part : " << max_mem_bin_part << "\n";
	cout << "Stage 2  : " << max_mem_stage2 << "\n";
#endif

	//if(max_mem_storer < (1 << 28))
	//	return false;
	//if(max_mem_fastq < 3 * (uint32)fastq_buffer_size)
	//	return false;

	if(max_mem_fastq < 2 * (uint32)fastq_buffer_size)
		return false;

	//if(max_mem_bin_part < bin_part_size * sizeof(uint64) * 2)
		//return false;

	return true;
}

bool CKMclass::New_kmer_len(CWKmersReader *w_kmers)
{
	kmer_len = w_kmers->get_kmer_len();
	return true;
}


void CKMclass::OpenOutFile()
{

	outMatch	= fopen(output_file_nameMatch.c_str(), "wt");
	outMisMatch	= fopen(output_file_nameMisMatch.c_str(), "wt");

	if(!outMatch)
		cout << "Error: Cannot create " << output_file_nameMatch << "\n";

	if(!outMisMatch)
		cout << "Error: Cannot create " << output_file_nameMisMatch << "\n";

}


void CKMclass::GetStats(double &time1, double &time2, uint64 &_num_alig_reads, uint64 &_num_REValig_reads, uint64 &_num_NOalig_reads, uint64 &_n_reads, uchar &_kmer_len, double &_matchcutoff, uint32 &_n_cores, uint32 &_n_splitters, uint32 &_n_classification, uint32 &_max_mem_size_in_gb)
{
	fclose(outMatch);
	fclose(outMisMatch);

	time1 = w1.getElapsedTime();
	time2 = w2.getElapsedTime();

	_num_alig_reads	  = num_alig_reads;
	_num_REValig_reads = num_REValig_reads;
	_num_NOalig_reads      = num_NOalig_reads;
	_n_reads      = n_reads;
	_kmer_len = kmer_len;
	_matchcutoff = matchcutoff;
	_n_cores = n_cores;
	_n_splitters = n_splitters;
	_n_classification = n_classification;
	_max_mem_size_in_gb = max_mem_size >> 30; 

}


//----------------------------------------------------------------------------------
bool CKMclass::Process()
{
	if(!initialized)
		return false;

	if(!AdjustMemoryLimits())
		return false;
		
	w1.startTimer();

	//// Create monitors
	mm = new CMemoryMonitor(max_mem_size);
	//dm = new CDiskMonitor(n_disks);
	//bo = new CBinOrdering(n_bins);

	//// Create queues
	int n_splitters_read = 1;
	pkmer  = new CPartKMERQueue(max_mem_seqKMER); // loaded file with kmers plik z kmerami	
	pq  = new CPartQueue(max_mem_fastq); // all loaded sequences	
	pread = new CReadsQueue(n_splitters_read, max_mem_reads); // loaded reads

		

	// ---------------------------------- kmers -------------
	//new -----------------
	bin_kmers = new CBinKmers_all_new(n_splitters, kmer_len, num_kmer, max_mem_binsKMER, b_sortKmer);
	if (b_sortKmer)
	{
		bin_kmers->ReadSortKmersFile(data_file_name);
	}
	else
	{		
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
	}
	w1.stopTimer();
	cout << "\n";

	w2.startTimer();
	OpenOutFile();


	// ------------------------------------ reads ----------------------
	pq  = new CPartQueue(max_mem_fastq);
	w_fastq = new CWFastqReader(mm, read_file_name, fastq_buffer_size, pq, is_fasta, true, kmer_len);
	w_splitters2class = new CWSplitter2class(mm, pq, pread, is_fasta);// rozdzielenie readow
		
	gr3.create_thread(boost::ref(*w_fastq)); // read reads
	gr4.create_thread(boost::ref(*w_splitters2class));	

	// ------------------------------------ classifier ----------------------
	w_classifier.resize(n_classification); // split reads
	
	//old ---------------------------------
	//for(int i = 0; i < n_classification; ++i)
	//{
	//		w_classifier[i] = new CWClassifier(mm, pread, bin_kmers, kmer_len, matchcutoff, outMatch, outMisMatch);
	//		gr5.create_thread(boost::ref(*w_classifier[i]));
	//}
	
	//new ----------------------
	for(int i = 0; i < n_classification; ++i)
	{
			w_classifier[i] = new CWClassifier_new(mm, pread, bin_kmers, kmer_len, step_k, matchcutoff, outMatch, outMisMatch);
			gr5.create_thread(boost::ref(*w_classifier[i]));
	}


	//--------------
	gr3.join_all(); //  read 'reads'		
	gr4.join_all();	//	split 'reads'
	gr5.join_all();	//	classification


	// ---------------- Delete  ---------------------------------------------
	delete w_fastq;
	w_splitters2class->GetTotal(n_reads);
	delete w_splitters2class;
	
	//......	
	num_alig_reads = num_REValig_reads = num_NOalig_reads = 0;
	for(int i = 0; i < n_classification; ++i)
	{
		uint64 _num_alig_reads, _num_REValig_reads, _num_NOalig_reads;
		w_classifier[i]->GetTotal(_num_alig_reads, _num_REValig_reads, _num_NOalig_reads);
		num_alig_reads += _num_alig_reads, 
		num_REValig_reads += _num_REValig_reads; 
		num_NOalig_reads += _num_NOalig_reads;
		delete w_classifier[i];
	}
	
	
	// ----------------------------------------------------
	delete pq;
	delete pread;
	delete pkmer;
	delete bin_kmers;

	//delete bpq;
	//delete bd;
	//delete bq;	
	//delete kq;

	w2.stopTimer();

	return true;
}


