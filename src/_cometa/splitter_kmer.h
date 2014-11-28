/*
	This file is a part of CoMeta software.
	The homepage of the CoMeta project is http://sun.aei.polsl.pl/cometa
 
	Authors: Jolanta Kawulok and Sebastian Deorowicz
	
	Version: 0.3
	Date   : 2014-11-26

*/

#ifndef _SPLITTER_kmer_H
#define _SPLITTER_kmer_H

#include "defs.h"
#include "kmer_bin.h"
#include "queues.h"
#include <stdio.h>
#include <iostream>
#include <vector>
using namespace std;


//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------
class CKmerBin_Sing {

	uchar prefix; // number bin [bin]
	bool empty; 
	uint64 act_number;		
	uint64 num_elements; //elements number
	uint64* start_next_poz; // number of occurrences
	uint64* buff;
	//uchar* buff_num;

public:
	CKmerBin_Sing()
		: empty(true)
		, act_number(0)
		, num_elements(0)
		, prefix(0)
		, start_next_poz(NULL)
		, buff(NULL)
		//, buff_num(NULL) // 		

{};
		

	CKmerBin_Sing(uchar _prefix, uint64 _num_elements);	
	~CKmerBin_Sing();
	

	
	void AddKmerSecNum(uchar second, uint64 xSecNum) // 'Second' separately, and 'suffix' and 'number' together
	{
		start_next_poz[second]++;
		buff[act_number++] = xSecNum;		
		empty = false;

	}

		
	


	void Complete()
	{
		if (!empty)
			for (int i=1; i<256; i++)
				start_next_poz[i] += start_next_poz[i-1];
	}

	uint64 get_val_next_pozSin(uchar second){
		return start_next_poz[second];
	}

	bool is_empty(){
		return empty;
	}

	uint64* get_buff(){
		return buff;
	}


};


//************************************************************************************************************
//************************************************************************************************************
//************************************************************************************************************
class CBinKmers_all_new{

	bool is_completed;
	int n_writers;
	uint64 cur_buffer_size;
	uint64 max_buffer_size;

	uint64 num_kmer;

	bool* empty; 
	uint64* num_elements; //includes elements
	uint64* start_next_poz; // number occurences
	uint64* buff;
	uint64* start_prefix; // where a new prefix starts in buffor

	uint64 mask_second;
	uint64 mask_ending;
	uchar len_ending;

	uint64 mask_second2;
	uint64 mask_ending2;
	uchar len_ending2;


	boost::mutex mn_writers;
	boost::mutex m_mutex;								// The mutex to synchronise on
	boost::condition_variable m_cond_queue_full;		// The condition to wait for
	boost::condition_variable m_cond_queue_empty;

	uchar prefix_len;
	int num_bins;


public:
	CBinKmers_all_new(int _n_writers, uint32 _kmer_len, uint64 _num_kmer, uint64 _max_buffer_size, bool b_sortKmer=0);
	~CBinKmers_all_new();
	
	void save_outKMERS(	FILE *outDatabase);

	
	uint64* fun_start_prefix(){
		return start_prefix;
	}


	bool ReadSortKmersFile(string data_file_name);

	void Complete()
	{
		for (int i=1; i<num_bins*256; i++)
		{
			start_next_poz[i] += start_next_poz[i-1];
		}
	}


	bool completed() {
		return (n_writers == 0 && is_completed) ;
	}

	void mark_completed() {

		mn_writers.lock();
		n_writers--;
		mn_writers.unlock();

		if(!n_writers)		{			
			
			m_cond_queue_empty.notify_all();
			Complete();
			is_completed = true;

		}
	}

	void initial(uchar _prefix){
		empty[_prefix]=false;	
	}


	void pushSecNum(uchar bin, uchar second, uint64 sing_kmer_suffixSecNum, uint64& act_number){ // bin - prefix, sing_kmer_suffix - suffix 
		start_next_poz[bin*256+second]++;
		buff[start_prefix[bin]+act_number++] = sing_kmer_suffixSecNum;		
	}
	

	bool is_empty(uchar prefix){
		return empty[prefix];
	}

	
	uint64 get_val_next_pozSin(uchar prefix, uchar second){
		return start_next_poz[prefix*256+second];
	}

	uint64* get_buff(){
		return buff;
	}


};

//************************************************************************************************************
//************************************************************************************************************
//************************************************************************************************************
class CSplitter_KMER_new {
	CMemoryMonitor *mm;

	uchar *part;
	uint64 part_size, part_pos;
	list<uint64>* list_num_suff;

	CBinKmers_all_new *bin_KMERnew; // -----------
	CPartKMERQueue *pkmer; // a string of chars containing a loaded kmers

	uchar kmer_len;
	uint32 prefix_len;
	uint32 second_len;
	uint32 n_bins;

	uint32 suffix_len;
	uint64 act_number;		//currently occupied cell in the buffer


public:
	CSplitter_KMER_new(CMemoryMonitor *_mm, CPartKMERQueue *_pkmer, CBinKmers_all_new *_bin_KMERnew, uchar _kmer_len);
	~CSplitter_KMER_new();

	bool ProcessReads(uchar *_part, uint64 _part_size, list<uint64>* _list_num_suff);

};


//----------------------------------------------------------------------------------



#endif
