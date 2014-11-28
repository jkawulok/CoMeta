/*
  This file is a part of CoMeta software distributed under GNU GPL 2 licence.
  The homepage of the CoMeta project is http://sun.aei.polsl.pl/cometa

	Authors: Jolanta Kawulok and Sebastian Deorowicz
	
	Version: 0.3
	Date   : 2014-11-26
*/

#include "stdafx.h"
#include <algorithm>
#include "defs.h"
#include "splitter_kmer.h"



//*****************************************************************
//                      CKmerBin_Sing
//*****************************************************************
CKmerBin_Sing::CKmerBin_Sing(uchar _prefix, uint64 _num_elements)
{
	empty = true;
	prefix = _prefix;
	num_elements = _num_elements;
	act_number = 0;

	start_next_poz = new uint64[256];
	for (int i=0; i<256; i++)
		start_next_poz[i] = 0;

	buff = new uint64[num_elements];
	//buff_num = new uchar[num_elements];
}



CKmerBin_Sing::~CKmerBin_Sing()
{
	delete[] start_next_poz;
	delete[] buff;
	//delete[] buff_num;

	empty = true;
}


//*****************************************************************
//                     CBinKmers_all - new
//*****************************************************************
CBinKmers_all_new::CBinKmers_all_new(int _n_writers, uint32 _kmer_len, uint64 _num_kmer, uint64 _max_buffer_size, bool b_sortKmer) {
		
	prefix_len = 4;
	int num_bin = 1<<(2*prefix_len);
		
	n_writers = _n_writers;
	is_completed = false;
	num_bins = 256;
	num_kmer = _num_kmer;
		
	empty = new bool[num_bins]; 
	num_elements= new uint64[num_bins]; 	 //number of elements	
	start_next_poz= new uint64[num_bins*256]; 	 // the number of occurrences
	buff = new uint64[num_kmer]; 
	start_prefix = new uint64[num_bins];
	start_prefix[0]=0;

	if (!b_sortKmer)
	{
		for (int i=0; i<256; i++)
		{
			empty[i]=true;
			num_elements[i]=0;
		
			for (int s=0; s<256; s++)
				start_next_poz[i*256+s] = 0;
		}
	}

	mask_second = (1<<8)-1;
	mask_second <<= ((_kmer_len-8)*2 + 8);

	mask_second2 = (mask_second >> 8);
		
	len_ending = (_kmer_len-8)*2 + 8;
	len_ending2 = (_kmer_len-8)*2;

	mask_ending = (uint64(1)<<len_ending) - 1;
	mask_ending2 = (uint64(1)<<len_ending2) - 1;
			
	cur_buffer_size = 0;
	max_buffer_size = _max_buffer_size;
}

CBinKmers_all_new::~CBinKmers_all_new() {

	delete[] empty;
	delete[] num_elements;
	delete[] start_next_poz;
	delete[] buff;
	//cur_buffer_size = 0;
}


void CBinKmers_all_new::save_outKMERS(FILE *outDatabase)
{
	uint64 _num_kmer = num_kmer;


	fwrite(&_num_kmer, 8, 1, outDatabase);
	fwrite(empty, 1, 256, outDatabase);
	fwrite(num_elements, 8, 256, outDatabase);
	fwrite(start_next_poz, 8, 256*256, outDatabase);
	fwrite(buff, 8, _num_kmer, outDatabase);
	fwrite(start_prefix, 8, 256, outDatabase);


}



bool CBinKmers_all_new::ReadSortKmersFile(string data_file_name)
{
	prefix_len = 4;
	int num_bin = 1<<(2*prefix_len);
		
	is_completed = false;
	num_bins = 256;
		
	empty = new bool[num_bins]; 
	num_elements= new uint64[num_bins]; 	 //how much includes elements
	start_next_poz= new uint64[num_bins*256]; 	 // the number of occurrences
	buff = new uint64[num_kmer]; 
	start_prefix = new uint64[num_bins];



	FILE* outDatabase = NULL;
	if((outDatabase = fopen(data_file_name.c_str(), "rb")) == NULL)
	{
		perror("fopen");
		cout << data_file_name.c_str() << " does not exist \n";
		return false;
	}


	int kmer_len;
	fread(&kmer_len, 1, 1, outDatabase);

	fread(&num_kmer, 8, 1, outDatabase);
	fread(empty, 1, 256, outDatabase);
	fread(num_elements, 8, 256, outDatabase);
	fread(start_next_poz, 8, 256*256, outDatabase);
	fread(buff, 8, num_kmer, outDatabase);
	fread(start_prefix, 8, 256, outDatabase);

	is_completed = 1;
	fclose(outDatabase);
	return true;
}


//*****************************************************************
//                      CSplitter_KMER_NEW
//*****************************************************************
CSplitter_KMER_new::CSplitter_KMER_new(CMemoryMonitor *_mm, CPartKMERQueue *_pkmer, CBinKmers_all_new *_bin_KMERnew, uchar _kmer_len)
{
	mm = _mm;
	pkmer  = _pkmer;
	kmer_len   = _kmer_len;
	prefix_len = 4;
	second_len = 4;
	act_number = 0;		// currently occupied cell in the buffer

	n_bins        = (1 << (2 * prefix_len)) * 4;
	bin_KMERnew = _bin_KMERnew;
	part = NULL;

}

//----------------------------------------------------------------------------------
CSplitter_KMER_new::~CSplitter_KMER_new()
{

}

//----------------------------------------------------------------------------------
bool CSplitter_KMER_new::ProcessReads(uchar *_part, uint64 _part_size, list<uint64>*_list_num_suff)
{
	part = _part;
	part_size = _part_size;
	list_num_suff = _list_num_suff;
	
	uchar count_kmer;
	uchar prefix;
	uint64 suff_a_num; //suffix and num occurences
	uint64 number_suff;
	
	uint64 suffix;
	uchar number;

	uchar second;
	int num_byte_ending =  ceil( (kmer_len - 4)/4.0) -1; // number of bytes for the suffix (without prefix and second)

	int k=0;
	while (k<part_size)
	{
		number_suff = (*list_num_suff).front(); //get number of suffixes (occurence)
		(*list_num_suff).pop_front();		
		k += 8;

		prefix = part[k++]; //get preffix
		bin_KMERnew->initial(prefix);
		act_number = 0;
		
		for(int i = 0; i < number_suff; i++)
		{
			suff_a_num = 0;
			second = (part[k++] & 0xFF);
			for (int j = 0; j <= num_byte_ending; j++)
			{			
				suff_a_num <<= 8;
				suff_a_num += (part[k++] & 0xFF);			
			}
			bin_KMERnew->pushSecNum(prefix, second, suff_a_num, act_number);
		}
	}

	putchar('*');
	fflush(stdout);

	return true;
}


// ***** EOF
