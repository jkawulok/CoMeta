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
#include "kmer_reader.h"

uint64 CKmerReader::OVERHEAD_SIZE = 1 << 16;

//----------------------------------------------------------------------------------
CKmerReader::CKmerReader(CMemoryMonitor *_mm,  uchar &_kmer, uint64* _start_prefix, uint64 _part_size)
{
	mm = _mm;
	kmer = _kmer;
	start_prefix = _start_prefix;

	in = NULL;
	part_size = _part_size;
	part = NULL;	
	pref_old = 0;

}

//----------------------------------------------------------------------------------
CKmerReader::~CKmerReader()
{
	if(in)
	{
		fclose(in);
		in = NULL;
	}
	if(part)
		delete[] part;
}

//----------------------------------------------------------------------------------
bool CKmerReader::SetNames(string _input_file_name)
{
	input_file_name = _input_file_name;

	return true;
}

//----------------------------------------------------------------------------------
bool CKmerReader::SetPartSize(uint64 _part_size)
{
	if(in)
		return false;

	//if(_part_size < (1 << 20) || _part_size > (1 << 30))
	//	return false;

	part_size = _part_size;

	return true;
}

//----------------------------------------------------------------------------------
bool CKmerReader::OpenFiles(uchar &_kmer_len)
{
	if(in)
		return false;

	if((in = fopen(input_file_name.c_str(), "rb")) == NULL)
	{
		perror("fopen");
		cout << input_file_name.c_str() << " does not exist \n";
		return false;
	}
	uchar x;
	fread(&x, 1, 1, in);
	_kmer_len = x;
	fread(&number_suff_char, 8, 1, in);

	number_suff = 0;
	for (int j = 0; j<8; j++)
	{			
		number_suff <<= 8;
		number_suff += (number_suff_char[j] & 0xFF);			
	}

	part = new uchar[part_size + 5];
	list_num_suff = new list<uint64>;
	kmer_len = _kmer_len;

	num_byte = ceil( (kmer_len - 4)/4.0) + 1; // number of bytes for the suffix and l occurrences
	return true;
}

//----------------------------------------------------------------------------------
bool CKmerReader::GetPart(uchar *&_part, uint64 &_size, list<uint64> *&_list_num_suff, int &_size_list)
{
	if(!in)
		return false;

	if(feof(in))
		return false;

	uint64 w_size = 0;
	int size_list = 0;
	int pref = -1;

	while(!feof(in)) 
	{
		uint64 addlen = 1 + number_suff * num_byte;

		if (addlen + 8 + w_size <= part_size)
		{
			(*list_num_suff).push_back(number_suff);
			size_list++;

			memcpy(part + w_size, number_suff_char, 8);
			w_size += 8;
			fread(part + w_size, 1, addlen, in); //download prefix and the rest 
		
			pref = part[w_size];
			

			if (pref<256 && pref > pref_old)
			{
				for (int pp = pref_old + 1; pp <= pref; pp++)
					start_prefix[pp] = start_prefix[pp-1];
				pref_old = pref;
			}

			if (pref<256 && pref == pref_old)
					start_prefix[pref+1] = number_suff + start_prefix[pref];

			
			pref_old = pref+1;
			w_size += addlen;
			

			fread(&number_suff_char, 8, 1, in);
			number_suff = 0;
			for (int j = 0; j<8; j++)
			{			
				number_suff <<= 8;
				number_suff += (number_suff_char[j] & 0xFF);			
			}

		}
		else
		{
			//if (w_size == 0) // added
			//{
			//	cout << "\npref: " << pref;  // --------------------------------------- 
			//	cout << "\naddlen: " << addlen;  // --------------------------------------- 	
			//}

			break;

		}

	}
	
	_part = part;
	_size = w_size;

	_list_num_suff = list_num_suff;
	_size_list = size_list;

	part = new uchar[part_size + 5];
	list_num_suff = new list<uint64>;

	return true;

}


//
// ***** EOF
