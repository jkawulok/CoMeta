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
#include "fastq_reader.h"


#ifdef _DEBUG
   #ifndef DBG_NEW
      #define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
      #define new DBG_NEW
   #endif
#endif  // _DEBUG


uint64 CFastqReader::OVERHEAD_SIZE = 1 << 16;

//----------------------------------------------------------------------------------
CFastqReader::CFastqReader(CMemoryMonitor *_mm, bool _is_fasta, bool _is_read, int _kmer)
{
	mm = _mm;
	is_fasta = _is_fasta;
	is_read = _is_read;
	kmer = _kmer;

	in = NULL;

	part_size = 1 << 20;

	part = NULL;
	part_seq = NULL;
}

//----------------------------------------------------------------------------------
CFastqReader::~CFastqReader()
{
	if(in)
		fclose(in);
	if(part)
		delete[] part;
	if(part_seq)
		delete[] part_seq;
}

//----------------------------------------------------------------------------------
bool CFastqReader::SetNames(list<string> _input_file_name)
{
	input_file_name = _input_file_name;

	return true;
}

//----------------------------------------------------------------------------------
bool CFastqReader::SetPartSize(uint64 _part_size)
{
	if(in || part)
		return false;

	//if(_part_size < (1 << 20) || _part_size > (1 << 30))
	//	return false;

	part_size = _part_size;

	return true;
}

//----------------------------------------------------------------------------------
bool CFastqReader::OpenFiles()
{
	if(in)
	{
		fclose(in);
		in = NULL;
		if(input_file_name.empty() )
			return false;
	}

	string input_file_nameSING = input_file_name.front();
	if((in = fopen(input_file_nameSING.c_str(), "rb")) == NULL)
	{
		perror("fopen");
		cout << input_file_nameSING.c_str() << " does not exist \n";
		return false;
	}
	
	if (!part)
		part = new uchar[part_size + OVERHEAD_SIZE];
	
	part_filled = 0;
	input_file_name.pop_front();

	return true;
}

//----------------------------------------------------------------------------------
bool CFastqReader::GetPart(uchar *&_part, uint64 &_size, char *&_part_seq, int &_kmer)
{
	_kmer = kmer;


	if(!in)
		return false;

	if(feof(in) && part_seq==NULL)
		return false;

	uint64 readed = fread(part+part_filled, 1, part_size, in);


	if (part_seq!=NULL)
	{			
		int kk=0;
		int kc=0;

		if (strlen(part_seq) == kmer-1) //  
		{
			while (kk<kmer-1 && kk<readed)
			{
				char cc = part[kc];
				kc++;

				if (cc>=32)
				{
					part_seq[kmer-1+kk] = cc;
					kk++;
				}
			}
			part_seq[kmer-1+kk] = 0;
		}

		if ((strlen(part_seq)) > kmer) // 
		{
			_part_seq = part_seq; 
			part_seq = NULL;
		}
	}


	int64 total_filled = part_filled + readed;
	int64 i;

	if(part_filled >= OVERHEAD_SIZE)
		cout << "Overhead 1!!!\n";

	if(feof(in))
	{
		_part = part;
		_size = total_filled;

		part = NULL;
		return true;
	}
	

	// Looking for a FASTA record at the end of the area
	int64 line_start[3];
	int32 j;

	if (is_read)
		i = total_filled - OVERHEAD_SIZE / 2;
	else
		i = total_filled - OVERHEAD_SIZE;

	if (i<0)
	{
		cout << "ii \n";
		i=0;
	}

	for(j = 0; j < 1; ++j)
	{
		if(!SkipNextEOL(part, i, total_filled, part_seq, kmer))
			break;
		line_start[j] = i;
	}

	
	_part = part;
	_size = line_start[0];

		
	if (_size > total_filled)
	{
		cout << "\n size: " << _size << "; total_filled: " << total_filled << "  lin st[0] " << line_start[0] <<  "i  " << i << "\n";
		_size = 0;
	}


	part = new uchar[part_size + OVERHEAD_SIZE];
	copy(_part+_size, _part+total_filled, part);
	part_filled = total_filled - _size;

	return true;
}

//----------------------------------------------------------------------------------
bool CFastqReader::SkipNextEOL(uchar *part, int64 &pos, int64 max_pos, char *&part_seq, int kmer)
{
	int64 i;
	for(i = pos; i < max_pos-2; ++i)
	{
		if (part[i+1] == '>')
			break;

	}



	if(i >= max_pos-2) //Debug!!!
	{
		part_seq = new char[2*kmer+10];

		int kk = kmer-2;
		int64 kc = max_pos-2;

		while (kk>=0)
		{
			char cc = part[kc];
			kc--;

			if (cc>=32)
			{
				if(kk >= 2*kmer+10 || kk < 0) 
					printf("Error: part_seq exceeded");
				part_seq[kk] = cc;
				kk--;
			}
		}

		part_seq[kmer-1] = 0;
	}

	pos = i+1;

	return true;
}

// ***** EOF
