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
#include <numeric>
#include <iostream>
#include "kmer_bin.h"
#include "radix.h"

#if _DEBUG
#include <crtdbg.h> 
#define DEBUG_NEW new(_NORMAL_BLOCK, __FILE__, __LINE__) 
#define ASSERT(expr) _ASSERT(expr) 
#endif



using namespace std;

extern uint64 total_reads;

//************************************************************************************************************
// CKmerBinCollector
//************************************************************************************************************
CKmerBinCollector::CKmerBinCollector(CMemoryMonitor *_mm, uint32 _prefix, uint32 _buffer_size, uint32 _prefix_len, uint32 _kmer_len)
{
	mm				  = _mm;
	prefix            = _prefix;

	buffer_size		  = min(_buffer_size, (uint32) 1 << 24);
	prefix_len		  = _prefix_len;
	kmer_len		  = _kmer_len;

	suffix_len		  = kmer_len - prefix_len;
	kmer_mask         = (((uint64) 1) << (2*suffix_len)) - 1;

	// Determine the best value of packing len
	uint64 storage = buffer_size * 8;
	count_len = 0;
	for(uint32 i = 0; i < min(suffix_len, (uint32) 12); ++i)
	{
		uint64 tmp_storage = 1;
		tmp_storage += buffer_size * ((suffix_len - i + 3) / 4);		// truncated kmers
		tmp_storage += ((uint64) 1 << (2*i)) * 3;						// counters

		if(tmp_storage < storage)
		{
			storage   = tmp_storage;
			count_len = i;
		}
	}

	buffer = NULL;
	raw_buffer = NULL;
	file_buffer = NULL;
	pref_pack_cnts = NULL;

	Start();

}

//----------------------------------------------------------------------------------
CKmerBinCollector::~CKmerBinCollector()
{
	Release();
}

//----------------------------------------------------------------------------------
void CKmerBinCollector::Start()
{
	if(raw_buffer)
	{
		delete[] raw_buffer;
		raw_buffer = NULL;
	}
	else if(buffer)
		delete[] buffer;
	
	buffer = new uint64[buffer_size];
	buffer_counter = 0;
	total_counter  = 0;

	uint32 n_buckets = 1 << (count_len*2);
	pref_pack_cnts = new uint32[n_buckets];
	bucket_shift = (suffix_len - count_len) * 2;

	file_buffer = NULL;
}

//----------------------------------------------------------------------------------
void CKmerBinCollector::Release()
{
	if(buffer)
	{
		if(raw_buffer)
		{
			delete[] raw_buffer;
			raw_buffer = NULL;
		}
		else
			delete[] buffer;
		buffer = NULL;
		if(file_buffer)
			delete[] file_buffer;
		file_buffer = NULL;
		delete[] pref_pack_cnts;
		pref_pack_cnts = NULL;
	}
}

//----------------------------------------------------------------------------------
bool CKmerBinCollector::GetBinPart(uchar *&_bin_part, uint64 &_bin_part_size, uint64 &_n_rec)
{
	uint32 n_buckets = 1 << (count_len*2);
	uint32 i;

	if(!buffer_counter)
		return false;

	fill_n(pref_pack_cnts, n_buckets, 0);

	for(i = 0; i < buffer_counter; ++i)
		pref_pack_cnts[buffer[i] >> bucket_shift]++;

	partial_sum(pref_pack_cnts, pref_pack_cnts+n_buckets, pref_pack_cnts);

	uint32 rec_len = (kmer_len - prefix_len + 3 - count_len) / 4;
	uint32 file_buffer_size = buffer_counter * rec_len + 3 * n_buckets + sizeof(uint32);
	file_buffer = new uchar[file_buffer_size];
	uchar *file_buffer_pos = file_buffer;

	*file_buffer_pos++ = buffer_counter >> 24;
	*file_buffer_pos++ = (buffer_counter >> 16) & 0xFF;
	*file_buffer_pos++ = (buffer_counter >>  8) & 0xFF;
	*file_buffer_pos++ = (buffer_counter      ) & 0xFF;

	// Counter of packed prefixes
	for(uint32 i = 0; i < n_buckets; ++i)
	{
		*file_buffer_pos++ = pref_pack_cnts[i] >> 16;
		*file_buffer_pos++ = (pref_pack_cnts[i] >> 8) & 0xFF;
		*file_buffer_pos++ = pref_pack_cnts[i] & 0xFF;
	}
	for(uint32 i = n_buckets-1; i > 0; --i)
		pref_pack_cnts[i] = pref_pack_cnts[i-1];
	pref_pack_cnts[0] = 0;

	// Truncated kmers
	if(rec_len == 6)
	{
		for(uint32 i = 0; i < buffer_counter; ++i)
		{
			uint32 x = 6 * pref_pack_cnts[buffer[i] >> bucket_shift]++;
			file_buffer_pos[x++] = (buffer[i] >> 40) & 0xFF;
			file_buffer_pos[x++] = (buffer[i] >> 32) & 0xFF;
			file_buffer_pos[x++] = (buffer[i] >> 24) & 0xFF;
			file_buffer_pos[x++] = (buffer[i] >> 16) & 0xFF;
			file_buffer_pos[x++] = (buffer[i] >>  8) & 0xFF;
			file_buffer_pos[x++] = (buffer[i]      ) & 0xFF;
		}
	}
	else if(rec_len == 5)
	{
		for(uint32 i = 0; i < buffer_counter; ++i)
		{
			uint32 x = 5 * pref_pack_cnts[buffer[i] >> bucket_shift]++;
			file_buffer_pos[x++] = (buffer[i] >> 32) & 0xFF;
			file_buffer_pos[x++] = (buffer[i] >> 24) & 0xFF;
			file_buffer_pos[x++] = (buffer[i] >> 16) & 0xFF;
			file_buffer_pos[x++] = (buffer[i] >>  8) & 0xFF;
			file_buffer_pos[x++] = (buffer[i]      ) & 0xFF;
		}
	}
	else if(rec_len == 4)
	{
		for(uint32 i = 0; i < buffer_counter; ++i)
		{
			uint32 x = 4 * pref_pack_cnts[buffer[i] >> bucket_shift]++;
			uint64 bi = buffer[i];
			file_buffer_pos[x++] = (bi >> 24) & 0xFF;
			file_buffer_pos[x++] = (bi >> 16) & 0xFF;
			file_buffer_pos[x++] = (bi >>  8) & 0xFF;
			file_buffer_pos[x]   = (bi      ) & 0xFF;
		}
	}
	else if(rec_len == 3)
	{
		for(uint32 i = 0; i < buffer_counter; ++i)
		{
			uint32 x = 3 * pref_pack_cnts[buffer[i] >> bucket_shift]++;
			file_buffer_pos[x++] = (buffer[i] >> 16) & 0xFF;
			file_buffer_pos[x++] = (buffer[i] >>  8) & 0xFF;
			file_buffer_pos[x]   = (buffer[i]      ) & 0xFF;
		}
	}
	else if(rec_len == 2)
	{
		for(uint32 i = 0; i < buffer_counter; ++i)
		{
			uint32 x = 2 * pref_pack_cnts[buffer[i] >> bucket_shift]++;
			file_buffer_pos[x++] = (buffer[i] >>  8) & 0xFF;
			file_buffer_pos[x]   = (buffer[i]      ) & 0xFF;
		}
	}
	else
	{
		for(uint32 i = 0; i < buffer_counter; ++i)
		{
			uint32 x = pref_pack_cnts[buffer[i] >> bucket_shift]++;
			file_buffer_pos[x] = (buffer[i]      ) & 0xFF;
		}
	}

	_bin_part	   = file_buffer;
	_bin_part_size = file_buffer_size;
	_n_rec         = buffer_counter;

	total_counter += buffer_counter;
	buffer_counter = 0;

	file_buffer = NULL;

	return true;
}

//************************************************************************************************************
// CKmerBinStorer
//************************************************************************************************************
CKmerBinStorer::CKmerBinStorer(CMemoryMonitor *_mm, int _n_bins, CBinPartQueue *_q, CBinDesc *_bd, vector<string> _working_directories, uint32 _prefix_len, uint32 _kmer_len, uint64 _max_mem_buffer)
{
	mm					= _mm;
	n_bins			    = _n_bins;
	q				    = _q;
	bd                  = _bd;
	working_directories = _working_directories;

	files = NULL;
	buffer_size_bytes = 0;
	max_mem_buffer    = _max_mem_buffer;

	buffer.resize(n_bins);

	write_buffer_size = 1 << 25;
	write_buffer = new uchar[write_buffer_size];
}


//----------------------------------------------------------------------------------
CKmerBinStorer::~CKmerBinStorer()
{
	Release();
}

//----------------------------------------------------------------------------------
void CKmerBinStorer::Release() // zmiany bins
{
	if(!files)
		return;
	for(int i = 0; i < n_bins; ++i)
	{
		//if(files[i])
		//{
			//files[i] = fopen(files_name[i].data(), "ab");	// 
			WriteBin(i);
			fclose(files[i]);
		//}
	}

	delete[] files;
	delete[] files_name;
	files = NULL;

	delete[] write_buffer;
}

//----------------------------------------------------------------------------------
string CKmerBinStorer::GetName(int n, uint32 &c_disk)
{
	char tmp[6];
	sprintf(tmp, "%d", n);
	string s_tmp(tmp);
	while(s_tmp.length() < 5)
		s_tmp = string("0") + s_tmp;
	
	if(working_directories.size() == 1)
	{
		c_disk = 0;
		if(*working_directories[0].rbegin() != '/' && *working_directories[0].rbegin() != '\\')
			working_directories[0] += "/";
		return working_directories[0]  + "cometa_" + s_tmp + ".bin";
	}

	uint32 n_disks = (uint32) working_directories.size();
	//if(n >= n_bins / 4)  // 
	//	n /= 4;

	c_disk = n % (2 * n_disks);
	c_disk = MIN(c_disk, 2*n_disks-1-c_disk);

	if(*working_directories[c_disk].rbegin() != '/' && *working_directories[c_disk].rbegin() != '\\')
		working_directories[c_disk] += "/";
	return working_directories[c_disk]  + "cometa_" + s_tmp + ".bin";
}

//----------------------------------------------------------------------------------
void CKmerBinStorer::CheckBuffer()
{
	int32 i, m;

	if(buffer_size_bytes < max_mem_buffer)
		return;

	for(i = m = 0; i < (int32) buffer.size(); ++i)
		if(buffer[i].size() > buffer[m].size())
			m = i;

	WriteBin(m);
}

//----------------------------------------------------------------------------------
void CKmerBinStorer::WriteBin(uint32 n)
{
	uint32 pos = 0;

	for(uint64 i = 0; i < buffer[n].size(); ++i)
	{
		if(pos + buffer[n][i].second > write_buffer_size)
		{
			uint32 to_copy = write_buffer_size - pos;
			memcpy(write_buffer+pos, buffer[n][i].first, to_copy);
			fwrite(write_buffer, 1, write_buffer_size, files[n]);
			pos = (uint32) (buffer[n][i].second-to_copy);
			memcpy(write_buffer, buffer[n][i].first+to_copy, pos);
		}
		else
		{
			memcpy(write_buffer+pos, buffer[n][i].first, buffer[n][i].second);
			pos += (uint32) buffer[n][i].second;
		}
		buffer_size_bytes -= buffer[n][i].second;
		delete[] buffer[n][i].first;
	}

	if(pos)
		fwrite(write_buffer, 1, pos, files[n]);

	buffer[n].clear();
}

//----------------------------------------------------------------------------------
bool CKmerBinStorer::OpenFiles()   // 
{
	string f_name;
	uint32 c_disk;

	files = new FILE*[n_bins];
	files_name = new string[n_bins];

	for(int i = 0; i < n_bins; ++i)
	{
		f_name = GetName(i, c_disk);

		//if(i < n_bins / 4)
		//{
			(files[i] = fopen(f_name.c_str(), "wb"));	// 
			bd->insert(i, f_name, c_disk, 0, 0);
			files_name[i] = f_name.c_str(); // 
			//fclose(files[i]);


		//}
		//else if(i % 4 == 0)
		//{
		//	files[i] = fopen(f_name.c_str(), "wb");	
		//	bd->insert(i, f_name, c_disk, 0, 0);
		//}
		//else
		//	files[i] = NULL;
	}

	return true;
}

//----------------------------------------------------------------------------------
bool CKmerBinStorer::CloseFiles()
{
	Release();
	cout << "\n";
	
	return true;
}

//----------------------------------------------------------------------------------
void CKmerBinStorer::ProcessQueue()
{
	while(!q->completed())
	{
		int32 bin_id;
		uchar *part;
		uint64 size;

		if(q->pop(bin_id, part, size))
		{
			buffer[bin_id].push_back(make_pair(part, size));
			buffer_size_bytes += size;
			CheckBuffer();
		}
	}
}


//************************************************************************************************************
// CKmerBinReader
//************************************************************************************************************
CKmerBinReader::CKmerBinReader(CMemoryMonitor *_mm, CDiskMonitor *_dm, CBinOrdering *_bo, CBinDesc *_bd, CBinQueue *_bq)
{
	mm = _mm;
	dm = _dm;
	bo = _bo;

	bd = _bd;
	bq = _bq;
}

//----------------------------------------------------------------------------------
CKmerBinReader::~CKmerBinReader()
{

}

//----------------------------------------------------------------------------------
void CKmerBinReader::ProcessBins()
{
	FILE *in;
	uchar *data;
	uint64 readed;

	int32 bin_id;
	string name;
	uint32 c_disk;
	uint64 size;
	uint64 n_rec;
//	boost::filesystem::path p_name;

	while((bin_id = bd->get_next_bin()) >= 0)
	{
		bd->read(bin_id, name, c_disk, size, n_rec);
//		p_name.assign(name.begin(), name.end());
#ifdef DEBUG_MODE
		cout << bin_id << ":  " << name << "  " << c_disk << "  " << size << "  " << n_rec << "\n";
#else
		cout << "*";
#endif
		fflush(stdout);
		bo->block(bin_id);
		
		mm->force_increase(2*n_rec*sizeof(uint64) + 20);  // Memory reservation
		
		dm->block(c_disk);
		bo->unblock();
		
		if(size > 0)
		{
			if((in = fopen(name.c_str(), "rb")) == NULL)
				continue;

			data = new uchar[size];
			readed = fread(data, 1, size, in);
			
			if(readed != size)
				cout << "Error in bin: " << bin_id << "\n";

			bq->push(bin_id, data, size, n_rec);

			fclose(in);
		}
		else
			bq->push(bin_id, NULL, 0, 0);
//		boost::filesystem::remove(p_name);

		dm->unblock(c_disk);
	}
	bq->mark_completed();
}

//************************************************************************************************************
// CKmerBinSorter
//************************************************************************************************************
CKmerBinSorter::CKmerBinSorter(CMemoryMonitor *_mm, int32 _n_bins, CBinDesc *_bd, CBinQueue *_bq, CKmerQueue *_kq, int _n_omp_threads, uint32 _kmer_len)
{
	mm     = _mm;
	n_bins = _n_bins;
	bd	   = _bd;
	bq	   = _bq;
	kq     = _kq;
	kmer_len = _kmer_len;

	n_omp_threads = _n_omp_threads;
}

//----------------------------------------------------------------------------------
CKmerBinSorter::~CKmerBinSorter()
{

}

//----------------------------------------------------------------------------------
void CKmerBinSorter::ProcessBins()
{
	uint64 tmp_size;
	uint64 tmp_n_rec;
	uint32 c_disk;

	omp_set_num_threads(n_omp_threads);
	SetMemcpyCacheLimit(1);
	
	while(!bq->completed())
	{
		if(!bq->pop(bin_id, data, size, n_rec))
			continue;

		bd->read(bin_id, desc, c_disk, tmp_size, tmp_n_rec, buffer_size, prefix_len, kmer_len, count_len);
		//cout << "*" ;

		Expand();

		delete[] data;

		Sort(); 
		//printf("sort, "); fflush(stdout); fflush(stdout);
		
		Compact();
		//printf("compact"); fflush(stdout); fflush(stdout);
	}
	kq->mark_completed();	
}

//----------------------------------------------------------------------------------
void CKmerBinSorter::Expand()
{
	uchar *data_p = data;
	raw_buffer = new uint64[n_rec + ALIGNMENT];  //  n_rec - number in bin

	buffer = raw_buffer;
    while(((unsigned long long) buffer) % ALIGNMENT)
        buffer++;
	
	uint32 n_buckets = 1 << (count_len*2);
	uint32 *pref_pack_cnts = new uint32[n_buckets];

	uint32 rec_len = (kmer_len - prefix_len + 3 - count_len) / 4;
	uchar *part_buffer;
	
	uint64 to_read = 0;
	uint64 bin_rec = 0;
	for(bin_rec = 0; bin_rec < n_rec; bin_rec += to_read)
	{
		to_read = *data_p++;
		to_read = (to_read << 8) + *data_p++;
		to_read = (to_read << 8) + *data_p++;
		to_read = (to_read << 8) + *data_p++;

		part_buffer = data_p;
		data_p += to_read * rec_len + 3 * n_buckets;

		uchar *file_buffer_pos = part_buffer;
		for(uint32 i = 0; i < n_buckets; ++i)
		{
			pref_pack_cnts[i]  = ((uint32) (*file_buffer_pos++)) << 16;
			pref_pack_cnts[i] += ((uint32) (*file_buffer_pos++)) << 8;
			pref_pack_cnts[i] += ((uint32) (*file_buffer_pos++));
		}

		uint32 pref_id = 0;
		uint32 pref_shift = (2 * (kmer_len - prefix_len - count_len));
		uint64 x;

		if(rec_len == 6)
		{
			for(uint32 i = 0; i < to_read; ++i)
			{
				x  = ((uint64) *file_buffer_pos++) << 40;
				x += ((uint64) *file_buffer_pos++) << 32;
				x += ((uint64) *file_buffer_pos++) << 24;
				x += ((uint64) *file_buffer_pos++) << 16;
				x += ((uint64) *file_buffer_pos++) <<  8;
				x += ((uint64) *file_buffer_pos++);

				while(i >= pref_pack_cnts[pref_id])
					pref_id++;
				x += ((uint64) pref_id) << pref_shift;

				buffer[bin_rec + i] = x;
			}
		}
		else if(rec_len == 5)
		{
			for(uint32 i = 0; i < to_read; ++i)
			{
				x  = ((uint64) *file_buffer_pos++) << 32;
				x += ((uint64) *file_buffer_pos++) << 24;
				x += ((uint64) *file_buffer_pos++) << 16;
				x += ((uint64) *file_buffer_pos++) <<  8;
				x += ((uint64) *file_buffer_pos++);

				while(i >= pref_pack_cnts[pref_id])
					pref_id++;
				x += ((uint64) pref_id) << pref_shift;

				buffer[bin_rec + i] = x;
			}
		}
		else if(rec_len == 4)
		{
			for(uint32 i = 0; i < to_read; ++i)
			{
				x  = ((uint64) *file_buffer_pos++) << 24;
				x += ((uint64) *file_buffer_pos++) << 16;
				x += ((uint64) *file_buffer_pos++) <<  8;
				x += ((uint64) *file_buffer_pos++);

				while(i >= pref_pack_cnts[pref_id])
					pref_id++;
				x += ((uint64) pref_id) << pref_shift;

				buffer[bin_rec + i] = x;
			}
		}
		else if(rec_len == 3)
		{
			for(uint32 i = 0; i < to_read; ++i)
			{
				x  = ((uint64) *file_buffer_pos++) << 16;
				x += ((uint64) *file_buffer_pos++) <<  8;
				x += ((uint64) *file_buffer_pos++);

				while(i >= pref_pack_cnts[pref_id])
					pref_id++;
				x += ((uint64) pref_id) << pref_shift;

				buffer[bin_rec + i] = x;
			}
		}
		else if(rec_len == 2)
		{
			for(uint32 i = 0; i < to_read; ++i)
			{
				x  = ((uint64) *file_buffer_pos++) <<  8;
				x += ((uint64) *file_buffer_pos++);

				while(i >= pref_pack_cnts[pref_id])
					pref_id++;
				x += ((uint64) pref_id) << pref_shift;

				buffer[bin_rec + i] = x;
			}
		}
		else if(rec_len == 1)
		{
			for(uint32 i = 0; i < to_read; ++i)
			{
				x  = ((uint64) *file_buffer_pos++);

				while(i >= pref_pack_cnts[pref_id])
					pref_id++;
				x += ((uint64) pref_id) << pref_shift;

				buffer[bin_rec + i] = x;
			}
		}
	}
	if (pref_pack_cnts)
		delete[] pref_pack_cnts;
}

//----------------------------------------------------------------------------------
void CKmerBinSorter::Sort()
{
	uint32 rec_len = (kmer_len - prefix_len + 3) / 4;

	RadixSort(raw_buffer, buffer, n_rec, rec_len, n_omp_threads);
}

//----------------------------------------------------------------------------------
void CKmerBinSorter::Compact()
{
	bool write_sing = 1;

	uint64 i;
	uint64 out_buffer_size;
	uint32 kmer_shift = (2 * (kmer_len - prefix_len));
	uint32 ending_shift = (2 * (kmer_len - prefix_len - 4));

	uint64 mask_ending = (1ULL << ending_shift) - 1;
	int start_shift = ceil( (kmer_len - 8)/4.0) * 8 - 8;
	

	out_buffer_size = (n_rec)*ceil(kmer_len/4.0) + 9 + 1; //(8 for number of n_total, 1 for prefix)

	
	uchar *out_buffer = new uchar[out_buffer_size];
	out_buffer[0]=0;

	uint32 out_pos = 0;
	uint32 count;
	uint64 act_kmer, complete_kmer;
	uint64 prefix = bin_id;

	n_unique     = 0;
	n_singletons = 0;
	n_total      = 0;

	mm->decrease(n_rec*8 + 10 - out_buffer_size);
	
	//if(prefix >= (uint32) n_bins / 4)   //
	//	prefix /= 4;

	act_kmer = buffer[0];
	count = 1;
	

	n_total = n_rec;

	if (n_rec>0)
		out_pos+=9;

	for(i = 1; i < n_rec; ++i)
	{
		if(act_kmer == buffer[i])
			count++;
		else
		{
			if(count == 1 && !write_sing)
			{
				act_kmer = buffer[i];
				n_singletons++;
				n_unique++;
			}
			else
			{
				if (count==1)
					n_singletons++;

				complete_kmer = act_kmer; // + (((uint64) prefix) << kmer_shift);

				out_buffer[out_pos++] = (complete_kmer >> ending_shift) & 0xFF; // second
				act_kmer &= mask_ending;
				for (int sh = start_shift; sh >= 0; sh-=8)
					out_buffer[out_pos++] = (act_kmer >> sh) & 0xFF; //the rest (without prefix and second)


				//for (int sh = 48; sh >= 0; sh-=8)
				//	out_buffer[out_pos++] = (act_kmer >> sh) & 0xFF;



				out_buffer[out_pos++] = min<uint32>(count, 255);// & 0xFF;
				act_kmer = buffer[i];
				count = 1;
				n_unique++;
			}
		}
	}


	// zapis ost kmera z bufforu
	if( (count > 1 ||  write_sing) && n_rec>0 )
	{
		if (count==1)
			n_singletons++;

		complete_kmer = act_kmer; //+ (((uint64) prefix) << kmer_shift);

		out_buffer[out_pos++] = (complete_kmer >> ending_shift) & 0xFF; // second
		act_kmer &= mask_ending;
		for (int sh = start_shift; sh >= 0; sh-=8)
			out_buffer[out_pos++] = (act_kmer >> sh) & 0xFF;


		//for (int sh = 48; sh >= 0; sh-=8)
		//	out_buffer[out_pos++] = (act_kmer >> sh) & 0xFF;

		out_buffer[out_pos++] = min<uint32>(count, 255);// & 0xFF;
	}
	else
		if (n_rec>0)
			n_singletons++;
	

	// Write the number of unique kmers of the prefix + prefix notation (total 9 bytes \ F3w)
	if (n_rec>0)
	{
		n_unique++;
		
		uint64 n_write = n_unique;

		if (!write_sing)
			n_write -= n_singletons;

		// write the number of unique kmers of the prefix
		out_buffer[0] = (n_write >> 56) & 0xFF;
		out_buffer[1] = (n_write >> 48) & 0xFF;
		out_buffer[2] = (n_write >> 40) & 0xFF;
		out_buffer[3] = (n_write >> 32) & 0xFF;
		out_buffer[4] = (n_write >> 24) & 0xFF;
		out_buffer[5] = (n_write >> 16) & 0xFF;
		out_buffer[6] = (n_write >>  8) & 0xFF;
		out_buffer[7] = (n_write      ) & 0xFF;


		out_buffer[8] = (prefix      ) & 0xFF; // write prefix
		
		//cout << "\nbin: " << prefix << "  num write: " << n_write;
	}

	//printf("; out_pos=%d,  out_buffer_size=%d \n", out_pos, out_buffer_size);
	if (out_pos > out_buffer_size)
		printf("\n!!!!!!!!!!!! out_buffer is too small");


	kq->push(bin_id, out_buffer, out_pos, n_unique, n_singletons, n_total);

	

	if(raw_buffer)
	{
		delete[] raw_buffer;
		raw_buffer = NULL;
	}
	else
		delete[] buffer;
	buffer = NULL;
}

//************************************************************************************************************
// CKmerBinCompleter
//************************************************************************************************************
CKmerBinCompleter::CKmerBinCompleter(CMemoryMonitor *_mm, string _file_name, CBinDesc *_bd, CKmerQueue *_kq, uint32 _kmer_len)
{
	mm		  = _mm;
	file_name = _file_name;
	kq        = _kq;
	bd		  = _bd;
	kmer_len  = _kmer_len;
}

//----------------------------------------------------------------------------------
CKmerBinCompleter::~CKmerBinCompleter()
{

}

//----------------------------------------------------------------------------------
void CKmerBinCompleter::ProcessBins()
{
	int32 bin_id = -1;
	uchar *data = NULL;
	uint64 size = 0;

	FILE *out = fopen(file_name.c_str(), "wb");
	if(!out)
	{
		cout << "Error: Cannot create " << file_name << "\n";
		return;
	}
	fwrite(&kmer_len, 1, 1, out);

	uint64 _n_unique, _n_singletons, _n_total;

	n_unique = n_singletons = n_total = 0;

	while(!kq->empty())
	{
		if(!kq->pop(bin_id, data, size, _n_unique, _n_singletons, _n_total))
			continue;

		fwrite(data, 1, size, out);
		

		//cout << "\nwrite nr. bin:" << bin_id;
		//cout << ", size: " << size << "\n";

		
		if (data);
			delete[] data;

		n_unique	 += _n_unique;
		n_singletons += _n_singletons;
		n_total      += _n_total;

		string name;
		uint32 c_disk;
		uint64 n_rec;
		uint64 raw_size;


		bd->read(bin_id, name, c_disk, raw_size, n_rec);
		mm->decrease( (n_rec)*ceil(kmer_len/4.0) + 10 + n_rec*8 + 10);
	}

	//number of total unique kmer
	fwrite(&n_unique, 1, sizeof(uint64), out);
	

	fclose(out);
	cout << "\n";
}

//----------------------------------------------------------------------------------
void CKmerBinCompleter::GetTotal(uint64 &_n_unique, uint64 &_n_singletons, uint64 &_n_total)
{
	_n_unique	  = n_unique;
	_n_singletons = n_singletons;
	_n_total      = n_total;
}

// ***** EOF
