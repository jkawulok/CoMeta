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

#ifndef _QUEUES_H
#define _QUEUES_H

#include "defs.h"
#include <stdio.h>
#include <iostream>
#include <queue>
//#include <tuple>
#include <list>
#include <map>
#include <string>

#include "boost/tuple/tuple.hpp"
//using boost::tuple;
//using boost::tuples::make_tuple;
//using boost::tuples::get;


using namespace std;
using namespace boost;

#define MAX_DISKS		16


//************************************************************************************************************
class CPartQueue {
	typedef queue<pair<uchar *, uint64> > queue_t;

	queue_t q;
	bool is_completed;
	uint64 cur_buffer_size;
	uint64 max_buffer_size;

	boost::mutex m_mutex;								// The mutex to synchronise on
	boost::condition_variable m_cond_queue_full;		// The condition to wait for
	boost::condition_variable m_cond_queue_empty;

public:
	CPartQueue(uint64 _max_buffer_size) {
		is_completed = false;
		cur_buffer_size = 0;
		max_buffer_size = _max_buffer_size;
	};
	~CPartQueue() {};

	bool empty() {
		return q.empty();
	}
	bool completed() {
		return q.empty() && is_completed;
	}
	void mark_completed() {
		is_completed = true;
		m_cond_queue_empty.notify_all();
	}
	void push(uchar *part, uint64 size) {
		boost::unique_lock<boost::mutex> lock(m_mutex);
		while ((cur_buffer_size + size) > max_buffer_size) m_cond_queue_full.wait(lock);
		q.push(make_pair(part, size));
		cur_buffer_size += size;
		m_cond_queue_empty.notify_one();
	}
bool pop(uchar *&part, uint64 &size) {
		boost::unique_lock<boost::mutex> lock(m_mutex);
		while ((q.size() == 0) && !is_completed) m_cond_queue_empty.wait(lock);

		if(q.size() != 0)
		{
			part = q.front().first;
			size = q.front().second;
			cur_buffer_size -= size;
			q.pop();
			m_cond_queue_full.notify_one();
			return true;
		}
		return false;
	}
};


//************************************************************************************************************
class CPartKMERQueue {
	typedef queue<tuple<uchar *, uint64, list<uint64>*, int> > queue_t;

	queue_t q;
	bool is_completed;
	uint64 cur_buffer_size;
	uint64 max_buffer_size;

	boost::mutex m_mutex;								// The mutex to synchronise on
	boost::condition_variable m_cond_queue_full;		// The condition to wait for
	boost::condition_variable m_cond_queue_empty;

public:
	CPartKMERQueue(uint64 _max_buffer_size) {
		is_completed = false;
		cur_buffer_size = 0;
		max_buffer_size = _max_buffer_size;
	};
	~CPartKMERQueue() {};

	bool empty() {
		return q.empty();
	}
	bool completed() {
		return q.empty() && is_completed;
	}
	void mark_completed() {
		is_completed = true;
		m_cond_queue_empty.notify_all();
	}
	void push(uchar *part, uint64 size, list<uint64> *list_num_suff, int size_list) {
		boost::unique_lock<boost::mutex> lock(m_mutex);
		while ((cur_buffer_size + size) > max_buffer_size) m_cond_queue_full.wait(lock);
		
		q.push(make_tuple(part, size, list_num_suff, size_list));
		cur_buffer_size += (size + 8 * (size_list+1)) ;
		m_cond_queue_empty.notify_one();
	}
	bool pop(uchar *&part, uint64 &size, list<uint64> *&list_num_suff, int &size_list ) {
		boost::unique_lock<boost::mutex> lock(m_mutex);
		while ((q.size() == 0) && !is_completed) 
			m_cond_queue_empty.wait(lock);

		if(q.size() != 0)
		{
			part			= get<0>(q.front());
			size			= get<1>(q.front());
			list_num_suff	= get<2>(q.front());
			size_list		= get<3>(q.front());			

			cur_buffer_size -= (size + 8 * (size_list+1)) ;
			q.pop();
			m_cond_queue_full.notify_one();
			return true;
		}
		return false;
	}
};

//************************** Ready **************************************** 
class CReadsQueue {
	typedef queue<tuple<char *, char *, uint64, uint32> > queue_t;
	
	queue_t q;
	bool is_completed;
	int n_writers;

	uint64 cur_buffer_size;
	uint64 max_buffer_size;

	boost::mutex mn_writers;
	boost::mutex m_mutex;								// The mutex to synchronise on
	boost::condition_variable m_cond_queue_full;		// The condition to wait for
	boost::condition_variable m_cond_queue_empty;

public:
	CReadsQueue(int _n_writers, uint64 _max_buffer_size) {
		n_writers = _n_writers;
		is_completed = false;
		cur_buffer_size = 0;
		max_buffer_size = _max_buffer_size;
	};
	~CReadsQueue() {};

	bool empty() {
		return q.empty();
	}
	bool completed() {
		return q.empty() && is_completed;
	}
	void mark_completed() {
		mn_writers.lock();
		n_writers--;
		mn_writers.unlock();
		if(!n_writers)
		{
			is_completed = true;
			m_cond_queue_empty.notify_all();			
		}
	}
	void push(char *name, char *seq, uint64 size, uint32 size_name) {
		boost::unique_lock<boost::mutex> lock(m_mutex);

		while ((cur_buffer_size + size + size_name) > max_buffer_size) 
			m_cond_queue_full.wait(lock);
		
		q.push(make_tuple(name, seq, size, size_name));
		cur_buffer_size += (size + size_name);
		m_cond_queue_empty.notify_one();
	}
	bool pop(char *&name, char *&seq, uint64 &size, uint32 &size_name) {
		boost::unique_lock<boost::mutex> lock(m_mutex);
		
		while (  ((q.size() == 0) && !is_completed) ) 
			m_cond_queue_empty.wait(lock);

		if(q.size() != 0)
		{
			name		= get<0>(q.front()); // single read's name
			seq			= get<1>(q.front()); // single read's content
			size		= get<2>(q.front()); // single read's length
			size_name	= get<3>(q.front()); // length name of read

			q.pop();
			cur_buffer_size -= (size + size_name);			
			m_cond_queue_full.notify_one();
			return true;
		}
		return false;
	}
};


//************************************************************************************************************
class CBinPartQueue {
	typedef queue<tuple<int32, uchar *, uint64> > queue_t;
	
	queue_t q;
	int n_writers;
	bool is_completed;
	uint64 cur_buffer_size;
	uint64 max_buffer_size;

	boost::mutex mn_writers;
	boost::mutex m_mutex;								// The mutex to synchronise on
	boost::condition_variable m_cond_queue_full;		// The condition to wait for
	boost::condition_variable m_cond_queue_empty;

public:
	CBinPartQueue(int _n_writers, uint64 _max_buffer_size) {
		n_writers = _n_writers;
		is_completed = false;
		cur_buffer_size = 0;
		max_buffer_size = _max_buffer_size;
	}
	~CBinPartQueue() {}

	bool empty() {
		return q.empty();
	}
	bool completed() {
		return q.empty() && n_writers == 0;
	}
	void mark_completed() {
		mn_writers.lock();
		n_writers--;
		mn_writers.unlock();
		if(!n_writers)
			m_cond_queue_empty.notify_all();
	}
	void push(int32 bin_id, uchar *part, uint64 size) {
		boost::unique_lock<boost::mutex> lock(m_mutex);
		
		while ((cur_buffer_size + size) > max_buffer_size) 
			m_cond_queue_full.wait(lock);

		q.push(make_tuple(bin_id, part, size)); // part: bin 
		cur_buffer_size += size;
		m_cond_queue_empty.notify_one();
	}
	bool pop(int32 &bin_id, uchar *&part, uint64 &size) {
		boost::unique_lock<boost::mutex> lock(m_mutex);
		while ((q.size() == 0) && (n_writers > 0)) 
			m_cond_queue_empty.wait(lock);

		if(q.size() != 0)
		{
			bin_id = get<0>(q.front());
			part   = get<1>(q.front()); // content of bin
			size   = get<2>(q.front()); // size of bin
			q.pop();
			cur_buffer_size -= size;
			m_cond_queue_full.notify_one();
			return true;
		}
		return false;
	}
};

//************************************************************************************************************
class CBinDesc {
	typedef tuple<string, uint32, uint64, uint64, uint32, uint32, uint32, uint32> desc_t;
	typedef map<int32, desc_t> map_t;

	map_t m;
	int32 bin_id;

	boost::mutex m_mutex;

public:
	CBinDesc() {
		bin_id = -1;
	}
	~CBinDesc() {}

	void reset_reading() {
		boost::unique_lock<boost::mutex> lock(m_mutex);
		bin_id = -1;
	}

	bool empty() {
		return m.empty();
	}

	int32 get_next_bin()
	{
		boost::unique_lock<boost::mutex> lock(m_mutex);
		map_t::iterator p;
		if(bin_id == -1)
			p = m.begin();
		else
		{
			p = m.find(bin_id);
			if(p != m.end())
				++p;
		}

		if(p == m.end())
			bin_id = -1000;
		else
			bin_id = p->first;		

		return bin_id;
	}
	void insert(int32 bin_id, string desc, uint32 c_disk, uint64 size, uint64 n_rec, uint32 buffer_size = 0, uint32 prefix_len = 0, uint32 kmer_len = 0, uint32 count_len = 0) {
		boost::unique_lock<boost::mutex> lock(m_mutex);
		map_t::iterator p = m.find(bin_id);
		if(p != m.end())
		{
			if(desc != "")
			{
				get<0>(m[bin_id]) = desc;
				get<1>(m[bin_id]) = c_disk;
			}
			get<2>(m[bin_id]) += size; // bin size
			get<3>(m[bin_id]) += n_rec; //number in bin

			if(buffer_size)
			{
				get<4>(m[bin_id]) = buffer_size;
				get<5>(m[bin_id]) = prefix_len;
				get<6>(m[bin_id]) = kmer_len;
				get<7>(m[bin_id]) = count_len;
			}
		}
		else
			m[bin_id] = make_tuple(desc, c_disk, size, n_rec, buffer_size, prefix_len, kmer_len, count_len);
	}
	void read(int32 bin_id, string &desc, uint32 &c_disk, uint64 &size, uint64 &n_rec, uint32 &buffer_size, uint32 &prefix_len, uint32 &kmer_len, uint32 &count_len) {
		boost::unique_lock<boost::mutex> lock(m_mutex);
		desc		= get<0>(m[bin_id]);
		c_disk      = get<1>(m[bin_id]);
		size		= get<2>(m[bin_id]);
		n_rec		= get<3>(m[bin_id]);
		buffer_size = get<4>(m[bin_id]);
		prefix_len  = get<5>(m[bin_id]);
		kmer_len    = get<6>(m[bin_id]);
		count_len   = get<7>(m[bin_id]);
	}
	void read(int32 &bin_id, string &desc, uint32 &c_disk, uint64 &size, uint64 &n_rec) {
		boost::unique_lock<boost::mutex> lock(m_mutex);
		desc		= get<0>(m[bin_id]);
		c_disk		= get<1>(m[bin_id]);
		size		= get<2>(m[bin_id]);
		n_rec		= get<3>(m[bin_id]);
	}
};

//************************************************************************************************************
class CBinQueue {
	typedef queue<tuple<int32, uchar *, uint64, uint64> > queue_t;
	queue_t q;

	int n_writers;

	boost::mutex mn_writers;
	boost::mutex m_mutex;								// The mutex to synchronise on
	boost::condition_variable m_cond_queue_full;		// The condition to wait for
	boost::condition_variable m_cond_queue_empty;

public:
	CBinQueue(int _n_writers) {
		n_writers = _n_writers;
	}
	~CBinQueue() {}

	bool empty() {
		return q.empty();
	}
	bool completed() {
		return q.empty() && n_writers == 0;
	}
	void mark_completed() {
		mn_writers.lock();
		n_writers--;
		mn_writers.unlock();
		if(n_writers == 0)
			m_cond_queue_empty.notify_all();
	}
	void push(int32 bin_id, uchar *part, uint64 size, uint64 n_rec) {
		boost::unique_lock<boost::mutex> lock(m_mutex);
		q.push(make_tuple(bin_id, part, size, n_rec));
		m_cond_queue_empty.notify_one();
	}
	bool pop(int32 &bin_id, uchar *&part, uint64 &size, uint64 &n_rec) {
		boost::unique_lock<boost::mutex> lock(m_mutex);
		while ((q.size() == 0) && (n_writers > 0)) 
			m_cond_queue_empty.wait(lock);

		if(q.size() != 0)
		{
			bin_id = get<0>(q.front());
			part   = get<1>(q.front());
			size   = get<2>(q.front());
			n_rec  = get<3>(q.front());
			q.pop();

			return true;
		}
		return false;
	}
};

//************************************************************************************************************
class CKmerQueue {
	typedef tuple<int32, uchar*, uint64, uint64, uint64, uint64> data_t;
	typedef list<data_t> list_t;

	int n_writers;

	boost::mutex mn_writers;
	boost::mutex m_mutex;								// The mutex to synchronise on
	boost::condition_variable m_cond_queue_full;		// The condition to wait for
	boost::condition_variable m_cond_queue_empty;

	list_t l;
	int32 n_bins;
	int32 cur_id;

public:
	CKmerQueue(int32 _n_bins, int _n_writers) {
		n_bins = _n_bins;
		cur_id = 0;

		n_writers = _n_writers;
	}
	~CKmerQueue() {
	}


	bool empty() {
		boost::unique_lock<boost::mutex> lock(m_mutex);
		return cur_id >= n_bins;
	}
	bool ready() {
		boost::unique_lock<boost::mutex> lock(m_mutex);
		if(l.empty())
			return false;
		return get<0>(l.front()) == cur_id;
	}
	bool completed() {
		return l.empty() && n_writers == 0;
	}
	void mark_completed() {
		mn_writers.lock();
		n_writers--;
		mn_writers.unlock();
		if(!n_writers)
			m_cond_queue_empty.notify_all();
	}
	void push(int32 bin_id, uchar *data, uint64 size, uint64 n_unique, uint64 n_singletons, uint64 n_total) {
		boost::unique_lock<boost::mutex> lock(m_mutex);

		list_t::iterator p;
		for(p = l.begin(); p != l.end(); ++p)
			if(bin_id < get<0>(*p))
			{
				l.insert(p, make_tuple(bin_id, data, size, n_unique, n_singletons, n_total));
				break;
			}

		if(p == l.end())
			l.insert(p, make_tuple(bin_id, data, size, n_unique, n_singletons, n_total));
		m_cond_queue_empty.notify_one();
	}
	bool pop(int32 &bin_id, uchar *&data, uint64 &size, uint64 &n_unique, uint64 &n_singletons, uint64 &n_total) {
		boost::unique_lock<boost::mutex> lock(m_mutex);
		while ((l.size() == 0) && (n_writers > 0)) 
			m_cond_queue_empty.wait(lock);
		
		if(l.empty() || get<0>(l.front()) > cur_id)
			return false;

		bin_id       = get<0>(l.front());
		data         = get<1>(l.front());
		size         = get<2>(l.front());
		n_unique     = get<3>(l.front());
		n_singletons = get<4>(l.front());
		n_total      = get<5>(l.front());
				
		l.pop_front();

		//if(cur_id < n_bins / 4) //
			cur_id++;
		//else
		//	cur_id += 4;

		if(cur_id >= n_bins)
			m_cond_queue_empty.notify_all();

		return true;
	}
};

//************************************************************************************************************
class CMemoryMonitor {

	uint64 max_memory;
	uint64 memory_in_use;

	boost::mutex m_mutex;								// The mutex to synchronise on
	boost::condition_variable m_cond_full;				// The condition to wait for

public:
	CMemoryMonitor(uint64 _max_memory) {
		max_memory    = _max_memory;
		memory_in_use = 0;
	}
	~CMemoryMonitor() {
	}

	void increase(uint64 n) {
		boost::unique_lock<boost::mutex> lock(m_mutex);
		while (memory_in_use + n > max_memory) 
			m_cond_full.wait(lock);
		
		memory_in_use += n;
	}
	void force_increase(uint64 n) {
		boost::unique_lock<boost::mutex> lock(m_mutex);
		
		if (n > max_memory)
		{
			printf("\n!!!!!!!!!!Not enough memory for bin (max_memory = %lu, n = %lu \n", max_memory, n);
		}
		
		
		while (memory_in_use + n > max_memory && memory_in_use > 0)
			m_cond_full.wait(lock);
		
		memory_in_use += n;
	}
	void decrease(uint64 n) {
		boost::unique_lock<boost::mutex> lock(m_mutex);
		memory_in_use -= n;
		m_cond_full.notify_one();
	}
	void info(uint64 &_max_memory, uint64 &_memory_in_use)
	{
		boost::unique_lock<boost::mutex> lock(m_mutex);
		_max_memory    = max_memory;
		_memory_in_use = memory_in_use;
	}
};

//************************************************************************************************************
class CDiskMonitor {
	uint32 n_disks;
	boost::mutex v_mutex[MAX_DISKS];
	boost::condition_variable v_cond[MAX_DISKS];
	bool v_disk[MAX_DISKS];

public:
	CDiskMonitor(uint32 _n_disks) {
		n_disks = NORM(_n_disks, 1, MAX_DISKS);
		for(uint32 i = 0; i < n_disks; ++i)
			v_disk[i] = true;
	}
	~CDiskMonitor() {
	}

	void block(uint32 n) {
		boost::unique_lock<boost::mutex> lock(v_mutex[n]);
		while(!v_disk[n]) v_cond[n].wait(lock);
		v_disk[n] = false;
	}
	void unblock(uint32 n) {
		boost::unique_lock<boost::mutex> lock(v_mutex[n]);
		v_disk[n] = true;
		v_cond[n].notify_one();
	}
};

//************************************************************************************************************
class CBinOrdering {
	int32 n_bins;
	int32 cur_bin;

	boost::mutex m_mutex;							// The mutex to synchronise on
	boost::condition_variable m_cond;				// The condition to wait for

public:
	CBinOrdering(int32 _n_bins) {
		n_bins  = _n_bins;
		cur_bin = -1;
	}
	~CBinOrdering() {
	}
	void block(int32 n) {
		boost::unique_lock<boost::mutex> lock(m_mutex);
		while(cur_bin + (cur_bin <  n_bins/4 ? 1 : 4) < n) m_cond.wait(lock);
	}
	void unblock() {
		boost::unique_lock<boost::mutex> lock(m_mutex);
		
		cur_bin += cur_bin < n_bins/4 ? 1 : 4;
		m_cond.notify_all();
	}
};

#endif
