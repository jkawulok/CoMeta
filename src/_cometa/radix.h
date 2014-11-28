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

#include <cassert>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <omp.h>
#include <algorithm>
#include "asmlib.h"
#include "defs.h"
#include <boost/static_assert.hpp>

#ifdef WIN32
typedef unsigned __int8 uint8_t;
#else
#include <stdint.h>
#endif

#define MAX_NUM_THREADS 32
#define BUFFER_WIDTH 32
#define ALIGNMENT 0x100
#define WIN_ALIGNMENT 64

void RadixOMP(uint64 *Source, uint64 *Dest, const int64 SourceSize, const unsigned n_phases, const unsigned n_threads);
void RadixSort(uint64 *&raw_data, uint64 *&data, const int64 size, const unsigned n_phases, const unsigned n_threads);

