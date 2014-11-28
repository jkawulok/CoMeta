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
#include <stdio.h>
#include "radix.h"
#include "defs.h"

//----------------------------------------------------------------------------------
void RadixOMP(uint64 *Source, uint64 *Dest, const int64 SourceSize, const unsigned n_phases, const unsigned n_threads)
{
    
#ifdef WIN32
	__declspec( align( WIN_ALIGNMENT ) ) unsigned ByteCounter[MAX_NUM_THREADS][256];  //WIN_ALIGNMENT 64
#else
	unsigned ByteCounter[MAX_NUM_THREADS][256] __attribute__((aligned(ALIGNMENT)));
#endif

#ifdef WIN32
	__declspec( align( WIN_ALIGNMENT ) ) unsigned globalHisto[256];
#else
	unsigned globalHisto[256] __attribute__((aligned(ALIGNMENT)));
#endif
	
	
	#pragma omp parallel
	{ 
		int myID = omp_get_thread_num();
		uint8_t ByteIndex = 0;
		long long i;					
		unsigned prevSum;
		unsigned temp;
		uint32 n;
		
		int index_x;
		int private_i;
		int byteValue;
		uint64 *tempSource = Source;
		uint64 *tempDest = Dest;
		uint64 *tempPtr;
					
		uint64 *raw_Buffer = new uint64[256 * BUFFER_WIDTH + ALIGNMENT]; //BUFFER_WIDTH 32; ALIGNMENT 0x100
		uint64 *Buffer = raw_Buffer;     
		while(((unsigned long long) Buffer) % ALIGNMENT)
			Buffer++; 

#ifdef WIN32
	__declspec( align( WIN_ALIGNMENT ) ) unsigned privateByteCounter[256] = {0};
#else
	__attribute__((aligned(ALIGNMENT)))  unsigned privateByteCounter[256] = {0};
#endif
		
   		for(uint32 privatePhaseCounter = 0; privatePhaseCounter < n_phases; privatePhaseCounter++)
		{
			#pragma omp for private(i) schedule(static) 
			for(i = 0; i < SourceSize; ++i)
			{
				byteValue = *(reinterpret_cast<const uint8_t*>(&tempSource[i]) + ByteIndex);
				++privateByteCounter[byteValue];
			}	
			memcpy(&ByteCounter[myID][0], privateByteCounter, sizeof(privateByteCounter)); //copying between buffers

			#pragma omp barrier

			#pragma omp for schedule(static)
			for(i = 0; i < 256; ++i)
			{
				prevSum = 0; 
				for(n = 0; n < n_threads; n++)
				{
					temp = ByteCounter[n][i];
					ByteCounter[n][i] = prevSum;
					prevSum += temp; 
				}
				globalHisto[i] = prevSum;
			}	

			#pragma omp single
			{
				prevSum = 0; 
				for(i = 0; i < 256; ++i)
				{
					temp = globalHisto[i];
					globalHisto[i] = prevSum;
					prevSum += temp; 
				}	
			}


			for (private_i = 0; private_i < 256; private_i++)
				ByteCounter[myID][private_i] += globalHisto[private_i];

			memcpy(privateByteCounter, &ByteCounter[myID][0], sizeof(privateByteCounter));
		

			#pragma omp for schedule(static)
			for(i = 0; i < SourceSize; ++i)
			{
				byteValue = *(reinterpret_cast<const uint8_t*>(&tempSource[i]) + ByteIndex);

				index_x = privateByteCounter[byteValue] % BUFFER_WIDTH;

				Buffer[byteValue * BUFFER_WIDTH + index_x] = tempSource[i];
			
				privateByteCounter[byteValue]++;

				if(index_x == (BUFFER_WIDTH -1))
					memcpy ( &tempDest[privateByteCounter[byteValue] - (BUFFER_WIDTH)], &Buffer[byteValue * BUFFER_WIDTH], BUFFER_WIDTH *sizeof(uint64) );
			} //end_for

			//#pragma omp barrier 
			int elemInBuffer;
			int index_stop;
			int index_start;
			int elemWrittenIntoBuffer;
		
			for(private_i = 0; private_i < 256; private_i++)
			{
				index_stop = privateByteCounter[private_i] % BUFFER_WIDTH;
				index_start = ByteCounter[myID][private_i] % BUFFER_WIDTH;
				elemWrittenIntoBuffer = privateByteCounter[private_i] - ByteCounter[myID][private_i];

				if((index_stop - elemWrittenIntoBuffer) <= 0)
					elemInBuffer = index_stop;
				else
					elemInBuffer = index_stop - index_start;

				if(elemInBuffer != 0)
					memcpy ( &tempDest[privateByteCounter[private_i] - elemInBuffer], &Buffer[private_i * BUFFER_WIDTH + (privateByteCounter[private_i] - elemInBuffer)%BUFFER_WIDTH], (elemInBuffer)*sizeof(uint64) );
			
			}
			#pragma omp barrier

			tempPtr = tempDest;
			tempDest = tempSource;
			tempSource = tempPtr;
			ByteIndex++;
			std::memset(privateByteCounter, 0, sizeof(privateByteCounter));
		}
		delete [] raw_Buffer;
	}
}

//----------------------------------------------------------------------------------
void RadixSort(uint64 *&raw_data, uint64 *&data, int64 size, const unsigned n_phases, const unsigned n_threads)
{
    uint64 *raw_tempData = new uint64[size+ALIGNMENT]; //ALIGNMENT 0x100
    uint64 *tempData = raw_tempData;     

    while(((uint64) tempData) % ALIGNMENT)
        tempData++; 	

    RadixOMP(data, tempData, size, n_phases, n_threads);
	//data - sorted

	if(n_phases % 2)
	{
		if (raw_data)
			delete[] raw_data;
		raw_data = raw_tempData;
		data = tempData;
	}
	else	
		delete [] raw_tempData;
}

// ***** EOF
