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

#pragma once

#ifndef _DEFS_H
#define _DEFS_H

#define _CRT_SECURE_NO_WARNINGS

#define MIN(x,y)				((x) < (y) ? (x) : (y))
#define NORM(x, lower, upper)	((x) < (lower) ? (lower) : (x) > (upper) ? (upper) : (x))
#define round(a)				((a) < (0) ? ceil( a - 0.5 ) : floor( a + 0.5 ) )

#define frac_eq(x,y)			(x-y<=0.000001 && y-x<=0.000001)
#define frac_eq_grat(x,y)		(x>y || (frac_eq(x,y)))



#define uchar	unsigned char
#define kmer_t	unsigned long long

#include <time.h>

//#define DEBUG_MODE

#ifdef _DEBUG
#define A_memcpy	memcpy
#define A_memset	memset
#endif

#ifdef WIN32
#define my_fopen	fopen
#define my_fseek	_fseeki64
#define my_ftell	_ftelli64
typedef int int32;
typedef unsigned int uint32;
typedef long long int64;
typedef unsigned long long uint64;
#else
#define my_fopen	fopen
#define my_fseek	fseek
#define my_ftell	ftell
#define _TCHAR	char
#define _tmain	main

typedef int int32;
typedef unsigned int uint32;
typedef long long int64;
typedef unsigned long long uint64;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#endif

const int32 FILE_BUFFER_SIZE = 1 << 20;
const int32 MAX_STR_LEN = 10240; 

void fun_STRES(char a);

#endif
