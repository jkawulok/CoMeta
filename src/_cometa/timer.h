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

#ifndef _TIMER_H
#define _TIMER_H

#ifdef WIN32
#include <windows.h>

typedef struct {
    LARGE_INTEGER start;
    LARGE_INTEGER stop;
} stopWatch;

class CStopWatch {

private:
	stopWatch timer;
	LARGE_INTEGER frequency;
	double LIToSecs( LARGE_INTEGER & L);
public:
	CStopWatch();
	void startTimer( );
	void stopTimer( );
	double getElapsedTime();
};

#else
#include <sys/time.h>

typedef struct {
	timeval start;
	timeval stop;
} stopWatch;

class CStopWatch {

private:
	stopWatch timer;
public:
	CStopWatch() {};
	void startTimer( );
	void stopTimer( );
	double getElapsedTime();
};

#endif

#endif
