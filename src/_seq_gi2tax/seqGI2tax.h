/*
	This file is a part of CoMeta software.
	The homepage of the CoMeta project is http://sun.aei.polsl.pl/cometa
 
	Authors: Jolanta Kawulok and Sebastian Deorowicz
	
	Version: 0.3
	Date   : 2014-11-26

*/
#ifndef _SEQ_GI_2_TAX_H
#define _SEQ_GI_2_TAX_H


#include <map>

class comp
{
public:
	bool operator() (const unsigned int& n1, const unsigned int& n2)
	{
		return n1 < n2;
	}
};

typedef std::map<unsigned int, unsigned int, comp> TIDMap;

#endif