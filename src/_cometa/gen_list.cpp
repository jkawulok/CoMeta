/*
	This file is a part of CoMeta software.
	The homepage of the CoMeta project is http://sun.aei.polsl.pl/cometa
 
	Authors: Jolanta Kawulok and Sebastian Deorowicz
	
	Version: 0.3
	Date   : 2014-11-26
*/

#include "stdafx.h"
#include <iostream>
#include <stdio.h>
#include <string>
#include <cstring>
#include <cstdlib> 
#include <vector>
#include <list>
#include <time.h>
#include <sys/time.h>
#include <map>

using namespace std;

int main(int argc, char* argv[])
{
	string s_nameSeq="";
	string s_nameGroupALL="";
	string s_nameMETA="";
	int i_kmer=24;
	int i_MC=5;
	bool b_newFile=1;
	string s_pathWyniki="";
	bool b_MM = 0;	
	string file_name_FILEscore ="name_files.txt";
	string file_name_FILEgroup ="name_groups.txt";

	
	int i;
	for(i = 1 ; i < argc; ++i)
	{
		if(argv[i][0] != '-')
			break;
		
		
		if(strncmp(argv[i], "-pw", 3) == 0)
			s_pathWyniki = (char*)(&argv[i][3]);
		else if(strncmp(argv[i], "-ns", 3) == 0)
			s_nameSeq = (char*)(&argv[i][3]);
		else if(strncmp(argv[i], "-ng", 3) == 0)
			s_nameGroupALL= (char*)(&argv[i][3]);
		else if(strncmp(argv[i], "-nm", 3) == 0)
			s_nameMETA = (char*)(&argv[i][3]);
		else if(strncmp(argv[i], "-NF", 3) == 0)
			file_name_FILEscore = (char*)(&argv[i][3]);
		else if(strncmp(argv[i], "-NG", 3) == 0)
			file_name_FILEgroup = (char*)(&argv[i][3]);
		else if(strncmp(argv[i], "-mc", 3) == 0)
			i_MC = atoi(&argv[i][3]);
		else if(strncmp(argv[i], "-k", 2) == 0)
			i_kmer = atoi(&argv[i][2]);
		else if(strncmp(argv[i], "-nw", 3) == 0)
			b_newFile = atoi(&argv[i][3]);
		else if(strncmp(argv[i], "-MM", 3) == 0)
			b_MM = atoi(&argv[i][3]);
	}


	FILE* pFile_outGROUP;
	FILE* pFile_outFILE;

	if (b_newFile)
	{
		pFile_outGROUP = fopen((s_pathWyniki + file_name_FILEgroup).c_str(), "wt");
		pFile_outFILE = fopen((s_pathWyniki + file_name_FILEscore).c_str(), "wt");

		if (!(pFile_outGROUP && pFile_outFILE))
			return 0;


		if (s_nameSeq == "")
		{
			fclose(pFile_outGROUP);
			fclose(pFile_outFILE);
			return 0;

		}

	}
	else
	{
		pFile_outGROUP = fopen((s_pathWyniki + file_name_FILEgroup).c_str(), "at");
		pFile_outFILE = fopen((s_pathWyniki + file_name_FILEscore).c_str(), "at");
	}



/////////
	int start = -1;
	string s_name_group="0";


	if (s_nameGroupALL == "")
	{

		int id_g = s_nameSeq.find("gi");

		if (id_g >= 0)
			start = id_g+2;
	
		
		if (start >= 0)
		{
			char c_char = s_nameSeq[start];
			while (c_char < 48 || c_char > 57)
				c_char = s_nameSeq[++start];
			
			s_name_group = "";
			while (c_char >= 48 && c_char <= 57)
			{
				s_name_group += c_char;
				c_char = s_nameSeq[++start];
			}
		}
	}
	else
		s_name_group = s_nameGroupALL;


	//fprintf(pFile_outFILE, "%s%s%d%s%s%s%d%s", s_nameSeq.c_str(), "_k", i_kmer, "_", s_nameMETA.c_str(), "_", i_MC, "_M.out\n"); //linux
	//fprintf(pFile_outFILE, "%s%s%d%s%s%s%d%s", s_nameSeq.c_str(), "", i_kmer, "_", s_nameMETA.c_str(), "_", i_MC, "_M.out\n");  //windows
	fprintf(pFile_outFILE, "%s%s", s_nameSeq.c_str(), "M.out\n");
	
	int n_end = 1;
	if (b_MM)
	{
		//fprintf(pFile_outFILE, "%s%s%d%s%s%s%d%s", s_nameSeq.c_str(), "_k", i_kmer, "_", s_nameMETA.c_str(), "_", i_MC, "_MM.out\n"); //linux
		//fprintf(pFile_outFILE, "%s%s%d%s%s%s%d%s", s_nameSeq.c_str(), "", i_kmer, "_", s_nameMETA.c_str(), "_", i_MC, "_MM.out\n"); //windows
		fprintf(pFile_outFILE, "%s%s", s_nameSeq.c_str(), "MM.out\n");
		n_end = 2;
	}

	
	for (int j = 0; j < n_end; j++)
		fprintf(pFile_outGROUP, "%s\n", s_name_group.c_str());


		
	fclose(pFile_outGROUP);
	fclose(pFile_outFILE);

	return 0;
}

