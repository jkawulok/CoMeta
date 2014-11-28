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

unsigned GetTickCount()
{
        struct timeval tv;
        if(gettimeofday(&tv, NULL) != 0)
                return 0;

        return (tv.tv_sec * 1000) + (tv.tv_usec / 1000);
}




struct str_seq
{
	char *sequence;
	char *name;	
};

struct CGroupMatch
{
	string GroupName;
	string SequenceName;
	int Match;
	int Min_err;
};

typedef list<CGroupMatch> L_CGrM;



bool Fread_Sin_score(string s_file_name_input, map<string, L_CGrM>& M_Group)
{
	FILE* pFile_in = NULL;
	pFile_in = fopen(s_file_name_input.c_str(), "rt");
	int nMaxLineSize = 1000000;
	char* name_line = new char[nMaxLineSize];
	if (pFile_in == NULL)
	{
		perror("fopen");
		printf("%s \n", s_file_name_input.c_str());
		cout << s_file_name_input;
		return 0;
	}
	
	int id_fs, id_fe;
	string seq_name;
	int match;
	int err;
	string name_lineS;

	while(fscanf(pFile_in, "%s", name_line) != EOF)
	{ 
		CGroupMatch grmat;
		fscanf(pFile_in, "%d", &(grmat.Min_err));  //err
		fscanf(pFile_in, "%s", name_line);
		fscanf(pFile_in, "%d", &(grmat.Match));  //match
		grmat.GroupName = "";

		fscanf(pFile_in, "%[^\n]", name_line);
		name_lineS = name_line;
		id_fe = name_lineS.find(9, 3);		
		int id_st = 0;
		while (name_line[id_st] > 8 && name_line[id_st] < 33)
		{
			id_st++;						
		}		
		seq_name = name_lineS.substr(id_st, id_fe);
		
		grmat.SequenceName = seq_name;

		char *cstr = new char[seq_name.length() + 5];
		strcpy(cstr, seq_name.c_str());
		
		int len_name = seq_name.find("<name|");						
		if (len_name < 0)
			len_name = strlen(cstr);
		else
			cstr[len_name]=0;
				
		while (cstr[len_name-1] > 8 && cstr[len_name-1] < 33)
		{
			len_name--;
			cstr[len_name]=0;						
		}	
		M_Group[string(cstr)].push_back(grmat);
		delete[] cstr;
	}

	delete[] name_line;
	fclose(pFile_in);

	return 1;
}


bool SetFileNamesList(string file_name_FILEsequence, list<string>& file_name_sequence, string path_seq)
{
	FILE *in;
	if((in = fopen(file_name_FILEsequence.c_str(), "rt")) == NULL)
	{
		perror("fopen");
		cout << file_name_FILEsequence.c_str() << "\n";		
		return 0;
	}

	char* pszCursor = new char[10000];
	char* pszCursorADD = new char[10000];
	string pszCursorS;

	while(fscanf(in,  "%[^\n]", pszCursor) != EOF) //whole line
	{
		pszCursorS = string(pszCursor);
		file_name_sequence.push_back( path_seq+pszCursorS);
		fscanf(in,  "%[\n]", pszCursorADD);
	}

	delete[] pszCursor;
	delete[] pszCursorADD;
	return 1;
}



bool compare_match(const CGroupMatch& a, const CGroupMatch& b)
{
		if (a.Match>b.Match)
			return true;
		else if (a.Match<b.Match)
			return false;
		else
			return a.Min_err<b.Min_err;
};


void funGetTimeSEK(double& sekD, int start)
{	
	int T=(GetTickCount()-start);
	sekD = T/1000.0;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
	unsigned long  T_start = GetTickCount();
	double sekTOT;

	map<string, L_CGrM> M_Group;
	int  i = 0;
	int MC = 0;

	
	string meta = "";
	string s_path_sum ="";
	string s_file_name_scoreclass = "scores";
	string s_name_group = "";
	

	string kmer = "0";
	string path_out = "";
	string path_score_IN = "";  

	string file_name_FILEscore ="name_files.txt";
	string file_name_FILEgroup ="name_groups.txt";
	string file_name_out = "Num_class.txt";


	list<string> Lfile_name_sequence;
	list<string> Ls_name_group;



	for(i = 1 ; i < argc; ++i)
	{
		if(argv[i][0] != '-')
			break;		
		//printf("\n %s\n", argv[i]);	
		

		if(strncmp(argv[i], "-FO", 3) == 0)
			file_name_out = (char*)(&argv[i][3]);
		else if(strncmp(argv[i], "-K", 2) == 0)
			kmer = (char*)(&argv[i][2]);
		else if(strncmp(argv[i], "-mc", 3) == 0)
			MC = atoi(&argv[i][3]);
		else if(strncmp(argv[i], "-NF", 3) == 0)
			file_name_FILEscore = (char*)(&argv[i][3]);
		else if(strncmp(argv[i], "-NG", 3) == 0)
			file_name_FILEgroup = (char*)(&argv[i][3]);
		else if(strncmp(argv[i], "-NO", 3) == 0)
			s_file_name_scoreclass = (char*)(&argv[i][3]);
		else if (strncmp(argv[i], "-WI", 3) == 0)
			path_score_IN =  (char*)(&argv[i][3]);
		else if (strncmp(argv[i], "-WO", 3) == 0)
			path_out =  (char*)(&argv[i][3]);
		else if (strncmp(argv[i], "-WS", 3) == 0)
			s_path_sum =  (char*)(&argv[i][3]);







	}

	if (path_out =="")
		path_out = path_score_IN;

	if (s_path_sum =="")
		s_path_sum = path_out;
	


		 
	if (!SetFileNamesList(path_score_IN+file_name_FILEscore, Lfile_name_sequence, path_score_IN))
		return 0;
	




	
	for( list<string>::iterator iter=Lfile_name_sequence.begin(); iter != Lfile_name_sequence.end(); )
	{
		printf("*");
		//cout << *iter;

		if (!Fread_Sin_score(*iter, M_Group))
		{
			cout << "\n ERROR";			
			return 0;
		}

		iter++;



	}

	printf("\nData is read \n");
	int num_class_seq = -1;
	num_class_seq = M_Group.size();
	printf("\nNum_class =%d \n", num_class_seq);

	
	


	printf("Save summation to file\n");
	FILE* pFile_out;
	string s_file_name_outScore = (s_path_sum + "/" + file_name_out);
	pFile_out = fopen(s_file_name_outScore.c_str(), "at");
	fprintf(pFile_out, "%s", s_file_name_scoreclass.c_str());

	fprintf(pFile_out, "\t%s", kmer.c_str());
	fprintf(pFile_out, "\tmc=%d", MC);
	fprintf(pFile_out, "\tnum=%d", num_class_seq);	
	fprintf(pFile_out, "\n");
	fclose(pFile_out);	



	printf("END\n");
}