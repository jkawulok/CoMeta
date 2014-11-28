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


// Program for adding tax number using gi number and file gi_taxid_nucl.dmp

//#include "targetver.h"
//#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <string>
#include <cstring>
#include <errno.h>
#include <cstdlib> 
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <assert.h>
#include <functional>  
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include "seqGI2tax.h"

#include <dirent.h>
#include <boost/lexical_cast.hpp>

using namespace std;


void usage()
{

	cout << "  -in<path_name> - path and file name for input file with sequence which include GI number - from I.2.1 point(e.g., ./NT/sequences_nt00.fa) \n";
	cout << "  -out<path_name> - path and file name for output file with sequence, where tax number is added (e.g., ./NT/sequences_nt00_TAX.fa) \n";
	cout << "  -filGT<name> - file name with relation between gi number and tax number (default: gi_taxid_nucl.dmp) \n";
	cout << "  -pGT<path> - path where is file with relation between gi number and tax number (path for gi_taxid_nucl.dmp file) \n";
	cout << "  -div<0/1> - dividing gi_taxid_nucl.dmp file (default 0). \n";

}





unsigned GetTickCount()
{
        struct timeval tv;
        if(gettimeofday(&tv, NULL) != 0)
                return 0;

        return (tv.tv_sec * 1000) + (tv.tv_usec / 1000);
}

//  --------------------------------------  counting time  -------------------------------------- 
void funGetTime(int& h, int& min, int& sek, int& milisek, int start)
{	
	int T=(GetTickCount()-start);

	h = floor(T/double(1000*60*60));
	min = floor((T-h*1000*60*60)/double(1000*60));
	sek = floor((T-h*1000*60*60-min*1000*60)/double(1000));
	milisek = T-h*1000*60*60-min*1000*60 - sek*1000;
}

//  -------------------------------------- 
void fun_write(char* name, char* seq, FILE*  pfile)
{
	int len = strlen(seq);
	
	// name
	fprintf(pfile, "%s\n", name);

	// sequence
	int step = 65;

	char * frag = new char[step+3];

	for (int jj=0; jj<len; jj+= step)
	{
		
		strncpy(frag, seq+jj, step);
		frag[step] = 0;
		fprintf(pfile, "%s\n", frag);
		
	}

	delete[] frag;
}


//  -------------------------------------- 
int ftakeGI(char* name)
{
	int len = strlen(name);
	char* c_gi = new char[len];

	for (int i=0; i<len; i++)
	{
		//if (name[i] == 'G' && name[i+1] == 'I' && name[i+2] == '=') 
		if (name[i] == 'g' && name[i+1] == 'i' && name[i+2] == '|') 
		{
			int c = 0;

			while (name[i+3+c] != '|') //			while (name[i+3+c] != ',')
			{
				c_gi[c] = name[i+3+c];
				c++;
			}
			c_gi[c] = 0;
			break;
		}
	}
	int gi = atoi(c_gi);
	delete[] c_gi;
	return gi;
}

//  -------------------------------------- 
bool fun_split_dmp(char* s_file_gi2tax, char* s_file_gi2tax_path)
{
	printf("Divide map to several smaller maps \n");
	int cut_lin = 15000000;
	//int cut_lin = 5500;
	FILE* pFile = NULL;
	string all_path =( string(s_file_gi2tax_path) + "/" + string(s_file_gi2tax));
	pFile = fopen(all_path.c_str(), "rt");
	if (pFile == NULL)
	{
		printf("\n Not exist file: '%s'", all_path.c_str());
		return 0;
	}

	else
	{
		fseek(pFile, 0L, SEEK_END);
		unsigned long long int nFileSize2 = ftell(pFile);
		fseek(pFile, 0L, SEEK_SET);
		int nL;
		int nSum=0;
		
		char* pszTextin = new char[nFileSize2+5];
		char* pszCursor = pszTextin;
		int num_seq=65;
		int len_p_name_seq = strlen(s_file_gi2tax);
		int len_path= strlen(s_file_gi2tax_path);
		char* p_name_seq_frag = new char[len_p_name_seq+25+len_path];



		int i_num = sprintf(p_name_seq_frag, "%s", s_file_gi2tax_path);
		i_num += sprintf(p_name_seq_frag + i_num, "%s", "/");
		i_num += sprintf(p_name_seq_frag + i_num, "%s", s_file_gi2tax);
		i_num += sprintf(p_name_seq_frag + i_num, "%s", "_part_");
		int i_num_stop=i_num;
		i_num += sprintf(p_name_seq_frag + i_num, "%c%s", num_seq, ".dmp");

		printf("\nMap: %c", num_seq);
		FILE* pfileTextBW = fopen(p_name_seq_frag, "wt"); 
		int num = 0;

		unsigned int id = 0;
		while(fscanf(pFile, "%ud", &id) != EOF)
		{
			bool nf = 0;
			
			unsigned int idTax = 0;
			fscanf(pFile, "%ud", &idTax);

			fprintf(pfileTextBW, "%u\t", id);
			fprintf(pfileTextBW, "%u\n", idTax);

			num++;
			if (num == cut_lin)
			{
				fclose(pfileTextBW);
				
				num_seq++;
				sprintf(p_name_seq_frag + i_num_stop, "%c%s", num_seq, ".dmp");
				printf("\nMap: %c", num_seq);
				FILE* pfileTextBW = fopen(p_name_seq_frag, "wt"); 
				num = 0;
			}
						
			
		
		}
		fclose(pfileTextBW);
		printf("\nMain file is dived into %d files \n\n", num_seq-65+1);
		string s_file_num=(string(s_file_gi2tax_path) + "/" + string(s_file_gi2tax)+"_NUMmap.txt");
		FILE* pFilenum = fopen(s_file_num.c_str(), "wt");
		fprintf(pFilenum, "%u\n", num_seq-65+1);
		fclose(pFilenum);

		delete []pszTextin;
		delete []p_name_seq_frag;
	}
	return 1;
}


//  -------------------------------------- 
bool fun_read_GI(char* s_file_gi2tax, TIDMap &mapGI)
{
	FILE* pFileGI = NULL;
	pFileGI = fopen(s_file_gi2tax, "rt");

	if (pFileGI == NULL)
	{
		perror("fopen");
		printf("%s \n", s_file_gi2tax);
		printf("Probably maps have not split. If you first start program, use the -div1 parameter\n");
		return 0;
	}
	
	fseek(pFileGI, 0L, SEEK_END);
	unsigned long long int nFileSizeGI = ftell(pFileGI);
	fseek(pFileGI, 0L, SEEK_SET);
	
	unsigned int id = 0;
	while(fscanf(pFileGI, "%ud", &id) != EOF)
	{
		unsigned int idTax = 0;
		fscanf(pFileGI, "%ud", &idTax);

		mapGI[id] = idTax;		

		if(id%4000000 == 0)
			printf("ID: %d \n", id);
	}
	fclose(pFileGI);
	return 1;
}


//  -------------------------------------- 
bool fun_read_write_seq(char* s_file_seq, char* s_file_out, TIDMap &mapGI)
{
	int nMaxLineSize = 1000000;
	char* pszCursor = new char[nMaxLineSize];
	FILE* pFile = NULL;
	pFile = fopen(s_file_seq, "rt");
	FILE* pFile_out = fopen(s_file_out, "wt");

	if (pFile == NULL)
	{
		perror("fopen");
		printf("%s \n", s_file_seq);
		return 0;
	}

	if (pFile_out == NULL)
	{
		perror("fopen");
		printf("%s \n", s_file_seq);
		return 0;
	}
	
	fseek(pFile, 0L, SEEK_END);
	unsigned long long int nFileSize = ftell(pFile);
	fseek(pFile, 0L, SEEK_SET);
	int nMAXseq = 100000000;

	bool first = true;

	char* name = new char[nMaxLineSize];
	char* sequence = new char[nMAXseq];
	char* name_TAX = new char[nMaxLineSize];
	int gi = -1;
	int nL = 0;
	
	int nuc_used = 0;	
	int jseq=0;
	long int num_pr=200000;

	int nDbg = 0;
	int nCursorRead = -1;
	while((nCursorRead = fscanf(pFile, "%s", pszCursor)) != EOF)
	{
		if(nCursorRead >= nMaxLineSize)
			printf("Error - line size exceeded (A): %d\n", nCursorRead);
		
		if (pszCursor[0]=='>')
		{
			nDbg++;
			if(nDbg%num_pr == 0)
				printf("Seq: %d \n", nDbg);

			if (!first)
			{
				int idTax = -1;
				TIDMap::iterator it = mapGI.find(gi);
				if(it != mapGI.end())
					idTax = (*it).second;

				if (name[strlen(name)-1] == '\r')
					name[strlen(name)-1] = 0;
				
	
					
					
				if (idTax > -1)
				{
					
					int dod = sprintf(name_TAX, "%s", name);
					dod += sprintf(name_TAX + dod, "%s", " <idTAX|");
					dod += sprintf(name_TAX + dod, "%d", idTax);
					dod += sprintf(name_TAX + dod, "%s", "| ");
					fun_write(name_TAX, sequence, pFile_out);

					if(dod >= nMaxLineSize)
						printf("Error - line size exceeded (B): %d\n", dod);
					
				}
				else
					fun_write(name, sequence, pFile_out);
			}
			else
				first = false;
	
			jseq=0;

			int j = sprintf(name, "%s", pszCursor );
			if(j >= nMaxLineSize)
				printf("Error - line size exceeded (Ba): %d\n", j);


			
			int ir = fscanf(pFile, "%[^\n]", pszCursor);
			if(ir >= nMaxLineSize)
				printf("Error - line size exceeded (C): %d\n", ir);
			if (ir!=0)
			{
				int lenCur = strlen(pszCursor);
				if(lenCur >= nMaxLineSize)
					printf("Error - line size exceeded (D): %d\n %s", lenCur, pszCursor);
				
				pszCursor[lenCur] = 0;
				j += sprintf(name + j, "%s", pszCursor); 

				if(j >= nMaxLineSize)
					printf("Error - line size exceeded (E): %d\n", j);
			}
			
			gi = ftakeGI(name);
				
			//if(nDbg >= 104788)
				//printf("B: %d \n", nDbg);

			
		}
		else
		{							
			nL = strlen(pszCursor);
			nuc_used += nL;
			pszCursor[nL] = 0;
			jseq+=sprintf(sequence + jseq, "%s", pszCursor ); 

			if(jseq >= nMAXseq)
				printf("Error - line size exceeded (F): %d\n", jseq);

		}

	}

	//printf("D ");

	int idTax = -1;
	TIDMap::iterator it = mapGI.find(gi);
	if(it != mapGI.end())
		idTax = (*it).second;

	//printf("E ");

	if (name[strlen(name)-1] == '\r')
		name[strlen(name)-1] = 0;

	if (idTax > -1)
	{
		
		int dod = sprintf(name_TAX, "%s", name);
		dod += sprintf(name_TAX + dod, "%s", " <idTAX|");
		dod += sprintf(name_TAX + dod, "%d", idTax);
		dod += sprintf(name_TAX + dod, "%s", "| ");
		
		if(dod >= nMaxLineSize)
			printf("Error - line size exceeded (A): %d\n", dod);
		fun_write(name_TAX, sequence, pFile_out);
	}
	else
		fun_write(name, sequence, pFile_out);
				
	//printf("F ");

	fclose(pFile);
	fclose(pFile_out);

	//printf("G ");

	delete [] sequence;
	delete [] name;
	delete [] pszCursor;
	delete [] name_TAX;

	//printf("H ");

	return 1;
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
//#ifdef _DEBUG
	//getchar();
//#endif
	
	bool b_divide = 0;
	long int  Start = GetTickCount();		
	long int  StartALL = GetTickCount();		
	int  h, min, sek, milisek;

	char* s_file_seq = (char*)"";
	char* s_file_outEND  = (char*)"";

	char* s_file_gi2tax = (char*)"gi_taxid_nucl.dmp";
	char* s_file_gi2tax_path = (char*)"";

	if(argc > 1)
	{

		for (int i=1; i<argc; i++)
		{

			if(strncmp(argv[i], "-in", 3) == 0)
				 s_file_seq = (char*)(&argv[i][3]);
			else if (strncmp(argv[i], "-out", 4) == 0)
				 s_file_outEND = (char*)(&argv[i][4]);
			else if (strncmp(argv[i], "-pGT", 4) == 0)
				 s_file_gi2tax_path = (char*)(&argv[i][4]);
			else if (strncmp(argv[i], "-filGT", 6) == 0)  // the file  gi_taxid_nucl.dmp
				 s_file_gi2tax = (char*)(&argv[i][4]);
			else if (strncmp(argv[i], "-div", 4) == 0)
				b_divide = atoi(&argv[i][4]);
			else if(  (strncmp(argv[i], "-h", 2) == 0) || (strncmp(argv[i], "help", 4) == 0) || (strncmp(argv[i], "-help", 5) == 0 ) )
				usage();
		}
	}


	if (b_divide)
	{
		fun_split_dmp(s_file_gi2tax, s_file_gi2tax_path);
		return 0;
	}


	printf("\nIN %s \n", s_file_seq);

	if (s_file_outEND == "")
	{
		printf("No parameter \"-out\" \n");
		return 0;
	}
	printf("OUT %s \n", s_file_outEND);

	FILE* pFile_IN = NULL;
	pFile_IN = fopen(s_file_seq, "rt");
	if (pFile_IN == NULL)
	{
		perror("fopen");
		printf("%s \n", s_file_seq);
		return 0;
	}
	fclose(pFile_IN);


	char* s_file_in = s_file_seq;
	char* s_file_out_tempA = new char[strlen(s_file_seq) + 20];
	strcpy(s_file_out_tempA, s_file_seq);
	strcat(s_file_out_tempA, "_seq_prob_XYZ_A");

	char* s_file_out_tempB = new char[strlen(s_file_seq) + 20];
	strcpy(s_file_out_tempB, s_file_seq);
	strcat(s_file_out_tempB, "_seq_prob_XYZ_B");

	
	
	char* s_file_out = s_file_out_tempA;
	int nMaxLineSize_namefile = 1000;
	char* s_file_gi2tax_nr = new char[nMaxLineSize_namefile];

	char c_end = 0;
	{
		string s_file_num=string(string(s_file_gi2tax_path) + "/" + string(s_file_gi2tax));
		FILE* pFilenum = fopen(s_file_num.c_str(), "rt");
		if (pFilenum == NULL)
		{
			perror("fopen");
			printf("%s \n", s_file_seq);
			return 0;
		}
		fclose(pFilenum);
		s_file_num=(string(s_file_gi2tax_path) + "/" + string(s_file_gi2tax)+"_NUMmap.txt");
		pFilenum = fopen(s_file_num.c_str(), "rt");
		if (pFilenum == NULL)
		{
			perror("fopen");
			printf("%s \n", s_file_gi2tax);
			printf("Probably maps have not split. If you first start program, use the -div1 parameter\n");
			return 0;
		}
	
		int num_map;
		fscanf(pFilenum, "%ud", &num_map);
		fclose(pFilenum);
	
		c_end='A'+num_map-1;
	}

	for (char nr_TAX = 'A'; nr_TAX <= c_end; nr_TAX++)
	{
		Start = GetTickCount();	
		printf("Map %c \n", nr_TAX);
		
		int i_num = sprintf(s_file_gi2tax_nr, "%s", s_file_gi2tax_path);
		i_num += sprintf(s_file_gi2tax_nr + i_num, "%s", "/");
		i_num += sprintf(s_file_gi2tax_nr + i_num, "%s", s_file_gi2tax);
		i_num += sprintf(s_file_gi2tax_nr + i_num, "%s", "_part_");
		int i_num_stop=i_num;
		i_num += sprintf(s_file_gi2tax_nr + i_num, "%c%s", nr_TAX, ".dmp");



		if(i_num >= nMaxLineSize_namefile)
			printf("Error - line size exceeded (X): %d\n", i_num);
					
		
		// reading file with taxid
		TIDMap mapGI;		
		if (!fun_read_GI(s_file_gi2tax_nr, mapGI))
			return 0;

		printf("is readed\n");
	
		// reading file with sequence
		if (!fun_read_write_seq(s_file_in, s_file_out, mapGI))
			return 0;

		if (nr_TAX != 'A' && (nr_TAX+1) != c_end)
		{
			char* s_file_in_old = s_file_in;
			s_file_in = s_file_out;
			s_file_out = s_file_in_old;

		}
		else if (nr_TAX == 'A')
		{
			s_file_in = s_file_out_tempA;
			s_file_out =  s_file_out_tempB;
		}
		else if (nr_TAX+1 == c_end)
		{
			char* s_file_in_old = s_file_in;
			s_file_in = s_file_out;
			s_file_out = s_file_outEND;
		}
		else if (nr_TAX >= c_end)
			printf("Error \n");

		funGetTime( h, min,sek, milisek, Start);
		printf("Counting for single map: \t %dh %02dmin %02dsek %03dmilisek \n\n", h, min, sek, milisek);
	}	
	
	
	remove(s_file_out_tempA);
	remove(s_file_out_tempB);
	
	delete [] s_file_gi2tax_nr;
	delete [] s_file_out_tempA;
	delete [] s_file_out_tempB;

	s_file_gi2tax_nr = NULL;
	s_file_out_tempA=NULL;
	s_file_out_tempB=NULL;


	funGetTime( h, min,sek, milisek, StartALL);
	printf("\nAll time\t %dh %02dmin %02dsek %03dmilisek \n\n", h, min, sek, milisek);

	printf("\nEND \n");
	
	return 0;
}


