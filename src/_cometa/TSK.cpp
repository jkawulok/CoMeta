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
#include <vector>
#include <time.h>
#include "Tr_KMclass.h"

#ifdef WIN32
#include <windows.h>
#endif

#ifdef _DEBUG
   #ifndef DBG_NEW
      #define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
      #define new DBG_NEW
   #endif
#endif  // _DEBUG


using namespace std;
string file_name_database;
string path_data;

TrKMclass *kmcTr;




void usage()
{
	cout << "Usage:\n TSK [options] \n";
	cout << "Options: \n";
	cout << "  -mr<size> - max amount of RAM in GB; default: 1.5\n";
	cout << "  -t<number> - total number of computation threads (default: equal to no. of system cores\n";
	cout << "  -NS<name> - name of file for out of kmers, score, \n";
	cout << "  -D<name> - name of database file \n";
	cout << "  -WK<name> - working directory for kmer database\n\n";
	cout << "  -k<len> - k-mer length\n";
	 

}






bool parse_parameters(TrKMclass *kmcTr, int argc, char *argv[])
{
	int i;
	
	int p_m = 2;
	int p_k = 0;
	int p_t = 0;

	bool p_f = 1; // fasta
	bool p_r = 0; // read
	string name_file = "";

	
	cout << "\nRead parameters\n";

	if(argc < 1)
		return false;

	for(i = 1 ; i < argc; ++i)
	{
		if(argv[i][0] != '-')
			break;
			
		
		if(strncmp(argv[i], "-mr", 3) == 0)
			p_m = atoi(&argv[i][3]);
		else if(strncmp(argv[i], "-t", 2) == 0)
			p_t = atoi(&argv[i][2]);
		else if (strncmp(argv[i], "-D", 2) == 0)
			file_name_database =  (string)(&argv[i][2]);
		else if (strncmp(argv[i], "-WK", 3) == 0)
			path_data =  (string)(&argv[i][3]);
		else if (strncmp(argv[i], "-NS", 3) == 0)
			name_file =  (string)(&argv[i][3]);
		else if(strncmp(argv[i], "-k", 2) == 0)
			p_k = atoi(&argv[i][2]);
		else if(  (strncmp(argv[i], "-h", 2) == 0) || (strncmp(argv[i], "help", 4) == 0) || (strncmp(argv[i], "-help", 5) == 0 ) )
			usage();
		
	}
	
	if (path_data != "")
		path_data+="/";

	if (file_name_database == "")	
	{

		char s_kC[20];
		sprintf(s_kC, "%d", p_k);
		string s_k = s_kC;
		file_name_database = "out_";
		file_name_database += (name_file + "_k" + s_k.c_str() + ".res");
	}


	p_r = 1;

	string sOutprefix = "Tr_";
	kmcTr->SetFileNames((path_data+file_name_database), path_data+sOutprefix+file_name_database);
	
	kmcTr->SetParams(p_t, p_m); // l threads,  mem_size

	return true;
}



int _tmain(int argc, _TCHAR* argv[])
{
	#ifdef WIN32
	_CrtSetDbgFlag( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);   //for Windows
	#endif

	//CStopWatch w0, w1;
	kmcTr = new TrKMclass;
	
	if(!parse_parameters(kmcTr, argc, argv))
	{		
		usage();
		return 0;
	}


	if(!kmcTr->Process())
	{
		cout << "Not enough memory or some other error\n";
		return 0;
	}
		
	kmcTr->GetStats();  	
	delete kmcTr;	



	return 0;
}

