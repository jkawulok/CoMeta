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
#include "timer.h"
#include "cometa.h"
#include "kmclass.h"


#ifdef _DEBUG
   #ifndef DBG_NEW
      #define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
      #define new DBG_NEW
   #endif
#endif  // _DEBUG


using namespace std;

uint64 total_reads, total_fastq_size;

void usage();
bool parse_parameters(CCOMETA *cometa, CKMclass *kmclass, int argc, char *argv[]);

list<string> file_name_sequence;
string file_name_FILEsequence = "";
string file_name_database;
string file_name_scoreclassM;
string file_name_scoreclassMM;
string working_directories;
string path_seq;
string path_data;
string path_outs;



CCOMETA *cometa;
CKMclass *kmclass;


void usage()
{
	cout << "\nClassifiication of metagenomes (CoMeta)\n";
	cout << "Usage:\n cometa [options] \n";
	cout << "Options: (R! - requested)\n";
	cout << "  -go<B/C> - if program build database (B) or classificate (C) R!\n\n"; 
	cout << "For both option (B and C)\n"; 
	cout << "  -mr<size> - max amount of RAM in GB; \n";
	cout << "  -t<number> - total number of computation threads (default: equal to no. of system cores\n";
	cout << "  -sp<number> - number of splitting threads\n";
	cout << "  -NS<name> - name of file for out of kmers, score, \n";
	cout << "  -S<name> - name of input file with sequence or reads in FASTA format \n";
	cout << "  -D<name> - name of database file (output in 'bulding', input in 'classification') \n";
	cout << "  -OS<name> - name of file summarize (default score_build.txt for <B> or score_class.txt for <C> \n";
	cout << "  -WS<name> - path for data witch includes sequences\n";
	cout << "  -WK<name> - path for kmer databases\n\n";

	cout << "For building\n"; 
	cout << "  -k<len> - k-mer length (k from 1 to 32); default: 24\n";
	cout << "  -L<name> - file name with list of input files with sequence in FASTA format\n"; 
	cout << "  -so<number> - number of sorting threads\n";
	cout << "  -sr<number> - number of threads per single sorter\n";
	cout << "  -WD<name> - working directory \n";
	cout << "  For building parameter \"L\" or \"S\" is requested! \n\n";

	cout << "For classification\n"; 
	cout << "  -sc<number> - number of classification threads\n";
	cout << "  -stepk<len> - step (stepk from 1 to k-mer length)\n";
	cout << "  -sK<0/1> - 1 if kmer database is sorted (after TSK application), 0 if not\n";
	cout << "  -mc<number> - the percent identity to classify a match [%] \n";
	cout << "  -OC<name> - name of output after classification (matching)\n";
	cout << "  -ON<name> - name of output after classification (mismatching)\n";
	cout << "  -WO<name> - working directory for outputs\n\n";
	

	cout << "Example:\n";
	cout << "cometa -goB -k20  -Schr01.fa -Ddata_chr01.res -WD.\\_BINS\n";
	cout << "cometa -goB -k20  -LlistofNames.txt -Ddata_mix.res -WD.\\_BINS\n";
	cout << "cometa -goC -mc10  -Sreads.fa -Ddata_chr01.res -WS.\\data\\sequence\\ -WK.\\data\\data_kmer\\ -WO.\\data\\output\\ \n";
}

//----------------------------------------------------------------------------------
bool parse_parameters(CCOMETA *cometa, CKMclass *kmclass, int argc, char *argv[], char &p_go, string &file_name_outSumarize)
{
	int i;
	
	bool p_sortKmer = 1;
	int p_m = 2;
	int p_k = 0;
	int p_stepk = 0;
	int p_t = 0;
	int p_sp = 0;
	int p_so = 0;
	int p_sr = 0;
	int p_sc = 0;
	double p_mc = 0.3;

	bool p_f = 1; // fasta
	bool p_r = 0; // read
	string name_file = "";

	
	cout << "\nRead parameters\n";

	if(argc < 4)
		return false;

	for(i = 1 ; i < argc; ++i)
	{
		if(argv[i][0] != '-')
			break;
		
		
		
		if(strncmp(argv[i], "-t", 2) == 0)
			p_t = atoi(&argv[i][2]);
		else if(strncmp(argv[i], "-k", 2) == 0)
			p_k = atoi(&argv[i][2]);
		else if (strncmp(argv[i], "-stepk", 6) == 0)
			p_stepk = atoi(&argv[i][6]);
		else if(strncmp(argv[i], "-mr", 3) == 0)
			p_m = atoi(&argv[i][3]);
		else if(strncmp(argv[i], "-sp", 3) == 0)
			p_sp = atoi(&argv[i][3]);
		else if(strncmp(argv[i], "-so", 3) == 0)
			p_so = atoi(&argv[i][3]);
		else if(strncmp(argv[i], "-sr", 3) == 0)
			p_sr = atoi(&argv[i][3]);
		else if(strncmp(argv[i], "-sc", 3) == 0)
			p_sc = atoi(&argv[i][3]);
		else if(strncmp(argv[i], "-go", 3) == 0)
			p_go = char(argv[i][3]);
		else if(strncmp(argv[i], "-sK", 3) == 0)
			p_sortKmer = atoi(&argv[i][3]);
		else if(strncmp(argv[i], "-mc", 3) == 0)
			p_mc = (atoi(&argv[i][3]))/100.0;
		else if (strncmp(argv[i], "-S", 2) == 0)
			file_name_sequence.push_back((string)(&argv[i][2]));
		else if (strncmp(argv[i], "-L", 2) == 0)
			file_name_FILEsequence =  (string)(&argv[i][2]);
		else if (strncmp(argv[i], "-D", 2) == 0)
			file_name_database =  (string)(&argv[i][2]);
		else if (strncmp(argv[i], "-OC", 3) == 0)	
			file_name_scoreclassM = (string)(&argv[i][3]);
		else if (strncmp(argv[i], "-ON", 3) == 0)
			file_name_scoreclassMM =  (string)(&argv[i][3]);
		else if (strncmp(argv[i], "-OS", 3) == 0)
			file_name_outSumarize =  (string)(&argv[i][3]);		
		else if (strncmp(argv[i], "-WD", 3) == 0)
			working_directories =  (string)(&argv[i][3]) + "/";
		else if (strncmp(argv[i], "-WS", 3) == 0)
			path_seq =  (string)(&argv[i][3]) + "/";
		else if (strncmp(argv[i], "-WK", 3) == 0)
			path_data =  (string)(&argv[i][3]) + "/";
		else if (strncmp(argv[i], "-WO", 3) == 0)
			path_outs =  (string)(&argv[i][3]) + "/";
		else if (strncmp(argv[i], "-NS", 3) == 0)
			name_file =  (string)(&argv[i][3]);
		else if(  (strncmp(argv[i], "-h", 2) == 0) || (strncmp(argv[i], "help", 4) == 0) || (strncmp(argv[i], "-help", 5) == 0 ) )
			usage();

		//else if(strncmp(argv[i], "-f", 2) == 0)
		//	p_f = atoi(&argv[i][2]);
		
	}
	

	if (p_go == 0)
		return 0;

	if (p_go == 'B' && p_k == 0)
		p_k = 24;

	if (file_name_database == "")	
	{
		
		char s_kC[20];
		sprintf(s_kC, "%d", p_k);
		string s_k = s_kC;
		
		
		if (p_sortKmer && p_go =='C')
			file_name_database = "Tr_out_";
		else
			file_name_database = "out_";

		
		file_name_database += (name_file + "_k" + s_k.c_str() + ".res");
	}




	if (p_go == 'B')
	{

		if ( (file_name_sequence.empty() == 1 && file_name_FILEsequence == "") || file_name_database == "" || working_directories == "")
		{

			if (file_name_database == "")
				printf("You must give name of database file (-D)\n");
			if (working_directories == "")
				printf("You must give working directory  (-WD)\n");
			if (file_name_sequence.empty() == 1 || file_name_FILEsequence == "")
				printf("You must give name of input file with sequence -S) or name file of list with names of input file with sequence (-L) \n");
			
			return 0;
		}

		if (file_name_FILEsequence != "")
			cometa->SetFileNamesList(file_name_FILEsequence, file_name_sequence, path_seq);
		else
			file_name_sequence.front() = (path_seq + file_name_sequence.front());

		

		cometa->SetFileNames(file_name_sequence, (path_data+file_name_database), working_directories, p_f, p_r);				// true = FASTA file ////// zmiana na zmienna p_f
	

		if(p_sp && p_so && p_sr)
			cometa->SetParams(p_k, p_sp, p_so, p_sr, p_m);		// k, n_splitters, n_sorters, omp_threads, mem_size
		else
			cometa->SetParams(p_k, p_t, p_m);
	}
	else if (p_go == 'C')
	{

		if (file_name_sequence.empty() == 1)
		{
			printf("Name of sequence (-S) is not known \n");
			return 0;
		}

		if (file_name_scoreclassM == "" || file_name_scoreclassMM == "")
		{			
			if (name_file == "" && p_k == -1)
			{
				if (file_name_scoreclassM == "")
					file_name_scoreclassM = "out_classM.txt";
				if (file_name_scoreclassMM == "")
					file_name_scoreclassMM = "out_classMM.txt";

			}
			else
			{
				
				char s_mcC[20];
				int p_mcINT =  round(p_mc*100);
				sprintf(s_mcC, "%d", p_mcINT);
				string s_mc = s_mcC;
			
				string name_seq = file_name_sequence.front();
				basic_string <char>::size_type indexCh = 0;
				indexCh	= name_seq.find_first_of ( ".");
				static const basic_string <char>::size_type npos = -1;
				if ( indexCh == npos )
					indexCh = name_seq.size();

				string s_add =  "_";
				s_add += s_mc.c_str();
				name_seq.replace(indexCh,  name_seq.size(), s_add);  // place start, place end , name file 


				char s_kC[20];
				sprintf(s_kC, "%d", p_k);
				string s_k = s_kC;

				cout << "\n s_k1 =  ";	
				cout << s_k.c_str() << "\n";

				//cout << "\n name_seq =  ";	
				//cout << name_seq.c_str() << "\n";
				cout << "\n name_file =  ";	
				cout << name_file.c_str() << "\n";

				
				if (file_name_scoreclassM == "")
					file_name_scoreclassM = (name_file + "_k" + s_k.c_str() + "_" + name_seq + "_M.txt");

				if (file_name_scoreclassMM == "")
					file_name_scoreclassMM = (name_file + "_k" + s_k.c_str() + "_" + name_seq + "_MM.txt");		
			
			}
		}
		
		
		
		p_r = 1;

		list<string> file_name_sequenceP;
		file_name_sequenceP.push_back(path_seq + file_name_sequence.front());

		kmclass->SetFileNames(file_name_sequenceP, (path_data+file_name_database), path_outs+file_name_scoreclassM, path_outs+file_name_scoreclassMM, p_mc, p_sortKmer);
	
		if(p_sp && p_sc)
			kmclass->SetParams(p_stepk, p_sp, p_sc, p_m);		//n_splitters, n_classification
		else
			kmclass->SetParams(p_stepk, p_t, p_m); // l threads,  mem_size

	}
	else
		return 0;

	return true;
}


//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
typedef unsigned char uint56[5];

int _tmain(int argc, _TCHAR* argv[])
{

	#ifdef WIN32
	_CrtSetDbgFlag( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);  //for windows
	#endif
	


	CStopWatch w0, w1;
	bool build_seq = 1;
	bool build_classification = 0;
	kmclass = new CKMclass;
	cometa = new CCOMETA;
	
	char p_go = 'x'; // B - build, C - classification

	string file_name_outSumarize = "";
	if(!parse_parameters(cometa, kmclass, argc, argv, p_go, file_name_outSumarize))
	{		
		usage();
		return 0;
	}



	if (p_go == 'B')
	{	
		delete kmclass;
		double time1, time2;
		uint64 n_unique, n_singletons, n_total, n_reads, tmp_size;
		uint32 n_cores;
		uchar kmer_len;
		uint32 n_splitters,n_sorters, n_omp_th, max_mem_size_in_gb;

		if(!cometa->Process())
			cout << "Not enough memory or some other error\n";

		cometa->GetStats(time1, time2, n_unique, n_singletons, n_total, n_reads, tmp_size, kmer_len, n_cores, n_splitters, n_sorters, n_omp_th, max_mem_size_in_gb);

		cout << "1st stage: " << time1 << "s\n";
		cout << "2nd stage: " << time2  << "s\n";
		cout << "Total    : " << (time1+time2) << "s\n";
		cout << "Tmp size : " << tmp_size / (1 << 20) << "MB\n";
		cout << "\nStats (for kmer length = " << uint32(kmer_len)  << "):\n";
		cout << "   No. of singleton k-mers            : " << setw(12) << n_singletons << "\n";
		cout << "   No. of unique k-mers               : " << setw(12) << n_unique << "\n";
		cout << "   No. of unique non-singleton k-mers : " << setw(12) << n_unique-n_singletons << "\n";
		cout << "   Total no. of k-mers                : " << setw(12) << n_total << "\n";
		cout << "   Total no. of reads                 : " << setw(12) << n_reads << "\n";



		string path_ALL_outs;		
		if (file_name_outSumarize == "")
			path_ALL_outs = path_outs + "/" + "score_build.txt";
		else
			path_ALL_outs = path_outs + "/" + file_name_outSumarize;
		
		FILE* pfile_out_class = fopen(path_ALL_outs.c_str(), "at");
		
		fprintf(pfile_out_class, "\n%s\t", (file_name_database.c_str()) );
		fprintf(pfile_out_class, "%.2f\t", time1+time2); //time all
		fprintf(pfile_out_class, "%.2f\t%.2f\t", time1, time2);
		fprintf(pfile_out_class, "k=%d\t", kmer_len);
		fprintf(pfile_out_class, "%llu\t", n_total); //Total no. of k-mers
		fprintf(pfile_out_class, "%llu\t", n_unique); //No. of unique k-mers 
		fprintf(pfile_out_class, "%llu\t", (tmp_size / (1 << 20)));  //Tmp size
		fprintf(pfile_out_class, "%d\t", n_cores);  //number of cores
		fprintf(pfile_out_class, "%llu\t", n_reads);  //number of reads
		fprintf(pfile_out_class, "%d\t", n_splitters);  
		fprintf(pfile_out_class, "%d\t", n_sorters);  
		fprintf(pfile_out_class, "%d\t", n_omp_th);  
		fprintf(pfile_out_class, "%d\t", max_mem_size_in_gb);  
		fclose(pfile_out_class);


		delete cometa;

	}


	if (p_go == 'C')
	{
		delete cometa;
		double time1, time2;
		uint64 num_alig_reads, num_REValig_reads, num_NOalig_reads, n_reads;
		uchar kmer_len;
		double matchcutoff;
		uint32 n_cores, n_splitters, n_classification, max_mem_size_in_gb;


		if(!kmclass->Process())
		{
			cout << "Not enough memory or some other error\n";
			return 0;
		}
		
		kmclass->GetStats(time1, time2, num_alig_reads, num_REValig_reads, num_NOalig_reads, n_reads, kmer_len, matchcutoff, n_cores, n_splitters, n_classification, max_mem_size_in_gb);  

		int mcp =  round(matchcutoff*100); 
		
		cout << "\n\nRead KMERS: " << time1 << "s\n";
		cout << "Read 'reads' and classification: " << time2  << "s\n";
		cout << "Total    : " << (time1+time2) << "s\n";

		printf("\nStats (for kmer length = %d, and match cut off = %d%% ) \n", kmer_len, mcp) ;
		cout << "   No. of alignment reads                : " << setw(12) << num_alig_reads << "\n";
		cout << "   No. of complementary alignment reads  : " << setw(12) << num_REValig_reads << "\n\n";

		cout << "   No. of all alignment reads            : " << setw(12) << (num_alig_reads+num_REValig_reads) << "\n";
		cout << "   No. of no alignment reads             : " << setw(12) << num_NOalig_reads << "\n";
		cout << "   Total no. of reads                    : " << setw(12) << n_reads << "\n";
		cout << "   No. of cores	                      : " << setw(12) << n_cores << "\n";


		string path_ALL_outs;
		if (file_name_outSumarize == "")
			path_ALL_outs = path_outs + "/" + "score_class.txt";
		else
			path_ALL_outs = path_outs + "/" + file_name_outSumarize;

		FILE* pfile_out_class = fopen(path_ALL_outs.c_str(), "at");

		fprintf(pfile_out_class, "\n%s\t %s\t", (file_name_sequence.front()).c_str(), file_name_database.c_str()  );
		fprintf(pfile_out_class, "mc=%.2f\t",matchcutoff);
		fprintf(pfile_out_class, "k=%d\t", kmer_len);
		fprintf(pfile_out_class, "%.2f\t", time1+time2); //time all
		fprintf(pfile_out_class, "%.2f\t%.2f\t", time1, time2);
		fprintf(pfile_out_class, "%lu\t", (num_alig_reads+num_REValig_reads));
		fprintf(pfile_out_class, "%lu\t", num_NOalig_reads);
		fprintf(pfile_out_class, "%d\t", n_cores);  //number of cores
		fprintf(pfile_out_class, "%d\t", n_splitters);  
		fprintf(pfile_out_class, "%d\t", n_classification);  
		fprintf(pfile_out_class, "%d\t", max_mem_size_in_gb);  

		fclose(pfile_out_class);

		
		delete kmclass;		
	}

	//_CrtDumpMemoryLeaks();

	return 0;
}

// ***** EOF

