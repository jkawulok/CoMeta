/*
	This file is a part of CoMeta software.
	The homepage of the CoMeta project is http://sun.aei.polsl.pl/cometa
 
	Authors: Jolanta Kawulok and Sebastian Deorowicz
	
	Version: 0.3
	Date   : 2014-11-26

*/

// class2best
// After comparison with low value of MC, algorithm chose the best match

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



void usage()
{
	cout << "  -NO<name> - prefix name of output after classification\n";
	cout << "  -NR<name> - file name, which contains read\n";
	cout << "  -NF<name> - file name which includes names of input files (which inlude scores)	\n";
	cout << "  -NG<name> - file name, which contains names of groups where reads were classified\n";
	cout << "  -WI<name> - path where there are files with the results of the comparison step (match, mistamtch)\n";
	cout << "  -WO<name> - path for output data \n";
	cout << "  -WS<name> - path for summarization file\n";
	cout << "  -cl<-1/0/1> - when reads are classified to a few groups, then reads are assignment to: -1 - any group; 0 - random group; 1 -   all of these groups \n" ;
	cout << "  -proc<number> - similarity the best results, which are taken into account; default 100%] \n";	
	
	cout << "\nParameters only for display them in the results file\n";
	cout << "  -k<len> - k-mer length\n";
	cout << "  -mc<number> - the percent identity to classify a match [%] \n";
	cout << "  -MM<0/1> - taking into the account the mismatch file; 0-NO, 1-YES\n";
	


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



bool Fread_Sin_score(string s_file_name_input, string s_GroupName, 	map<string, L_CGrM>& M_Group)
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
		grmat.GroupName = s_GroupName;

		fscanf(pFile_in, "%[^\n]", name_line);
		name_lineS = name_line;
		id_fe = name_lineS.find(9, 3);		
		int id_st = 0;
		while (name_line[id_st] > 8 && name_line[id_st] < 33)
		{
			id_st++;						
		}		
		grmat.SequenceName = name_lineS.substr(id_st, name_lineS.length());		
		seq_name = name_lineS.substr(id_st, id_fe);
		//grmat.SequenceName = seq_name;
		

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

void fPrintScore(FILE* pFile_out, 	map<string, L_CGrM>& M_Group, string& ind_name_seq, string& name_seq)
{
		fprintf(pFile_out, "%s\t", "Min_err:"); 
		fprintf(pFile_out, "%d\t",  (M_Group[ind_name_seq].begin())->Min_err ); //min_err
		fprintf(pFile_out, "%s\t", "Match:"); 
		fprintf(pFile_out, "%d\t",  (M_Group[ind_name_seq].begin())->Match ); //match
		fprintf(pFile_out, "%s\t", name_seq.c_str()); 
}

void fPrintScoreDUP(FILE* pFile_out, CGroupMatch data)
{
		fprintf(pFile_out, "%s\t", "Min_err:"); 
		fprintf(pFile_out, "%d\t",  data.Min_err);  //(M_Group[ind_name_seq].begin())->Min_err ); //min_err
		fprintf(pFile_out, "%s\t", "Match:"); 
		fprintf(pFile_out, "%d\t",   data.Match ); //match
		fprintf(pFile_out, "%s\t", data.SequenceName.c_str()); 
		fprintf(pFile_out, "%s%s%s\n", "<class|", data.GroupName.c_str() , "|") ; 


}

void funGetTimeSEK(double& sekD, int start)
{	
	int T=(GetTickCount()-start);
	sekD = T/1000.0;
}

bool read_seq_query(string file_query, map<string, str_seq>& M_Reads)
{
	printf("\nREADING 'reads'\n");
	int nL;
	FILE* pFile = NULL;
	pFile = fopen(file_query.c_str(), "rt");
	int nr_seq = 0;

	if (pFile == NULL)
	{
		perror("fopen");
		cout << file_query << "\n";
		return 0;
	}


	//FILE* pFile_outNOTclass = fopen("path/repeat_reads.txt", "wt");  //// check if it is name reads repeatet

	
	fseek(pFile, 0L, SEEK_END);
	unsigned long long int nFileSize = ftell(pFile);
	fseek(pFile, 0L, SEEK_SET);
	
	//sQuerys->SEQ = new str_seq[nFileSize];
		
	char* pszCursor = new char[1000];
	bool first = true;

	str_seq sinSEQ;
	sinSEQ.name = new char[1000];
	sinSEQ.sequence = new char[nFileSize+2];	
	
	unsigned int nuc_used = 0;	
	unsigned int jseq=0;

	while(fscanf(pFile, "%s", pszCursor) != EOF)
	{

		if (pszCursor[0]=='>')
		{
			nr_seq++;
			if (!first)
			{
				// remove ending name
				string name_seq = string(sinSEQ.name);
				int len_name =  name_seq.find("<name|");	
				
				if (len_name < 0)
					len_name = strlen(sinSEQ.name);
				else
					sinSEQ.name[len_name]=0;
				
				while (sinSEQ.name[len_name-1] > 8 && sinSEQ.name[len_name-1] < 33)
				{
					len_name--;
					sinSEQ.name[len_name]=0;						
				}
				

				//if (!(M_Reads.find(sinSEQ.name) == M_Reads.end()))   // check if it is name reads repeatet
				//	fprintf(pFile_outNOTclass, "%s\n", sinSEQ.name ) ; 

				M_Reads[sinSEQ.name].name = new char[len_name + 2];
				strcpy( M_Reads[sinSEQ.name].name, sinSEQ.name);

				M_Reads[sinSEQ.name].sequence = new char[jseq + 10];
				strcpy(M_Reads[sinSEQ.name].sequence, sinSEQ.sequence);
				nuc_used = 0;
			}
			else
				first = false;
	
			jseq=0;			
			int j = sprintf(sinSEQ.name, "%s", pszCursor );
			int ir = fscanf(pFile, "%[^\n]", pszCursor);
			if (ir!=0)
			{
				int lenCur = strlen(pszCursor);
				pszCursor[lenCur]=0;
				j += sprintf(sinSEQ.name+j, "%s", pszCursor ); 
			}

		}
		else
		{							
			nL = strlen(pszCursor);
			nuc_used += nL;
			pszCursor[nL] = 0;
			if (jseq == 0)
			{
				jseq+=sprintf(sinSEQ.sequence + jseq, "%s", pszCursor ); 
			}
			else
			{
				jseq+=sprintf(sinSEQ.sequence + jseq, "%c%s", '\n', pszCursor ); 
			}
		}
	}

	
	string name_seq = string(sinSEQ.name);
	int len_name =  name_seq.find("<name|");	
				
	if (len_name < 0)
		len_name = strlen(sinSEQ.name);
	else
		sinSEQ.name[len_name]=0;
				
	while (sinSEQ.name[len_name-1] > 8 && sinSEQ.name[len_name-1] < 33)
	{
		len_name--;
		sinSEQ.name[len_name]=0;						
	}
	

	//if (!(M_Reads.find(sinSEQ.name) == M_Reads.end()))  // check if it is name reads repeatet
	//	fprintf(pFile_outNOTclass, "%s\n", sinSEQ.name ) ; 

	M_Reads[sinSEQ.name].name = new char[strlen(sinSEQ.name) + 2];
	strcpy( M_Reads[sinSEQ.name].name, sinSEQ.name);

	M_Reads[sinSEQ.name].sequence = new char[jseq + 10];
	strcpy(M_Reads[sinSEQ.name].sequence, sinSEQ.sequence);
	

	delete[] sinSEQ.name;
	delete[] sinSEQ.sequence;
	delete[] pszCursor;

	printf("%d sequnces are in the metagenomic set", nr_seq);
	return 1;
	
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
	map<string, string> M_Trans;
	map<string, str_seq> M_Reads;
	map<string, FILE*> M_File_class;
	map<int, int> M_count_dup;

	double i_procMatch = 1;
	int step = 80; //
	int  i = 0;
	int MC = 0;
	int mismatch = -1;
	
	string meta = "";
	string s_path_sum ="";
	string s_file_name_scoreclass = "scores";
	string s_name_group = "";
	string s_file_name_reads = ""; 

	string kmer = "0";
	string path_out = "";
	string path_score_IN = "";  

	
	string file_name_FILEscore ="name_files.txt";
	string file_name_FILEgroup ="name_groups.txt";
	int classALL = 0;
	bool check_class = 0;

	list<string> Lfile_name_sequence;
	list<string> Ls_name_group;

	vector<int> v_num_score(4,0); // TP FP NotKnow  Multi


	for(i = 1 ; i < argc; ++i)
	{
		if(argv[i][0] != '-')
			break;		
		//printf("\n %s\n", argv[i]);	
		
		if(  (strncmp(argv[i], "-h", 2) == 0) || (strncmp(argv[i], "help", 4) == 0) || (strncmp(argv[i], "-help", 5) == 0 ) )
			usage();
		else if(strncmp(argv[i], "-NO", 3) == 0)
			s_file_name_scoreclass = (char*)(&argv[i][3]);
		else if(strncmp(argv[i], "-NR", 3) == 0)
			s_file_name_reads = (char*)(&argv[i][3]);
		else if(strncmp(argv[i], "-NF", 3) == 0)
			file_name_FILEscore = (char*)(&argv[i][3]);
		else if(strncmp(argv[i], "-NG", 3) == 0)
			file_name_FILEgroup = (char*)(&argv[i][3]);
		else if (strncmp(argv[i], "-WI", 3) == 0)
			path_score_IN =  (char*)(&argv[i][3]);
		else if (strncmp(argv[i], "-WO", 3) == 0)
			path_out =  (char*)(&argv[i][3]);
		else if (strncmp(argv[i], "-WS", 3) == 0)
			s_path_sum =  (char*)(&argv[i][3]);
		else if(strncmp(argv[i], "-k", 2) == 0)
			kmer = (char*)(&argv[i][2]);
		else if(strncmp(argv[i], "-mc", 3) == 0)
			MC = atoi(&argv[i][3]);
		else if(strncmp(argv[i], "-cl", 3) == 0)
			classALL = atoi(&argv[i][3]);
		else if(strncmp(argv[i], "-ch", 3) == 0)
			check_class = atoi(&argv[i][3]);
		else if(strncmp(argv[i], "-proc", 5) == 0)
			i_procMatch = (atoi(&argv[i][5])) * 0.01;
		else if(strncmp(argv[i], "-MM", 3) == 0)
			mismatch = (atoi(&argv[i][3]));



	}

	if (path_out =="")
		path_out = path_score_IN;

	if (s_path_sum =="")
		s_path_sum = path_out;
	
	path_out+="/";
	s_path_sum+="/";
	path_score_IN+="/";
	

	if (!read_seq_query(s_file_name_reads, M_Reads))
		return 0;
		 
	if (!SetFileNamesList(path_score_IN+file_name_FILEscore, Lfile_name_sequence, path_score_IN))
		return 0;
	
	if (!SetFileNamesList(path_score_IN+file_name_FILEgroup, Ls_name_group, ""))
		return 0;

	if (Lfile_name_sequence.size() != Ls_name_group.size())
	{
		printf("Length of list of sequence names is to equals  length of list of group names\n");
		return 0;
	}

	list<string>::iterator iterG = Ls_name_group.begin();
	for( list<string>::iterator iter=Lfile_name_sequence.begin(); iter != Lfile_name_sequence.end(); )
	{
		printf("*");
		//cout << *iter;

		if (!Fread_Sin_score(*iter, *iterG, 	M_Group))
		{
			cout << "\n ERROR";			
			return 0;
		}

		iter++;
		iterG++;


	}

	printf("\nData is read \n");
	int icounter=0;
	
	FILE* pFile_outTP = NULL;
	FILE* pFile_outFP = NULL;
	FILE* pFile_outNK = NULL;
	FILE* pFile_outMultiClass = NULL;
	if (check_class)
	{
		pFile_outTP = fopen((path_out+s_file_name_scoreclass + "_TP.txt").c_str(), "wt");
		pFile_outFP = fopen((path_out+s_file_name_scoreclass + "_FP.txt").c_str(), "wt");
		pFile_outNK = fopen((path_out+s_file_name_scoreclass + "_NK.txt").c_str(), "wt");
	}
	
	pFile_outMultiClass = fopen((path_out+s_file_name_scoreclass + "_MultiClass.txt").c_str(), "wt");
	M_File_class["notknow"]=fopen((path_out+s_file_name_scoreclass + "_notknow.txt").c_str(), "wt");

	int id_f;
	string name_seq;
	string line_out;
	list<CGroupMatch>::iterator iter2;
	L_CGrM l_grups;
	int i_Match, num, num_size;

	for (map<string, L_CGrM>::iterator iter = M_Group.begin(); iter !=  M_Group.end(); iter++)
	{		
		//printf(".");
		string ind_name_seq = iter->first;
		name_seq = (iter->second.begin())->SequenceName;
		//printf("\nind_name_seq %s", ind_name_seq.c_str());

		
		M_Group[ind_name_seq].sort(compare_match);
		l_grups = iter->second; 
		vector<string> v_grups_MATCH;
		v_grups_MATCH.reserve(15);

		iter2 = l_grups.begin();
		v_grups_MATCH.push_back(iter2->GroupName);
		i_Match = i_procMatch * iter2->Match;
		iter2++;
		
		
		for (; iter2 !=  l_grups.end(); iter2++)
		{
			if (iter2->Match < i_Match)
				break;

			v_grups_MATCH.push_back(iter2->GroupName);
		}

		num_size = v_grups_MATCH.size();	


		if (M_count_dup.find(num_size) == M_count_dup.end())			
			M_count_dup[num_size] = 1;// new
		else
			M_count_dup[num_size]++;		
		
		
		num = 0;
		bool b_NOmulti=0;
		//printf("nazwa %s\t", (iter->first).c_str());
		//printf("num size %d\n", num_size);
		
		
		if ( num_size > 1)
		{
			if (classALL == 0)
			{
				srand (time(NULL));
				num = rand() % num_size;
				num_size = num+1;
			}
			else if (classALL == -1) 
			{
				b_NOmulti = 1;

			}							
		}

		icounter++;
		bool b_NotKnow = 1;

		if (!(M_Group[iter->first].empty()))
		{		

							
			if ( (M_Group[iter->first].begin())->Match > 0 )
			{	
				b_NotKnow = 0;
				
				int id_f_nam =  name_seq.find("name|");	
				//num = num_size+1; // check

				
				if (num_size > 1)
				{
					int icounter=0;
					for (iter2 = l_grups.begin(); iter2 !=  l_grups.end(); iter2++)
					{
						icounter++;
						fPrintScoreDUP(pFile_outMultiClass, *iter2);
						if (icounter>=num_size)
							break;
					}
					fprintf(pFile_outMultiClass, "\n") ; 	
				}				
				
				if (b_NOmulti)
				{
					v_num_score[2]++; //NotKnow because it was classified to few group
					v_num_score[3]++; //
					continue;
				}
						
				for (; num < num_size; num++)
				{
					s_name_group = v_grups_MATCH[num];
					if (M_File_class.find(s_name_group) == M_File_class.end())
						M_File_class[s_name_group] = fopen((path_out+s_file_name_scoreclass + "_" + s_name_group + ".txt").c_str(), "wt");

						
					if (M_Reads.find(ind_name_seq) == M_Reads.end())
					{
						printf("A%sA\n", ind_name_seq.c_str());
					}
					else
					{
						fprintf(M_File_class[s_name_group], "%s\n%s\n", ind_name_seq.c_str(),  M_Reads[ind_name_seq].sequence );
						//fprintf(M_File_class[s_name_group], "%s\n%s\n", name_seq.c_str(),  M_Reads[ind_name_seq].sequence );
					}
				
					if (!check_class)
						continue;

					if (id_f_nam > 0)
					{
						unsigned int IDstart = name_seq.find("<name|");
						IDstart += 6;
						unsigned int IDend = name_seq.rfind("|");
						string IDTAX=name_seq.substr(IDstart, (IDend-IDstart));
						//printf ("%s\n", IDTAX.c_str());
						//id_f = name_seq.find(s_name_group);
						id_f = s_name_group.find(IDTAX);

						if (id_f >= 0)
						{
							fPrintScore(pFile_outTP, 	M_Group, ind_name_seq, name_seq);
							fprintf(pFile_outTP, "%s%s%s\n", "<class|", s_name_group.c_str() , "|") ; 
							v_num_score[0]++; //TP
						}
						else
						{
							fPrintScore(pFile_outFP, 	M_Group, ind_name_seq, name_seq);
							fprintf(pFile_outFP, "%s%s%s\n", "<class|", s_name_group.c_str() , "|") ; 
							v_num_score[1]++; //FP		

						}
					}	
				}
				if (id_f_nam <= 0 && check_class) // if check_class and notknow name
				{
					fPrintScore(pFile_outNK, 	M_Group, ind_name_seq, name_seq);
					fprintf(pFile_outNK, "%s%s%s\n", "<class|", s_name_group.c_str() , "|") ; 
					v_num_score[2]++; //NotKnow
				}				
			}
		}				
		if (b_NotKnow)
		{
			if (check_class)
			{
				fPrintScore(pFile_outNK, 	M_Group, ind_name_seq, name_seq);
				fprintf(pFile_outNK, "%s\n", "<class|0|");
			}
			fprintf(M_File_class["notknow"], "%s\n%s\n", M_Reads[string(ind_name_seq)].name,  M_Reads[string(ind_name_seq)].sequence );
			v_num_score[2]++; //NotKnow
		}
	}

	cout << "END\n";


	funGetTimeSEK(sekTOT, T_start);
	
	if (check_class)
	{
		fclose(pFile_outTP);
		fclose(pFile_outFP);
		fclose(pFile_outNK);
		
		printf("Save summation to file\n");
		FILE* pFile_out;
		string s_file_name_outScore = (s_path_sum + "/Summarize_Class.txt");
		pFile_out = fopen(s_file_name_outScore.c_str(), "at");
		fprintf(pFile_out, "%s", s_file_name_scoreclass.c_str());
		for (int i = 0; i < 4; i++)
			fprintf(pFile_out, "\t%d", v_num_score[i]); // TP FP 0 0 NotKnow

		fprintf(pFile_out, "\t%s", kmer.c_str());
		//fprintf(pFile_out, "\t%s", meta.c_str());
		fprintf(pFile_out, "\t%.3f", sekTOT);
		fprintf(pFile_out, "\tmc=%d", MC);
		fprintf(pFile_out, "\tsim=%d", int(i_procMatch*100));
		fprintf(pFile_out, "\tclALL=%d", classALL);
		fprintf(pFile_out, "\tMM=%d", mismatch);	
		fprintf(pFile_out, "\n");
		fclose(pFile_out);
	}
	//return 0;
	for (map<string, FILE*>::iterator itf =  M_File_class.begin(); itf != M_File_class.end(); itf++)
		fclose(itf->second);
	
	string s_file_name_outDUP = (s_path_sum + "/Summarize_Class_dup.txt");
	FILE* pFile_outDUP = fopen(s_file_name_outDUP.c_str(), "at");

	fprintf(pFile_outDUP, "%s\n", s_file_name_scoreclass.c_str());
	for (map<int, int>::iterator iter = M_count_dup.begin(); iter !=  M_count_dup.end(); iter++)
	{
		fprintf(pFile_outDUP, "%d\t%d\n", iter->first, iter->second );
	}
	fclose(pFile_outDUP);
	
	printf("END\n");
}