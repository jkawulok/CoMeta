/*
	This file is a part of CoMeta software.
	The homepage of the CoMeta project is http://sun.aei.polsl.pl/cometa
 
	Authors: Jolanta Kawulok and Sebastian Deorowicz
	
	Version: 0.3
	Date   : 2014-11-26
*/

#include "stdafx.h"
#include <algorithm>
#include "defs.h"
#include "classifier.h"


static pthread_mutex_t CriticalSectionM = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t CriticalSectionMM = PTHREAD_MUTEX_INITIALIZER;

CClassifier_new::CClassifier_new(CMemoryMonitor *_mm, CBinKmers_all_new *_bin_kmers, CReadsQueue *_pread,  uchar _kmer_len, uchar _step_k, double _matchcutoff, FILE* _outMatch, FILE* _outMisMatch)
{	
	mm = _mm;
	bin_kmers = _bin_kmers;
	pread = _pread;
	kmer_len = _kmer_len;
	step_k = _step_k;
	matchcutoff =  _matchcutoff;
	MISmatchcutoff = 1-matchcutoff;

	name = NULL;
	seq = NULL;
	num_alig_reads = 0;
	num_REValig_reads = 0;
	num_NOalig_reads = 0;


	kmer_mask = (uint64(1) << 2*kmer_len)-1;
	second_mask = (uint32(1) << 8)-1;
	ending_mask = (uint64(1) << ((kmer_len-8)*2))-1;

	outMatch		= _outMatch;
	outMisMatch		= _outMisMatch;


};

CClassifier_new::~CClassifier_new()
{

}



void CClassifier_new::Write_score(int s)
{
	// for SIMHC short names
	//int dl_n = strlen(name);
	//for (int u=0; u<dl_n; u++)
	//{
	//	if (name[u] == ' ')
	//	{
	//			name[u]=0;
	//			break;
	//	}
	//}

	//cout << name << "\n";
	
	if (s==1)
	{
		pthread_mutex_lock( &CriticalSectionM );
		fprintf(outMatch, "Min_err:\t%d\tMatch:\t%d\t", min_error, match_nucl);
		fprintf(outMatch, ">%s \tM_ORG\n", name);
		pthread_mutex_unlock( &CriticalSectionM );
	}
	else if (s==2)
	{
		pthread_mutex_lock( &CriticalSectionM );
		fprintf(outMatch, "Min_err:\t%d\tMatch:\t%d\t", min_error, match_nucl);
		fprintf(outMatch, ">%s \tM_COMP\n", name);
		pthread_mutex_unlock( &CriticalSectionM );
	}
	else
	{
		pthread_mutex_lock( &CriticalSectionMM );
		fprintf(outMisMatch, "Min_err:\t%d\tMatch:\t%d\t", min_error, match_nucl);
		fprintf(outMisMatch, ">%s \tMM\n", name);
		pthread_mutex_unlock( &CriticalSectionMM );
	}
	
}

bool CClassifier_new::ProcessClass(char* _name, char* _seq, uint64 _size)
{
	name = _name;
	seq = _seq;
	size = _size;

	match_nucl = 0;
	trackmatch = 0;
	trackposition = 0;
	min_error = 0;

	Check_seqence();
	//if (double(match_nucl)/size >= matchcutoff) //	
	if (frac_eq_grat(double(match_nucl)/size, matchcutoff))
	{
		num_alig_reads++;
		Write_score(1);
		return true;
	}

		
	Reverse_Comp();
	match_nucl = 0;
	trackmatch = 0;
	trackposition = 0;
	min_error = 0;

	Check_seqence();
	//if (double(match_nucl)/size >= matchcutoff) //	
	if (frac_eq_grat(double(match_nucl)/size, matchcutoff))
	{
		num_REValig_reads++;
		Write_score(2);
		return true;
	}

	num_NOalig_reads++;
	Write_score(0);
	return false; 
}


// --------------------------------------------------------------------------------------------------------
void CClassifier_new::Check_seqence()
{
	uint64 kmer_str = 0;
	uint64 kmer_strQP = 0;

	int i;
	uchar prefix;
	uchar second;
	uint64 ending;

	int omit_next_n_kmers = 0;
	uint64 num_finded = 0;
	uint64 num_finded_BIN = 0;
	uint64 num_finded_INTER = 0;
	int kmermatchcountqp = 0;
	int nrqp = 1;


	for(i = 0; i<kmer_len-1; ++i)
	{
		if(seq[i] < 0)
		{
			omit_next_n_kmers = i+1;
			kmer_str = (kmer_str << 2) + 0;
		}
		else
			kmer_str = (kmer_str << 2) + seq[i];
	}
	
	int k = i;
	if (omit_next_n_kmers==0)
	{
		kmer_strQP = (kmer_str << 2) + seq[k];

		prefix = (uchar) (kmer_strQP >> ( (kmer_len-4)*2) );
		second = (uchar) ((kmer_strQP >> ( (kmer_len-8)*2) )  & second_mask);
		ending = (kmer_strQP & ending_mask);
		if (Check_kmer_BIN(prefix, second, ending) > 0) 		//if (Check_kmer_INTER(prefix, second, ending) > 0)	
			kmermatchcountqp++;
	}



	// #Quick pass of query
	//k++;	
	k = step_k;	//

	if (step_k > 1)
	{
		while (kmermatchcountqp < nrqp && (k+kmer_len) < size)
		{
			int j;
			kmer_strQP = 0;
			for (j=0; j<kmer_len; j++)
			{
				if(seq[j+k] < 0)
					break;
				else
					kmer_strQP = (kmer_strQP << 2) + seq[j+k];
			}

			if(seq[j+k] >= 0)
			{
				prefix = (uchar) (kmer_strQP >> ( (kmer_len-4)*2) );
				second = (uchar) ((kmer_strQP >> ( (kmer_len-8)*2) )  & second_mask);
				ending = (kmer_strQP & ending_mask);
				if (Check_kmer_BIN(prefix, second, ending) > 0)//if (Check_kmer_INTER(prefix, second, ending)>0)
					kmermatchcountqp++;
			}

			k += step_k;
		}


		if (kmermatchcountqp == 0)
			return;
	}
	



	// -------------
	for(; i < size; ++i)
	{
		if(seq[i] < 0)
		{
			omit_next_n_kmers = kmer_len;
			trackmatch = 0;
			kmer_str = ((kmer_str << 2) + 0) & kmer_mask;	
		}
		else
			kmer_str = ((kmer_str << 2) + seq[i]) & kmer_mask;		
			
		
		
		if(omit_next_n_kmers > 0)
		{
			omit_next_n_kmers--;
			num_finded = 0;
		}
		else
		{
			prefix = (uchar) (kmer_str >> ( (kmer_len-4)*2) );
			second = (uchar) ((kmer_str >> ( (kmer_len-8)*2) )  & second_mask);
			ending = (kmer_str & ending_mask);

			num_finded = Check_kmer_BIN(prefix, second, ending);
			//num_finded = Check_kmer_INTER(prefix, second, ending);

		}
		


		if (num_finded>0)
		{	 
			if (trackmatch == 0) 
			{		
				if (trackposition >= kmer_len || match_nucl==0 )
					match_nucl += kmer_len; // Adds to match score 
				else
					match_nucl += trackposition;

				trackmatch = 1;
				trackposition = 0;
			}
			else 
			{
				match_nucl++;
			}
		}    
		else 
		{	

			if ( (double(i-trackposition+1) - match_nucl) / size - MISmatchcutoff > 0.000001)
			{
				return;
			}


			if (trackposition%kmer_len == 0)
			{
				min_error++;
			}

			trackposition++;
			trackmatch = 0;
		}

	}
}

uint64 CClassifier_new::Check_kmer_BIN(uchar prefix, uchar second, uint64 searching)
{
	uint64 num_finded = 0;
	uint64 indS = 0;
	uint64 indE;
	uchar mask_count = 255;

	if (bin_kmers->is_empty(prefix))
		return 0;

	if (second>0)
		indS = bin_kmers->get_val_next_pozSin(prefix, second-1);


	indE = (bin_kmers->get_val_next_pozSin(prefix, second));
	if (indE > 0)
		indE -= 1;
	else
		return 0;


	if (indS > indE)
		return 0;	

	uint64 *buffor = bin_kmers->get_buff();
	
	
	if (indS == indE)
	{
		if ( (buffor[indS]) >> 8 == searching)
			return (buffor[indS] & mask_count );		//if ( (buffor[indS]) == searching)		//	return (1);
		else
			return 0;
	}
	else
	{
		if ( ((buffor[indS]) >> 8 > searching) || ((buffor[indE]) >> 8 < searching) )		//if ( ((buffor[indS]) > searching) || ((buffor[indE]) < searching) )
			return 0;


		uint64 indM;
		while (indS <= indE)
		{
			indM = (indS+indE)/2;
			
			if ( (buffor[indM]) >> 8 == searching)
				return (buffor[indM] & mask_count );		//if ( (buffor[indM]) == searching)		//	return (1);
			else
			{
				if ( (buffor[indM]) >> 8 > searching)			//if ( (buffor[indM]) > searching)
					indE = indM - 1;
				else
					indS = indM +1 ;
			}
		}
	}
	return 0;
}

//---------------------------------------------------------------------------------------------

uint64 CClassifier_new::Check_kmer_INTER(uchar prefix, uchar second, uint64 searching)
{
	uint64 num_finded = 0;
	uint64 indS = 0;
	uint64 indE;
	uchar mc = 255; //mask count

	if (bin_kmers->is_empty(prefix))
		return 0;

	if (second>0)
		indS = bin_kmers->get_val_next_pozSin(prefix, second-1);
	
	indE = (bin_kmers->get_val_next_pozSin(prefix, second));

	if (indS > indE)
		return 0;	

	uint64 *buffor = bin_kmers->get_buff();
	
	if (indS == indE)
	{
		if ( (buffor[indS]) >> 8 == searching)
			return (buffor[indS] & mc );
		else
			return 0;
	}
	else
	{
		uint64 indM;
		while ( ((buffor[indE] >> 8) >= searching) && ((buffor[indS] >> 8) <= searching) )
		{
			indM = indS + (((searching - (buffor[indS] >> 8)) * (indE-indS)) / ((buffor[indE] >> 8) - (buffor[indS] >> 8)));
			
			if ( (buffor[indM]) >> 8 == searching)
				return (buffor[indM] & mc );
			else
			{
				if ( (buffor[indM]) >> 8 > searching)
					indE = indM - 1;
				else
					indS = indM +1 ;
			}
		}
	}
	return 0;
}


void CClassifier_new::Reverse_Comp()
{
	char* seq_new = new char[size+2];
	char nuc;
	char new_nuc;

	for (int i=0; i<size; ++i)
	{
		nuc = seq[i];

		if (nuc == 0) //a
			new_nuc = 3;
		else
			if (nuc == 1) //c
				new_nuc = 2;
			else
				if (nuc == 2) //g
					new_nuc = 1;
				else
					if (nuc == 3) // t
						new_nuc = 0;
					else
						new_nuc = nuc;


		seq_new[size-i-1] = new_nuc;
	}
	seq = seq_new;

}



void CClassifier_new::GetTotal(uint64 &_num_alig_reads, uint64 &_num_REValig_reads, uint64 &_num_NOalig_reads)
{
	_num_alig_reads = num_alig_reads;
	_num_REValig_reads = num_REValig_reads;
	_num_NOalig_reads = num_NOalig_reads;
}