CoMeta
======

CoMeta (<b>C</b>lassification <b>o</b>f <b>meta</b>genomes) is a tool used to assign a query read (a DNA fragment) from metagenomic sample into one of the groups previously prepared by the user. Typically, this is one of the taxonomic rank (e.g., phylum, genus), however prepared groups may contain sequences having various functions.

Licence: CoMeta software distributed under GNU GPL 2 licence.

This repository contains the current version of the program.
The initial version of the Cometa is described in paper: <u><i>Jolanta Kawulok, Sebastian Deorowicz: CoMeta: Classification of metagenomes using k-mers. submitted to journal </i> </u>
and is available at http://sun.aei.polsl.pl/cometa/.

Contact: 
<a href="mailto:jolanta.kawulok@polsl.pl">Jolanta Kawulok</a>



Contents
------------
<ol type="A">
<li> Preparation of the working environment </li>
<li> Description of the use of files that automate the databases building and classification </li>
<li> Description of the various stages of building k-mer databases and read classification </li>
<li> Selection of CoMeta parameters </li>
</ol>




A.	Preparation of the working environment:
------------
Installation, download files and programs, adding the taxonomic id (tax number) to the reference sequence

	I ***** Necessary files and programs *****
	1.	Programs:
		1.1 BLAST+ - Useful for reference sequences extraction from NCBI database
		1.2 bioperl - For dividing reference sequence to the taxonomic groups

	2.	Files:
		2.1	Set of reference sequences:
			Prepare set of reference sequences which include number gi (in FASTA
			format). In order classify to taxon, data includes nucleotide sequences could be downloaded from
			ftp://ftp.ncbi.nlm.nih.gov/blast/db/ website (nt.00.tar.gz, nt.01.tar.gz, nt.02.tar.gz,...) and extracted.
		
			For extracting sequences (from NCBI), unzip the nt.xx.tar.gz file and use command "blastdbcmd" from blast package from NCBI. E.g.: 
				$ ./blastdbcmd -entry all -db nt.xx -out sequences_ntxx.fa
			where xx is number of nt file.
		
		2.2	Taxonomic data:
			Download and unzip the two taxonomic data from NCBI website: ftp://ftp.ncbi.nih.gov/pub/taxonomy
				i.	file: taxdump.tar.gz, which include:
					- names.dmp 	– Taxonomy names
					- nodes.dmp 	– Taxonomy nodes (hierarchy)
						$ wget -c  -P ./NCBI_tree_tax  ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
						$ cd  ./NCBI_tree_tax/
						$ gunzip taxdump.tar.gz
						$ tar -xvf taxdump.tar
					
				ii.	file: gi_taxid_nucl.dmp.gz, which include
					- gi_taxid_nucl.dmp
						$ wget -c -P ./NCBI_tree_tax  ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz
						$ gunzip -c ./NCBI_tree_tax/gi_taxid_nucl.dmp.gz > ./NCBI_tree_tax/gi_taxid_nucl.dmp
				
		2.3	Module for bioperl:
			Bio::LITE::Taxonomy::NCBI - For dividing reference sequence to the taxonomic groups
			
		2.4 Boost version 1.51 or higher (for Boost/filesystem and Boost/thread libraries) - for instalation CoMeta program
			change BOOST_LIB and BOOST_H in makefile to the directories where Boost is installed
			
			
	II ***** CoMeta program *****
	1.	Directory structure:
			bin		- main directory of CoMeta (programs after compilation will be stored here), also it includes two perl scipts
					(gen_list_names.pl and class_seq_to_taxon_all.pl)
			src		- source codes
			example - folder with sample data
		
		
	2.	Binaries:
		After compilation you will obtain six binaries:
			- tsk
			- seq_gi2tax
			- cometa
			- class2best
			- genlist
			- num_class
			
			
	III ***** Preparing reference sequences for taxonomic classification *****
	In order to build the k-mer database for taxonomic classification, the taxonomic id (tax number) have to be added to the 
	single-line description of each reference sequence based on gi number. 
	In this purpose, use the program seq_gi2tax with file gi_taxid_nucl.dmp. Due to the huge file size, we suggest to divide it 
	into smaller parts. Before the first start of the program, use the argument -div1 to split the file:
		$ ./seq_gi2tax -filGT<name> -pGT<path> 	-div1
	And then attributing tax number can be started:
		$ ./seq_gi2tax -filGT<name> -pGT<path> -in<path_name> -out<path_name> 		
	where:
		-in<path_name> - path and file name for input file with sequence which include GI number - from I.2.1 point(e.g., ./NT/sequences_nt00.fa) 
		-out<path_name> - path and file name for output file with sequence, where tax number is added (e.g., ./NT/sequences_nt00_TAX.fa)
		-filGT<name> - file name with relation between gi number and tax number (default: gi_taxid_nucl.dmp)
		-pGT<path> - path where is file with relation between gi number and tax number (path for gi_taxid_nucl.dmp file)
		-div<0/1> - dividing gi_taxid_nucl.dmp file (default 0).



B.	Description of the use of files that automate the databases building and classification 
------------
(a simple example of usage is in "Example_auto_cometa.txt" file)

	There are two scripts to automate the building of databases and reads classification:
	* CoMeta_bud_class_first - The script compares ready with each group. This is a classification used to start the taxonomic classification 
								or classification groups created by the user.
	* CoMeta_bud_class_sing - The script for reads classification, which takes into account where reads have been classified to the higher level 
								(for the taxonomic classification).
								
	I ***** CoMeta_bud_class_first *****
		$ CoMeta_bud_class_first -maindir <DIR_MAIN_BIN> -S <FILE_META_SET>  -WS <PATH_META_SET>  -WD <DIR_TEMP_BIN>  -tall <NR_THREADS>  -t <NR_THREADCOMETA> \
		-mr <NR_MEM>  -proc <PROC_SIM>  -k <LEN_KMER>  -stepk <LEN_STEP>  -mc <MATCH_CUT_OFF>  -mm <0/1>  -cl <-1/0/1>  -suffDESCRP <SUFFIX_DESCRP> \
		-listgr <LIST_GROUP>  -OTU <TAX_SING_LEV>  -dirin <DIR_MAIN_SCORE_CLASS>  -WK <DIR_KMER_DATABASE>  -diroutref <PATH_SEQ_REF>  -OTUPREV <TAX_PREV_LEV> \
		-dirncbitax <DIR_TAX_NCBI>  -dirinref <PATH_IN_START_SEQ_REF>  -divseq <0/1>  
	
	Description of the script parameters:
		-maindir <DIR_MAIN_BIN> -the path where are all scripts 
		-S <FILE_META_SET> - file name of metagenomic set
		-WS <PATH_META_SET> - path where is FILE_META_SET file
		-WD <DIR_TEMP_BIN> - working directory for temporary files
		-t <NR_THREADCOMETA> - total number of computional threads for single k-mer database and classification  (default: 4)
		-tall <NR_THREADS> - total number of computation threads (default: equal to no. of system cores). At the same time, "NUMJOBS" databases 
							are built and reads are compared with "NUMJOBS" groups, where NUMJOBS=NR_THREADS/NR_THREADCOMETA. Therefore, the multiple 
							of "NR_THREADCOMETA" is recommended. 
		-mr <NR_MEM> - max amount of RAM in GB
		-proc <PROC_SIM> - similarity the best results, which are taken into account; default 100[%];
		-k <LEN_KMER> - k-mer length (max 32); default: 24
		-stepk <LEN_STEP> - k' - length of offset sliding window (default: length of k-mer)
		-mc <MATCH_CUT_OFF> - the percent identity to classify a match (default: 5);
		-mm <0/1> - taking into account the mismatch files; 0 - NO; 1 - YES
		-cl<-1/0/1> - when reads are classified to a few groups, then reads are assignment to: -1 - any group; 0 - random group; 1 -   all of these groups;
						default: -1 \n"
		-suffDESCRP <SUFFIX_DESCRP> - additional description of results (suffix)
		-listgr <LIST_GROUP> - the list of substrings, which must be included in the name of the group to which the query reads are compared. For example, 
								in the folder of reference sequences/k-mer databases, there are data for bacterium, viruses, and eukaryotes, and we only 
								want to classify to bacteria and virus. Then, this command would be <-listgr "Bacteria Viruses"> (assuming that these names
								appear in the file names).
		-OTU <TAX_SING_LEV> - the generic name of the group to which reads are classified. For example, for the taxonomic classification: "phylum", 
								our: "mikedb".
		-dirin <DIR_MAIN_SCORE_CLASS> - the path, where "TAX_SING_LEV" folder is created with the results
		-WK <DIR_KMER_DATABASE> - the path, where "TAX_SING_LEV" folder is created with the k-mer databases
		-diroutref <PATH_SEQ_REF> - the path, where "DIR_KMER_DATABASE/TAX_SING_LEV" folder contains reference sequences. For taxonomic classification, 
									the script creates "TAX_SING_LEV" folder with reference sequences (see the following parameters).		
	
	Parameters only for taxonomic classification:
		-OTUPREV <TAX_PREV_LEV> - the name of a higher level to which reads are classified (eg, if TAX_SING_LEV = "phylum", "TAX_PREV_LEV" = "superkingom".
		-dirncbitax <DIR_TAX_NCBI> - the path to the files of names.dmp and nodes.dmp	
		-dirinref <PATH_IN_START_SEQ_REF> - the path to folder with reference sequences that have not yet been divided by "TAX_SING_LEV" but containing 
											the taxonomy id
		-divseq <0/1> - if reference sequences shall be divided by TAX_SING_LEV. 1 - yes, for the taxonomic classification, 0 - no, for their own 
						classification, or if it was earlier done
		
		
	II ***** CoMeta_bud_class_sing *****
	Reads classification is based on a higher level to which reads were classified.	This script is executed for the target and each intermediate 
	classification levels.
	
		$ CoMeta_bud_class_sing -maindir <DIR_MAIN_BIN> -S <FILE_META_SET>  -WS <PATH_META_SET>  -WD <DIR_TEMP_BIN>  -tall <NR_THREADS>  -t <NR_THREADCOMETA> \
		-mr <NR_MEM>  -proc <PROC_SIM>  -k <LEN_KMER>  -stepk <LEN_STEP>  -mc <MATCH_CUT_OFF>  -mm <0/1>  -cl <-1/0/1>  -suffDESCRP <SUFFIX_DESCRP> \
		-OTU <TAX_SING_LEV>  -dirin <DIR_MAIN_SCORE_CLASS>  -WK <DIR_KMER_DATABASE>  -diroutref <PATH_SEQ_REF>  -OTUPREV <TAX_PREV_LEV> \
		-dirncbitax <DIR_TAX_NCBI> 
	
	Description of the parameters is the same as for the script: CoMeta_bud_class_first.
	Compared to the script CoMeta_bud_class_first, CoMeta_bud_class_sing does not contain the following parameters: -listgr <LIST_GROUP>, 
	-dirinref <PATH_IN_START_SEQ_REF>, and -divseq <0/1>

	The lower level are, the less memory is needed, however it is recommended to use a larger number of threads (-tall <NR_THREADS>). Because greater number
	of groups is used.

	
*******************************************************************************************************************************************************
C.	Description of the various stages of building k-mer databases and read classification (
------------
a simple example of usage is in "Example_cometa.txt" file)


	I ***** Data preparation for the building of k-mer databases (used only in the taxonomic classification) ***** 
		1. Make sure that the reference sequences contain a tax numbers (A III subsection).
			
		2. Divide reference sequences into groups on taxomic rank.
			Use the program class_seq_to_taxon_all.pl in perl.
				$ perl class_seq_to_taxon_all.pl –i <fileseq> –os <suffix> -op <prefix> –Pin <pathIn> –Pout <pathOut> –LS <OTU> -LF <OTUH> -Pncbi <pathNCBItax>
			where:		
				-i <fileseq> – file include reference sequences (e.g. sequences_nt00_TAX.fa, output from II.1 point)
				-os <suffix> - suffix after divided – it is not necessary if you have only one fileseq. E.g. _NT_00 for sequences_nt00_TAX.fa
				-op <prefix> - prefix after divided
				-Pin <pathIN> – path for input data (<fileseq>)
				-Pout <pathOUT> – path for output data 
				-LS <OTU> - taxonomic level of sampling selected for classification e.g., 'phylum'
				-LF <OTUH> - higher taxonomic level than OTU e.g.,'superkingdom'
				-Pncbi <pathNCBItax> - path to the files of names.dmp and nodes.dmp		

	II ***** Build k-mers database *****
		Single file "group file" includes all reference sequences belonging to the group (eg. 'Proteobacteria') to which the query read
		is compared. The user can create own group, which contains a set of sequences of specific attribution.
		For taxonomic classification, beginning the classification from the phylum is recommended, however, it is possible to start at the 
		lower rank. 
		
		1.	Generate a file containg file names (for taxonomic classifcation):
			The sequences from the same kind of rank can be in several files, e.g.,proteobacteria sequnces are in nt.00, nt.01, ... 
			files, therefore before building k-mers database, file with list of file names have to be created. This can be made as follows:
				$ perl gen_list_names.pl -nf <filenames> –pf <pathIn> –ps <pathOut>
			where:
				-nf <filenames> – file that include list of all file names after dived sequnces into OTUs (file have to be in <pathIN> folder)
				-pf <pathIN> – path for input data (filenames)
				-ps <pathOUT> – path for output data
				
		2.	Build kmer databases – main program
			To build the database use program cometa with parameter <B> for '-go' parameter
				$ cometa -goB -mr<num_mr> -t<num_cor> -NS<name_sub> -k<numK> -L<nameLIST>* -D<nameDB> -WD<pathBIN> -WS<pathSEQ> -WK<pathKMER> -WO<pathOUT> \
						-OS<nameSum>
			where:
				-go<B/C> - use parameter “B” for building kmer database
				-mr<num_mr> - max amount of RAM in GB
				-t<num_cor> - total number of computation threads (default: equal to no. of system cores)
				-NS<name_sub> - name of file (subfix) for out of kmers, score;
				-k<numK> - k-mer length (max 32); default: 24
				-L<nameLIST> - if the sequences belonging to the same group are in several files, 
								a file containing a list of the names of these files
				-S<nameFILEin> - file name of input file with sequences ("group file")
				-D<nameDB> - name of database file (output)
				-WK<pathKMER> -path for k-mer database (nameDB)
				-WD<pathBIN> - working directory for temporary files
				-WS<pathSEQ> - path for inputdata with sorted sequences in OTU
				-WO<pathOUT> - path for save summary score 
				-OS<nameSum> - file name for save summary score (default: score_build.txt)

				* You can use single file with input sequences (parameter: -S) or list with names of input file with sequences (parameter –L)
				
		3.	Decompressed build k-mer databases
			In order to faster loaded databases into the program during the classification, the file with k-mer database could be decompressed.
				$ tsk -mr<num_mr> -t<num_cor> -WK<pathKMER> -D<nameDB>
			where:
				-mr<num_mr>-max amount of RAM in GB
				-t<num_cor> - total number of computation threads (default: equal to no. of system cores
				-WK<pathKMER> -path for k-mer database
				-D<nameDB> - name of database file

			

	III ***** Read classification *****
		1.	Comparisons:
			To compare the reads with every reference groups use program cometa with parameter <C> for '-go' parameter
			
				$ cometa -goC -mr<num_mr> -t<num_cor> -NS<name_sub> -mc<num_MC> -S<file_meta> -D<nameDB> -sK<0/1> -stepk <LEN_STEP>  -WS<pathSEQ> \
				-WK<pathKMER> -WO<pathOUT> 
			where:
				-go<B/C> - use parameter “C” for comparison
				-mr<num_mr>-max amount of RAM in GB;
				-t<num_cor> - total number of computation threads (default: equal to no. of system cores)
				-NS<name_sub> - name of file (subfix) for out of kmers, score;
				-mc<num_MC> - the percent identity to classify a match (default: 5);
				-S<file_meta> - name of input file with metagenomic sequences
				-D<nameDB> - name of database file (input)
				-SK<0/1> - 0 – if kmer database is compressed; 1 – if kmer database is decompressed (output after tsk program)
				-stepk <LEN_STEP> - k' - length of offset sliding window (default: length of k-mer)
				-WK<pathKMER> -path for k-mer database (nameDB)
				-WS<pathSEQ> - path for metagenomic file 
				-WO<pathOUT> - path for save match and mismatch file and summary score
				
		2.	Assignment to the best group
			
				$ class2best -NO<namePrefix> -WI<path_IN> -WO<pathOUT> -NF<inputfiles> -NG<inputgroup>  -NR<pathNAMEreads>  -proc<numP> -cl<-1/0/1> \
				-k<numK> -mc<MC>		
			where:
				-NO<namePrefix> - prefix name e.g., metagenomic sample name
				-WO<pathOUT> – path for output data
				-NF<inputfiles> - file name which includes names of input files (which inlude scores)	
				-NG<inputgroup> - file name which includes groups of input files, where reads were classified		
				-NR<pathNAMEreads> - path and name file of reads which was classified 
				-cl<-1/0/1> - when reads are classified to a few groups, then reads are assignment to: -1 - any group; 0 - random group; 
								1 - all of these groups; ; default: -1 \n"
				-proc<numP> - similarity the best results, which are taken into account; default 100[%];			
				-WI<path_IN> - path where there are files with the results of the comparison step (match, mistamtch)
				-WS<path_sum> - path for summarization file
				-mc<numMC> - the percent identity to classify a match [%] 
				-k<numK> - k-mer length
	




