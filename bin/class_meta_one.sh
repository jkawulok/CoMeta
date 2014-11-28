#!/bin/env bash 
# classification 

#	This file is a part of CoMeta software.
#	The homepage of the CoMeta project is http://sun.aei.polsl.pl/cometa

echo -e "\n\n####################### Classify reads part1 ###########################" >&2


ERROR_PROC=1;
tempout=TEMPOUT_"$PBS_ARRAYID"_"$RANDOM"_classmetaone_



pathOut=$DIR_MAIN_SCORE_CLASS/tax_"$TAX_SING_LEV"  # where save scores after classification


if [ ! -d "$pathOut" ]; then
	mkdir "$pathOut" ;
fi 


suffixOUT="$PREF_META_SET"_"$MATCH_CUT_OFF"; 
scoreNO_TR=$DIR_SUMMARY/File_Tr_notknow.txt ;  # if there is error, becaouse TR doesn't exist

echo "NR_THREADS=$NR_THREADS size of memeber=$NR_MEM	   MC=$MATCH_CUT_OFF  $FILE_META_SET  classification to: $TAX_SING_LEV "  #

pathKmerpre=$DIR_KMER_DATABASE/"$TAX_SING_LEV"/  ; #
pathKmer="$pathKmerpre"/kmer_"$LEN_KMER"/;
nameMetaSING="$FILE_META_SET"



# +++++++++++++++++++++++++
limit=$NR_THREADS
if [ $limit -gt "$NR_THREADS_SERVER" ]  ||  [ $limit -eq 0 ]
then
	limit="$NR_THREADS_SERVER"
fi

echo Num thread: $limit

# NRtempTHREADS=$limit
NRtempTHREADS="$NR_THREADCOMETA"
NRjobs=$[limit/NRtempTHREADS]
if [ $NRjobs -eq 0 ]
then
	NRjobs=1;
fi


function job_wait {
 
	while [ $(jobs | wc -l) -ge $NRjobs ] ; do
		sleep 0.1;
		jobs &> /dev/null
	done
}
# +++++++++++++++++++++++++++

function process {

	FILE_TR=$1;

	if [[ -e  "$FILE_TR" ]] 
	then
		echo -e "\nFile Tr:  $FILE_TR ";
		file_TRpref=${FILE_TR%\"*};					
		file_TRpref=${file_TRpref#\"*};				
		file_TRpref=${file_TRpref##*/};				
		file_TRpref=${file_TRpref%_k*.*};
		file_TRpref=${file_TRpref#Tr_out_NT_};

		# echo file_TRpref  $file_TRpref;
		pathOutEND="$pathOut"/OTU_ ;

		if [ ! -d "$pathOutEND" ] 
		then
			mkdir "$pathOutEND" 
		fi 
	
		nameMATCH="$file_TRpref"_k"$LEN_KMER"_"$suffixOUT"_M.out
		nameMISMATCH="$file_TRpref"_k"$LEN_KMER"_"$suffixOUT"_MM.out


		if [[ -e  "$pathOutEND"/"$nameMATCH" &&  -e  "$pathOutEND"/"$nameMISMATCH" ]] 
		then
			echo "FOR $file_TRpref was classified ("$nameMATCH" exist)" ;
		else
			# WK  is empty, because  database name include path 
			echo $DIR_MAIN_BIN/cometa 	-goC -t"$NRtempTHREADS" -k"$LEN_KMER" -sK1 -mc$MATCH_CUT_OFF -mr$NR_MEM	-S"$nameMetaSING"   -NS"$file_TRpref"   -D"$FILE_TR"  -WS"$PATH_META_SET"	-WK  -WO"$pathOutEND"   -OSscore_class_mc"$MATCH_CUT_OFF"_k"$LEN_KMER"_"$PBS_ARRAYID".txt  -OC"$nameMATCH"  -ON"$nameMISMATCH" -stepk"$LEN_STEP"  
			
			$DIR_MAIN_BIN/cometa 	-goC -t"$NRtempTHREADS" -k"$LEN_KMER" -sK1 -mc$MATCH_CUT_OFF -mr$NR_MEM	-S"$nameMetaSING"   \
			-NS"$file_TRpref"   -D"$FILE_TR"  -WS"$PATH_META_SET"	-WK  -WO"$pathOutEND"   \
			-OSscore_class_mc"$MATCH_CUT_OFF"_k"$LEN_KMER"_"$PBS_ARRAYID".txt  -OC"$nameMATCH"  -ON"$nameMISMATCH" \
			-stepk"$LEN_STEP"  

		fi
		
	else
		echo "$FILE_TR" >> $scoreNO_TR ;
	fi # -e  "$FILE_TR"

}

# +++++++++++++++++++++++++++

it=0;
# for group in $(cat $FILE_LIST_GROUP)
for group in $LIST_GROUP
do
	# echo $group 
	
	for PAR in "$pathKmer"/Tr_*"$group"*"$numK".res 
	do 
		it=$[it+1]
		process "$PAR" > $DIR_MAIN_BIN/"$tempout""$it".o  2> $DIR_MAIN_BIN/"$tempout""$it".e &
		job_wait
		
		
	done # PAR
done # group

wait


for (( i=1; i<=it; i++ ))
do
	echo -e "\n\niiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii $i iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii" 
	# echo -e "\n\niiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii $i iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii" >&2

	cat $DIR_MAIN_BIN/"$tempout""$i".e >&2
	cat $DIR_MAIN_BIN/"$tempout""$i".o
	rm $DIR_MAIN_BIN/"$tempout""$i".e
	rm $DIR_MAIN_BIN/"$tempout""$i".o
done

ERROR_PROC=0;

