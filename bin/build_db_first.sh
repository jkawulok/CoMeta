#!/bin/bash
# Build kmer database using all reference sequences

#	This file is a part of CoMeta software.
#	The homepage of the CoMeta project is http://sun.aei.polsl.pl/cometa

echo -e "\n\n####################### Build kmer database ###########################" >&2

ERROR_PROC=1;

tempout=TEMPOUT_"$PBS_ARRAYID"_"$RANDOM"_buildFIR_



path_dataIN=$PATH_SEQ_REF/"$TAX_SING_LEV"/; # 
pathKmer=$DIR_KMER_DATABASE/"$TAX_SING_LEV"/  ; #

if [ ! -d $pathKmer ]; then
	mkdir $pathKmer;
fi





pathKmerALL="$pathKmer"/kmer_"$LEN_KMER"/;	
if [ ! -d $pathKmerALL ]; then
	mkdir $pathKmerALL;
fi 





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

	numrand=$RANDOM
	nameBin=bin_"$LEN_KMER"_"$TAX_SING_LEV"_"$numrand"

	temp_bin=$DIR_TEMP_BIN/"$nameBin" 
	if [ ! -d $temp_bin ]; then
	  mkdir $temp_bin;
	fi 


	fileNW=$1;
	echo -e "\n\n ===== file: $fileNW ===="
	
	nameSeq=${fileNW%\"*}	# 
	nameSeq=${nameSeq#\"*}	# 
	nameSeq=${nameSeq##*/};	#
	nameSeq=${nameSeq%.*}; 	# 
				
	NAME_OUT=out_NT_NW_"$nameSeq"_k"$LEN_KMER".res;  #
	if [ -e $pathKmerALL/"$NAME_OUT" ] 
	then
		echo " $NAME_OUT exists ";
	else
		$DIR_MAIN_BIN/cometa    -goB -t$NRtempTHREADS -mr$NR_MEM -k$LEN_KMER -S"$fileNW"  -D"$NAME_OUT" \
		-WD"$temp_bin" 	-WS  -WK$pathKmerALL   -WO$DIR_SUMMARY  -OS"$PBS_ARRAYID"_score_build.txt ;  
	fi

	echo "###"
	FILE_TR="$pathKmerALL"/Tr_"$NAME_OUT"
	if [ -e "$FILE_TR" ] 
	then
		echo " Tr_$NAME_OUT exist";
	else
		echo 	$DIR_MAIN_BIN/tsk -mr$NR_MEM -WK"$pathKmerALL"  -D"$NAME_OUT"  -t$NRtempTHREADS;
		$DIR_MAIN_BIN/tsk -mr$NR_MEM -WK"$pathKmerALL"  -D"$NAME_OUT"  -t$NRtempTHREADS;
	fi


	rm -r "$temp_bin"

}


echo "*******************************************************************************************************************"
it=0
for group in $LIST_GROUP
do
	# echo $group 
	for ITER in "$path_dataIN"/*"$group"*.fa
	do

		if [ ! -e "$ITER" ] 
		then
			continue;	
		fi
	
	
		it=$[it+1]
		process "$ITER" > $DIR_MAIN_BIN/"$tempout""$it".o  2> $DIR_MAIN_BIN/"$tempout""$it".e &
		job_wait

	done # ITER

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

