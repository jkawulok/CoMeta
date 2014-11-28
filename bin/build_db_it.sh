#!/bin/bash
# Build kmer database using scores after previous classification

#	This file is a part of CoMeta software.
#	The homepage of the CoMeta project is http://sun.aei.polsl.pl/cometa

echo -e "\n\n####################### Build kmer database ###########################" >&2

ERROR_PROC=1;
tempout=TEMPOUT_"$PBS_ARRAYID"_"$TAX_SING_LEV"_"$RANDOM"_build_



path_dataIN=$PATH_SEQ_REF/"$TAX_SING_LEV"/; # folder with reference sequences  np. 2class
pathSCOREafterClassify=$DIR_MAIN_SCORE_CLASS/tax_"$TAX_PREV_LEV"/  # scores after previous classification

pathKmer=$DIR_KMER_DATABASE/"$TAX_SING_LEV"/  ; #

if [ ! -d "$pathKmer" ]; then
	mkdir "$pathKmer";
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

	nameBin=bin_"$LEN_KMER"_"$TAX_SING_LEV"_$RANDOM
	temp_bin=$DIR_TEMP_BIN/"$nameBin" 
	if [ ! -d $temp_bin ]; then
		mkdir $temp_bin;
	fi 

	y=$1;

	fileTr_Pre=${y%\"*};	#
	#echo "1 $fileTr_Pre";
	fileTr_Pre=${fileTr_Pre#\"*};	
	#echo "2 $fileTr_Pre";
	fileTr_Pre=${fileTr_Pre##*/};	
	echo "$fileTr_Pre";
	fileTr_Pre=${fileTr_Pre#$PREF_META_SET*};
	#echo FILE: $fileTr_Pre		
	#fileTr_Pre=${fileTr_Pre##*\_};	# 
	fileTr_Pre=${fileTr_Pre%.*}; # 
	#echo "fileTr_Pre:  $fileTr_Pre "
	fileTr_NAME=${fileTr_Pre%_idTAX*};	# 
	fileTr_NAME=${fileTr_NAME##*\_};	# idNAME
	
	fileTr_Pre=${fileTr_Pre##*\_};	# idTAX
	echo "fileTr_Pre:  $fileTr_Pre  fileTr_NAME  $fileTr_NAME"
			  
	if [[ $fileTr_Pre = *"notknow"* || $fileTr_Pre = *"FP"*  || $fileTr_Pre = *"TP"*  \
			|| $fileTr_Pre = *"NK"* || $fileTr_Pre = *"MultiClass"* ]]		  
	then 
		continue;
	else

		echo ""
		echo "++++++++++++++++++++++++"
		fileTr_prefix="$fileTr_NAME" ;
		b_NIEist=1;
		echo "$TAX_PREV_LEV": $fileTr_prefix  
		
		subfix1="."
		subfix2=""
		if [[ $fileTr_NAME = *"undef"* ]]
		then
		
			idprev=${fileTr_Pre##idTAX}
			subfix="$idprev"
			echo UNDEF: $subfix
			subfix1=idTAX"$subfix".
			subfix2=idPREV"$subfix"_
			
		fi
		
		
		for sub2 in $subfix1 $subfix2
		do
			echo ". $sub2 ."
			for par in "$path_dataIN"/*"$fileTr_prefix"*"$sub2"*"fa"
			do

				if [ ! -e "$par" ]; then
					echo "FILE for $fileTr_prefix ($sub2) do not exist in $path_dataIN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "

				else		
				
					file=${par%\"*};	
					file=${file#\"*};	#
					file=${file##*/};	# 

					echo " "
					echo "#### file: $file  ####";  # file includes sequences to build Kmer
					
					nameSeq=${file/$'\n'/};
					nameSeq=${nameSeq/$'\r'/};


					nameSeq=${nameSeq%.*}; 	#
					
					NAME_OUT=out_NT_NW_"$nameSeq"_k"$LEN_KMER".res;  #

					
					#echo  "NAME_OUT: $NAME_OUT"	
					if [ -e $pathKmerALL/"$NAME_OUT" ] 
					then
						echo "$NAME_OUT exists ";
					else
						$DIR_MAIN_BIN/cometa    -goB -t$NRtempTHREADS -mr$NR_MEM -k$LEN_KMER -S"$file"  -D"$NAME_OUT"  \
						 -WD"$temp_bin" -WS"$path_dataIN"  -WK"$pathKmerALL"   -WO"$DIR_SUMMARY"  -OS"$PBS_ARRAYID"_score_build.txt 
						
					fi

					echo "###"
					FILE_TR="$pathKmerALL"/Tr_"$NAME_OUT"
					if [ -e "$FILE_TR" ] 
					then
						echo "Tr_$NAME_OUT exist";
					else
						echo 	$DIR_MAIN_BIN/tsk -mr$NR_MEM -WK"$pathKmerALL"  -D"$NAME_OUT"  -t$NRtempTHREADS;
						$DIR_MAIN_BIN/tsk -mr$NR_MEM -WK"$pathKmerALL"  -D"$NAME_OUT"  -t$NRtempTHREADS;
					fi
				fi
			done  # par 
			
		done # sub2
	fi
	
	rm -r "$temp_bin"
}


echo "*******************************************************************************************************************"
it=0
for folderTAX_SING_LEV in "$pathSCOREafterClassify"/OTU_*  # 
do
	
	echo -e "\n\nffffffffffffffffffffffffffffffffffffffff"
	echo folder $folderTAX_SING_LEV;
	echo -e "ffffffffffffffffffffffffffffffffffffffff"
	folderTAX_SING_LEVEL2=$folderTAX_SING_LEV/NEW/
	if [[ ! -d $folderTAX_SING_LEVEL2 ]]; then
		folderTAX_SING_LEVEL2=$folderTAX_SING_LEV/;
	fi 

	if [ $TAKE_MM -eq 1 ]; then
		folderTAX_SING_LEVEL2="$folderTAX_SING_LEVEL2"/MisMatch_all"$CLASS_ALL"_sim"$PROC_SIM"
	else
		folderTAX_SING_LEVEL2="$folderTAX_SING_LEVEL2"/Match_all"$CLASS_ALL"_sim"$PROC_SIM"
	fi

	
	if [[ ! -d $folderTAX_SING_LEVEL2 ]]; then
		folderTAX_SING_LEVEL2=$(echo $folderTAX_SING_LEVEL2 | perl -nle '$_=~ s/all/rand/; print $_;')
	fi



	for ITER in $folderTAX_SING_LEVEL2/"$PREF_META_SET"*_k"$LEN_KMER"_*
	do
		# echo ITER $ITER
		if test ! -s "$ITER"
		then
			continue ;
		fi

		it=$[it+1]
		process "$ITER" > $DIR_MAIN_BIN/"$tempout""$it".o  2> $DIR_MAIN_BIN/"$tempout""$it".e &
		job_wait
	done  # ITER for iterating through results
done # for iterating through  TAX_SING_LEV_* folders

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

