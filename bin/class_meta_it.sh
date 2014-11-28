#!/bin/env bash 
# Classification using a higher level,

#	This file is a part of CoMeta software.
#	The homepage of the CoMeta project is http://sun.aei.polsl.pl/cometa


echo -e "\n\n####################### Classify reads part1 ###########################" >&2



ERROR_PROC=1;
tempout=TEMPOUT_"$PBS_ARRAYID"_"$TAX_SING_LEV"_"$RANDOM"_classmetait_

pathOut=$DIR_MAIN_SCORE_CLASS/tax_"$TAX_SING_LEV"/  # where save scores after classification

pathSCOREafterClassify=$DIR_MAIN_SCORE_CLASS/tax_"$TAX_PREV_LEV"/  # scores after previous classification


if [ ! -d "$pathOut" ]; then
	mkdir "$pathOut" ;
fi 




suffixOUT="$PREF_META_SET"_"$MATCH_CUT_OFF"; 
scoreNO_TR="$DIR_SUMMARY"/File_Tr_notknow.txt ;  

pathKmerpre=$DIR_KMER_DATABASE/"$TAX_SING_LEV"/  ; #
pathKmer="$pathKmerpre"/kmer_"$LEN_KMER"/;



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

	nameMetaSING=$1
	# echo -e "\n\nwwwwwwwwwwwwwwwwww $1 wwwwwwwwwwwwwwwwwwwww"		
	fileTr_prefix=${nameMetaSING%\"*};	#
	fileTr_prefix=${fileTr_prefix#\"*};	# 
	fileTr_prefix=${fileTr_prefix##*/};	# 
	fileTr_prefix=${fileTr_prefix#*_k"$LEN_KMER"*"$SUFFIX_DESCRP"_}; #
	fileTr_prefix=${fileTr_prefix%.*}; # 


	#echo fileTr_prefix  $fileTr_prefix ;


	nameMetaSING=${nameMetaSING/$'\n'/}; #
	nameMetaSING=${nameMetaSING/$'\r'/};

	if [[ $nameMetaSING = "" ]]
	then
		continue ;
	fi

	nameMetaSING=${nameMetaSING%\"*};	# 
	nameMetaSING=${nameMetaSING#\"*};	# 
	path_Meta=${nameMetaSING%/*}/;   	#

	nameMetaSING=${nameMetaSING##*/};	# 
	file_prefix=${nameMetaSING%.*};	# 
	# echo "3 $file_prefix"





	file_prefix=${file_prefix#*_k"$LEN_KMER"*"$SUFFIX_DESCRP"_}; # 
	echo "";
	echo -e "\n=================== $file_prefix =================";
	echo "META: $nameMetaSING  LEN_KMER: $LEN_KMER" ;
	# file_prefix=${file_prefix#*\_};	# name without TAX_PREV_LEV name
	fileTr_NAME=${file_prefix%_idTAX*};	# 
	fileTr_NAME=${fileTr_NAME##*\_};	# idNAME
	fileTr_ID=${file_prefix##*\_};	# idTAX
	
	echo "== NAME: $fileTr_NAME   ID: $fileTr_ID ==";
	
	subfix=""
	if [[ $fileTr_NAME = *"undef"* ]]
	then
		#suffix="$fileTr_ID"_
		#echo UNDEF: $suffix
		idprev=${fileTr_ID##idTAX}
		subfix="$idprev"
		echo UNDEF: $subfix

	fi


	EXIST_Tr=1 #for FILE_TR in "$pathKmer"/Tr_*NW_"$fileTr_NAME"*"$suffix"*
	for FILE_TR in "$pathKmer"/Tr_*NW_"$fileTr_NAME"*"$subfix"*
	do 
		if [[ -e "$FILE_TR" ]]
		then
		
		echo -e "\nFile Tr:  $FILE_TR ";
		file_TRpref=${FILE_TR%\"*};					
		file_TRpref=${file_TRpref#\"*};				
		file_TRpref=${file_TRpref##*/};				
		file_TRpref=${file_TRpref%_k*.*};
		file_TRpref=${file_TRpref#Tr_out_NT_};

		# echo file_TRpref  $file_TRpref;
			pathOutEND="$pathOut"/OTU_"$file_prefix"/ ;
			if [ ! -d "$pathOutEND" ]; then
				mkdir "$pathOutEND" ;
			fi 


			nameMATCH="$file_TRpref"_k"$LEN_KMER"_"$suffixOUT"_M.out
			nameMISMATCH="$file_TRpref"_k"$LEN_KMER"_"$suffixOUT"_MM.out

			if [[ -e  "$pathOutEND"/"$nameMATCH" &&  -e  "$pathOutEND"/"$nameMISMATCH" ]] 
			then
				echo "FOR $file_TRpref was classified" ;
			else
				# WK bez niczego bo w nazwie bazy jest juz sciezka
				echo $DIR_MAIN_BIN/cometa   -goC -t"$NRtempTHREADS" -k"$LEN_KMER" -sK1 -mc$MATCH_CUT_OFF -mr$NR_MEM	 -S"$nameMetaSING"   -NS"$file_TRpref"   -D"$FILE_TR"  -WS"$path_Meta"	-WK  -WO"$pathOutEND"   -OSscore_class_mc"$MATCH_CUT_OFF"_k"$LEN_KMER"_"$PBS_ARRAYID".txt  -OC"$nameMATCH"  -ON"$nameMISMATCH" 


				$DIR_MAIN_BIN/cometa   -goC -t"$NRtempTHREADS" -k"$LEN_KMER" -sK1 -mc$MATCH_CUT_OFF -mr$NR_MEM	 -S"$nameMetaSING" \
				-NS"$file_TRpref"   -D"$FILE_TR"  -WS"$path_Meta"	-WK  -WO"$pathOutEND"  \
				-OSscore_class_mc"$MATCH_CUT_OFF"_k"$LEN_KMER"_"$PBS_ARRAYID".txt  -OC"$nameMATCH"  -ON"$nameMISMATCH" \
				-stepk"$LEN_STEP"  

			
			fi
		else
			if [[ $fileTr_prefix = *"notknow"* ]] 
			then
				continue;
				# echo "it is notknow or undef";
			else
				echo NOT KNOW: TR "$FILE_TR" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				echo "$FILE_TR" >> "$scoreNO_TR" ;
			fi
		fi

	done

}


# +++++++++++++++++++++++++++
echo "NR_THREADS=$NR_THREADS size of memeber=$NR_MEM	   MC=$MATCH_CUT_OFF  $FILE_META_SET    $TAX_PREV_LEV -> $TAX_SING_LEV "  #

it=0;
for folderOTU in "$pathSCOREafterClassify"/OTU_*  # 	
do

	folderOTU_END=$folderOTU/NEW/
	if [[ ! -d $folderOTU_END ]]; then
		folderOTU_END=$folderOTU/;
	fi 

	if [ $TAKE_MM -eq 1 ]; then
		folderOTU_END="$folderOTU_END"/MisMatch_all"$CLASS_ALL"_sim"$PROC_SIM"
	else
		folderOTU_END="$folderOTU_END"/Match_all"$CLASS_ALL"_sim"$PROC_SIM"
	fi
	
	if [[ ! -d $folderOTU_END ]]; then
		folderOTU_END=$(echo $folderOTU_END | perl -nle '$_=~ s/all/rand/; print $_;')
	fi


	
	# echo folder $folderOTU_END ;
	for nameMetaSINGpre in $folderOTU_END/*"$PREF_META_SET"*_k"$LEN_KMER"*
	do

		if [ ! -s "$nameMetaSINGpre" ]
		then
			continue ;
		fi
				
		if [[ $nameMetaSINGpre = *"notknow"*   || $nameMetaSINGpre = *"FP"*  || $nameMetaSINGpre = *"TP"* \
			|| $nameMetaSINGpre = *"NK"* || $nameMetaSINGpre = *"MultiClass"* ]]
		then
			continue ;
		fi

		it=$[it+1]
		process "$nameMetaSINGpre" > $DIR_MAIN_BIN/"$tempout""$it".o  2> $DIR_MAIN_BIN/"$tempout""$it".e &
		job_wait
		
	done
done


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

