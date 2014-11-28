#!/bin/bash
#!/usr/bin/perl

#	This file is a part of CoMeta software.
#	The homepage of the CoMeta project is http://sun.aei.polsl.pl/cometa

echo -e "\n\n####################### Divide sequence ###########################" >&2

ERROR_PROC=1;
tempout=TEMPOUT_"$PBS_ARRAYID"_"$TAX_SING_LEV"_"$RANDOM"_divseq_


# echo NR_THREADS $NR_THREADS
# echo NR_THREADS_SERVER "$NR_THREADS_SERVER"

limit=$NR_THREADS
if [ $limit -gt "$NR_THREADS_SERVER" ]  ||  [ $limit -eq 0 ]
then
	limit="$NR_THREADS_SERVER"
fi

echo Num thread: $limit



function job_wait {
 
	while [ $(jobs | wc -l) -ge $limit ] ; do
		sleep 0.1;
		jobs &> /dev/null
	done
}



pathIn="$PATH_SEQ_REF"/"$TAX_PREV_LEV"/
pathOut="$PATH_SEQ_REF"/"$TAX_SING_LEV"/

pathSCOREafterClassify=$DIR_MAIN_SCORE_CLASS/tax_"$TAX_PREV_LEV"/  # scores after previous classification
echo -e "**** pathSCOREafterClassify $pathSCOREafterClassify **** \n\n "

if [[ ! -d $pathIn ]]; then
	echo "$pathIn do not exist."
	exit;
fi  


if [[ ! -d $pathOut ]]; then
	mkdir $pathOut;
fi  

it=0;
for folderOTU in "$pathSCOREafterClassify"/OTU_*  
do



	pathSCORE="$folderOTU"/NEW/
	if [[ ! -d "$pathSCORE" ]]; then
		pathSCORE="$folderOTU"/;
	fi 

	if [ $TAKE_MM -eq 1 ]; then
		pathSCORE="$pathSCORE"/MisMatch_all"$CLASS_ALL"_sim"$PROC_SIM"
	else
		pathSCORE="$pathSCORE"/Match_all"$CLASS_ALL"_sim"$PROC_SIM"
	fi
	
	if [[ ! -d $pathSCORE ]]; then
		pathSCORE=$(echo $pathSCORE | perl -nle '$_=~ s/all/rand/; print $_;')
	fi

 
	echo -e "\n\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	echo "pathSCORE $pathSCORE"
	

	str_test="`ls "$pathSCORE"`"
	if [[ -z $str_test ]] # empty folder
	then
		print "Empty folder"
		echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		continue;
	fi
	echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	


	for y in "$pathSCORE"/"$PREF_META_SET"*_k"$LEN_KMER"*
	do

		if [[ ! -e $y ]]
		then
			continue;
		fi

		fileTr_prefix=${y%\"*};	# 
		# echo "1 $fileTr_prefix";
		fileTr_prefix=${fileTr_prefix#\"*};	#
		# echo "2 $fileTr_prefix";
		fileTr_prefix=${fileTr_prefix##*/};	# name without path
		echo FILE: $fileTr_prefix
		fileTr_prefix=${fileTr_prefix#$PREF_META_SET*}; #  name without group names  
		# echo FILE: $fileTr_prefix
		fileTr_Pre=${fileTr_prefix%.*}; # 
		# echo "fileTr_Pre:  $fileTr_Pre "
		fileTr_NAME=${fileTr_Pre%_idTAX*};	# 
		fileTr_NAME=${fileTr_NAME##*\_};	# idNAME
		
		fileTr_Pre=${fileTr_Pre##*\_};	# idTAX
		echo "fileTr_Pre:  $fileTr_Pre  fileTr_NAME  $fileTr_NAME"

		if [[ $fileTr_Pre = *"notknow"* || $fileTr_Pre = *"FP"*  || $fileTr_Pre = *"TP"*  \
			|| $fileTr_Pre = *"NK"* || $fileTr_Pre = *"MultiClass"* ]]
		then
			continue;				
		else

			echo "..."
			fileTr_prefix="$fileTr_Pre"
			b_NIEist=1;
			echo "$TAX_PREV_LEV": $fileTr_prefix 

			suffix=".fa";
			if [[ $fileTr_NAME = *"undef"* ]]
			then
				suffix="$fileTr_Pre".fa
				# echo UNDEF: $y
			fi			
			
			for x in "$pathOut"/"$fileTr_NAME"*"$suffix"	
			do 
				# echo x: $x
				if [[ -e $x ]]
				then
					echo "$fileTr_NAME ($fileTr_prefix) has been divided"
					b_NIEist=0;
				fi
				break;
			done


			if [[ $b_NIEist -eq 1 ]]
			then
				for files in "$pathIn"/*_"$fileTr_prefix".fa
				do
					if [[ ! -e $files  ]]
					then
						echo "$files do not exist"
						continue;
					fi
			
					file=${files%\"*};	# 
					file=${file#\"*};	# 
					file=${file##*/};	# 
				
					echo "Start divide "$file" file"

					it=$[it+1]
					perl $DIR_MAIN_BIN/class_seq_to_taxon_all.pl -i "$file" -Pin "$pathIn" -Pout "$pathOut"   -LF "$TAX_PREV_LEV"  -LS "$TAX_SING_LEV" -Pncbi "$DIR_TAX_NCBI" \
						 > $DIR_MAIN_BIN/"$tempout""$it".o  2> $DIR_MAIN_BIN/"$tempout""$it".e &
					job_wait

				done
			fi			
		fi

	done #y


done # folderOTU

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






