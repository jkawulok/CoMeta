#!/bin/env bash 

#	This file is a part of CoMeta software.
#	The homepage of the CoMeta project is http://sun.aei.polsl.pl/cometa


ERROR_PROC=1;
tempout=TEMPOUT_"$PBS_ARRAYID"_"$TAX_SING_LEV"_"$RANDOM"_class2one_

echo -e "\n\n####################### Classify reads part2 ###########################" >&2

# +++++++++++++++++++++++++
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
# +++++++++++++++++++++++++
Trans=0 

FOLDER_List=$DIR_MAIN_SCORE_CLASS/tax_"$TAX_SING_LEV"/ ; #folder where folders are classified


echo "################### numK $LEN_KMER  $FILE_META_SET  OTU: $TAX_SING_LEV  MM=$TAKE_MM  classALL=$CLASS_ALL  proc_sim=$PROC_SIM    ########################" ;
echo ""


function process {


	pFolder="$1"


	echo FOLDER $pFolder;

	pathSCORES=$pFolder/NEW/
	if [[ ! -d $pathSCORES ]]; then
		pathSCORES=$pFolder/;
	fi 	
	

	folderOTU=${pFolder%\"*};	# 
	folderOTU=${folderOTU#\"*};	# 
	folderOTU=${folderOTU##*/};	# 
	OTUH=${folderOTU#OTU_*};  #

	rand_num=$RANDOM
	name_file=name_files_"$PBS_ARRAYID"_"$TAKE_MM"_"$CLASS_ALL"_"$LEN_KMER"_sim"$PROC_SIM"_"$rand_num".txt
	name_group=name_groups_"$PBS_ARRAYID"_"$TAKE_MM"_"$CLASS_ALL"_"$LEN_KMER"_sim"$PROC_SIM"_"$rand_num".txt
	
	"$DIR_MAIN_BIN"/genlist -nw1 -pw"$pathSCORES" -ns"" -nm"" -mc$MATCH_CUT_OFF -k$LEN_KMER    -NF$name_file -NG$name_group  

	sth_is=0;
	echo "Create group files"	
	for x in "$pathSCORES"/*_k"$LEN_KMER"*"$MATCH_CUT_OFF"_M.out
	do 
		
		if test -e "$x"  # file exist
		then
			if test ! -s "$x"  # empty file
			then
				if [ $TAKE_MM -eq 1 ]; then
					sth_is=1;
				else
					continue;  #					
				fi
			else
				sth_is=1; #not empty file
			fi
		else
			continue; # not exist
		fi


		#echo x $x
		fileTr_prefix=${x%\"*};	# 
		# echo "1 $fileTr_prefix";
		fileTr_prefix=${fileTr_prefix#\"*};	# 
		# echo "2 $fileTr_prefix";
		fileTr_prefix=${fileTr_prefix##*/};	# 
		# fileTr_prefix=${fileTr_prefix%.*};	# 
		fileTr_prefix=${fileTr_prefix%%M.*};	#

		# echo "3 $fileTr_prefix";

		nameGroupS=${fileTr_prefix%%_k*};  #group
		nameGroup=${nameGroupS#*"$OTUH"_*};
		# echo "groupS: $nameGroupS";
		#echo "group: $nameGroup";

		
		#if [ "$TAX_SING_LEV" = "phylum" ]
		#then
		#	nameGroup=${nameGroup#*_}
		#fi


		if [ "$nameGroup" = "undef" ]
		then

				nameGroup=${nameGroupS#NW_*};  #group np. proteobacteria_undef_IDTAX345
				echo "grupaINNE: $nameGroup";
		fi

		# Create a list of groups and files
		"$DIR_MAIN_BIN"/genlist -nw0 -pw"$pathSCORES" -ns"$fileTr_prefix"    -ng"$nameGroup"   -nm"$FILE_META_SET" -mc -k$LEN_KMER  -MM$TAKE_MM   -NF$name_file -NG$name_group 
		
	done # x


	# proper grouping-------------------------------------------------------------------------------
	if [ $TAKE_MM -eq 1 ]; then
		path_OU=$pathSCORES/MisMatch_all"$CLASS_ALL"_sim"$PROC_SIM"
	else
		path_OU=$pathSCORES/Match_all"$CLASS_ALL"_sim"$PROC_SIM"
	fi
	


	if [ ! -d "$path_OU" ]; then
		mkdir "$path_OU";
	fi 

	name_group_class="$PREF_META_SET"_k"$LEN_KMER"_"$OTUH"_"$SUFFIX_DESCRP"

	if [ $sth_is -eq 1 ]
	then
		echo "Start class to best score"
		#echo "$DIR_MAIN_BIN"/class2best -NO"$name_group_class" -WO"$path_OU"  -WI"$pathSCORES" -WS"$FOLDER_List" -P1 -K"$LEN_KMER"  -Tr$Trans  -mc"$MATCH_CUT_OFF"  -cl"$CLASS_ALL" -proc"$PROC_SIM"   -NR"$PATH_META_SET"/"$FILE_META_SET"    -ch"$CHECK_CLASS"    -NF"$name_file" -NG"$name_group"  -MM"$TAKE_MM"
		
		
		"$DIR_MAIN_BIN"/class2best -NO"$name_group_class" -WO"$path_OU"  -WI"$pathSCORES" -WS"$FOLDER_List" \
		-P1 -K"$LEN_KMER"  -Tr$Trans  -mc"$MATCH_CUT_OFF"  -cl"$CLASS_ALL" -proc"$PROC_SIM"   \
		-NR"$PATH_META_SET"/"$FILE_META_SET"    -ch"$CHECK_CLASS"    -NF"$name_file" -NG"$name_group"  -MM"$TAKE_MM"



		if [ $CLASS_ALL -eq 1 ]
		then
			echo Check number of classification reads 
			nameSUM="$name_group_class"_mm"$TAKE_MM"_cl"$CLASS_ALL"
			"$DIR_MAIN_BIN"/num_class  -NO"$nameSUM"  -WO"$FOLDER_List" -WI"$pathSCORES"  -K$LEN_KMER -mc"$MATCH_CUT_OFF" -FONum_class_"$PBS_ARRAYID".txt -NF$name_file -NG$name_group
		fi

	else
		echo "Not exist for" "$PREF_META_SET"_k"$LEN_KMER"_"$OTUH";
	fi

	echo ""
	echo ""
	
	rm "$pathSCORES"/"$name_file"
	rm "$pathSCORES"/"$name_group"
}
	
# +++++++++++++++++++++++++++	



it=0;	
for PAR in $FOLDER_List/OTU_*   # for each folder (each TAX_SING_LEV) at a given level
do
	
	it=$[it+1]
	process "$PAR" > $DIR_MAIN_BIN/"$tempout""$it".o  2> $DIR_MAIN_BIN/"$tempout""$it".e &
	
done

wait




for (( i=1; i<=it; i++ ))
do
	echo -e "\n\niiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii $i iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii"
	cat $DIR_MAIN_BIN/"$tempout""$i".e >&2
	cat $DIR_MAIN_BIN/"$tempout""$i".o
	rm $DIR_MAIN_BIN/"$tempout""$i".e
	rm $DIR_MAIN_BIN/"$tempout""$i".o
done

ERROR_PROC=0;