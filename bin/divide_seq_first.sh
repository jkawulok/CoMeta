#!/bin/bash
#!/usr/bin/perl

#	This file is a part of CoMeta software.
#	The homepage of the CoMeta project is http://sun.aei.polsl.pl/cometa




echo -e "\n\n####################### Divide sequence ###########################" >&2
ERROR_PROC=1;
pathOut="$PATH_SEQ_REF"/"$TAX_SING_LEV"/
fileout="$PREFIX_SEQ_REF"_"$SUFFIX_SEQ_REF"_outtemp.txt


if [[ ! -d $pathOut ]]; then
	mkdir $pathOut;
fi  


# echo pathOut $pathOut
# echo fileout $fileout

touch "$pathOut""$fileout"


for fileNW in "$PATH_IN_START_SEQ_REF"/"$PREFIX_SEQ_REF"*"$SUFFIX_SEQ_REF"
do

	echo $fileNW
	cat  $fileNW >> "$pathOut""$fileout"

done



perl $DIR_MAIN_BIN/class_seq_to_taxon_all.pl -i $fileout -Pin $pathOut -Pout $pathOut   -LF $TAX_PREV_LEV  -LS $TAX_SING_LEV -Pncbi $DIR_TAX_NCBI;

rm "$pathOut""$fileout"

ERROR_PROC=0;