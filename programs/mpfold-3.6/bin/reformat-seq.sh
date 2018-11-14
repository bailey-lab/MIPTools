#!/bin/bash

export _POSIX2_VERSION=0

# reformat-seq.sh
# Input: Sequence file $1
# FILE_PREFIX: The sequence file name is stripped of .seq, .fasta, etc.
# Output file is $FILE_PREFIX-local.seq
# Cases:

# IG (Intelligenetics) is left as is an the output file is just a
# symbolic link to the input file.

# FASTA is left as is but the header line is trucated at 80 characters
# and the sequence itself has all white space removed and is wrapped
# to 75 nt per record. 

# Both GenBank and EMBL are parsed properly. White spaces are
# compressed to single spaces. The output is is FASTA format.

# GCG is no longer accepted.

# If the file does not appear to be in IG, FASTA, GB or EMBL format,
# then it is treated as raw sequence data. The internal name for the
# sequence is $FILE_PREFIX and the output is in FASTA format.

if [ $# -lt 1 ] ; then
	echo -e ' *** Usage: reformat-seq.sh in_fil,\nwhere in_fil is the name of the input sequence file.'
	exit 1
elif [ ! -s $1 ] ; then
    echo -e "The file, $1, does not exist or is empty"
    exit 2
fi

FILE_PREFIX=`basename $1 .seq|sed -e 's/\.gb$//'|sed -e 's/\.embl$//'|\
            sed -e 's/\.fasta$//'|sed -e 's/\.SEQ$//'|sed -e 's/\.GB$//'|\
            sed -e 's/\.EMBL$//'|sed -e 's/\.FASTA$//'`

#Test for GenBank
GB=`grep -n '^LOCUS ' $1`
if [ $? = 0 ] ; then
    NAME_INDEX=`echo -e $GB|cut -d: -f1`
    NAME=`tail -n+$NAME_INDEX $1|head -n1|tr -s ' ' ' '|sed -e 's/^ //'|cut -c1-72`
    START=`grep -n '^ORIGIN ' $1`
    if [ $? = 0 ] ; then
	START=`echo -e $START|cut -d: -f1`
	START=`expr $START + 1`
    else
	echo -e 'Corrupted GenBank file. No ORIGIN line found.'
	exit 3
    fi
    STOP=`grep -n '^//' $1`
    if [ ! $? = 0 ] ; then
	STOP=`wc $i|tr -s ' ' ' '|sed -e 's/^ //'`
    else
	STOP=`echo -e $STOP|cut -d: -f1`
    fi

    echo -e ">$NAME" > ${FILE_PREFIX}-local.seq
    STOP=`expr $STOP - 1`
    NLINES=`expr $STOP - $START + 1`
    tail -n+$START $1|tr -cd '\012 A-Za-z'|tr -s ' ' ' '|head -n$NLINES \
	>> ${FILE_PREFIX}-local.seq
    echo -e "$FILE_PREFIX"
    exit 0
fi

#Test for EMBL
EMBL=`grep -n '^ID ' $1`
if [ $? = 0 ] ; then
    NAME_INDEX=`echo -e $EMBL|cut -d: -f1`
    NAME=`tail -n+$NAME_INDEX $1|head -n1|tr -s ' ' ' '|sed -e 's/^ //'|cut -c1-72`
    START=`grep -n '^SQ ' $1`
    if [ $? = 0 ] ; then
	START=`echo -e $START|cut -d: -f1`
	START=`expr $START + 1`
    else
	echo -e 'Corrupted EMBL file. No SQ line found.'
	exit 3
    fi
    STOP=`grep -n '^//' $1`
    if [ ! $? = 0 ] ; then
	STOP=`wc $i|tr -s ' ' ' '|sed 's/^ //'`
    else
	STOP=`echo -e $STOP|cut -d: -f1`
    fi

    echo -e ">$NAME" > ${FILE_PREFIX}-local.seq
    STOP=`expr $STOP - 1`
    NLINES=`expr $STOP - $START + 1`
    tail -n+$START $1|tr -cd '\012 A-Za-z'|tr -s ' ' ' '|head -n$NLINES \
	>> ${FILE_PREFIX}-local.seq
    echo -e "$FILE_PREFIX"
    exit 0
fi

#Test for FASTA
FASTA=`head -n1 $1|grep '^>'`
if [ $? = 0 ] ; then
    NAME=`echo -e $FASTA|cut -c2-80`
    echo -e ">$NAME" > ${FILE_PREFIX}-local.seq
    tail -n+2 $1|tr -cd 'A-Za-z'|fold -w 75 >> ${FILE_PREFIX}-local.seq
    echo -e '' >> ${FILE_PREFIX}-local.seq
    echo -e "$FILE_PREFIX"
    exit 0
fi

#Test for IG
IG=`head -n1 $1|grep '^;'`
if [ $? = 0 ] ; then
    rm -f ${FILE_PREFIX}-local.seq
    ln -s $1 ${FILE_PREFIX}-local.seq
    echo -e "$FILE_PREFIX"
    exit 0
fi

#At this point, assume raw sequence data and hope for the best
echo -e ">$FILE_PREFIX" > ${FILE_PREFIX}-local.seq
cat $1|tr -cd 'A-Za-z'|fold -w 75 >> ${FILE_PREFIX}-local.seq
echo -e '' >> ${FILE_PREFIX}-local.seq
echo -e "$FILE_PREFIX"
exit 0

