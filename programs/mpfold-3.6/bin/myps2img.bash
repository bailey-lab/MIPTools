#!/bin/bash

if [ $# = 0 ] ; then
    echo -e "Please enter the name of a PostScript file."
    exit 1
fi

NAME=$1
CALL=`basename $0 .sh`
CALL=`basename $CALL .bash`
DIRNAME=`dirname $NAME`
PREFIX=`basename $NAME .ps`
PREFIX=`basename $PREFIX .eps`

if [ "$CALL" = "myps2pdf" ] ; then
    epstopdf $NAME || convert $NAME $DIRNAME/${PREFIX}.pdf
elif [ "$CALL" = "myps2jpg" ] ; then
    convert $NAME ${DIRNAME}/${PREFIX}.jpg
elif [ "$CALL" = "myps2png" ] ; then
    convert $NAME ${DIRNAME}/${PREFIX}.png
else
    echo -e "This script must be called as myps2pdf.bash, myps2jpg.bash or as myps2png.bash."
    exit 2
fi
