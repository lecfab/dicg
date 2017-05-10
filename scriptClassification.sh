#!/bin/zsh
zmodload -i zsh/mathfunc

##YOUR Classification progname
CLASSIFPROG="python3 src/classificationOld.py"

#nb of source images to test
NBIMGTESTS=20

#nb of noisyfied copies
NBTESTS=2

#No noise here
MAXNOISE=1.0
sum=0.0
variance=0.0
cpt=0.0

for ((i=0; i < $NBIMGTESTS; i++)); do
    ##Pick a random class
    CLASSID=`expr $RANDOM % 61 + 1`
    CLASSNAME=`head -$CLASSID classes.csv | tail -1 | sed -e 's/,//'`

    ##Pick a random image in this class
    IMGID=`expr $RANDOM % 10 + 1`
    IMGNAME=database/$CLASSNAME-$IMGID.pgm
    echo "Classname: "$CLASSNAME "\t\t"$IMGNAME

    if [[ -e $IMGNAME ]]; then
        for ((j=0; j < $NBTESTS; j+=1)); do
            ##Random scale+rotation
            ANGLE=$((rand48()*3.1415))
            SCALE=$((rand48()*3))
            NOISE=$((rand48()*MAXNOISE))
            ./imgRotate -i $IMGNAME -o tmp.pgm -a $ANGLE 2> /dev/null
            ./imgScale -i tmp.pgm -o tmp2.pgm -s $SCALE 2> /dev/null
            ./imgAddNoise -i tmp2.pgm -o tmp.pgm -n $NOISE  2> /dev/null
            #mv tmp2.pgm tmp.pgm

            ##Running the retrieval
            eval `echo $CLASSIFPROG tmp.pgm` >! scores_tmp.txt
            RANK=` ./getRank scores_tmp.txt $CLASSID`
            echo "Rank=$RANK\t\t$ANGLE\t$SCALE\t$NOISE"
            
            ##Number of correct results in the first 1
            if [[ $RANK -le 1 ]]; then
               correct1=$(($correct1 + 1.0))
            fi
            ##Number of correct results in the first 10 
            if [[ $RANK -le 10 ]]; then
               correct=$(($correct + 1.0))
            fi
            cpt=$(($cpt + 1.0))
        done
    fi
done
echo
echo "Total number of tests= "$cpt
echo "Number of 'perfect' classification (first classe)= "$correct1
echo "Perfect score= " $(( $correct1 / $cpt ))
echo "Number of 'correct' classification (first 10 classes)= "$correct
echo "Final score= " $(( $correct / $cpt ))
rm tmp.pgm
#rm tmp2.pgm
rm scores_tmp.txt
