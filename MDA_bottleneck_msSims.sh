#!/bin/sh
##read from stdin args ./MDA_bottleneck_msSims.sh 8.438 0.0000342 5 6 .30
THETA=$1    #8.438
MU=$2       #.0000342  #per gene of 450bp
GENS_YEAR=$3    #5
STAT=$4        #6 TAJD, X Theta
##population sizes##"scale=9; $x *1.423" | bc
MUS=$( echo "scale=9; $MU * 4" | bc )
NOW_N=$( echo "scale=1; $THETA / $MUS"| bc )  #effective size for time scaling
MDA_N=$5      #% reduced by the MDA; population was 40% of current size; 30,50,75,90
ANCESTRAL_N=5  #population is 5 times larger than current size
##time##
YEAR=$( echo "scale=12; $GENS_YEAR / $NOW_N"|bc )

#call_ms number_inds reps flag theta change_size time size_to_Now pip samplestats cut_only_column sstats

#before MDA; #1990-1993

ms 100 5000 -t $THETA -eN 0 $ANCESTRAL_N | sample_stats | cut -f $STAT | stats 0.025 0.975 > ms_MDAsims_$MDA_N.txt

i=0
while [ $i -le 5 ]
do
#ms 100 100 -t $THETA -eN 0 $ANCESTRAL_N | msstats | cut -f $STAT | python sstats.py > ms_MDAsims_$MDA_N.txt
    ms 100 5000 -t $THETA -eN 0 $ANCESTRAL_N | sample_stats | cut -f $STAT | stats 0.025 0.975 >> ms_MDAsims_$MDA_N.txt
    i=$(( $i + 1 ))
done

##during MDA; #1994 -1998
for s in 1 2 3 4 5
do
#ms 100 100 -t $THETA -eN 0 $MDA_N -eN $( echo "scale=12; $YEAR * $s"|bc ) $ANCESTRAL_N | msstats | cut -f $STAT | python sstats.py >> ms_MDAsims_$MDA_N.txt
ms 100 5000 -t $THETA -eN 0 $MDA_N -eN $( echo "scale=12; $YEAR * $s"|bc ) $ANCESTRAL_N | sample_stats | cut -f $STAT | stats 0.025 0.975 >> ms_MDAsims_$MDA_N.txt
done

##after MDA round 1; #1999-2013
for s in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
do
#ms 100 100 -t $THETA -eN $( echo "scale=12; $YEAR * $s"|bc ) $MDA_N -eN $( echo "scale=12; $YEAR * $[ $s + 5 ]"|bc ) $ANCESTRAL_N | msstats | cut -f $STAT | python sstats.py >> ms_MDAsims_$MDA_N.txt
ms 100 5000 -t $THETA -eN $( echo "scale=12; $YEAR * $s"|bc ) $MDA_N -eN $( echo "scale=12; $YEAR * $[ $s + 5 ]"|bc ) $ANCESTRAL_N | sample_stats | cut -f $STAT | stats 0.025 0.975 >> ms_MDAsims_$MDA_N.txt
done

##during MDA round 2; #2013-2018
for s in 1 2 3 4 5
do
#ms 100 100 -t $THETA -eN 0 $MDA_N -eN $( echo "scale=12; $YEAR * $s"|bc ) 1 -eN $( echo "scale=12; $YEAR * $[ $s + 15 ]"|bc ) $MDA_N -eN $( echo "scale=12; $YEAR * $[ $s + 20 ]"|bc ) $ANCESTRAL_N | msstats | cut -f $STAT | python sstats.py > ms_MDAsims_$MDA_N.txt
ms 100 5000 -t $THETA -eN 0 $MDA_N -eN $( echo "scale=12; $YEAR * $s"|bc ) 1 -eN $( echo "scale=12; $YEAR * $[ $s + 15 ]"|bc ) $MDA_N -eN $( echo "scale=12; $YEAR * $[ $s + 20 ]"|bc ) $ANCESTRAL_N | sample_stats | cut -f $STAT | stats 0.025 0.975 >> ms_MDAsims_$MDA_N.txt
done

##after MDA round 2; #2018-2023
#for s in 1 2 3 4 5
#do
#./ms 100 1000 -t $THETA -eN $[ $1_YEAR * $s ] $MDA_N -eN $ [ $1_YEAR * $[ $s + 5 ] ] 1 -eN $[ $1_YEAR * $[ $s + 15 ] ] $MDA_N -eN $[ $1_YEAR * $[ $s + 20 ] ] $ANCESTRAL_N | msstats | cut -f $STAT | python sstats.py > ms_MDAsims_$MDA_N.txt
#done
