#!/bin/bash
#sims
#seq .01 .05 1
for m in $(seq 1 .05 100);do
    n=0
    while [ $n -lt 1000 ]; do
        ./ms 40 1 -T -I 2 20 20 $m | tail -n +4 | grep -v // >treefile
        ./seq-gen -mHKY -l 1000 -s .01 <treefile >seqfile
        head -1 seqfile >seqfile1;tail -n +2 seqfile | sort -n >>seqfile1
        ./dnadist < myoptionfile_dnadist.txt
        python Anderson_sim.py
        let n=n+1 
    done
done

