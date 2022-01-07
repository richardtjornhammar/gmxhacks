#!/bin/bash
a1=10; a2=100; a3=1000;

rm energydump energy_new newest_energy tmp

choice=SS6
touch qmcstd.log

for k in `seq 0 2 100`;
do 
if [ $k -lt $a1 ] 
then
b="l0.0$k"
lam="0.0$k"
else
    if [ $k -lt $a2 ]
    then
        b="l0.$k"
        lam="0.$k"
    else
        if [ $k -lt $a3 ]
            then
               b="l1.00"
	       lam="1.00"
               fi
    fi
fi

cat $b/qmcsigma.dat >> qmcstd.log
done
