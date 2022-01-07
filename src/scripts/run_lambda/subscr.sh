#!/bin/bash
a1=10; a2=100; a3=1000;

for k in `seq 0 2 100`;
do 
if [ $k -lt $a1 ] 
then
b="l0.0$k"
else
    if [ $k -lt $a2 ]
    then
        b="l0.$k"
        
    else
        if [ $k -lt $a3 ]
            then
               b="l1.00"
               fi
    fi
fi

cd $b
qsub runscript.sh
cd ..

done

