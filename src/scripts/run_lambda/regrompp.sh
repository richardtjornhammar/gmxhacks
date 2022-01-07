#!/bin/bash

ext="mdp"
dirpref="lambda_qm2mm"
filnamn1="qm"
filnamn2="wrk-$filnamn1"
dirnamn="l"
base="sch3cluster"
optional="-maxwarn 1"

a1=10; a2=100; a3=1000;
ftext="QMlambda"

#mkdir $dirpref$base

for k in `seq 0 2 100`;
do 
if [ $k -lt $a1 ] 
then
b="0.0$k"
else
    if [ $k -lt $a2 ]
    then
        b="0.$k"
        
    else
        if [ $k -lt $a3 ]
            then
               b="1.00"
               fi
    fi
fi

cd $dirnamn$b
grompp -f $filnamn2.$ext -c qm.gro -p qm.top -n qm.ndx -o $filnamn2.tpr $optional
chmod ug+rwx *
rm \#* step.* job.*
cd ..
done

