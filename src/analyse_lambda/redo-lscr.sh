#!/bin/bash

ext="mdp"
dirpref="lambda-"
filnamn1="qm"
filnamn2="wrk-$filnamn1"
dirnamn="l"
base="zinc_in_cccc"

a1=10; a2=100; a3=1000;
ftext="init_lambda"

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

intext="$ftext = $b"
echo "$ftext and $intext"
cat dummy.$ext | grep -v $ftext >  $dirnamn$b/$filnamn2.$ext
echo "bQMlambda = yes" >> $dirnamn$b/$filnamn2.$ext
echo $intext >> $dirnamn$b/$filnamn2.$ext
cd $dirnamn$b
cp ../qm.* .
grompp -f $filnamn2.$ext -c qm.gro -p qm.top -n qm.ndx -o $filnamn2.tpr -maxwarn 1
chmod ug+rwx *
rm \#* step.* job.*
cd ..
done

