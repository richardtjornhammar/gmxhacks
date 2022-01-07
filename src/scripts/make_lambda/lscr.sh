#!/bin/bash

ext="mdp"
dirpref="lambda-"
filnamn1="accc"
filnamn2="wrk-$filnamn1"
dirnamn="l"
base="."

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

#mkdir $base/$dirnamn$b
cp $base/qm* $dirpref$base/$dirnamn$b/.
cp $base/*itp $dirpref$base/$dirnamn$b/.
cp $base/ff* $dirpref$base/$dirnamn$b/.
cp runscript.sh $dirpref$base/$dirnamn$b/.
intext="$ftext = $b"
echo "$ftext and $intext"
cat $base/$filnamn1.$ext | grep -v $ftext >  $base/$dirnamn$b/$filnamn2.$ext
echo "bQMlambda = yes" >> $base/$dirnamn$b/$filnamn2.$ext
echo $intext >> $base/$dirnamn$b/$filnamn2.$ext
cd $base/$dirnamn$b
grompp -f $filnamn2.$ext -c $filnamn1.gro -p $filnamn1.top -n $filnamn1.ndx -o $filnamn2.tpr -maxwarn 1
chmod ug+rwx *
rm \#*
cd ../..
done

