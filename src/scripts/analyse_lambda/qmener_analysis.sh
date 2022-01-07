#!/bin/bash
a1=10; a2=100; a3=1000;

rm energydump energy_new newest_energy tmp

choice=SS1

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

cp $b/qmener.dat $b/qmener.xvg
g_analyze -f $b/qmener.xvg > tmp
cat tmp | grep $choice > energydump

cat energydump | sed "s/$choice/$lam/g" > energy_new
sed '$!N;s/\n/ /' energy_new >> newest_energy

done

#cut -c1-77 newest_energy > tmp
#cat tmp > newest_energy

mv newest_energy qmene_av.dat
