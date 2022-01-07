#!/bin/bash
a1=10; a2=100; a3=1000;

rm distdump dist_new newest_dist tmp

FILE="mindist.xvg"
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

g_analyze -f $b/$FILE > tmp 
cat tmp | grep $choice > distdump

cat distdump | sed "s/$choice/$lam/g" > dist_new
cat tmp | grep err.est | sed "s/Set   1:  err.est.//g" >> dist_new
sed '$!N;s/\n/ /' dist_new >> dist.dat

done

#cut -c1-77 dist_energy > tmp
#cat tmp > dist_energy
