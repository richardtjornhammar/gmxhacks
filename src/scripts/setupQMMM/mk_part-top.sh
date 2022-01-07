#!/bin/bash

#echo "; Include forcefield parameters" > made.top
#echo "#include \"ffG43a1.itp\"">> made.top
#echo " " >> made.top
echo "[ defaults ]" > made.top
echo "; nbfunc	comb-rule	gen-pairs	fudgeLJ	fudgeQQ" >> made.top
echo "  1		1		no		1.0	1.0" >> made.top
echo " " >> made.top
echo "[ atomtypes ]" >> made.top
echo "; full atom spec found in ffG43a1nb.itp" >> made.top
echo "; name  at.num    mass      charge ptype        c6           c12" >> made.top
echo "  ORT    8     15.9994    -0.8476    A   0.0022619536  1.505529e-06" >> made.top
echo "  HRT    1      1.0080     0.4238    A           0           0" >> made.top
echo " " >> made.top
echo "#define   OQ   ORT" >> made.top
echo "#define   HQ   HRT" >> made.top
echo "#define   hq    0.4238     1.0080" >> made.top
echo "#define   oq   -0.8476    15.9994" >> made.top
echo " " >> made.top

for j in $(seq 1 1 64)
do
for i in $(seq 1 1 3)
do
echo " " >> made.top
echo " [ moleculetype ]" >> made.top
echo "; Name            nrexcl" >> made.top
   if [ $i = 1 ]; then
	echo "QMO$(($i + 3*($j - 1)))              3" >> made.top
   else 
        echo "QMH$(($i + 3*($j - 1)))              3" >> made.top
   fi
echo " " >> made.top
echo "[ atoms ]" >> made.top
echo ";   nr       type  resnr residue  atom   cgnr     charge       mass  typeB    chargeB	   massB" >> made.top
   if [ $i = 1 ]; then
	echo "    1          OQ      1    QMW     OW      1     oq ;" >> made.top
   else 
        echo "    1          HQ      1    QMW     HW      1     hq ;" >> made.top
   fi
done
done

echo " " >> made.top
echo "; Include Position restraint file"  >> made.top
echo "#ifdef POSRES " >> made.top
echo "#include \"posre.itp\" " >> made.top
echo "#endif " >> made.top
echo " " >> made.top
echo "; Include water topology " >> made.top
#echo "#include \"spce.itp\" " >> made.top
echo " " >> made.top
echo "#ifdef POSRES_WATER " >> made.top
echo "; Position restraint for each water oxygen " >> made.top
echo "[ position_restraints ] " >> made.top
echo ";  i funct       fcx        fcy        fcz " >> made.top
echo "   1    1       1000       1000       1000 " >> made.top
echo "#endif " >> made.top
echo " " >> made.top
echo "; Include generic topology for ions " >> made.top
echo "#include \"ions.itp\" " >> made.top
echo " " >> made.top
echo "[ system ] " >> made.top
echo "; Name " >> made.top
echo "WATER " >> made.top
echo " " >> made.top
echo "[ molecules ] " >> made.top
echo "; Compound        #mols " >> made.top
for j in $(seq 1 1 64)
do
for i in $(seq 1 1 3)
do
   if [ $i = 1 ]; then
	echo "QMO$(($i + 3*($j - 1)))              1" >> made.top
   else 
        echo "QMH$(($i + 3*($j - 1)))              1" >> made.top
   fi
done
done

echo " " >> made.top
