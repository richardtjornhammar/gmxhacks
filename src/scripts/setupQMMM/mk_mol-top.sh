#!/bin/bash
touch generictop
touch generictail

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
echo "    O    8 	0.000      0.000     A  0.0022619536  7.4149321e-07 " >> made.top
echo "   OM    8 	0.000      0.000     A  0.0022619536  7.4149321e-07 " >> made.top
echo "   OA    8 	0.000      0.000     A  0.0022619536  1.505529e-06 " >> made.top
echo "   OW    8 	0.000      0.000     A  0.0026173456  2.634129e-06 " >> made.top
echo "   HC    1 	0.000      0.000     A   8.464e-05  1.5129e-08 " >> made.top
echo "    H    1 	0.000      0.000     A           0           0 " >> made.top
echo " HCHL    1  	0.000      0.000     A  3.76996e-05  4.2999495e-09 " >> made.top
echo "ODMSO    8 	0.000      0.000     A  0.0022707131  7.5144626e-07 " >> made.top
echo " OWT3    8     15.9994     0.000     A  0.24889E-02   0.24352E-05" >> made.top
echo " OWT4    8     15.9994     0.000     A  0.25519e-02   0.25104e-05 " >> made.top
echo "  ORT    8     15.9994    -0.8476    A  0.0024889    2.4352e-06 " >> made.top
#0.0022619536  1.505529e-06     too loose  double CHECK
# 0.0024494  1.7856e-06  UNSURE BUT SEEMS BAD
#0.0024396496  1.123511105e-06 :: too dense
#0.0024396496  1.323511105e-06 :: looks promising for the liquid 
#0.0022619536    1.323511105e-06  !! BEST SO FAR !! might be because of QM .. too dense
#0.001897432  1.07150318350742e-06 uses pure rdf... not that bad ended on 960
# put in tip5p which SUCKS
# using TIP3P  0.0024889   2.4352e-6 
# sigma=0.3227078
#TRYING HC AND OA ::way to loose
#USING SPCe TYPE VDW INTERACTION IS WAY TOO LOOSE ( H AND OW )
echo "  HRT    1      1.0080     0.4238    A           0           0" >> made.top
echo " " >> made.top
echo "#define   OQ   OW" >> made.top
echo "#define   HQ   HRT" >> made.top
echo "#define   hq    0.4238     1.0080" >> made.top
echo "#define   oq   -0.8476    15.9994" >> made.top
echo " " >> made.top

for k in `seq 1 1 128`;
do 
echo " ">>generictop
echo " [ moleculetype ]" >> generictop
echo "; Name            nrexcl" >> generictop
echo "QMW$k              3" >> generictop
echo " ">> generictop
echo "[ atoms ]" >> generictop
echo ";   nr       type  resnr residue  atom   cgnr     charge       mass  typeB    chargeB	   massB" >> generictop
echo "    1          OQ      1    QMW     OW      1     oq ;" >> generictop
echo "    2          HQ      1    QMW    HW1      1     hq ;" >> generictop
echo "    3          HQ      1    QMW    HW2      1     hq ;" >> generictop

echo "QMW$k           1" >> generictail 
done

cat generictop >> made.top
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
cat generictail >> made.top

rm generic*


 

