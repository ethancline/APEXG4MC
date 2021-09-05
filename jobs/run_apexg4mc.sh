#!/bin/sh

# --------------------------------------------------------------------------------------------
# set-up
# --------------------------------------------------------------------------------------------
export JOBDIR=`pwd`/$1
mkdir -p $JOBDIR
cd $JOBDIR
ln -sf ../../build/APEXG4MC APEXG4MC
ln -sf ../../build/Septa-JB_map.table Septa-JB_map.table
ln -sf ../../build/macros macros
# --------------------------------------------------------------------------------------------
# run geant4 and convert root tree for SIMC
# --------------------------------------------------------------------------------------------
./APEXG4MC macros/$2.mac
cp ../ConvertToSIMC.C .
root -b -q ConvertToSIMC.C
cp ../ext_part1 $SIMC/infiles
cp ../ext_part2 $SIMC/infiles
cp ../ext_part3 $SIMC/infiles
cp ../inpL.txt  $SIMC
cp ../inpR.txt  $SIMC
cp ../ext_part3 $SIMC/infiles
cp LQ1_TCS.dat $SIMC/infiles
cp RQ1_TCS.dat $SIMC/infiles
# --------------------------------------------------------------------------------------------
# run SIMC for LHRS and convert output to root
# --------------------------------------------------------------------------------------------
cd $SIMC/infiles
ln -sf LQ1_TCS.dat generated_events.dat
nn=`cat generated_events.dat | wc -l`
rm -f extendedL.inp
cat ext_part1 > extendedL.inp
echo "     $nn	Monte-Carlo trials" >> extendedL.inp
cat ext_part2 >> extendedL.inp
scp $SIMC/infiles/* david@npc69.physics.gla.ac.uk:/home/david/geant4/aAPEX_G4MC/SIMC/infiles
cd ..
./mc_hrs_single < inpL.txt
cd worksim
h2root extendedl.rzdat 
cp extendedl.root $JOBDIR
# --------------------------------------------------------------------------------------------
# run SIMC for RHRS and convert output to root
# --------------------------------------------------------------------------------------------
cd $JOBDIR
cd $SIMC/infiles
ln -sf RQ1_TCS.dat generated_events.dat
nn=`cat generated_events.dat | wc -l`
rm -f extendedR.inp
cat ext_part1 > extendedR.inp
echo "     $nn	Monte-Carlo trials" >> extendedR.inp
cat ext_part3 >> extendedR.inp
scp $SIMC/infiles/* david@npc69.physics.gla.ac.uk:/home/david/geant4/aAPEX_G4MC/SIMC/infiles
cd ..
./mc_hrs_single < inpR.txt
cd worksim
h2root extendedr.rzdat 
cp extendedr.root $JOBDIR
# --------------------------------------------------------------------------------------------
# merge original G4 tree with simc trees
# --------------------------------------------------------------------------------------------
cd $JOBDIR
cp ../MergeG4SIMC.C .
root -b -q MergeG4SIMC.C
