#! /bin/sh

# reading 4 files from a file list for each run

read -p " First Set Index  : " s2
read -p " Total Set(each set with 4 files)   : " s3

for (( i=0; i<$s3; i++ ))
do

s1=$(($s2*4-3))
r=$(($s2+$i))
t=$(($s1+$i*4))

RUNFOLDER="QCD_$r"

# create python file
PATOUT="qcd_fall08_$r.root"
JETOUT="qcd_JetEtAnalysis_$r.root"
MUONOUT="qcd_IsoMuAnalysis_$r.root"
echo $PATOUT
RUNPY="qcd_$r.py" 
echo $RUNPY 
cat qcd4.py | sed s/QCDPAT1/QCD_PATSkim1_$t/ | sed s/QCDPAT2/QCD_PATSkim1_$(($t+1))/ | sed s/QCDPAT3/QCD_PATSkim1_$(($t+2))/ | sed s/QCDPAT4/QCD_PATSkim1_$(($t+3))/ |  sed s/ANAFILE/$PATOUT/ |  sed s/JETFILE/$JETOUT/ |  sed s/MUFILE/$MUONOUT/> $RUNPY

# create run script
RUNSH="qcd_$r.sh" 
echo $RUNSH 
cat submit.sh | sed s/RUNSCRIPT/$RUNPY/ | sed s/RUNDIR/$RUNFOLDER/ > $RUNSH
chmod a+x $RUNSH

# create submit file
RUNSUB="qcd_$r.sub" 
echo $RUNSUB 
cat condor_0.sub | sed s/CONDOR_SH/$RUNSH/ | sed s/RUNDIR1/$RUNFOLDER/ > $RUNSUB

# create run folder and move all stuff in this folder
mkdir $RUNFOLDER
mv $RUNSUB $RUNFOLDER
mv $RUNSH $RUNFOLDER
mv $RUNPY $RUNFOLDER
cd $RUNFOLDER
condor_submit $RUNSUB
sleep 1s
cd ../

done
