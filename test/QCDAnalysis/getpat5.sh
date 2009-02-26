#! /bin/sh

# reading 10 files from a file list for each run

read -p " Run Lumi Number(each for 2/pb) : " s

echo " Total Lumi= "$(($s*2))"/pb   Total Event= "$(($s*56830)) 

for (( i=11; i<$s; i++ ))
do

SS=$(($i*56830))  

r=$(($i+1))

RUNFOLDER="QCD_$r"

# create python file
PATOUT="qcd_fall08_$r.root"
echo $PATOUT
RUNPY="qcd_$r.py" 
echo $RUNPY 
cat qcd.py | sed s/SKIPNUM/$SS/ |  sed s/ANAFILE/$PATOUT/ > $RUNPY

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
