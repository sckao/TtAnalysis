#! /bin/sh

input=$1
output=$2

if [ -z $1 ]; then
   echo "Please specify the source file you want to compile"
   exit
fi

if [ -z $2 ]; then
   output=a.out
fi

cflags=`sh $ROOTSYS/bin/root-config --cflags`
libs=`sh $ROOTSYS/bin/root-config --libs`
glibs=`sh $ROOTSYS/bin/root-config --glibs`
cxxflags="-g -fPIC -Wno-deprecated -O -ansi -D_GNU_SOURCE -g -O2"
cxx="-m32"

g++ MassFitFunction.cc -c -o MassFitFunction.o $libs $cflags $cxx $cxxflags
g++ MassAnaInput.cc -c  -o MassAnaInput.o $libs $cflags $cxx $cxxflags
g++ MassAna.cc -c -o MassAna.o $libs $cflags $cxx $cxxflags $glibs
g++ AlgoZero.cc -c -o AlgoZero.o $libs $cflags $cxx $cxxflags $glibs
g++ AlgoKcon.cc -c -o AlgoKcon.o $libs $cflags $cxx $cxxflags $glibs
g++ JES.cc -c -o JES.o $libs $cflags $cxx $cxxflags $glibs
g++ MassAnaOutput.cc -c -o MassAnaOutput.o $libs $cflags $cxx $cxxflags $glibs
g++ HadWMassFitter.cc -c -o HadWMassFitter.o $libs $cflags $cxx $cxxflags $glibs -lMinuit
g++ LepTopMassFitter.cc -c -o LepTopMassFitter.o $libs $cflags $cxx $cxxflags $glibs -lMinuit
g++ PseudoExp.cc -c -o PseudoExp.o $libs $cflags $cxx $cxxflags $glibs
g++ WAnalysis.cc -c -o WAnalysis.o $libs -lMinuit $cflags $cxx $cxxflags $glibs 
g++ JetSpectrum.cc -c -o JetSpectrum.o $libs -lMinuit $cflags $cxx $cxxflags $glibs 
g++ BgEstimation.cc -c -o BgEstimation.o $libs -lMinuit $cflags $cxx $cxxflags $glibs 

g++ $input -o $output MassFitFunction.o MassAnaInput.o MassAna.o AlgoZero.o AlgoKcon.o JES.o MassAnaOutput.o HadWMassFitter.o LepTopMassFitter.o PseudoExp.o WAnalysis.o JetSpectrum.o BgEstimation.o $libs $cflags $cxx $cxxflags $glibs -lMinuit
#g++ $input -o $output MassFitFunction.o MassAnaInput.o MassAna.o AlgoZero.o AlgoKcon.o JES.o MassAnaOutput.o PseudoExp.o $libs $cflags $cxx $cxxflags
