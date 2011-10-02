#ifndef PseudoExp_H
#define PseudoExp_H

#include <TObject.h>
#include <TMinuit.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TRandom.h>
#include <TRandom2.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <string>

#include "MassAnaInput.h"

class PseudoExp : public TObject {

private:


   MassAnaInput*    fitInput;
   vector<double> jetCuts;
   vector<double> muonCuts;

public:

   PseudoExp();     
   ~PseudoExp();     
 
   // random seed = 0 => using the system time for the seed
   vector< pair<int,int> > GetEnsemble( string fileName, double pMean, int RandomSeed = 0 );
   vector<int> GetEnsemble( string fileName, TString treeName, double pMean, int RandomSeed = 0 );
 
   vector<int> EventShuffle( int theSize, int RandomSeed );
 
   void PhaseSmearing( vector<TLorentzVector>& vs, int RandomSeed = 0, double jes = 0, bool ReMET = false ) ;
 
   void JetEtReSort( vector<TLorentzVector>& vs );
   
   //ClassDef(PseudoExp, 1);

};

//#if !defined(__CINT__)
//    ClassImp(PseudoExp);
#endif

