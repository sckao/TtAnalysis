#ifndef BgEstimation_H
#define BgEstimation_H

#include <TObject.h>
#include <TMinuit.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>

#include "MassFitFunction.h"
#include "MassAnaInput.h"
#include "HadWMassFitter.h"
#include "PseudoExp.h"

//#include "LepTopMassFitter.h"

class BgEstimation : public TObject {

private:

   MassFitFunction*   fitFunc;
   MassAnaInput*      fitInput;
   HadWMassFitter*    wmfitter;
   PseudoExp*         pseudoExp;

   string hfolder;
   TString theFolder;
   bool smearing ;
   vector<double> inputMean ;

public:

   BgEstimation();     
   ~BgEstimation();     
 
   vector<double> Ratio42( int statIdx = 0 );

   void MethodTest() ;

   vector<double> XSection( vector<double>& R42, double n2J, double n4J ) ;

   void EnsembleTest( int randomSeed );
   void EnsembleTest( int nRun, int randomSeed );

   vector<double> RunEnsembles( int treeSize, string fileName, double pMean, int nRun, int randomSeed, TTree* theTree = NULL, TH1D* hGen = NULL );

   vector<double> StatErr( double m );
   vector<double> ErrAovB( double A, double s_A, double B, double s_B );
   vector<double> ErrAxB( double A, double s_A, double B, double s_B );

   //ClassDef(BgEstimation, 1);

};

//#if !defined(__CINT__)
//    ClassImp(BgEstimation);
#endif

