#ifndef WAnalysis_H
#define WAnalysis_H

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

#include "MassAnaInput.h"
#include "MassAna.h"
#include "HadWMassFitter.h"
#include "LepTopMassFitter.h"
#include "PseudoExp.h"

class WAnalysis : public TObject {

private:

   double bTh;
   int    n_btag;
   int    n_Jets;

   MassAnaInput*      fitInput;
   HadWMassFitter*    wmfitter;
   LepTopMassFitter*  ltmfitter;
   PseudoExp*         pseudoExp;

   string hfolder;

public:

   WAnalysis( double massL, double massH );     
   ~WAnalysis();     
 
   void HadTopFitter( string mName, TString DrawOpt = "COLZ", bool isMCMatched = false );

   void LepTopFitter( string mName, TString DrawOpt = "COLZ", bool isMCMatched = false );

   void EnsembleTest( int randomSeed = 0 ,TString DrawOpt = "COLZ" );
   void LepTEnsembleTest( int randomSeed = 0 ,TString DrawOpt = "COLZ" );

   void MixBG( TString DrawOpt = "COLZ" );

   void MixAll( TString DrawOpt = "COLZ" );

   void M2M3Plotter( vector<TH2D*> h2Ds, string fileName, TString DrawOpt = "COLZ", bool isMCMatched = false );
   void LepTopPlotter( vector<TH2D*> h2Ds, string fileName, TString DrawOpt = "COLZ", bool isMCMatched = false );

   void BJetEff( string fileName ) ;

   //ClassDef(WAnalysis, 1);

};

//#if !defined(__CINT__)
//    ClassImp(WAnalysis);
#endif

