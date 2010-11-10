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
#include <TGaxis.h>

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

   bool isFolderCreate ;
   bool smearing ;
   vector<double> inputMean ;

   string plotType;
   string phaseSmear;
   string hfolder;

   MassAnaInput*      fitInput;
   HadWMassFitter*    wmfitter;
   LepTopMassFitter*  ltmfitter;
   PseudoExp*         pseudoExp;


public:

   WAnalysis();     
   ~WAnalysis();     

   void CreateFolders(); 

   void HadTopFitter( string mName, TString DrawOpt = "COLZ", bool isMCMatched = false );

   void LepTopFitter( string mName, TString DrawOpt = "COLZ", bool isMCMatched = false );

   void Had_SBRatio();
   void SBCEPlotter();

   void EnsembleTest( int randomSeed = 0 ,TString DrawOpt = "COLZ" );
   void LepTEnsembleTest( int randomSeed = 0 ,TString DrawOpt = "COLZ" );

   void MixBG( TString DrawOpt = "COLZ" );

   void MixAll( vector<string>& flist, TString DrawOpt = "COLZ" );

   void M2M3Plotter( vector<TH2D*> h2Ds, string fileName, TString DrawOpt = "COLZ", bool isMCMatched = false );
   void AN_M2M3Plotter( vector<TH2D*> h2Ds, string fileName, TString DrawOpt = "COLZ", bool isMCMatched = false );
   void LepTopPlotter( vector<TH2D*> h2Ds, string fileName, TString DrawOpt = "COLZ", bool isMCMatched = false );

   void BJetEff( string fileName ) ;

   //ClassDef(WAnalysis, 1);

};

//#if !defined(__CINT__)
//    ClassImp(WAnalysis);
#endif

