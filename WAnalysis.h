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
#include "PseudoExp.h"

class WAnalysis : public TObject {

private:

   double mL;
   double mH;

   double bTh;
   int    n_btag;
   MassAnaInput*    fitInput;
   HadWMassFitter*  wmfitter;
   PseudoExp*       pseudoExp;


   /*
   TCanvas* c1;
   TCanvas* c2;
   TCanvas* c3;
   TCanvas* c4;
   TCanvas* c5;
   TCanvas* c6;
   TCanvas* c7;
   */
   string hfolder;

public:

   WAnalysis( double massL, double massH );     
   ~WAnalysis();     
 
   void TWFitter( string mName, int type,  TString DrawOpt = "COLZ", bool isMCMatched = false );
   void TWFitter1( int njets, string mName, int type, TString DrawOpt = "COLZ", bool isMCMatched = false );

   void EnsembleTest( int type, int randomSeed = 0 ,TString DrawOpt = "COLZ" );
   void EnsembleTest1( int njets, int randomSeed = 0 ,TString DrawOpt = "COLZ" );

   void MixBG( int type, TString DrawOpt = "COLZ" );
   void MixBG1( int njets, int type, TString DrawOpt = "COLZ" );

   void MixAll( int type, TString DrawOpt = "COLZ" );
   void MixAll1( int njets, TString DrawOpt = "COLZ" );

   void M2M3Plotter( vector<TH2D*> h2Ds, string fileName, TString DrawOpt = "COLZ", bool isMCMatched = false );

   void BJetEff( string fileName ) ;

   //ClassDef(WAnalysis, 1);

};

//#if !defined(__CINT__)
//    ClassImp(WAnalysis);
#endif

