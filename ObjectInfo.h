#ifndef ObjectInfo_H
#define ObjectInfo_H

#include <TObject.h>
#include <TMinuit.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom2.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TGraphErrors.h>
#include <TSystem.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>

#include "MassAnaInput.h"
#include "WFormat.h"
#include "PseudoExp.h"
#include "MassFitFunction.h"

class ObjectInfo : public TObject {

private:

   // for hadronic permutation

   MassAnaInput*    fitInput;
   PseudoExp*       pseudoExp;
   MassFitFunction* fitFunc;

   string hfolder;
   vector<double> jetCuts;
   vector<double> muonCuts;
   int n_Jets ;
   string plotType ;
   string Inclusive ;
   bool inclu ;
   bool isWeight ;

   double bgMtCut;
   double bgMETCut;

public:

   ObjectInfo();     
   ~ObjectInfo();     

   void EvtSelector( string fileName, recoObj* histos, bool smearing = false, double scale = 1., vector<int>* evtlist = NULL, int evtSplit = 0 );

   double MtFitter( TH1D* hs, TF1& fjb1, TCanvas* c1a = NULL, double minMt = 40., double maxMt = 120., bool doPlot = false );

   void ObjHistoPlotter( string fileName, bool smearing = false, bool doWMtfit = false );

   void DataPlotter( string DataName, vector<string>& fakeData, bool doScale = false );
   void Data2DPlotter( string DataName, vector<string>& fakeData, bool doScale = false );

   void JacobTester();

   void CombinedMCPlotter( vector<string>& fakeData, bool doWMtFit = false, bool doScale = false );

   void MCPlotter1( vector<string>& fakeData, double norm = 1000. );

   vector<double> BinErr( double m ) ;

   void Reset(int idx, bool newValue );
   void Reset(int idx, int newValue );
   void Reset(int idx, double newValue );
   void ResetBGCuts(int idx, double newValue );

   double EvtScaling( double w_pt , string fileName ) ;
   double EvtScaling( int NJets , string fileName ) ;

   void QCDSelector( string fileName, recoObj* histos, bool smearing = false, double scale = 1., bool doQCD = true, int mode = 0, int evtSplit = 0 );
   double QCDScaling( double w_pt ) ;
   double QCDScaling( int NJets, string fileName ) ;
   void QCDBGPlotter( string DataName, vector<string>& fakeData, bool doQCD, double MtCut, double METCut, double scale = 1 );
   void QCDBG2DPlotter( string DataName, vector<string>& fakeData, bool doQCD, double MtCut, double METCut, double scale = 1 );

   void WScaleStudy( vector<string>& fakeData ) ;

   double PtRel (TLorentzVector v1, TLorentzVector v2 );
   //ClassDef(ObjectInfo, 1);

};

//#if !defined(__CINT__)
//    ClassImp(ObjectInfo);
#endif

