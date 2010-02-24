#ifndef HadWMassFitter_H
#define HadWMassFitter_H

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
#include "MassFitFunction.h"

class HadWMassFitter : public TObject {

private:

   double mL;
   double mH;

   double sigma1;
   double sigma2;
   double zigma1;
   double zigma2;

   double wmass;

   double bTh;
   int    n_btag;

   // for hadronic permutation
   int eventId ;
   std::vector<int> jd1;
   std::vector<int> jd2;
   std::vector<int> jd3;

   MassAnaInput*    fitInput;
   MassAna*         fitTools;
   MassFitFunction* fitFunc;

   TCanvas* c1;
   TCanvas* c2;

   string hfolder;

public:

   HadWMassFitter( double massL, double massH );     
   ~HadWMassFitter();     
 
   void FitW( TLorentzVector jetv1, TLorentzVector jetv2, Double_t* pars, Double_t* errs, bool isJES );

   static void WMFCN(Int_t &npar, Double_t *, Double_t &wChi2, Double_t *par, Int_t iflag);

   double ReFitWMass( TLorentzVector v1, TLorentzVector v2, Double_t *par, Double_t *ProbW );

   double ReFitHadTopMass( TLorentzVector v1, TLorentzVector v2, TLorentzVector v3, Double_t *par );
   double ReFitLepTopMass( TLorentzVector v1, TLorentzVector v2, TLorentzVector v3, Double_t *par );

   Double_t jetAngle( TLorentzVector v1, TLorentzVector v2 );

   void ReFitSolution( string mName, TH2D* hM2M3, TH2D* hM3M3, int type, bool isMCMatched = false );

   void FitMCSolution( string mName, TH2D* hM2M3, int type );

   bool HadPermutation( int i1, int i2, int i3, int evtId );

   bool JetPtFilter( vector<TLorentzVector> jpv );

   void TWFitter( string mName, int type, TString DrawOpt = "COLZ", bool isMCMatched = false );

   //TLorentzVector NeuLooper( TTree* tr, int evtIdx, int objIdx );

   ClassDef(HadWMassFitter, 1);

};

class FitVectors : public TObject {

public:

   FitVectors(){}
   ~FitVectors(){}

   TLorentzVector jv1;

   TLorentzVector jv2;

   ClassDef(FitVectors, 1);
};

//#if !defined(__CINT__)
//    ClassImp(HadWMassFitter);
#endif

