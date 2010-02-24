#ifndef JES_H
#define JES_H

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
#include "MassFitFunction.h"

class JES : public TObject {

private:

   double mL;
   double mH;
   double jetscale;

   // for hadronic permutation

   MassAnaInput*    fitInput;
   MassFitFunction* fitFunc;

   TCanvas* c1;
   TCanvas* c2;

   TString hfolder;

public:

   JES( double massL, double massH );     
   ~JES();     
 
   double WMSC( TLorentzVector v1, TLorentzVector v2, double jes = 1. );
   double deltaR(TLorentzVector v1, TLorentzVector v2 );

   void WMassSpectrum( TString mName, TH1D* hWM0,  TH1D* hWM, TH1D* hMCW, double jes = 1. );

   void GetJESFromW( TString mName, Double_t* para = NULL, double jes = 1. );
 
   void CalibJES( TString mName );

   void Matching( TString mName, double jes = 1.0, TH1D* hMatchW = NULL );
   void Smearing( TString mName, double sgm = 0.7, TH1D* hSmearW = NULL );
   void SmearAndMatch( TString mName );

   ClassDef(JES, 1);

};

//#if !defined(__CINT__)
//    ClassImp(JES);
#endif

