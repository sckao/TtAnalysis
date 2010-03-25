#ifndef MassAna_H
#define MassAna_H

#include "TObject.h"
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <THStack.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TString.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TFrame.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TGraphErrors.h>

#include "MassAnaInput.h"
#include "MassFitFunction.h"


class MassAna : public TObject {

private:

   FILE* logfile;
   FILE* parfile;
   FILE* errfile;
   FILE* Bgpara;
   FILE* Sgpara;

   TString fname;
   TString cname;
   TString hfolder;
   TString ptype;

   TString plot1;
   TString plot2;
   TString plot3;
   TString plot8;
   TString plot9;

   TCanvas* c3;
   TCanvas* c7;
   TCanvas* c8;
   TCanvas* c9;

   MassAnaInput*    fitInput;
   MassFitFunction* fitFunc;

   double mL;
   double mH;

public:

   MassAna( TString channel, double massL = 0, double massH = 480 );     
   ~MassAna();     
   
   //void MoreCombinedFitting( TString mName, int rbin, int lowBound, int upBound, Bool_t* comp, int NBTag );
   //void CombinedFitting( TString mName, int rbin, int lowBound, int upBound, int NBTag );
   //void ConstrainFitting( TString mName, int rbin, int lowBound, int upBound, int NBTag );

   double Chi2Test( TString mName, TH1D* theData, int lowBound, int upBound, int nPar, int NBTag = -1, Double_t* statErr=NULL, int rbin = 10, bool isWeight = false );
   double getChi2( TH1D* theData,  TF1* theFunc, TF1* fS, TF1* fB, TF1* fW, double lowBound , double upBound, int nPar );

   void FitSignal1( TString mName, int rbin );
   void FitSignal( TString mName, int rbin, Double_t *para = NULL, Double_t *perr = NULL );
   void FitTtbar( TString mName, int rbin, Double_t *para = NULL, Double_t *perr = NULL );
   void FitTtbar2( TH1D* h1, TH1D* h2, TString mName );
   void FitBackground( TString mName, int rbin, int lowBound, int upBound, Double_t *para = NULL, Double_t *perr = NULL );

   void GetAllCoeff( TString mName, int rbin, int lowBound, int upBound, Bool_t *comp = NULL );

   double MassDigi( TString mString );

   void SetFitParameters( double mass, Double_t* para, int nPara, int NBTag, int nbin );

   //ClassDef(MassAna, 1);

};

//#if !defined(__CINT__)
//    ClassImp(MassAna);
#endif

