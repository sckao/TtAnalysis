#ifndef AlgoZero_H
#define AlgoZero_H

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

#include "MassAna.h"


class AlgoZero : public TObject {

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

   TCanvas* c7;
   TCanvas* c3;
   TCanvas* c8;
   TCanvas* c9;

   MassAna*         fitTools;
   MassAnaInput*    fitInput;

   double mL;
   double mH;

public:

   AlgoZero();     
   ~AlgoZero();     
 
   void MoreCombinedFitting( TString mName, int rbin, int lowBound, int upBound, int NBTag );

   void CombinedFitting( TString mName, int rbin, int lowBound, int upBound, int NBTag );
 
   //void getFakeData( TString mName,  TH1D* ttadd, THStack* ttstk, TH1D* dth0, TH1D* dth1, 
   //                  TH1D* dth2 = NULL, TH1D* dth3 = NULL, TH1D* dth4 = NULL, TH1D* dth5 = NULL );

   void SetFitParameters( double mass, Double_t* para, int nPara, int NBTag, int nbin );

   //ClassDef(AlgoZero, 1);

};

#endif

