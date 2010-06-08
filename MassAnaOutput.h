#ifndef MassAnaOutput_H
#define MassAnaOutput_H

#include "TObject.h"
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TString.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TColor.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFormula.h>

#include "MassAna.h"
#include "AlgoZero.h"
#include "AlgoKcon.h"

class MassAnaOutput : public TObject {

private:

   FILE* pfile;
   FILE* efile;

   TString channel;
   TString ch_name;
   TString ptype;
   TString hfolder;

   MassAna*         fitter ;
   //MassFitFunction* fitfunc;
   MassAnaInput*    fitInput;
   AlgoZero*        algo0;
   AlgoKcon*        algok;

public:

   MassAnaOutput();     
   ~MassAnaOutput();     

   void test();

   void CoeffCalib( int rbin, int lowBound, int upBound, Bool_t* comp = NULL);
   void MassCalib( int rbin, int lowBound, int upBound, int NBTag, int NPara, bool isWeight = false ); 

   TString GiveParTitle( int id );

   //ClassDef(MassAnaOutput, 1);

};

//#if !defined(__CINT__)
//    ClassImp(MassAnaOutput);
#endif

