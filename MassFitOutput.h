#include "TObject.h"
#include "TemplateMassFit.h"
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

class MassFitOutput : public TObject {

private:

   FILE* pfile;
   FILE* efile;

   TString channel;
   TString ch_name;
   TString ptype;
   TString hfolder;

   TemplateMassFit* fitter ;
   MassFitFunction* fitfunc;
   MassFitInput*    fitInput;

public:

   MassFitOutput( TString channel, int NBTag );     
   ~MassFitOutput();     

   void CoeffCalib( int rbin, int lowBound, int upBound, Bool_t* comp );
   void MassCalib( int rbin, int lowBound, int upBound,  Bool_t* comp, int NBTag, int NPara ); 

   TString GiveParTitle( int id );

   ClassDef(MassFitOutput, 1);

};

#if !defined(__CINT__)
    ClassImp(MassFitOutput);
#endif

