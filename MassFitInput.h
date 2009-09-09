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
#include <TString.h>
#include <TSystem.h>

class MassFitInput : public TObject {

private:

   TString channel;
   TString sname;
   TString hname;
   TString ch_name;
   int n_btag;
   double luminosity ;
   double N_tt;
   double N_wj;
   double N_stt;
   double N_stw;
   double N_qcd;

   // file names for fake data and templates
   TString subtag;
   TString theSG;   // ttbar signal and wrong combinatorics
   TString theBG1;  // W+jets
   TString theBG2;  // single top t-channel
   TString theBG3;  // single top tW-channel
   TString theBG4;  // QCD

   // files for templates
   TFile* file1 ;
   TFile* file2 ;
   TH2D* ttmass;

public:

   MassFitInput( TString channel, int NBTag );     
   ~MassFitInput();     
 
   void Initialize( TString* hfolder );

   void getSignal( TH1D* h_Sg, int rbin, TString mName );
   void getBackground( TH1D* h_Bg, int type, int rbin, TString mName );


   void getFakeData(TString mName, int rbin, TH1D* h_data, THStack* ttstk, TH1D* dth0, TH1D* dth1, TH1D* dth2 = NULL, TH1D* dth3 = NULL,
                                                            TH1D* dth4 = NULL, TH1D* dth5 = NULL );
   void getFakeData(TH1D* h_data, TString thefileName, int rbin, double theScale = 1. );
   void getFakeData(TH1D* h_data, TString thefileName, bool isSignal, int rbin, double theScale = 1. );
 
   void GetFileName( TString mName, int type );
   void NormalizeComponents( double lumi, double nEvents, int channel, TH1D* tmp );
   void combineBG( TString mName, TH1D* allbg, int rbin );

   void getData(TH1D* h_data, TString thefileName, int rbin );

   ClassDef(MassFitInput, 1);

};

#if !defined(__CINT__)
    ClassImp(MassFitInput);
#endif

