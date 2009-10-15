#ifndef MassAnaInput_H
#define MassAnaInput_H

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

class MassAnaInput : public TObject {

private:

   TString channel;
   TString hname;
   TString ch_name;
   TString probName;

   int n_btag;
   double bTh;
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

   double mL;
   double mH;
   bool   weighting;

public:

   MassAnaInput(){ hname ="hadTM", ch_name ="had", n_btag = -1, mL =0, mH = 400; }     
   MassAnaInput( TString channel, int NBTag, double massL = 0 , double massH = 480 );     
   ~MassAnaInput();     
 
   void Initialize( TString* hfolder );
   void GetFileName( TString mName, int type, TString* fNameList = NULL );

   void get_h1Obj(TString fName, TString TName, TString BName, TH1D* h1, double theScale = 1. );
   void getMcMatching(TString mName, TString BName, TH1D* h1, double theScale = 1. );

   void getHadPermutation(TString fName, TString BName, TH1D* h1, double theScale = 1., std::vector<bool>* blist = NULL );
   void getLepPermutation(TString fName, TString BName, TH1D* h1, double theScale = 1., std::vector<bool>* blist = NULL );

   void ListWithBTag(TString fName, int nbtags, std::vector<bool>* blist, std::vector<double>* bProb = NULL );

   void getMostProb(TString mName, TString BName, TH1D* h1, double theScale = 1. );

   //void getFakeData(TString mName, int rbin, TH1D* h_data, THStack* ttstk, TH1D* dth0, TH1D* dth1, 
   //                 TH1D* dth2 = NULL, TH1D* dth3 = NULL, TH1D* dth4 = NULL, TH1D* dth5 = NULL );
   //void getFakeData(TH1D* h_data, TString thefileName, bool isSignal, int rbin, double theScale = 1. );

   void getSignal( TH1D* h_Sg, int type, TString mName );
   void getBackground( TH1D* h_Bg, int type, int nbin, TString mName );

   void NormalizeComponents( double lumi, double nEvents, int channel, TH1D* tmp );
   /*
 
   void combineBG( TString mName, TH1D* allbg, int rbin );

   void getData(TH1D* h_data, TString thefileName, int rbin );
   */


   ClassDef(MassAnaInput, 1);

};

//#if !defined(__CINT__)
//    ClassImp(MassAnaInput);
#endif

