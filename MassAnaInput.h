#ifndef MassAnaInput_H
#define MassAnaInput_H

#include "TObject.h"
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <math.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <THStack.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TString.h>
#include <TSystem.h>
#include <TLorentzVector.h>
#include "WFormat.h"

struct jlist {
    int w1 ;
    int w2 ;
    int bh ;
    int bl ;
};

class MassAnaInput : public TObject {

private:

   vector<TTree*> MuJForest ;
   vector<TTree*> mcTtForest ;

   TString channel;
   TString hname;
   TString ch_name;
   TString probName;

   int n_btag;
   double bTh;

   int neu_str;
   int jet_str;
   int mu_str;

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
   MassAnaInput( TString channel, double massL = 0 , double massH = 480 );     
   ~MassAnaInput();     
 
   void Initialize( TString* hfolder );

   vector<TTree*> GetForest( string DataSet, TString treeName );

   TTree* GetTree( string chName, TString treeName, TFile* file = NULL );

   void get_h1Obj(TString fName, TString TName, TString BName, TH1D* h1, double theScale = 1., bool weight = false );
   void get_h1Obj(TChain* tChain, TString BName, TH1D* h1, double theScale = 1., bool weight = false );

   void getMcMatching(TString mName, TString BName, TH1D* h1, double theScale = 1. );

   void getHadPermutation(TString fName, TString BName, TH1D* h1, double theScale = 1., std::vector<bool>* blist = NULL );
   void getLepPermutation(TString fName, TString BName, TH1D* h1, double theScale = 1., std::vector<bool>* blist = NULL );

   void ListWithBTag(TString fName, int nbtags, std::vector<bool>* blist, std::vector<double>* bProb = NULL );

   void getMostProb(TString mName, TString BName, TH1D* h1, double theScale = 1. );

   //void getFakeData(TString mName, int rbin, TH1D* h_data, THStack* ttstk, TH1D* dth0, TH1D* dth1, 
   //                 TH1D* dth2 = NULL, TH1D* dth3 = NULL, TH1D* dth4 = NULL, TH1D* dth5 = NULL );
   //void getFakeData(TH1D* h_data, TString thefileName, bool isSignal, int rbin, double theScale = 1. );

   void getTt( TH1D* h_Sg, double tmass, bool matched = true );
   void getBackground( TH1D* allBg, THStack* bgstk, vector<TH1D*>& hlist, int groupId = 0);

   void getFakeData( int rbin, TH1D* ttadd, THStack* ttstk, vector<TH1D*>& hlist );

   vector<TLorentzVector> GetNeutrinos( TTree* tr, int evtIdx, int& synId, bool nextEvt = true );
   vector<TLorentzVector> GetJets( TTree* tr, int evtIdx, vector<double>& bCuts, int& synId, bool nextEvt = true );
   vector<TLorentzVector> GetMuons( TTree* tr, int evtIdx, int& synId, bool nextEvt = true );

   void GetPermutes( int njets, vector<jlist>& jlistV );

   void NormalizeComponents( double lumi, double nEvents, int channel, TH1D* tmp );
   void NormalizeComponents( string theChhannel, TH1D* tmp );
   double NormalizeComponents( string theChhannel );
   
   void GetParameters( string paraName, int* thePara );
   void GetParameters( string paraName, double* thePara );
   void GetParameters( string paraName, string* thePara );
   void GetParameters( string paraName, vector<double>* thePara );
   void GetParameters( string paraName, vector<string>* thePara );
   void GetParameters( string paraName, vector<int>* thePara );

   /*
   void combineBG( TString mName, TH1D* allbg, int rbin );
   void getData(TH1D* h_data, TString thefileName, int rbin );
   */

   //ClassDef(MassAnaInput, 1);

};

//#if !defined(__CINT__)
//    ClassImp(MassAnaInput);
#endif

