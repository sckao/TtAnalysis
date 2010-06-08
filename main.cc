#include <iostream> 
#include <vector>
#include <fstream>
#include <string>
#include <algorithm>
#include <stdio.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TMinuit.h>

#include "WFormat.h"
#include "MassFitFunction.h"
#include "MassAnaInput.h"
#include "MassAna.h"
#include "AlgoZero.h"
#include "AlgoKcon.h"
#include "JES.h"
#include "MassAnaOutput.h"
#include "HadWMassFitter.h"
#include "LepTopMassFitter.h"
#include "PseudoExp.h"
#include "WAnalysis.h"
#include "JetSpectrum.h"
#include "BgEstimation.h"

using namespace std; 

int main() { 

  MassAnaInput   *massInput = new MassAnaInput();
  //MassAna        *theFitter = new MassAna();
  //MassAnaOutput  *FitResult = new MassAnaOutput();
  PseudoExp      *pExp      = new PseudoExp();
  WAnalysis      *WFitter   = new WAnalysis();
  JetSpectrum    *h_objs    = new JetSpectrum();
  BgEstimation   *bgEst     = new BgEstimation();

  //theFitter->FitSignal1("171",10) ;
  
  TString DrawOpt = "COLZ";

  vector<int> module;
  massInput->GetParameters( "Module", &module );

  string mcMatching ;  // seed = 0 -> using system time as seed
  massInput->GetParameters( "MCMatching", &mcMatching );

  int randomSeed = 0 ;  // seed = 0 -> using system time as seed
  massInput->GetParameters( "RandomSeed", &randomSeed );

  vector<string> theFiles ;
  massInput->GetParameters( "FakeData", &theFiles );

  vector<string> the2JFiles ;
  massInput->GetParameters( "2JSamples", &the2JFiles );

  int nPseudoExp = 0 ;  // seed = 0 -> using system time as seed
  massInput->GetParameters( "nPseudoExp", &nPseudoExp );

  string phaseSmear ;  // seed = 0 -> using system time as seed
  massInput->GetParameters( "PhaseSmear", &phaseSmear );

  if ( module[0] == 1) {
     bool smearing =  ( phaseSmear == "ON" ) ? true : false ;
     h_objs->ObjHistoPlotter( theFiles[0], smearing );
     h_objs->ObjHistoPlotter( theFiles[1], smearing );
     h_objs->ObjHistoPlotter( theFiles[2], smearing );
     h_objs->ObjHistoPlotter( the2JFiles[0], smearing );
     h_objs->ObjHistoPlotter( the2JFiles[1], smearing );
     h_objs->ObjHistoPlotter( the2JFiles[2], smearing );
     //h_objs->ObjHistoPlotter( theFiles[3] );
  }

  if ( module[1] == 1 ) {
     if ( mcMatching == "ON" )   WFitter->HadTopFitter( theFiles[0],  DrawOpt, true );
     WFitter->HadTopFitter( theFiles[0],  DrawOpt );
     WFitter->HadTopFitter( theFiles[1],  DrawOpt );
     WFitter->HadTopFitter( theFiles[2],  DrawOpt );
     WFitter->HadTopFitter( the2JFiles[0],  DrawOpt );
     WFitter->HadTopFitter( the2JFiles[1],  DrawOpt );
     WFitter->HadTopFitter( the2JFiles[2],  DrawOpt );
     //WFitter->HadTopFitter( theFiles[3],  DrawOpt );
  }

  if ( module[2] == 1 ) {
     if ( mcMatching == "ON" )   WFitter->LepTopFitter( theFiles[0],  DrawOpt, true );
     WFitter->LepTopFitter( theFiles[0],  DrawOpt );
     WFitter->LepTopFitter( theFiles[1],  DrawOpt );
     WFitter->LepTopFitter( theFiles[2],  DrawOpt );
     //WFitter->LepTopFitter( theFiles[3],  DrawOpt );
  }

  if ( module[3] == 1 ) {
     WFitter->MixBG(  DrawOpt );
     WFitter->MixAll( DrawOpt );
     WFitter->EnsembleTest( randomSeed , DrawOpt );
  }

  if ( module[4] == 1)  {
     WFitter->Had_SBRatio();
  }
  if ( module[5] == 1)  {
     bgEst->MethodTest() ;
  }
  if ( module[6] == 1)  {
     bgEst->EnsembleTest( nPseudoExp , randomSeed ) ;
  }

  return 0;

}
