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

  MassAnaInput   *massInput = new MassAnaInput( "had", 0, 480 );
  //MassAna        *theFitter = new MassAna( "had", 0, 480 );
  //MassAnaOutput  *FitResult = new MassAnaOutput( "had", 0., 480. );
  PseudoExp      *pExp      = new PseudoExp( 0., 480. );
  WAnalysis      *WFitter   = new WAnalysis( 0., 480. );
  JetSpectrum    *h_objs    = new JetSpectrum(0., 480. );
  BgEstimation   *bgEst     = new BgEstimation(0., 480. );

  //theFitter->FitSignal1("171",10) ;
  
  TString DrawOpt = "COLZ";

  int module = 0;
  massInput->GetParameters( "Module", &module );

  string mcMatching ;  // seed = 0 -> using system time as seed
  massInput->GetParameters( "MCMatching", &mcMatching );

  int randomSeed = 0 ;  // seed = 0 -> using system time as seed
  massInput->GetParameters( "RandomSeed", &randomSeed );

  vector<string> theFiles ;
  massInput->GetParameters( "FakeData", &theFiles );

  int nPseudoExp = 0 ;  // seed = 0 -> using system time as seed
  massInput->GetParameters( "nPseudoExp", &nPseudoExp );

  if ( module == 0 || module == 1 ) {
     if ( mcMatching == "ON" )   WFitter->HadTopFitter( theFiles[0],  DrawOpt, true );
     WFitter->HadTopFitter( theFiles[0],  DrawOpt );
     WFitter->HadTopFitter( theFiles[1],  DrawOpt );
     WFitter->HadTopFitter( theFiles[2],  DrawOpt );
     //WFitter->HadTopFitter( theFiles[3],  DrawOpt );
  }

  if ( module == 0 || module == 2 ) {
     if ( mcMatching == "ON" )   WFitter->LepTopFitter( theFiles[0],  DrawOpt, true );
     WFitter->LepTopFitter( theFiles[0],  DrawOpt );
     WFitter->LepTopFitter( theFiles[1],  DrawOpt );
     WFitter->LepTopFitter( theFiles[2],  DrawOpt );
     //WFitter->LepTopFitter( theFiles[3],  DrawOpt );
  }

  if ( module == 0 || module == 3 ) {
     WFitter->MixBG(  DrawOpt );
     WFitter->MixAll( DrawOpt );
     WFitter->EnsembleTest( randomSeed , DrawOpt );
  }

  if ( module == 0 || module == 4 ) {
     h_objs->ObjHistoPlotter( theFiles[0] );
     h_objs->ObjHistoPlotter( theFiles[1] );
     h_objs->ObjHistoPlotter( theFiles[2] );
     //h_objs->ObjHistoPlotter( theFiles[3] );
  }

  if ( module == 0 || module == 5 ) {
     //bgEst->MethodTest() ;
     bgEst->EnsembleTest( nPseudoExp , randomSeed ) ;
  }

  return 0;

}
