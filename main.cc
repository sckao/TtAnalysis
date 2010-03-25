#include <iostream> 
#include <vector>
#include <stdio.h>
#include <fstream>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <algorithm>
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
#include "PseudoExp.h"
#include "WAnalysis.h"

using namespace std; 

int main() { 


  PseudoExp      *pExp      = new PseudoExp( 0., 480. );
  MassAna        *theFitter = new MassAna( "had", 0, 480 );
  MassAnaOutput  *FitResult = new MassAnaOutput( "had", 0., 480. );
  WAnalysis      *WFitter   = new WAnalysis( 0., 480. );

  
  int njets = 3 ;
  int randomSeed = 0 ;
  bool mcmatched =  true ;
  TString DrawOpt = "COLZ";

  WFitter->TWFitter1( njets, "tt171_336a", 0, DrawOpt, mcmatched );
  WFitter->TWFitter1( njets, "tt171_336a", 0, DrawOpt );
  WFitter->TWFitter1( njets, "wj_336a", 0, DrawOpt );
  WFitter->TWFitter1( njets, "qcd1_336a", 0, DrawOpt );
  WFitter->MixBG1( njets, 0, DrawOpt );
  WFitter->MixAll1( njets, DrawOpt );
  WFitter->EnsembleTest1( njets, randomSeed , DrawOpt );

  //pExp->GetEnsemble("tt171_336a", 10, 1);

  //theFitter->FitSignal1("171",10) ;

  return 0;

}
