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
#include "ObjectInfo.h"
#include "XSection.h"
#include "BgEstimation.h"

using namespace std; 

int main() { 

  MassAnaInput   *massInput = new MassAnaInput();
  //MassAna        *theFitter = new MassAna();
  //MassAnaOutput  *FitResult = new MassAnaOutput();
  PseudoExp      *pExp      = new PseudoExp();
  WAnalysis      *WFitter   = new WAnalysis();
  ObjectInfo     *h_objs    = new ObjectInfo();
  XSection       *xSec      = new XSection();
  BgEstimation   *bgEst     = new BgEstimation();

  //theFitter->FitSignal1("171",10) ;
  
  TString DrawOpt = "COLZ";

  vector<int> module;
  massInput->GetParameters( "Module", &module );

  string mcMatching ;  // seed = 0 -> using system time as seed
  massInput->GetParameters( "MCMatching", &mcMatching );

  int randomSeed = 0 ;  // seed = 0 -> using system time as seed
  massInput->GetParameters( "RandomSeed", &randomSeed );

  vector<string> RealData ;
  massInput->GetParameters( "TheData", &RealData );

  vector<string> theFiles ;
  massInput->GetParameters( "FakeData", &theFiles );

  vector<string> the2JFiles ;
  massInput->GetParameters( "2JSamples", &the2JFiles );

  vector<string> the4JFiles ;
  massInput->GetParameters( "4JSamples", &the4JFiles );

  int nPseudoExp = 0 ;  // seed = 0 -> using system time as seed
  massInput->GetParameters( "nPseudoExp", &nPseudoExp );

  string phaseSmear ;  // seed = 0 -> using system time as seed
  massInput->GetParameters( "PhaseSmear", &phaseSmear );

  massInput->LinkForests( "muJets" );

  if ( module[0] > 0) {
     bool smearing =  ( phaseSmear == "ON" ) ? true : false ;
     if ( module[0] == 2 ) h_objs->MCPlotter1( theFiles, 1000 );
     if ( module[0] == 3 ) h_objs->CombinedMCPlotter( theFiles, false );
     
     if ( module[0] == 4 ) {
        
        bool doQCD = true ;
        h_objs->DataPlotter( RealData[0], theFiles, false );
        h_objs->QCDBGPlotter( RealData[0], theFiles,  doQCD, 35, 20, 1. );
        h_objs->QCDBGPlotter( RealData[0], theFiles, !doQCD, 50, 30, 1. );
        
        //h_objs->DataPlotter( RealData[0], theFiles, true );
     }
     if ( module[0] == 5 ) {
        
        //bgEst->MtFitter( RealData[0], theFiles ) ;
        bool doQCD = true ;
        bool isVjNorm = false ;
        /*
        bgEst->MeasureScale1D( RealData[0], theFiles, 13, 35, 20, doQCD, isVjNorm );
        bgEst->MeasureScale1D( RealData[0], theFiles, 14, 35, 20, doQCD );
        bgEst->MeasureScale1D( RealData[0], theFiles, 15, 35, 20, doQCD );
        bgEst->MeasureScale1D( RealData[0], theFiles, 16, 35, 20, doQCD );
        bgEst->MeasureScale1D( RealData[0], theFiles, 17, 35, 20, doQCD );
        */
        bgEst->MeasureScale1D( RealData[0], theFiles, 13, 50, 30, !doQCD, isVjNorm );
        /*
        bgEst->MeasureScale1D( RealData[0], theFiles, 14, 50, 30, !doQCD );
        bgEst->MeasureScale1D( RealData[0], theFiles, 15, 50, 30, !doQCD );
        bgEst->MeasureScale1D( RealData[0], theFiles, 16, 50, 30, !doQCD );
        bgEst->MeasureScale1D( RealData[0], theFiles, 17, 50, 30, !doQCD );
        */

        // min chi2 method
        //bgEst->MeasureScale2D( RealData[0], theFiles, 35, 20, true );
        //bgEst->MeasureScale2D( RealData[0], theFiles, 50, 30, false );
        //bgEst->MeasureScale2D( RealData[0], theFiles, 120, 120, true, false );
        // normalize by total number of event
        //bgEst->MeasureScale( RealData[0], theFiles );
     }
     if ( module[0] == 6 ) {
        bool doPlots = true ;
        bool normMC  = false ;
        //bool withTt = false ;
        //bgEst->RatioXY( the4JFiles, the2JFiles, -1, withTt, doPlots, 20 ) ;
        
        //bgEst->CombinedRatio( 2, 1 ) ;
        //bgEst->CombinedRatio( 4, 1, true, false ) ;
        bgEst->RatioXY( 2, 1, theFiles,  1, doPlots, false, false, normMC ) ;
        //bgEst->RatioXY( 2, 1, theFiles,  2, doPlots, false, false, normMC ) ;
        bgEst->RatioXY( 2, 1, theFiles,  3, doPlots, false, false, normMC ) ;
        //bgEst->RatioXY( 2, 1, theFiles,  4, doPlots, false, false, normMC ) ;
        //bgEst->RatioXY( 2, 1, RealData,  3, doPlots, false, false, normMC ) ;
        //bgEst->RatioXY( 2, 1, RealData,  4, doPlots, false, false, normMC ) ;
          
        //bgEst->RatioXY( 1, 2, theFiles, -1, doPlots ) ;
        //bgEst->RatioXY( 2, 1, RealData, -1, doPlots ) ;
        //bgEst->RatioXY( 4, 1, RealData, -1, doPlots ) ;
        
     }
     
     if ( module[0] == 1 ) {
     /*
        h_objs->ObjHistoPlotter( the4JFiles[0], smearing );
	h_objs->ObjHistoPlotter( the4JFiles[1], smearing );
	h_objs->ObjHistoPlotter( the4JFiles[2], smearing );
	h_objs->ObjHistoPlotter( the4JFiles[3], smearing );
	h_objs->ObjHistoPlotter( the4JFiles[4], smearing );
	h_objs->ObjHistoPlotter( the4JFiles[5], smearing );
     */
	h_objs->ObjHistoPlotter( the2JFiles[0], smearing );
	h_objs->ObjHistoPlotter( the2JFiles[1], smearing );
	h_objs->ObjHistoPlotter( the2JFiles[2], smearing );
	h_objs->ObjHistoPlotter( the2JFiles[3], smearing );
	h_objs->ObjHistoPlotter( the2JFiles[4], smearing );
	h_objs->ObjHistoPlotter( the2JFiles[5], smearing );
        //h_objs->JacobTester();
     }

  }

  if ( module[1] == 1 ) {
     if ( mcMatching == "ON" )   WFitter->HadTopFitter( the4JFiles[0],  DrawOpt, true );
     
     WFitter->HadTopFitter( RealData[0],  DrawOpt );
     WFitter->HadTopFitter( the4JFiles[0],  DrawOpt );
     WFitter->HadTopFitter( the4JFiles[1],  DrawOpt );
     WFitter->HadTopFitter( the4JFiles[2],  DrawOpt );
     WFitter->HadTopFitter( the4JFiles[3],  DrawOpt );
     WFitter->HadTopFitter( the4JFiles[4],  DrawOpt );
     WFitter->HadTopFitter( the4JFiles[5],  DrawOpt );

     /*
     WFitter->HadTopFitter( the2JFiles[0],  DrawOpt );
     WFitter->HadTopFitter( the2JFiles[1],  DrawOpt );
     WFitter->HadTopFitter( the2JFiles[2],  DrawOpt );
     WFitter->HadTopFitter( the2JFiles[3],  DrawOpt );
     WFitter->HadTopFitter( the2JFiles[4],  DrawOpt );
     WFitter->HadTopFitter( the2JFiles[5],  DrawOpt );
     */
  }

  if ( module[2] == 1 ) {
     if ( mcMatching == "ON" )   WFitter->LepTopFitter( theFiles[0],  DrawOpt, true );
     WFitter->LepTopFitter( theFiles[0],  DrawOpt );
     WFitter->LepTopFitter( theFiles[1],  DrawOpt );
     WFitter->LepTopFitter( theFiles[2],  DrawOpt );
     WFitter->LepTopFitter( theFiles[3],  DrawOpt );
  }

  if ( module[3] == 1 ) {
     WFitter->MixAll( the4JFiles, DrawOpt );
     //WFitter->MixBG(  DrawOpt );
     //WFitter->EnsembleTest( randomSeed , DrawOpt );
  }

  if ( module[4] == 1)  {
     WFitter->Had_SBRatio();
  }
  if ( module[5] == 1)  {
     //xSec->BgClosureTest(4,1,true,false);
     xSec->BgClosureTest(2,1,false,false);
     //xSec->RealBackground(2,1,false,false);
     xSec->RealBackground(4,1,true,false);
     //xSec->MethodTest1() ;
  }
  if ( module[6] == 1)  {
     xSec->EnsembleTest1( nPseudoExp , randomSeed ) ;
  }

  return 0;
}

