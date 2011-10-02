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
#include "SystInfo.h"

using namespace std; 

int main() { 

  MassAnaInput   *massInput = new MassAnaInput();
  WAnalysis      *WFitter   = new WAnalysis();
  ObjectInfo     *h_objs    = new ObjectInfo();
  XSection       *xSec      = new XSection();
  BgEstimation   *bgEst     = new BgEstimation();
  SystInfo       *sysInfo   = new SystInfo();
  //MassAna        *theFitter = new MassAna();
  //MassAnaOutput  *FitResult = new MassAnaOutput();
  //PseudoExp      *pExp      = new PseudoExp();

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
     if ( module[0] == 2 ) h_objs->MCPlotter1( theFiles, 1500 );
     if ( module[0] == 3 ) h_objs->CombinedMCPlotter( theFiles, false, false );
     if ( module[0] == 4 ) {
        bool doScale = false ;
        h_objs->DataPlotter(  RealData[0], theFiles, doScale );
        h_objs->Data2DPlotter( RealData[0], theFiles, doScale );
     }
     if ( module[0] == 5 ) {
        bool doScale = true ;
        h_objs->DataPlotter( RealData[0], theFiles, doScale );
        //h_objs->Data2DPlotter( RealData[0], theFiles, doScale );
     }
     if ( module[0] == 14 ) {
        bool doQCD   = true ;
        // scaleMode = 0 : no scaling , scaleMode = 1 : only scale Vjets componenet , scaleMode = 2 : regular scaling
        int scaleMode = 2 ; 
        h_objs->QCDBGPlotter( RealData[0], theFiles,  doQCD, 35, 20, scaleMode );
        h_objs->QCDBGPlotter( RealData[0], theFiles, !doQCD, 50, 30, scaleMode );
     }
     // only use the scale-up value for JESType here, JESType =2,4,6
     if ( module[0] == 15 ) sysInfo->JESPlotter( theFiles );
     if ( module[0] == 16 ) sysInfo->SystPlotter();

     if ( module[0] == 6 ) {
        //bgEst->MtFitter( RealData[0], theFiles ) ;
        bool doQCD = true ;
        bool isVjNorm = false ;
        bgEst->MeasureScale1D( RealData[0], theFiles, 13, 50, 30, !doQCD, isVjNorm, 1 );
        bgEst->MeasureScale1D( RealData[0], theFiles, 13, 50, 30, !doQCD, isVjNorm, 2 );
        //bgEst->MeasureScale1D( RealData[0], theFiles, 13, 50, 30, !doQCD, isVjNorm, 3 );
        /*
        bgEst->MeasureScale1D( RealData[0], theFiles, 14, 50, 30, !doQCD );
        bgEst->MeasureScale1D( RealData[0], theFiles, 15, 50, 30, !doQCD );
        bgEst->MeasureScale1D( RealData[0], theFiles, 16, 50, 30, !doQCD );
        bgEst->MeasureScale1D( RealData[0], theFiles, 17, 50, 30, !doQCD );
        */
     }
     if ( module[0] == 7 ) {
        bool doQCD = true ;
        bool isVjNorm = true ;
        bgEst->MeasureScale1D( RealData[0], theFiles, 13, 35, 20, doQCD, isVjNorm, 1 );
        bgEst->MeasureScale1D( RealData[0], theFiles, 13, 35, 20, doQCD, isVjNorm, 2 );
        /*
        bgEst->MeasureScale1D( RealData[0], theFiles, 14, 35, 20, doQCD );
        bgEst->MeasureScale1D( RealData[0], theFiles, 15, 35, 20, doQCD );
        bgEst->MeasureScale1D( RealData[0], theFiles, 16, 35, 20, doQCD );
        bgEst->MeasureScale1D( RealData[0], theFiles, 17, 35, 20, doQCD );
        */
     }
     // module1: 8 & 9 are normalize MC -> MC_PseudoData                  doQCD   isVjNorm
     //if ( module[0] == 8 )  bgEst->MeasureScale1D( theFiles, 13, 50, 30,  false,  false );
     //if ( module[0] == 9 )  bgEst->MeasureScale1D( theFiles, 13, 35, 20,  true,   true );
     // module1: 8 & 9 are normalize system MC ->  Data                  doQCD   isVjNorm
     if ( module[0] == 8 )  {
        bgEst->MeasureScale1D( RealData[0], "VJetSystematic.txt", 13, 50, 30,  false,  false, 1 );
        bgEst->MeasureScale1D( RealData[0], "VJetSystematic.txt", 13, 50, 30,  false,  false, 2 );
     }
     if ( module[0] == 9 )  {
        bgEst->MeasureScale1D( RealData[0], "VJetSystematic.txt", 13, 35, 20,  true,   true,  1 );
        bgEst->MeasureScale1D( RealData[0], "VJetSystematic.txt", 13, 35, 20,  true,   true,  2 );
     }
     if ( module[0] == 10 )  {
        bgEst->MeasureScale1D( RealData[0], "OtherSystematic.txt", 13, 50, 30,  false,  false, 1 );
        bgEst->MeasureScale1D( RealData[0], "OtherSystematic.txt", 13, 50, 30,  false,  false, 2 );
     }
     if ( module[0] == 11 )  {
        bgEst->MeasureScale1D( RealData[0], "OtherSystematic.txt", 13, 35, 20,  true,   true,  1 );
        bgEst->MeasureScale1D( RealData[0], "OtherSystematic.txt", 13, 35, 20,  true,   true,  2 );
     }
     if ( module[0] == 12 )  {
        bgEst->MeasureScale1D( RealData[0], "TtSystematic.txt", 13, 50, 30,  false,  false, 1 );
        bgEst->MeasureScale1D( RealData[0], "TtSystematic.txt", 13, 50, 30,  false,  false, 2 );
     }
     if ( module[0] == 13 )  {
        bgEst->MeasureScale1D( RealData[0], "TtSystematic.txt", 13, 35, 20,  true,   true,  1 );
        bgEst->MeasureScale1D( RealData[0], "TtSystematic.txt", 13, 35, 20,  true,   true,  2 );
     }
     
     if ( module[0] == 1 ) {
     
        h_objs->Reset(1,true); // normalize the MC
        h_objs->ObjHistoPlotter( the4JFiles[0], smearing );
        h_objs->ObjHistoPlotter( the4JFiles[1], smearing );
        h_objs->ObjHistoPlotter( the4JFiles[2], smearing );
        h_objs->ObjHistoPlotter( the4JFiles[3], smearing );
        h_objs->ObjHistoPlotter( the4JFiles[4], smearing );
        h_objs->ObjHistoPlotter( the4JFiles[5], smearing );
        h_objs->ObjHistoPlotter( the4JFiles[6], smearing );

        h_objs->ObjHistoPlotter( the2JFiles[0], smearing );
        h_objs->ObjHistoPlotter( the2JFiles[1], smearing );
        h_objs->ObjHistoPlotter( the2JFiles[2], smearing );
        h_objs->ObjHistoPlotter( the2JFiles[3], smearing );
        h_objs->ObjHistoPlotter( the2JFiles[4], smearing );
        h_objs->ObjHistoPlotter( the2JFiles[5], smearing );
        h_objs->ObjHistoPlotter( the2JFiles[6], smearing );
        
        //h_objs->JacobTester();
     }

  }
  if ( module[1] == 1 ) {
        bool doPlots = false ;
        bool normMC  = true ;
        //bgEst->RatioXY( the4JFiles, the2JFiles, -1, withTt, doPlots, 20 ) ;
        //bgEst->CombinedRatio( 2, 1 ) ;
        //bgEst->CombinedRatio( 4, 1, true, false ) ;

        cout<<"  === R21 MC truth without Normalization === "<<endl;
        bgEst->RatioXY( 2, 1, theFiles,  1, doPlots, false, false, false ) ;
        bgEst->RatioXY( 2, 1, theFiles,  2, doPlots, false, false, false ) ;
        cout<<"  === R21 MC using cotrol region without Normalization === "<<endl;
        bgEst->RatioXY( 2, 1, theFiles,  3, doPlots, false, false, false ) ;
        bgEst->RatioXY( 2, 1, theFiles,  4, doPlots, false, false, false ) ;
        cout<<"  === R21 MC truth With  Normalization === "<<endl;
        bgEst->RatioXY( 2, 1, theFiles,  1, doPlots, false, false, normMC ) ;
        bgEst->RatioXY( 2, 1, theFiles,  2, doPlots, false, false, normMC ) ;
        cout<<"  === R21 MC using cotrol region with Normalization === "<<endl;
        bgEst->RatioXY( 2, 1, theFiles,  3, doPlots, false, false, normMC ) ;
        bgEst->RatioXY( 2, 1, theFiles,  4, doPlots, false, false, normMC ) ;
        
        cout<<"  === R21 Data in control region === "<<endl;
        bgEst->RatioXY( 2, 1, RealData,  3, doPlots, false, false, false ) ;
        bgEst->RatioXY( 2, 1, RealData,  4, doPlots, false, false, false ) ;
        /*
        cout<<"  === R21 using full MET MT region === "<<endl;
        bgEst->RatioXY( 2, 1, theFiles,  0, doPlots, false, false, normMC ) ;
        bgEst->RatioXY( 2, 1, RealData,  0, doPlots, false, false, false ) ;
        */

        /*
        cout<<"  === R42 === "<<endl;
        bgEst->RatioXY( 4, 2, theFiles,  1, doPlots, true, false, normMC ) ;
        bgEst->RatioXY( 4, 2, theFiles,  2, doPlots, true, false, normMC ) ;
        bgEst->RatioXY( 4, 2, theFiles,  0, doPlots, true, false, normMC ) ;
        cout<<"  === RData === "<<endl;
        cout<<"  === R21 using full MET MT region === "<<endl;
        bgEst->RatioXY( 2, 1, theFiles,  0, doPlots, false, false, normMC ) ;
        */
  }

  if ( module[1] == 2 ) {
        bool doPlots = false ;
        bool normMC  = true ;
        bgEst->RatioXY( 4, 2, theFiles, 0, doPlots, true, false, normMC ) ;
  }
  if ( module[1] == 3) bgEst->RatioForSystematic( "OtherSystematic.txt" ); 
  if ( module[1] == 4) bgEst->RatioForSystematic( "TtSystematic.txt" ); 
  if ( module[1] == 5) bgEst->RatioForSystematic( "VJetSystematic.txt" ); 
  if ( module[1] == 6 ) bgEst->RatioX2( theFiles,  RealData ) ;
  if ( module[1] == 7 ) {
        bool normMC  = false ;
        bgEst->RatioScan( RealData, theFiles, normMC ) ;
  }
  if ( module[1] == 8 ) bgEst->MeasureScale2D( RealData[0], theFiles, 130, 130, true, false );
  if ( module[1] == 9 ) bgEst->VQFraction( theFiles );

  if ( module[2] == 1 ) {
     if ( mcMatching == "ON" )   WFitter->HadTopFitter( the4JFiles[0],  DrawOpt, true );
     
     WFitter->HadTopFitter( RealData[0],  DrawOpt );
     WFitter->HadTopFitter( theFiles[0],  DrawOpt );
     WFitter->HadTopFitter( theFiles[1],  DrawOpt );
     /*
     WFitter->HadTopFitter( theFiles[2],  DrawOpt );
     WFitter->HadTopFitter( theFiles[3],  DrawOpt );
     WFitter->HadTopFitter( theFiles[4],  DrawOpt );
     WFitter->HadTopFitter( theFiles[5],  DrawOpt );
     */
     WFitter->HadTopFitter( theFiles[6],  DrawOpt );
  }

  if ( module[2] == 2 ) {
     if ( mcMatching == "ON" )   WFitter->LepTopFitter( theFiles[0],  DrawOpt, true );
     WFitter->LepTopFitter( theFiles[0],  DrawOpt );
     WFitter->LepTopFitter( theFiles[1],  DrawOpt );
     WFitter->LepTopFitter( theFiles[2],  DrawOpt );
     WFitter->LepTopFitter( theFiles[6],  DrawOpt );
  }
  if ( module[2] == 3 ) WFitter->M2M3_1DPlots( RealData[0], theFiles, false ); // no scale
  if ( module[2] == 4 ) WFitter->M2M3_1DPlots( RealData[0], theFiles, true );  // scaled

  if ( module[3] == 1 ) {
     WFitter->MixAll( theFiles, DrawOpt );
     //WFitter->EnsembleTest( randomSeed , DrawOpt );
  }

  if ( module[4] == 1)  WFitter->Had_SBRatio( theFiles );
  if ( module[4] == 2)  WFitter->SBCEPlotter();
  if ( module[4] == 3)  WFitter->SBPlotter();
  
  // check the master switch for MC-Data Normalization
  if ( module[5] == 1)  xSec->BgClosureTest(4,2,true,false);
  if ( module[5] == 2)  xSec->CutEff( 4 ) ;
  if ( module[5] == 3)  xSec->MethodTest1( 4 ) ;
  if ( module[5] == 4)  xSec->EnsembleTest1( nPseudoExp , randomSeed ) ;
  if ( module[5] == 5)  xSec->CutEff( "Systematic.txt", 4 ) ;
  if ( module[5] == 6)  xSec->MethodTest2( "Systematic.txt", 4 ) ;
  if ( module[5] == 7)  xSec->MethodTest3( 4 ) ;  // verify root working normally
  
  if ( module[6] == 1)   xSec->RealDataAnalysis( 4, 2 );
  if ( module[6] == 2)   xSec->DataForSystematic( "Systematic.txt", 4 );

  return 0;
}

