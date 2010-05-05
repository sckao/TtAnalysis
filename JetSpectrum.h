#ifndef JetSpectrum_H
#define JetSpectrum_H

#include <TObject.h>
#include <TMinuit.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom2.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TGraphErrors.h>
#include <TSystem.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>

#include "MassAnaInput.h"
#include "WFormat.h"

class JetSpectrum : public TObject {

private:

   double mL;
   double mH;

   // for hadronic permutation

   MassAnaInput*    fitInput;

   string hfolder;

public:

   JetSpectrum( double massL, double massH );     
   ~JetSpectrum();     
 

   void EtSpectrum( string fileName, recoObj* histos );

   void ObjHistoPlotter( string fileName );

   //ClassDef(JetSpectrum, 1);

};

//#if !defined(__CINT__)
//    ClassImp(JetSpectrum);
#endif

