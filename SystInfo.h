#ifndef SystInfo_H
#define SystInfo_H

#include <TObject.h>
#include <TMinuit.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TColor.h>
#include <TPaveText.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <TGaxis.h>

#include "MassAnaInput.h"
#include "ObjectInfo.h"
#include "WFormat.h"


class SystInfo : public TObject {

private:

   string hfolder ;
   string plotType ;
   int JESType ;

   MassAnaInput*    fitInput;
   ObjectInfo*      objInfo;

public:

   SystInfo();     
   ~SystInfo();     

   void SystPlotter();
   void JESPlotter( vector<string>& fakeData, double norm = 0);


   //ClassDef(SystInfo, 1);

};

//#if !defined(__CINT__)
//    ClassImp(SystInfo);
#endif

