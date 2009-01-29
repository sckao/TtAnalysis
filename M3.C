#include <vector>
#include <stdio.h>
#include <TFile.h>
#include <iostream>
#include <fstream>

void M3() {

 TFile *file  = TFile::Open("ttj_fall08.root");
 TFile *file1 = TFile::Open("wjets_fall08.root");
 TFile *file2 = TFile::Open("qcd_fall08.root");
 //TString hfolder = "wjets_test";
 TString hfolder = "tt_fall08";

 TString name1 = "m3_j3";

 TString plot3 = "m3_j3.gif";

 m3_j3   = (TH2F *) file->Get("Jets/"+name1); 
 m3_j3a  = (TH2F *) file1->Get("Jets/"+name1); 
 m3_j3b  = (TH2F *) file2->Get("Jets/"+name1); 

 gSystem->mkdir(hfolder);
 gSystem->cd(hfolder);

 TCanvas *c3 = new TCanvas("c3","",800,600);
 c3->SetFillColor(10);
 c3->SetFillColor(10);
 c3->Divide(2,2);

 gStyle->SetStatY(0.95); 
 gStyle->SetStatTextColor(1);

 c3->cd(1);
 gPad->SetGrid();
 m3_j3->RebinX(4);
 m3_j3->RebinY(4);
 m3_j3->SetAxisRange(1,150,"Y");
 m3_j3->DrawCopy("BOX");

 c3->cd(2);
 gPad->SetGrid();
 m3_j3a->RebinX(4);
 m3_j3a->RebinY(4);
 m3_j3a->SetAxisRange(1,150,"Y");
 m3_j3a->DrawCopy("BOX");

 c3->cd(3);
 gPad->SetGrid();
 m3_j3b->RebinX(4);
 m3_j3b->RebinY(4);
 m3_j3b->SetAxisRange(1,150,"Y");
 m3_j3b->DrawCopy("BOX");

 c3->cd(4);
 m3_j3->ProjectionX("m3_ttjets",1,500,"");
 m3_j3a->ProjectionX("m3_wjets",1,500,"");
 m3_j3b->ProjectionX("m3_qcd",1,500,"");

 c3->Update();
 gStyle->SetStatY(0.95); 
 gStyle->SetStatTextColor(1);
 m3_ttjets->SetTitle(" M3 distribution ");
 m3_ttjets->SetLineColor(1);
 m3_ttjets->DrawCopy();
 c3->Update();
 gStyle->SetStatY(0.70); 
 gStyle->SetStatTextColor(2);
 m3_wjets->SetLineColor(2);
 m3_wjets->DrawCopy("sames");
 c3->Update();
 gStyle->SetStatY(0.45); 
 gStyle->SetStatTextColor(4);
 m3_qcd->SetLineColor(4);
 m3_qcd->DrawCopy("sames");

 c3->Update();
 c3->Print(plot3);

 gStyle->SetStatY(0.95); 
 gStyle->SetStatTextColor(1);

 gSystem->cd("../");
 file2->Close();
 file1->Close();
 file->Close();

}


