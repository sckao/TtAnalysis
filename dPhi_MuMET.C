void dPhi_MuMET( TString name1 ) {

 TFile *file  = TFile::Open("ttj_fall08.root");
 TFile *file1 = TFile::Open("wjets_fall08.root");
 TFile *file2 = TFile::Open("qcd_fall08.root");
 TString hfolder = "tt_fall08";

 Int_t rb = 2;

 //TString name1 = "MET_dPhi";
 //TString plot1 = "met_dphi.gif";
 
 TString plot1 = name1+".gif";

 met_df   = (TH2F *) file->Get("METs/"+name1); 
 met_dfa  = (TH2F *) file1->Get("METs/"+name1); 
 met_dfb  = (TH2F *) file2->Get("METs/"+name1); 
 //met_dfb  = (TH2F *) file2->Get("METs/MET_dPhi"); 

 gSystem->mkdir(hfolder);
 gSystem->cd(hfolder);
 
 TCanvas *c1 = new TCanvas("c1","", 900, 700);
 c1->SetFillColor(10);
 c1->SetFillColor(10);
 c1->Divide(2,2);

 c1->cd(1);
 met_df->RebinX(rb);
 met_df->SetAxisRange(1,150,"X");
 met_df->DrawCopy("BOX");

 c1->cd(2);
 met_dfa->RebinX(rb);
 met_dfa->SetAxisRange(1,150,"X");
 met_dfa->DrawCopy("BOX");

 c1->cd(3);
 met_dfb->RebinX(rb);
 met_dfb->SetAxisRange(1,150,"X");
 met_dfb->DrawCopy("BOX");

 c1->Update();
 c1->Print(plot1);

 gSystem->cd("../");
 file2->Close();
 file1->Close();
 file->Close();

}

