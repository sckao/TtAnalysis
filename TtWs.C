void TtWs() {

 TFile *file  = TFile::Open("ttj_fall08_3j30.root");
 TFile *file1 = TFile::Open("wjets_fall08_3j30.root");
 TFile *file2 = TFile::Open("qcd_fall08_3j30.root");
 TString hfolder = "tt_fall08";

 //TString name1 = "hRecott0";

 TString plot1 = "recoLepW.gif";
 TString plot2 = "recoHadW.gif";

 tlepWa  = (TH1F *) file->Get("Tops/allRecolepW"); 
 wlepWa  = (TH1F *) file1->Get("Tops/allRecolepW"); 
 qlepWa  = (TH1F *) file2->Get("Tops/allRecolepW");

 tlepWp  = (TH1F *) file->Get("Tops/RecolepW"); 
 wlepWp  = (TH1F *) file1->Get("Tops/RecolepW"); 
 qlepWp  = (TH1F *) file2->Get("Tops/RecolepW");

 thadWa  = (TH1F *) file->Get("Tops/allRecohadW"); 
 whadWa  = (TH1F *) file1->Get("Tops/allRecohadW"); 
 qhadWa  = (TH1F *) file2->Get("Tops/allRecohadW");

 thadWp  = (TH1F *) file->Get("Tops/RecohadW"); 
 whadWp  = (TH1F *) file1->Get("Tops/RecohadW"); 
 qhadWp  = (TH1F *) file2->Get("Tops/RecohadW");

 tSelW  = (TH2F *) file->Get("Tops/selRecoW"); 
 wSelW  = (TH2F *) file1->Get("Tops/selRecoW"); 
 qSelW  = (TH2F *) file2->Get("Tops/selRecoW");

 tMCW  = (TH2F *) file->Get("Tops/selMCW"); 

 // re-bin value for top mass plot
 int rbin = 4;

 gSystem->mkdir(hfolder);
 gSystem->cd(hfolder);

 gStyle->SetOptStat(kTRUE);
 gStyle->SetOptStat(1000101);

 TCanvas *c1 = new TCanvas("c1","", 1000, 800);
 c1->SetFillColor(10);
 c1->SetFillColor(10);
 c1->Divide(2,2);

 c1->cd(1);
 c1_1->SetGridx();
 gStyle->SetStatY(0.95); 
 gStyle->SetStatX(0.80); 
 gStyle->SetStatTextColor(1);

 tlepWa->Rebin(rbin);
 tlepWa->SetLineColor(1);
 tlepWa->DrawCopy();
 c1->Update();

 gStyle->SetStatY(0.75); 
 gStyle->SetStatTextColor(2);
 tSelW->ProjectionX("tSelW_px",0,-1,"");
 tSelW_px->Rebin(rbin);
 tSelW_px->SetLineWidth(2);
 tSelW_px->SetLineColor(2);
 tSelW_px->DrawCopy("sames");
 c1->Update();

 gStyle->SetStatY(0.55); 
 gStyle->SetStatTextColor(4);
 tlepWp->Rebin(rbin);
 tlepWp->SetLineColor(4);
 tlepWp->DrawCopy("sames");
 c1->Update();

 c1->cd(2);
 c1_2->SetGridx();
 gStyle->SetStatY(0.95); 
 gStyle->SetStatX(0.8); 
 gStyle->SetStatTextColor(1);

 wlepWa->Rebin(rbin);
 wlepWa->SetLineColor(1);
 wlepWa->DrawCopy();
 c1->Update();

 gStyle->SetStatY(0.75); 
 gStyle->SetStatTextColor(2);
 wSelW->ProjectionX("wSelW_px",0,-1,"");
 wSelW_px->Rebin(rbin);
 wSelW_px->SetLineWidth(2);
 wSelW_px->SetLineColor(2);
 wSelW_px->DrawCopy("sames");
 c1->Update();

 gStyle->SetStatY(0.55); 
 gStyle->SetStatTextColor(4);
 wlepWp->Rebin(rbin);
 wlepWp->SetLineColor(4);
 wlepWp->DrawCopy("sames");
 c1->Update();

 c1->cd(3);
 c1_3->SetGridx();
 gStyle->SetStatY(0.95); 
 gStyle->SetStatX(0.8); 
 gStyle->SetStatTextColor(1);

 qlepWa->Rebin(rbin);
 qlepWa->SetLineColor(1);
 qlepWa->DrawCopy();
 c1->Update();

 gStyle->SetStatY(0.75); 
 gStyle->SetStatTextColor(2);
 qSelW->ProjectionX("qSelW_px",0,-1,"");
 qSelW_px->Rebin(rbin);
 qSelW_px->SetLineWidth(2);
 qSelW_px->SetLineColor(2);
 qSelW_px->DrawCopy("sames");
 c1->Update();

 gStyle->SetStatY(0.55); 
 gStyle->SetStatTextColor(4);
 qlepWp->Rebin(rbin);
 qlepWp->SetLineColor(4);
 qlepWp->DrawCopy("sames");
 c1->Update();

 c1->cd(4);
 c1_4->SetGridx();
 gStyle->SetStatY(0.95); 
 gStyle->SetStatX(0.80); 
 gStyle->SetStatTextColor(1);

 tMCW->ProjectionX("tMCW_px",0,-1,"");
 tMCW_px->Rebin(rbin);
 tMCW_px->SetLineColor(1);
 tMCW_px->SetAxisRange(0,180,"Y");
 tMCW_px->DrawCopy();
 c1->Update();

 c1->Print(plot1);

 TCanvas *c2 = new TCanvas("c2","", 1000, 800);
 c2->SetFillColor(10);
 c2->SetFillColor(10);
 c2->Divide(2,2);

 c2->cd(1);
 c2_1->SetGridx();
 gStyle->SetStatX(0.85); 
 gStyle->SetStatY(0.95); 
 gStyle->SetStatTextColor(1);
 c2_1->SetLogy();

 thadWa->Rebin(rbin);
 thadWa->SetLineColor(1);
 thadWa->DrawCopy();
 c2->Update();

 gStyle->SetStatX(0.65); 
 gStyle->SetStatTextColor(2);
 tSelW->ProjectionY("tSelW_py",0,-1,"");
 tSelW_py->Rebin(rbin);
 tSelW_py->SetLineWidth(2);
 tSelW_py->SetLineColor(2);
 tSelW_py->DrawCopy("sames");
 c2->Update();

 gStyle->SetStatX(0.45); 
 gStyle->SetStatTextColor(4);
 thadWp->Rebin(rbin);
 thadWp->SetLineColor(4);
 thadWp->DrawCopy("sames");
 c2->Update();

 c2->cd(2);
 c2_2->SetLogy();
 c2_2->SetGridx();
 gStyle->SetStatX(0.85); 
 gStyle->SetStatTextColor(1);

 whadWa->Rebin(rbin);
 whadWa->SetLineColor(1);
 whadWa->DrawCopy();
 c2->Update();

 gStyle->SetStatX(0.65); 
 gStyle->SetStatTextColor(2);
 wSelW->ProjectionY("wSelW_py",0,-1,"");
 wSelW_py->Rebin(rbin);
 wSelW_py->SetLineWidth(2);
 wSelW_py->SetLineColor(2);
 wSelW_py->DrawCopy("sames");
 c2->Update();

 gStyle->SetStatX(0.45); 
 gStyle->SetStatTextColor(4);
 whadWp->Rebin(rbin);
 whadWp->SetLineColor(4);
 whadWp->DrawCopy("sames");
 c2->Update();

 c2->cd(3);
 c2_3->SetLogy();
 c2_3->SetGridx();
 gStyle->SetStatX(0.85); 
 gStyle->SetStatTextColor(1);

 qhadWa->Rebin(rbin);
 qhadWa->SetLineColor(1);
 qhadWa->DrawCopy();
 c2->Update();

 gStyle->SetStatX(0.65); 
 gStyle->SetStatTextColor(2);
 qSelW->ProjectionY("qSelW_py",0,-1,"");
 qSelW_py->Rebin(rbin);
 qSelW_py->SetLineWidth(2);
 qSelW_py->SetLineColor(2);
 qSelW_py->DrawCopy("sames");
 c2->Update();

 gStyle->SetStatX(0.45); 
 gStyle->SetStatTextColor(4);
 qhadWp->Rebin(rbin);
 qhadWp->SetLineColor(4);
 qhadWp->DrawCopy("sames");
 c2->Update();

 c2->cd(4);
 c2_4->SetLogy();
 c2_4->SetGridx();
 gStyle->SetStatX(0.85); 
 gStyle->SetStatY(0.95); 
 gStyle->SetStatTextColor(1);

 tMCW->ProjectionY("tMCW_py",0,-1,"");
 tMCW_py->Rebin(rbin);
 tMCW_py->SetLineColor(1);
 tMCW_py->SetAxisRange(1,10000,"Y");
 tMCW_py->DrawCopy();
 c2->Update();

 c2->Print(plot2);

 // Reset gStyle
 gStyle->SetStatY(0.98); 
 gStyle->SetStatX(0.98); 
 gStyle->SetStatTextColor(1);
 gStyle->SetOptStat(1111);

 gSystem->cd("../");

 file2->Close();
 file1->Close();
 file->Close();

}
