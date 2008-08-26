#include <vector>
#include <stdio.h>
#include <TFile.h>
#include <iostream>
#include <fstream>

void PatTreeReader(Int_t evid ) {

  // n = (RUN_number*100000) + (event_number)

  // open files and set up names
  TFile *file = TFile::Open("tt_test_full210.root");
 
  TString hfolder = "tt_test_full";
  TString suffix = ".gif";
  gSystem->mkdir(hfolder);
  gSystem->cd(hfolder);

  FILE *dfile = fopen("TtInfoDump.txt","a");

  fprintf (dfile," ============================================================== \n" ); 
  fprintf (dfile," =======     Event(%i) Display in eta-phi plane     ======= \n", evid ); 
  fprintf (dfile," ============================================================== \n" ); 
  fprintf (dfile,"  \n" ); 

  // name the histogram
  char evtID[9];
  sprintf(evtID,"EVT%d",evid);
  TString evt_id = evtID;

  TString plot01 = "qPt_"+evt_id+suffix;
  TString plot02 = "qE_"+evt_id+suffix;

  TString plot11 = "Jpt_"+evt_id+suffix;
  TString plot12 = "Jem_"+evt_id+suffix;
  TString plot13 = "Jhd_"+evt_id+suffix;

  TString plot21 = "Ept_"+evt_id+suffix;
  TString plot22 = "Eem_"+evt_id+suffix;
  TString plot23 = "Ehd_"+evt_id+suffix;

  TString plot31 = "Mpt_"+evt_id+suffix;
  TString plot32 = "Mem_"+evt_id+suffix;
  TString plot33 = "Mhd_"+evt_id+suffix;

  TString plot41 = "Gpt_"+evt_id+suffix;
  TString plot42 = "Gem_"+evt_id+suffix;
  TString plot43 = "Ghd_"+evt_id+suffix;

  TString plot51 = "Npt_"+evt_id+suffix;
  TString plot52 = "Nem_"+evt_id+suffix;
  TString plot53 = "Nhd_"+evt_id+suffix;

  TString plot61 = "pt_"+evt_id+suffix;
  TString plot62 = "calo_"+evt_id+suffix;

  // read the tree and branch
  TTree *jetT = (TTree*) file->Get("jetT");
  TBranch *b0 = jetT->GetBranch("gen");
  TBranch *b1 = jetT->GetBranch("patJ");
  TBranch *b2 = jetT->GetBranch("patE");
  TBranch *b3 = jetT->GetBranch("patMu");
  TBranch *b4 = jetT->GetBranch("patGa");
  TBranch *b5 = jetT->GetBranch("patNu");

  // read the leaves from the branch
  TLeaf *h0  = b0->GetLeaf("eta");
  TLeaf *f0  = b0->GetLeaf("phi");
  TLeaf *lE  = b0->GetLeaf("energy");
  TLeaf *lpt = b0->GetLeaf("pt");
  TLeaf *lId = b0->GetLeaf("eventId");

  TLeaf *h1  = b1->GetLeaf("eta");
  TLeaf *f1  = b1->GetLeaf("phi");
  TLeaf *em1 = b1->GetLeaf("caloE");
  TLeaf *hd1 = b1->GetLeaf("caloH");
  TLeaf *pt1 = b1->GetLeaf("pt");
  TLeaf *Id1 = b1->GetLeaf("eventId");
  
  TLeaf *h2 = b2->GetLeaf("eta");
  TLeaf *f2 = b2->GetLeaf("phi");
  TLeaf *em2 = b2->GetLeaf("caloE");
  TLeaf *hd2 = b2->GetLeaf("caloH");
  TLeaf *pt2 = b2->GetLeaf("pt");
  TLeaf *Id2 = b2->GetLeaf("eventId");
  
  TLeaf *h3 = b3->GetLeaf("eta");
  TLeaf *f3 = b3->GetLeaf("phi");
  TLeaf *em3 = b3->GetLeaf("caloE");
  TLeaf *hd3 = b3->GetLeaf("caloH");
  TLeaf *pt3 = b3->GetLeaf("pt");
  TLeaf *Id3 = b3->GetLeaf("eventId");
  
  TLeaf *h4 = b4->GetLeaf("eta");
  TLeaf *f4 = b4->GetLeaf("phi");
  TLeaf *em4 = b4->GetLeaf("caloE");
  TLeaf *hd4 = b4->GetLeaf("caloH");
  TLeaf *pt4 = b4->GetLeaf("pt");
  TLeaf *Id4 = b4->GetLeaf("eventId");
  
  TLeaf *h5 = b5->GetLeaf("eta");
  TLeaf *f5 = b5->GetLeaf("phi");
  TLeaf *em5 = b5->GetLeaf("caloE");
  TLeaf *hd5 = b5->GetLeaf("caloH");
  TLeaf *pt5 = b5->GetLeaf("pt");
  TLeaf *Id5 = b5->GetLeaf("eventId");
  

  // Define the histograms
  h1_etaphi_pt  = new TH2D("h1_etaphi_pt","eta-phi  pt   for patJet",79,-3.95,3.95,63,-3.15,3.15);
  h1_etaphi_em  = new TH2D("h1_etaphi_em","eta-phi caloE for patJet",79,-3.95,3.95,63,-3.15,3.15);
  h1_etaphi_hd  = new TH2D("h1_etaphi_hd","eta-phi caloH for patJet",79,-3.95,3.95,63,-3.15,3.15);
  h2_etaphi_pt  = new TH2D("h2_etaphi_pt","eta-phi  pt   for patElc",79,-3.95,3.95,63,-3.15,3.15);
  h2_etaphi_em  = new TH2D("h2_etaphi_em","eta-phi caloE for patElc",79,-3.95,3.95,63,-3.15,3.15);
  h2_etaphi_hd  = new TH2D("h2_etaphi_hd","eta-phi caloH for patElc",79,-3.95,3.95,63,-3.15,3.15);
  h3_etaphi_pt  = new TH2D("h3_etaphi_pt","eta-phi  pt   for patMuon",79,-3.95,3.95,63,-3.15,3.15);
  h3_etaphi_em  = new TH2D("h3_etaphi_em","eta-phi caloE for patMuon",79,-3.95,3.95,63,-3.15,3.15);
  h3_etaphi_hd  = new TH2D("h3_etaphi_hd","eta-phi caloH for patMuon",79,-3.95,3.95,63,-3.15,3.15);
  h4_etaphi_pt  = new TH2D("h4_etaphi_pt","eta-phi  pt   for patGam",79,-3.95,3.95,63,-3.15,3.15);
  h4_etaphi_em  = new TH2D("h4_etaphi_em","eta-phi caloE for patGam",79,-3.95,3.95,63,-3.15,3.15);
  h4_etaphi_hd  = new TH2D("h4_etaphi_hd","eta-phi caloH for patGam",79,-3.95,3.95,63,-3.15,3.15);
  h5_etaphi_pt  = new TH2D("h5_etaphi_pt","eta-phi  pt   for patMET",79,-3.95,3.95,63,-3.15,3.15);
  h5_etaphi_em  = new TH2D("h5_etaphi_em","eta-phi caloE for patMET",79,-3.95,3.95,63,-3.15,3.15);
  h5_etaphi_hd  = new TH2D("h5_etaphi_hd","eta-phi caloH for patMET",79,-3.95,3.95,63,-3.15,3.15);
  h0_etaphi_pt  = new TH2D("h0_etaphi_pt","eta-phi  pt   for patMET",79,-3.95,3.95,63,-3.15,3.15);
  h0_etaphi_E   = new TH2D("h0_etaphi_E" ,"eta-phi E for qq from W " ,79,-3.95,3.95,63,-3.15,3.15);

  // get entries for branch
  Int_t b0sz = b0->GetEntries();
  Int_t b1sz = b1->GetEntries();
  Int_t b2sz = b2->GetEntries();
  Int_t b3sz = b3->GetEntries();
  Int_t b4sz = b4->GetEntries();
  Int_t b5sz = b5->GetEntries();

  fprintf (dfile," ---------------    PAT JET INFO      ------------------------- \n" ); 
  fprintf (dfile," -------------------------------------------------------------- \n" ); 
  fprintf (dfile,"      eta    phi         emE       hdE     caloE      pt \n" ); 
  fprintf (dfile," -------------------------------------------------------------- \n" ); 

  vector<Double_t> x1;
  vector<Double_t> y1;
  vector<Double_t> ptV1;
  vector<Double_t> emE1;
  vector<Double_t> hdE1;
  for (int k=0; k<b1sz; k++) {

      Int_t q  = b1->GetEntry(k);
      Int_t nu = Id1->GetValue();
      if ( evid != nu ) continue;
       
      Double_t h = h1->GetValue();
      Double_t f = f1->GetValue();
      Double_t pt = pt1->GetValue(); 
      Double_t em = em1->GetValue(); 
      Double_t hd = hd1->GetValue(); 
  
      x1.push_back( h );  
      y1.push_back( f );
      ptV1.push_back( pt );
      emE1.push_back( em );
      hdE1.push_back( hd );
      cout<<"ID:"<<nu<<" h:"<<h<<" f:"<<f<<" em:"<<em<<" hd:"<<hd<<endl; 
      fprintf (dfile,"   %6.2f,%6.2f, | %8.2f, %8.2f, %8.2f, %8.2f \n" , 
                             h, f, em, hd, em+hd, pt );
  }

  cout<<" number of jets "<< b1sz <<endl;
  fprintf (dfile," -------------------------------------------------------------- \n" ); 
  fprintf (dfile,"    (Jets in this Event)/(All selected jets) = %i/%i \n", x1.size(), b1sz );
  fprintf (dfile," ============================================================== \n" ); 
  fprintf (dfile,"  \n" ); 

  fprintf (dfile," ---------------    PAT ELECTRON INFO     --------------------- \n" ); 
  fprintf (dfile," -------------------------------------------------------------- \n" ); 
  fprintf (dfile,"      eta    phi         emE       hdE     caloE      pt \n" ); 
  fprintf (dfile," -------------------------------------------------------------- \n" ); 

  vector<Double_t> x2;
  vector<Double_t> y2;
  vector<Double_t> ptV2;
  vector<Double_t> emE2;
  vector<Double_t> hdE2;
  for (int k=0; k<b2sz; k++) {

      Int_t q  = b2->GetEntry(k);
      Int_t nu = Id2->GetValue();
      if ( evid != nu ) continue;
       
      Double_t h  = h2->GetValue();
      Double_t f  = f2->GetValue();
      Double_t pt = pt2->GetValue(); 
      Double_t em = em2->GetValue(); 
      Double_t hd = hd2->GetValue(); 
  
      x2.push_back( h );  
      y2.push_back( f );
      ptV2.push_back( pt );
      emE2.push_back( em );
      hdE2.push_back( hd );
      cout<<"ID:"<<nu<<" h:"<< h <<" f:"<< f <<" em:"<< em <<" hd:"<< hd <<endl; 
      fprintf (dfile,"   %6.2f,%6.2f, | %8.2f, %8.2f, %8.2f, %8.2f \n" , 
                             h ,f ,em ,hd ,em+hd ,pt );

  }
  cout<<" number of electrons "<< b2sz <<endl;
  fprintf (dfile," -------------------------------------------------------------- \n" ); 
  fprintf (dfile,"    (Electrons in this Event )/(All selected Electrons) = %i/%i \n", x2.size(), b2sz );
  fprintf (dfile," ============================================================== \n" ); 
  fprintf (dfile,"  \n" ); 


  fprintf (dfile," ---------------    PAT MUON INFO     ------------------------- \n" ); 
  fprintf (dfile," -------------------------------------------------------------- \n" ); 
  fprintf (dfile,"      eta    phi         emE       hdE     caloE      pt \n" ); 
  fprintf (dfile," -------------------------------------------------------------- \n" ); 

  vector<Double_t> x3;
  vector<Double_t> y3;
  vector<Double_t> ptV3;
  vector<Double_t> emE3;
  vector<Double_t> hdE3;
  for (int k=0; k<b3sz; k++) {

      Int_t q  = b3->GetEntry(k);
      Int_t nu = Id3->GetValue();
      if ( evid != nu ) continue;
       
      Double_t h  = h3->GetValue();
      Double_t f  = f3->GetValue();
      Double_t pt = pt3->GetValue(); 
      Double_t em = em3->GetValue(); 
      Double_t hd = hd3->GetValue(); 
  
      x3.push_back( h );  
      y3.push_back( f );
      ptV3.push_back( pt );
      emE3.push_back( em );
      hdE3.push_back( hd );
      cout<<"ID:"<<nu<<" h:"<< h <<" f:"<< f <<" em:"<< em <<" hd:"<< hd <<endl; 
      fprintf (dfile,"   %6.2f,%6.2f, | %8.2f, %8.2f, %8.2f, %8.2f \n" , 
                             h ,f ,em ,hd ,em+hd ,pt );

  }
  cout<<" number of Muons "<< b3sz <<endl;
  fprintf (dfile," -------------------------------------------------------------- \n" ); 
  fprintf (dfile,"    (Muons in this Event )/(All selected Muons) = %i/%i \n", x3.size(), b3sz );
  fprintf (dfile," ============================================================== \n" ); 
  fprintf (dfile,"  \n" ); 
  

  fprintf (dfile," ---------------    PAT Photon INFO     ----------------------- \n" ); 
  fprintf (dfile," -------------------------------------------------------------- \n" ); 
  fprintf (dfile,"      eta    phi         emE       hdE     caloE      pt \n" ); 
  fprintf (dfile," -------------------------------------------------------------- \n" ); 

  vector<Double_t> x4;
  vector<Double_t> y4;
  vector<Double_t> ptV4;
  vector<Double_t> emE4;
  vector<Double_t> hdE4;
  for (int k=0; k<b4sz; k++) {

      Int_t q  = b4->GetEntry(k);
      Int_t nu = Id4->GetValue();
      if ( evid != nu ) continue;
       
      Double_t h  = h4->GetValue();
      Double_t f  = f4->GetValue();
      Double_t pt = pt4->GetValue(); 
      Double_t em = em4->GetValue(); 
      Double_t hd = hd4->GetValue(); 
  
      x4.push_back( h );  
      y4.push_back( f );
      ptV4.push_back( pt );
      emE4.push_back( em );
      hdE4.push_back( hd );
      cout<<"ID:"<<nu<<" h:"<< h <<" f:"<< f <<" em:"<< em <<" hd:"<< hd <<endl; 
      fprintf (dfile,"   %6.2f,%6.2f, | %8.2f, %8.2f, %8.2f, %8.2f \n" , 
                             h ,f ,em ,hd ,em+hd ,pt );

  }
  cout<<" number of Photns "<< b4sz <<endl;
  fprintf (dfile," -------------------------------------------------------------- \n" ); 
  fprintf (dfile,"    (Photons in this Event )/(All selected Photons) = %i/%i \n", x4.size(), b4sz );
  fprintf (dfile," ============================================================== \n" ); 
  fprintf (dfile,"  \n" ); 

  fprintf (dfile," ------------------    PAT MET INFO     ----------------------- \n" ); 
  fprintf (dfile," -------------------------------------------------------------- \n" ); 
  fprintf (dfile,"      eta    phi         emE       hdE     caloE      pt \n" ); 
  fprintf (dfile," -------------------------------------------------------------- \n" ); 

  vector<Double_t> x5;
  vector<Double_t> y5;
  vector<Double_t> ptV5;
  vector<Double_t> emE5;
  vector<Double_t> hdE5;
  for (int k=0; k<b5sz; k++) {

      Int_t q  = b5->GetEntry(k);
      Int_t nu = Id5->GetValue();
      if ( evid != nu ) continue;
       
      Double_t h  = h5->GetValue();
      Double_t f  = f5->GetValue();
      Double_t pt = pt5->GetValue(); 
      Double_t em = em5->GetValue(); 
      Double_t hd = hd5->GetValue(); 
  
      x5.push_back( h );  
      y5.push_back( f );
      ptV5.push_back( pt );
      emE5.push_back( em );
      hdE5.push_back( hd );
      cout<<"ID:"<<nu<<" h:"<< h <<" f:"<< f <<" em:"<< em <<" hd:"<< hd <<endl; 
      fprintf (dfile,"   %6.2f,%6.2f, | %8.2f, %8.2f, %8.2f, %8.2f \n" , 
                             h ,f ,em ,hd ,em+hd ,pt );

  }
  cout<<" number of METs "<< b4sz <<endl;
  fprintf (dfile," -------------------------------------------------------------- \n" ); 
  fprintf (dfile,"    (with MET in this Event )/(All selected METs) = %i/%i \n", x5.size(), b5sz );
  fprintf (dfile," ============================================================== \n" ); 
  fprintf (dfile,"  \n" ); 

  fprintf (dfile," ---------------       GEN INFO       ------------------------- \n" ); 
  fprintf (dfile," -------------------------------------------------------------- \n" ); 
  fprintf (dfile,"      eta    phi                          total_E     pt \n" ); 
  fprintf (dfile," -------------------------------------------------------------- \n" ); 

  vector<Double_t> x0;
  vector<Double_t> y0;
  vector<Double_t> E0;
  vector<Double_t> ptV0;
  for (int k=0; k < b0sz; k++) {

      Int_t q  = b0->GetEntry(k);
      Int_t nu = lId->GetValue();
      if ( nu != evid ) continue;

      Double_t h = h0->GetValue();
      Double_t f = f0->GetValue();
      Double_t E = lE->GetValue(); 
      Double_t pt = lpt->GetValue(); 

      x0.push_back( h );  
      y0.push_back( f );
      E0.push_back( E );
      ptV0.push_back( pt );

      cout<<"ID:"<< nu <<" h:"<< h <<" f:"<< f <<" E:"<< E <<endl; 
      fprintf (dfile,"   %6.2f,%6.2f, |            %17.2f, %8.2f \n" , 
                           h, f, E, pt );

  }
  cout<<" number of qq "  << b0sz <<endl;
  fprintf (dfile," -------------------------------------------------------------- \n" ); 
  fprintf (dfile,"    (selected qq )/(total qq in Event) = %i/%i \n", x0.size(), b0sz );
  fprintf (dfile," ============================================================== \n" ); 
  fprintf (dfile,"  \n" ); 


  // *******************************
  // *    Fill the Histograms      *
  // *******************************
  Double_t xi,yi,qE0;
  Double_t ptV[6] = { 0. ,0. ,0. ,0. ,0. ,0. };
  Double_t emV[6] = { 0. ,0. ,0. ,0. ,0. ,0. };
  Double_t hdV[6] = { 0. ,0. ,0. ,0. ,0. ,0. };  
  for (int i=0; i<79; i++) {
      xi = -3.9 + (i*0.1) ;

      for (int j=0; j<63; j++) {
           yi  = -3.1 + (j*0.1);
           qE0 = 0.0;
           for (int r=0; r<6; r++) {
               ptV[r] = 0.0;
               emV[r] = 0.0;
               hdV[r] = 0.0;
           }

           // locate and fill pat jets
           for (int k=0; k < x1.size(); k++) {
              Double_t dx = fabs(x1[k]-xi);
              Double_t dy = fabs(y1[k]-yi);
              if ( dx <= 0.05 && dy <= 0.05) {
                ptV[1]  = ptV1[k];               
                emV[1]  = emE1[k];               
                hdV[1]  = hdE1[k];
                continue;
              } 
           }
           if ( ptV[1] != 0.0 ) { h1_etaphi_pt->Fill(xi,yi,ptV[1]); }
           if ( emV[1] != 0.0 ) { h1_etaphi_em->Fill(xi,yi,emV[1]); }
           if ( hdV[1] != 0.0 ) { h1_etaphi_hd->Fill(xi,yi,hdV[1]); }
  
           // locate and fill pat electrons
           for (int k=0; k < x2.size(); k++) {
              Double_t dx = fabs(x2[k]-xi);
              Double_t dy = fabs(y2[k]-yi);
              if ( dx <= 0.05 && dy <= 0.05) {
                ptV[2]  = ptV2[k];               
                emV[2]  = emE2[k];               
                hdV[2]  = hdE2[k];
                continue;
              } 
           }
           if ( ptV[2] != 0.0 ) { h2_etaphi_pt->Fill(xi,yi,ptV[2]); }
           if ( emV[2] != 0.0 ) { h2_etaphi_em->Fill(xi,yi,emV[2]); }
           if ( hdV[2] != 0.0 ) { h2_etaphi_hd->Fill(xi,yi,hdV[2]); }

           // locate and fill pat muons
           for (int k=0; k < x3.size(); k++) {
              Double_t dx = fabs(x3[k]-xi);
              Double_t dy = fabs(y3[k]-yi);
              if ( dx <= 0.05 && dy <= 0.05) {
                ptV[3]  = ptV3[k];               
                emV[3]  = emE3[k];               
                hdV[3]  = hdE3[k];
                continue;
              } 
           }
           if ( ptV[3] != 0.0 ) { h3_etaphi_pt->Fill(xi,yi,ptV[3]); }
           if ( emV[3] != 0.0 ) { h3_etaphi_em->Fill(xi,yi,emV[3]); }
           if ( hdV[3] != 0.0 ) { h3_etaphi_hd->Fill(xi,yi,hdV[3]); }

           // locate and fill pat photons
           for (int k=0; k < x4.size(); k++) {
              Double_t dx = fabs(x4[k]-xi);
              Double_t dy = fabs(y4[k]-yi);
              if ( dx <= 0.05 && dy <= 0.05) {
                ptV[4]  = ptV4[k];               
                emV[4]  = emE4[k];               
                hdV[4]  = hdE4[k];
                continue;
              } 
           }
           if ( ptV[4] != 0.0 ) { h4_etaphi_pt->Fill(xi,yi,ptV[4]); }
           if ( emV[4] != 0.0 ) { h4_etaphi_em->Fill(xi,yi,emV[4]); }
           if ( hdV[4] != 0.0 ) { h4_etaphi_hd->Fill(xi,yi,hdV[4]); }

           // locate and fill pat METs
           for (int k=0; k < x5.size(); k++) {
              Double_t dx = fabs(x5[k]-xi);
              Double_t dy = fabs(y5[k]-yi);
              if ( dx <= 0.05 && dy <= 0.05) {
                ptV[5]  = ptV5[k];               
                emV[5]  = emE5[k];               
                hdV[5]  = hdE5[k];
                continue;
              } 
           }
           if ( ptV[5] != 0.0 ) { h5_etaphi_pt->Fill(xi,yi,ptV[5]); }
           if ( emV[5] != 0.0 ) { h5_etaphi_em->Fill(xi,yi,emV[5]); }
           if ( hdV[5] != 0.0 ) { h5_etaphi_hd->Fill(xi,yi,hdV[5]); }

           // locate and fill generators
           for (int k=0; k < x0.size(); k++) {
              Double_t dx = fabs(x0[k]-xi);
              Double_t dy = fabs(y0[k]-yi);
              if ( dx <= 0.05 && dy <= 0.05) {
                qE0  = E0[k];               
                ptV[0]  = ptV0[k];               
                continue;
              } 
           }
           if ( ptV[0] != 0.0 ) { h0_etaphi_pt->Fill(xi,yi,ptV[0]); }
           if (    qE0 != 0.0 ) { h0_etaphi_E->Fill(xi,yi,qE0); }
      }
  }


 // plot the gen information
 gStyle->SetOptStat(kTRUE);
 gStyle->SetOptFit(0111);  
 c0 = new TCanvas("c0","",200,8,1000,700);
 c0->cd();

 h0_etaphi_pt->SetTitle("pt in eta-phi distribution");
 h0_etaphi_pt->GetXaxis()->SetTitle("  #eta  ");
 h0_etaphi_pt->GetYaxis()->SetTitle("  #phi  ");
 h0_etaphi_pt->GetZaxis()->SetTitle(" pt ");
 h0_etaphi_pt->Draw("LEGO2");
 c0->Update();
 c0->Print(plot01);

 h0_etaphi_E->SetTitle(" qq E in eta-phi distribution");
 h0_etaphi_E->GetXaxis()->SetTitle("  #eta  ");
 h0_etaphi_E->GetYaxis()->SetTitle("  #phi  ");
 h0_etaphi_E->GetZaxis()->SetTitle("  Energy ");
 h0_etaphi_E->Draw("LEGO2");
 c0->Update();
 c0->Print(plot02);

 // an empty histogram for base
 h_v  = new TH2D("h_v","",79,-3.95,3.95,63,-3.15,3.15);
 h_v->Draw("LEGO");

 // plot histograms for PAT Jets
 gStyle->SetOptStat(kTRUE);
 gStyle->SetOptFit(0111); 
 c1 = new TCanvas("c1","",200,8,1000,700);
 c1->cd();

 h1s_etaphi_pt = new THStack("h1s_etaphi_pt"," pt in eta-phi distribution ");
 h1s_etaphi_pt->Add(h_v);
 h1_etaphi_pt->Draw("LEGO");
 h1_etaphi_pt->SetFillColor(3);
 h1s_etaphi_pt->Add(h1_etaphi_pt);
 h1s_etaphi_pt->Draw();

 c1->Update();
 c1->Print(plot11);

 h1s_etaphi_em = new THStack("h1s_etaphi_em"," emE in eta-phi distribution ");
 h1s_etaphi_em->Add(h_v);
 h1_etaphi_em->Draw("LEGO");
 h1_etaphi_em->SetFillColor(3);
 h1s_etaphi_em->Add(h1_etaphi_em);
 h1s_etaphi_em->Draw();

 c1->Update();
 c1->Print(plot12);


 h1s_etaphi_hd = new THStack("h1s_etaphi_hd"," hdE in eta-phi distribution ");
 h1s_etaphi_hd->Add(h_v);
 h1_etaphi_hd->Draw("LEGO");
 h1_etaphi_hd->SetFillColor(8);
 h1s_etaphi_hd->Add(h1_etaphi_hd);
 h1s_etaphi_hd->Draw();
 c1->Update();
 c1->Print(plot13);

 // plot histograms for PAT Electrons
 gStyle->SetOptStat(kTRUE);
 gStyle->SetOptFit(0111);  
 c2 = new TCanvas("c2","",200,8,1000,700);
 c2->cd();

 h2s_etaphi_pt = new THStack("h2s_etaphi_pt"," pt in eta-phi distribution ");
 h2s_etaphi_pt->Add(h_v);
 h2_etaphi_pt->Draw("LEGO");
 h2_etaphi_pt->SetFillColor(7);
 h2s_etaphi_pt->Add(h2_etaphi_pt);
 h2s_etaphi_pt->Draw();

 c2->Update();
 c2->Print(plot21);

 h2s_etaphi_em = new THStack("h2s_etaphi_em"," emE in eta-phi distribution ");
 h2s_etaphi_em->Add(h_v);
 h2_etaphi_em->Draw("LEGO");
 h2_etaphi_em->SetFillColor(7);
 h2s_etaphi_em->Add(h2_etaphi_em);
 h2s_etaphi_em->Draw();

 c2->Update();
 c2->Print(plot22);


 h2s_etaphi_hd = new THStack("h2s_etaphi_hd"," hdE in eta-phi distribution ");
 h2s_etaphi_hd->Add(h_v);
 h2_etaphi_hd->Draw("LEGO");
 h2_etaphi_hd->SetFillColor(4);
 h2s_etaphi_hd->Add(h2_etaphi_hd);
 h2s_etaphi_hd->Draw();
 c2->Update();
 c2->Print(plot23);

 // plot histograms for PAT Muons
 gStyle->SetOptStat(kTRUE);
 gStyle->SetOptFit(0111);  
 c3 = new TCanvas("c3","",200,8,1000,700);
 c3->cd();

 h3s_etaphi_pt = new THStack("h3s_etaphi_pt"," pt in eta-phi distribution ");
 h3s_etaphi_pt->Add(h_v);
 h3_etaphi_pt->Draw("LEGO");
 h3_etaphi_pt->SetFillColor(2);
 h3s_etaphi_pt->Add(h3_etaphi_pt);
 h3s_etaphi_pt->Draw();

 c3->Update();
 c3->Print(plot31);

 h3s_etaphi_em = new THStack("h3s_etaphi_em"," emE in eta-phi distribution ");
 h3s_etaphi_em->Add(h_v);
 h3_etaphi_em->Draw("LEGO");
 h3_etaphi_em->SetFillColor(2);
 h3s_etaphi_em->Add(h3_etaphi_em);
 h3s_etaphi_em->Draw();

 c3->Update();
 c3->Print(plot32);


 h3s_etaphi_hd = new THStack("h3s_etaphi_hd"," hdE in eta-phi distribution ");
 h3s_etaphi_hd->Add(h_v);
 h3_etaphi_hd->Draw("LEGO");
 h3_etaphi_hd->SetFillColor(6);
 h3s_etaphi_hd->Add(h3_etaphi_hd);
 h3s_etaphi_hd->Draw();
 c3->Update();
 c3->Print(plot33);

 // plot histograms for PAT Photons
 gStyle->SetOptStat(kTRUE);
 gStyle->SetOptFit(0111);  
 c4 = new TCanvas("c4","",200,8,1000,700);
 c4->cd();

 h4s_etaphi_pt = new THStack("h4s_etaphi_pt"," pt in eta-phi distribution ");
 h4s_etaphi_pt->Add(h_v);
 h4_etaphi_pt->Draw("LEGO");
 h4_etaphi_pt->SetFillColor(5);
 h4s_etaphi_pt->Add(h4_etaphi_pt);
 h4s_etaphi_pt->Draw();

 c4->Update();
 c4->Print(plot41);

 h4s_etaphi_em = new THStack("h4s_etaphi_em"," emE in eta-phi distribution ");
 h4s_etaphi_em->Add(h_v);
 h4_etaphi_em->Draw("LEGO");
 h4_etaphi_em->SetFillColor(5);
 h4s_etaphi_em->Add(h4_etaphi_em);
 h4s_etaphi_em->Draw();

 c4->Update();
 c4->Print(plot42);


 h4s_etaphi_hd = new THStack("h4s_etaphi_hd"," hdE in eta-phi distribution ");
 h4s_etaphi_hd->Add(h_v);
 h4_etaphi_hd->Draw("LEGO");
 h4_etaphi_hd->SetFillColor(28);
 h4s_etaphi_hd->Add(h4_etaphi_hd);
 h4s_etaphi_hd->Draw();
 c4->Update();
 c4->Print(plot43);

 // plot histograms for PAT METs
 gStyle->SetOptStat(kTRUE);
 gStyle->SetOptFit(0111);  
 c5 = new TCanvas("c5","",200,8,1000,700);
 c5->cd();

 h5s_etaphi_pt = new THStack("h5s_etaphi_pt"," pt in eta-phi distribution ");
 h5s_etaphi_pt->Add(h_v);
 h5_etaphi_pt->Draw("LEGO");
 h5_etaphi_pt->SetFillColor(16);
 h5s_etaphi_pt->Add(h5_etaphi_pt);
 h5s_etaphi_pt->Draw();

 c5->Update();
 c5->Print(plot51);

 h5s_etaphi_em = new THStack("h5s_etaphi_em"," emE in eta-phi distribution ");
 h5s_etaphi_em->Add(h_v);
 h5_etaphi_em->Draw("LEGO");
 h5_etaphi_em->SetFillColor(16);
 h5s_etaphi_em->Add(h5_etaphi_em);
 h5s_etaphi_em->Draw();

 c5->Update();
 c5->Print(plot52);


 h5s_etaphi_hd = new THStack("h5s_etaphi_hd"," hdE in eta-phi distribution ");
 h5s_etaphi_hd->Add(h_v);
 h5_etaphi_hd->Draw("LEGO");
 h5_etaphi_hd->SetFillColor(13);
 h5s_etaphi_hd->Add(h5_etaphi_hd);
 h5s_etaphi_hd->Draw();
 c5->Update();
 c5->Print(plot53);


 // combine all pt information
 gStyle->SetOptStat(kTRUE);
 gStyle->SetOptFit(0111);  
 c6 = new TCanvas("c6","",200,8,1000,700);
 c6->cd();

 hs_pt = new THStack("hs_pt"," pt of all objects in eta-phi space ");
 hs_pt->Add(h_v);

 // add jet info
 h1_etaphi_pt->SetFillColor(3);
 hs_pt->Add(h1_etaphi_pt);

 // add electron info
 h2_etaphi_pt->SetFillColor(7);
 hs_pt->Add(h2_etaphi_pt);

 // add muon info
 h3_etaphi_pt->SetFillColor(2);
 hs_pt->Add(h3_etaphi_pt);

 // add photon info
 h4_etaphi_pt->SetFillColor(5);
 hs_pt->Add(h4_etaphi_pt);

 // add MET info
 h5_etaphi_pt->SetFillColor(16);
 hs_pt->Add(h5_etaphi_pt);

 hs_pt->Draw();
 c6->Update();
 c6->Print(plot61);

 // combine all calo information
 gStyle->SetOptStat(kTRUE);
 gStyle->SetOptFit(0111);  
 c7 = new TCanvas("c7","",200,8,1000,700);
 c7->cd();

 hs_caloE = new THStack("hs_caloE"," ECal + HCal energy deposit ");
 hs_caloE->Add(h_v);

 // add jet info
 h1_etaphi_em->SetFillColor(3);
 hs_caloE->Add(h1_etaphi_em);
 h1_etaphi_hd->SetFillColor(8);
 hs_caloE->Add(h1_etaphi_hd);

 // add electron info
 h2_etaphi_em->SetFillColor(7);
 hs_caloE->Add(h2_etaphi_em);
 h2_etaphi_hd->SetFillColor(4);
 hs_caloE->Add(h2_etaphi_hd);

 // add muon info
 h3_etaphi_em->SetFillColor(2);
 hs_caloE->Add(h3_etaphi_em);
 h3_etaphi_hd->SetFillColor(6);
 hs_caloE->Add(h3_etaphi_hd);

 // add photon info
 h4_etaphi_em->SetFillColor(5);
 hs_caloE->Add(h4_etaphi_em);
 h4_etaphi_hd->SetFillColor(28);
 hs_caloE->Add(h4_etaphi_hd);

 // add MET info
 h5_etaphi_em->SetFillColor(16);
 hs_caloE->Add(h5_etaphi_em);
 h5_etaphi_hd->SetFillColor(13);
 hs_caloE->Add(h5_etaphi_hd);

 hs_caloE->Draw();
 c7->Update();
 c7->Print(plot62);


 fclose(dfile);
 file->Close();
 gSystem->cd("../");
 //gROOT->Reset();
 //gROOT->ProcessLine(".q");

 
}
