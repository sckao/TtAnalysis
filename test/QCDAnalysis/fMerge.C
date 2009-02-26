void fMerge( int num ) {

 TString fname      = "qcd_fall08";
 TString mergeName  = "qcd_fall08.root" ; 
 TString fname1     = "qcd_JetEtAnalysis";
 TString mergeName1 = "qcd_JetEtAnalysis.root" ; 
 TString fname2     = "qcd_IsoMuAnalysis";
 TString mergeName2 = "qcd_IsoMuAnalysis.root" ; 

 TFileMerger *fileMerger  = new TFileMerger();
 TFileMerger *fileMerger1 = new TFileMerger();
 TFileMerger *fileMerger2 = new TFileMerger();

 // Loop all files to be merged 
 for (int i=1; i<=num; i++) {
     char fID[8];
     sprintf(fID,"_%d.root",i);
     TString subfile = "/data/top/sckao/Fall08QCDAna/" + fname + fID  ;
     TString subfile1 = "/data/top/sckao/Fall08QCDAna/" + fname1 + fID  ;
     TString subfile2 = "/data/top/sckao/Fall08QCDAna/" + fname2 + fID  ;
     cout<<"  Merge file"<<i<<" "<< subfile <<" - "<<subfile1<<" - "<<subfile2<<endl;
     fileMerger->AddFile(subfile);	
     fileMerger1->AddFile(subfile1);	
     fileMerger2->AddFile(subfile2);	
 }
 
 fileMerger->OutputFile(mergeName);
 bool work = fileMerger->Merge();
 fileMerger1->OutputFile(mergeName1);
 bool work1 = fileMerger1->Merge();
 fileMerger2->OutputFile(mergeName2);
 bool work2 = fileMerger2->Merge();

 if ( work && work1 && work2 ) {
    cout<<"  *** Merge Successful *** " <<endl;
 } else {
    cout<<"  !!! Fuck ! Merge Fail !!! " <<endl;
 }

}
