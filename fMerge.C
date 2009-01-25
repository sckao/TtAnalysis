void fMerge( TString fname, TString outName, int num ) {

 TString mergeName = outName+".root" ; 

 TFileMerger *fileMerger = new TFileMerger();

 // Loop all files to be merged 
 for (int i=1; i<=num; i++) {
     char fID[8];
     sprintf(fID,"_%d.root",i);
     TString subfile = fname + fID  ;
     cout<<"  Merge file"<<i<<" "<< subfile <<endl;
     fileMerger->AddFile(subfile);	
 }
 
 fileMerger->OutputFile(mergeName);
 bool work = fileMerger->Merge();

 if (work) {
    cout<<"  *** Merge Successful *** " <<endl;
 } else {
    cout<<"  !!! Fuck ! Merge Fail !!! " <<endl;
 }

}
