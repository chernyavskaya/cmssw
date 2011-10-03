#include "CMGTools/H2TauTau/interface/Sample.h"

//ClassImp(Sample);

Sample::Sample():
  TNamed("Sample","Sample"),
  crossection_(1.0),
  lumi_(0.0),
  genEvents_(0),
  nEvents_(0),
  chainNEvents_(0),
  processeff_(1.0),
  sampleChain_(NULL),
  firstrun_(0),
  lastrun_(0),
  effCorrFactor_(1.0),
  histFile_(NULL),
  dataType_(""),
  color_(0),
  lstyle_(0),
  plotOrder_(0),
  init_(0)
{}

Sample::Sample(const char * name, const char * path):
  TNamed(name,path),
  crossection_(1.0),
  lumi_(0.0),
  genEvents_(0),
  nEvents_(0),
  chainNEvents_(0),
  processeff_(1.0),
  sampleChain_(NULL),
  firstrun_(0),
  lastrun_(0),
  effCorrFactor_(1.0),
  histFile_(NULL),
  dataType_(""),
  color_(0),
  lstyle_(0),
  plotOrder_(0),
  init_(0)
{
  
}

Sample::~Sample(){
  
  if(sampleChain_!=NULL)delete sampleChain_;
  if(histFile_!=NULL)  delete histFile_;
  
}


bool Sample::init(){
  if(init_)return 1;

  //read in the processing efficiency
  struct stat stFileInfo;
  if(stat((const char*)(TString(GetTitle())+"/efficiency.txt"),&stFileInfo) != 0) {
    cout<<"unable to find "<<TString(GetTitle())+"/efficiency.txt"<<endl;
    return 0;
  }
  ifstream input;
  input.open((const char*)(TString(GetTitle())+"/efficiency.txt"));
  TString effname="";
  input>>effname>>processeff_;
  if(effname!="efficiency" || processeff_<0.1 || processeff_ > 1.001 ){
    cout<<" processing efficiency read is invalid"<<endl;
    return 0;
  }

  //print the trigger paths  
  cout<<"TriggerPaths:"<<endl;
  for(std::vector<std::string>::const_iterator trig=trigPaths_.begin(); trig!=trigPaths_.end(); trig++)
    cout<<trig->c_str()<<endl;
  

  return init_=1;
}


fwlite::ChainEvent*  Sample::getEvents(){
  
  if(!sampleChain_){
    //check if the summary file exists
    struct stat stFileInfo;
    if(stat((const char*)(TString(GetTitle())+"/summary.txt"),&stFileInfo) != 0) {
      cout<<"unable to find "<<TString(GetTitle())+"/summary.txt"<<endl;
      return NULL;
    }
    //read in the files 
    ifstream input;
    input.open((const char*)(TString(GetTitle())+"/summary.txt"));
    int good=0; string colln; int nevt=0; int nevtpass=0; int N=-1;
    input>>good>>colln>>nevt>>nevtpass;
    while(colln != "total" && N <10000){
      //test the files
      if(good==1){
	TFile ftest(colln.c_str(),"read");
	if(ftest.IsZombie()){
	  cout<<colln<<"  is declared good but is Zombie"<<endl;
	  good=0;
	}
	if(!ftest.GetListOfKeys()){
	  cout<<colln<<"  is declared good but does not have Keys"<<endl;
	  good=0;
	}
	TTree*tree=(TTree*)ftest.Get("Events");
	if(!tree){
	  cout<<colln<<"  is declared good but does not have Events"<<endl;
	  good=0;
	}
      
	//add the good files to the chain
	if(good==1){
	  sampleList_.push_back(colln);
	  nEvents_+=nevt;
	  chainNEvents_+=tree->GetEntriesFast();
	}
      }

      input>>good>>colln>>nevt>>nevtpass;
      N++;
    }
    if(N==10000){cout<<GetName()<<" Too many files"<<endl; return 0;}
    cout<<GetName()<<" "<<sampleList_.size()<<" files added for "<<nEvents_<<" events, with "<<chainNEvents_<<" chain events."<<endl;
    input.close();

    if(sampleList_.size()==0)return NULL;
    sampleChain_=new fwlite::ChainEvent(sampleList_);
  }

 
  if(!sampleChain_)cout<<"Error event chain not set"<<endl;
  return sampleChain_;
}


bool Sample::save(TFile* file){
  
  if(sampleHist_.size()==0){
    cout<<GetName()<<": no histos to save"<<endl; return 0;
  }

  TFile* fi;
  if(file)fi=file;
  else {
    system(TString("rm -f ")+GetTitle()+"/Sample_Histograms.root");
    fi=new TFile(TString(GetTitle())+"/Sample_Histograms.root","recreate");    
    if(fi->IsZombie()){
      cout<<"Sample "<<GetName()<<" Unable to create histograms file"<<endl;
      return 0;
    }
  }
    
  fi->cd();
  std::vector<TH1*>::const_iterator h=sampleHist_.begin();
  for( int i=0; h!=sampleHist_.end(); ++h, i++){
    cout<<"save: "<<i<<" "<<sampleHist_[i]->GetName()<<endl;
    sampleHist_[i]->Write();
  }

  cout<<"Sample "<<GetName()<<" saved histograms in "<<fi->GetName()<<endl;

  
  if(!file){
    fi->Close();
    delete fi;
  }
 
  return 1;
}



void Sample::print(){

  cout<<"Sample "<<GetName()<<" "<<this<<endl;
  cout<<"  path: "<<GetTitle()<<endl;
  cout<<"  NEvents: "<<nEvents_<<endl;

  std::vector<TH1*>::const_iterator h=sampleHist_.begin();
  for(int i=0; h!=sampleHist_.end(); ++h, i++)
    cout<<"  histo"<<i<<": "<<sampleHist_[i]->GetName()<<endl;


  
}
