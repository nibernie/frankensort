#ifndef __MATRICESHELPER_H__
#define __MATRICESHELPER_H__


#include <iostream>
#include <iomanip>
#include <utility>
#include <vector>
#include <cstdio>
#include <sys/stat.h>

#include "TROOT.h"
#include "TTree.h"
#include "TTreeIndex.h"
#include "TVirtualIndex.h"
#include "TChain.h"
#include "TFile.h"
#include "TList.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TCutG.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "TGRSIRunInfo.h"
#include "TGRSISortInfo.h"
#include "Globals.h"
#include "TVectorD.h"
#include "TPPG.h"

#include "TGRSIUtilities.h"




//global helpers.

TGRSIRunInfo *info   =0;
TPPG         *ppg    =0;
TChain       *gChain =0;
TFile       *outfile =0;
TList       *outlist =0;



TCutG *prompt=0;
TCutG *gb_prompt=0;
TCutG *timerandom=0;

///////////////////////////////////// SETUP ///////////////////////////////////////
//Histogram paramaters
Double_t low     = 0;
Double_t high    = 7000;
Double_t nbins   = 14000;


//Coincidence Parameters
Double_t ggTlow  = -350.;   //Times are in  ns
Double_t ggThigh = 250.;
Double_t ggBGlow  = 600.; //abs
Double_t ggBGhigh = 2000.; //abs
Double_t ggBGScale = (ggThigh - ggTlow)/(2*(ggBGhigh - ggBGlow));


Double_t gbTlow  =  -200.;
Double_t gbThigh = 300.;
Double_t gbBGlow  = 600.; //abs
Double_t gbBGhigh = 800.; //abs
Double_t gbBGScale = (gbThigh - gbTlow)/(2*(gbBGhigh - gbBGlow));


//Double_t beta_threshold = 800.;
Double_t beta_threshold = 200000.;

//this is in ms
//Double_t cycleBins = 14.5e8 ; 
//Int_t cycleBins = 13500;
Int_t cycleBins = 13500;
//if(ppg) {
//   cycleBins = ppg->GetCycleLength()/1e5;
//   printf("cycle length = %.02f/n",cycleBins);
//}

Int_t    cycleLength = 13.5e8; 
Double_t bgStart  =  1.0e8;  // 1 second tape move.
Double_t bgEnd    =  1.5e8;  // 0.5 sec background.
Double_t onStart  =  1.5e8;  //
Double_t onEnd    = 11.5e8;  // 10 sec beam on.
Double_t offStart = 11.5e8;  //
Double_t offEnd   = 13.5e8;  // 2 sec beam off.



double FindSceptarOffset(int id) {
  switch(id) {
    case 64:
      return 2.3;
    case 65:
      return -1.26;
    case 66:
      return 0.87;
    case 67:
      return 1.69;
    case 68:
      return 4.33;
    case 69:
      return 9.09;
    case 70:
      return 2.80;
    case 71:
      return 2.42;
    case 72:
      return -0.03;
    case 73:
      return 1.70;
    case 74:
      return 2.80;
    case 75:
      return -0.82;
    case 76:
      return 1.30;
    case 77:
      return 0.07;
    case 78:
      return -1.86;
    case 79:
      return -1.97;
  };
  return 0.0;
}


bool CheckBetaGammaTime(double time1,double time2, int type) {
  double timediff = time1-time2;
  if(type==0) { // prompt 
   if((timediff>=gbTlow) && (timediff<=gbThigh)) 
     return true;
  } else {
   if((timediff>=gbBGlow) &&(timediff<=gbBGhigh))
     return true;
 }
 return false;
}

//bool CheckGammaGammaTime(double time1,double time2, int type) {
bool CheckGammaGammaTime(TGriffinHit *hit1, TGriffinHit *hit2, int type) {
  double timediff;
  if(type==0) { //looking for a prompt
    if(hit1->GetEnergy()>hit2->GetEnergy()) {
      timediff = hit1->GetTime() - hit2->GetTime();
    } else {
      timediff = hit2->GetTime() - hit1->GetTime();
    }
    if((timediff>=ggTlow) && (timediff<=ggThigh)) 
      return true;
  } else {
    timediff = fabs(hit1->GetTime()-hit2->GetTime());
    if((timediff>=ggBGlow) &&(timediff<=ggBGhigh))
      return true;
  }
  return false;
}

int CheckCycle(TGriffinHit *hit) {

  // Double_t cycleLength = 13.5e8; 
  // Double_t bgStart  =  1.0e8;  // 1 second tape move.
  // Double_t bgEnd    =  1.5e8;  // 0.5 sec background.
  // Double_t onStart  =  1.5e8;  //
  // Double_t onEnd    = 11.5e8;  // 10 sec beam on.
  // Double_t offStart = 11.5e8;  //
  // Double_t offEnd   = 13.5e8;  // 2 sec beam off.

  double check = hit->GetTimeStamp() % cycleLength;
 
  if(check < bgStart) {
     return 1;  //tape move
  } 
  if(check < bgEnd) {
     return 2;  //background
  } 
  if(check < onEnd) {
     return 3; // beam on
  }
  if(check < offEnd) {
     return 4; // beam off
  }
 
  return 256; // unknown?
}

bool CheckCycleLess(TGriffinHit *hit) {

  double check = hit->GetTimeStamp() % cycleLength;

  double on  =  1.5e8;  //  start beam on
  double off  =  4.5e8;  //  3 sec beam on

  if((check>on) && (check<off))  {
     return true;
  }
  return false;
}


















TGriffin* griffin = 0;
TSceptar* sceptar = 0;

long GetEntry(long entry,TChain *chain) {
  griffin->Clear();
  sceptar->Clear();
  return chain->GetEntry(entry);
}









std::vector<std::string> RootFiles;
std::vector<std::string> CalFiles;
std::vector<std::string> InfoFiles;




void Subtract(const char *mat1,const char* mat2,const char *name,const char* title,double scale) {
   TH1 *m1 = (TH1*)outlist->FindObject(mat1);
   TH1 *m2 = (TH1*)outlist->FindObject(mat2);
   TH1* m3 = (TH1*)m1->Clone(name);
   m3->SetNameTitle(name,title);
   m3->Add(m2,-scale);
   outlist->Add(m3); 
}




void FillHist(const char *name,const char *title,int xbin,double xlow, double xhigh,double xvalue,
              int ybin=0,double ylow=0.0,double yhigh=0.0,double yvalue=0.0) {
   if(!outlist)
      outlist = new TList;
   TH1 *hist = (TH1*)outlist->FindObject(name);
   if(!hist) {
      if(ybin) 
         hist = new TH2D(name,title,xbin,xlow,xhigh,ybin,ylow,yhigh);
      else
         hist = new TH1D(name,title,xbin,xlow,xhigh);
      outlist->Add(hist);
   }
   if(ybin)
      hist->Fill(xvalue,yvalue);
   else
      hist->Fill(xvalue);
}











void FillHistSym(const char *name,const char *title,int xbin,double xlow, double xhigh,double xvalue,
      int ybin,double ylow,double yhigh,double yvalue) {
   if(!outlist)
      outlist = new TList;
   TH1 *hist = (TH1*)outlist->FindObject(name);
   if(!hist) {
      if(ybin) 
         hist = new TH2D(name,title,xbin,xlow,xhigh,ybin,ylow,yhigh);
      else
         hist = new TH1D(name,title,xbin,xlow,xhigh);
      outlist->Add(hist);
   }
   hist->Fill(xvalue,yvalue);
   hist->Fill(yvalue,xvalue);
}


int Trigger(TGriffin *griffin,TSceptar *sceptar) {
   int gsize = griffin->GetMultiplicity();
   int bsize = sceptar->GetMultiplicity();
   if(gsize==1 && bsize==0) {
      return 1;
   } else if (gsize>1 && bsize==0) {
      return 2;       
   } else if (gsize==1 && bsize>0) {
      return 4;
   } else if (gsize>1 && bsize>0) {
      return 6;
   }
   return 9;
}


std::vector<TGriffinHit> DoAddback(std::vector<TGriffinHit>* hits,TCutG *timegate,double distance) {
  // I am assuming the hits have been ordered in energy, largest to smallest.
  std::vector<TGriffinHit>::iterator it1;
  std::vector<TGriffinHit>::iterator it2;
  std::vector<TGriffinHit> addback;
  //std::vector<TGriffinHit> addback2;
  for(size_t i=0;i<hits->size();i++) {
    TGriffinHit hit(hits->at(i));
    addback.push_back(hit);
  }
  for(it1=addback.begin();it1!=addback.end();) {
    for(it2=it1+1;it2!=addback.end();) {
      //lets check the time gate.
      if(!timegate->IsInside(it1->GetTimeStamp()-it2->GetTimeStamp(),it2->GetEnergy())) {
        if(it1->GetTimeStamp()!=it2->GetTimeStamp()) {
          it2++;
          continue;
        }
      }
      double s = (it1->GetPosition() - it2->GetPosition()).Mag();
      if(s>distance) {
        it2++;
        continue;
      }
      it1->Add(&(*it2));
      addback.erase(it2);
      
      //addback2.push_back(*it1);
    }
    it1++;
  }
  return addback;
}


class Notifier : public TObject {
   public:
      Notifier() { }
      ~Notifier() { }

      void AddChain(TChain *chain)       { fChain = chain; }
      void AddRootFile(std::string name) { RootFiles.push_back(name); }
      void AddInfoFile(std::string name) { InfoFiles.push_back(name); }
      void AddCalFile(std::string name)  { CalFiles.push_back(name); }

      bool Notify() { 
         printf("%s loaded.\n",fChain->GetCurrentFile()->GetName());
         //fChain->GetCurrentFile()->ls();
         ppg = (TPPG*)fChain->GetCurrentFile()->Get("TPPG");

         //if(!TChannel::GetNumberOfChannels())   
         //  TChannel::ReadCalFile("Cd128.cal");
         //if(CalFiles.size()>0)   
         //   TChannel::ReadCalFile(CalFiles.at(0).c_str());
         for(int i=0;i<CalFiles.size();i++) 
            TChannel::ReadCalFile(CalFiles.at(i).c_str());

         if(InfoFiles.size()>0)
            TGRSIRunInfo::ReadInfoFile(InfoFiles.at(0).c_str());       

         return true;
      };
   private:
      TChain *fChain;
      std::vector<std::string> RootFiles;
      std::vector<std::string> CalFiles;
      std::vector<std::string> InfoFiles;
};


Notifier *notifier = new Notifier;
void OpenRootFile(std::string filename) {
   TFile f(filename.c_str());
   if(f.Get("AnalysisTree")) {
      if(!gChain) {
         gChain = new TChain("AnalysisTree");
         notifier->AddChain(gChain);
         gChain->SetNotify(notifier);
      }
      gChain->Add(filename.c_str());
      std::cout << "added: " << filename << std::endl;
   }
   if(f.Get("prompt")) {
      prompt = (TCutG*)f.Get("prompt");
      std::cout << "found prompt banana." << std::endl;
   }
   if(f.Get("gb_prompt_gate")) {
      gb_prompt = (TCutG*)f.Get("gb_prompt_gate");
      std::cout << "found gb_prompt banana." << std::endl;
   }
   if(f.Get("timerandom")) {
      timerandom = (TCutG*)f.Get("timerandom");
      std::cout << "found time random banana." << std::endl;
   }

}



void AutoFileDetect(std::string filename)  {
   size_t dot_pos = filename.find_last_of('.');
   std::string ext = filename.substr(dot_pos+1);
   if(ext=="root") {
      //RootFiles.push_back(filename);
      OpenRootFile(filename);
   } else if(ext=="cuts") {
     OpenRootFile(filename);
   } else if(ext=="cal") {
      notifier->AddCalFile(filename);
   } else if(ext=="info") {
      notifier->AddInfoFile(filename);
   } else {
      fprintf(stderr,"discarding unknown file: %s\n",filename.c_str());
   }
}


#endif
