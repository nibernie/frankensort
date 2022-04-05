//g++ doit11.cxx `grsi-config --cflags --all-libs --root` -O3 -lMathMore -ltbb -odo_it11

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
//#include "TRunInfo.h" // GRSISort4
#include "TGRSISortInfo.h"
#include "Globals.h"
#include "TVectorD.h"
#include "TPPG.h"

#include "TGRSIUtilities.h"

#ifndef __CINT__ 
#include "TGriffin.h"
#include "TSceptar.h"
#endif

//This code is an example of how to write an analysis script to analyse an analysis tree
//The compiled script works like this
//
//  1. Starts in the main function by finding the analysis tree, setting up input and output files
//  2. Calls the LeanMatrices function
//  3. LeanMatrices creates 1D and 2D histograms and adds them to a list
//  4. Some loops over event "packets" decides which histograms should be filled and fills them.
//  5. The list of now filled histograms is returned to the main function
//  6. The list is written (meaning all of the histograms are written) to the output root file
//  7. Papers are published, theses are granted, high-fives are made
//
/////////////////////////////////////////////////////////////////////////////////////////


//global helpers.

TGRSIRunInfo *info   =0;
TPPG         *ppg    =0;
TChain       *gChain =0;
TFile       *outfile =0;
TList       *outlist =0;

///////////////////////////////////// SETUP ///////////////////////////////////////
//Histogram paramaters
Double_t low     = 0;
Double_t high    = 5000;
Double_t nbins   = 5000;

//Coincidence Parameters
Double_t ggTlow  = 0.;   //Times are in 10 ns (timestamp)
Double_t ggThigh = 80.;
Double_t ggBGlow  = 100.; 
Double_t ggBGhigh = 180.; 
Double_t ggBGScale = 1.;	//(ggThigh - ggTlow)/(ggBGhigh - ggBGlow);


Double_t gbTlow  =  0.;
Double_t gbThigh = 80.;
Double_t gbBGlow  = 100.; 
Double_t gbBGhigh = 180.; 
Double_t gbBGScale = 1.;	//(gbThigh - gbTlow)/(gbBGhigh - gbBGlow);

//Double_t beta_threshold = 800.;
Double_t beta_threshold = 200000.;

//this is in ms
//Double_t cycleBins = 14.5e8 ; 
//Int_t cycleBins = 13500;
Int_t cycleBins = 13500.;
//if(ppg) {
//   cycleBins = ppg->GetCycleLength()/1e5;
//   printf("cycle length = %.02f/n",cycleBins);
//}

Long_t   cycleLength = 13.5e8; 
Double_t bgStart  =  1.0e8;  // 1 second tape move.
Double_t bgEnd    =  1.5e8;  // 0.5 sec background.
Double_t onStart  =  1.5e8;  //
Double_t onEnd    = 11.5e8;  // 10 sec beam on.
Double_t offStart = 11.5e8;  //
Double_t offEnd   = 13.5e8;  // 2 sec beam off.


void Subtract(const char *mat1,const char *mat2,const char *name,const char *title,double scale) {
  TH1 *m1 = (TH1*)outlist->FindObject(mat1);
  TH1 *m2 = (TH1*)outlist->FindObject(mat2);
  TH1 *m3 = (TH1*)m1->Clone(name);
  m3->SetNameTitle(name,title);
  m3->Add(m2,-scale);
  outlist->Add(m3); 
}


//TCutG *prompt=0;
//TCutG *gb_prompt=0;
//TCutG *timerandom=0;


std::map<int,unsigned long> glasttime;
std::map<int,unsigned long> blasttime;

std::map<int,double> glastenergy;

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

void FillHist(const char *name,const char *title,int xbin,double xlow, double xhigh,const char *xvalue,
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
    ((TH2*)hist)->Fill(xvalue,yvalue,1);
  else
    hist->Fill(xvalue,1);
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

/*
//From mymatrices_12_15.cxx
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

bool CheckGammaGammaTime(double time1,double time2, int type) {
  double timediff = fabs(time1-time2);
  if(type==0) { // prompt 
   if((timediff>=ggTlow) && (timediff<=ggThigh)) 
     return true;
  } else {
   if((timediff>=ggBGlow) &&(timediff<=ggBGhigh))
     return true;
 }
 return false;
}
*/






int Trigger(TGriffin *griffin,TSceptar *sceptar) {

  if(sceptar==0) return 1;
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


bool CheckCycleLess(TGriffinHit *hit) {

  double check = hit->GetTimeStamp() % cycleLength;

  double on  =  1.5e8;  //  start beam on
  double off  =  4.5e8;  //  3 sec beam on

  if((check>on) && (check<off))  {
    return true;
  }
  return false;
}


int CheckCycle(TGriffinHit *hit) {

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


/*
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
*/



bool UseSceptar = false;


//TList *LeanMatrices(TTree* tree, TPPG* ppg, TGRSIRunInfo* runInfo, long maxEntries = 0, TStopwatch* w = NULL) {
void MakeHistos(long maxEntries=0) {

  TGriffin* griffin = 0;
  TSceptar* sceptar = 0;
  gChain->SetBranchAddress("TGriffin", &griffin); //We assume we always have a Griffin
  if(gChain->FindBranch("TSceptar")) {
    UseSceptar = true;
    gChain->SetBranchAddress("TSceptar", &sceptar); //We assume we always have a Griffin
  }

  //store the last timestamp of each channel
  std::vector<long> lastTimeStamp(65,0);

  std::cout<<std::fixed<<std::setprecision(1); //This just make outputs not look terrible
  //size_t angIndex;
  if(maxEntries == 0 || maxEntries > gChain->GetEntries()) {
    maxEntries = gChain->GetEntries();
  }

  long entry=0;
  for(entry = 0; entry < maxEntries; entry++) { //Only loop over the set number of entries
    griffin->Clear();
    if(UseSceptar) sceptar->Clear();

    gChain->GetEntry(entry);

    int trigger = Trigger(griffin,sceptar);

    std::vector<double> betatimes;
    std::vector<int> betachans;
    double betatime=-1;
    int    betachan=-1;
    bool   GoodBeta=false;
    if(sceptar) {
      FillHist("gamma_beta_multi","#gamma multi vs #beta mutli",70,0,70,griffin->Size(),
          20,0,20,sceptar->Size());
      for(int bhit=0;bhit<sceptar->GetMultiplicity();bhit++) {
        int chan  =sceptar->GetSceptarHit(bhit)->GetChannel()->GetNumber();
        double chg=sceptar->GetSceptarHit(bhit)->Charge();
        //double tim=sceptar->GetSceptarHit(bhit)->GetTimeStamp() - sceptar->GetSceptarHit(bhit)->GetTime();
        FillHist("BetaSingles_All","BetaSingles",100,0,100,chan,
            8000,0,5120000,chg);
        //FillHist(Form("BetaSingles_%i",chan),Form("BetaSingles_%i",chan),8000,0,5120000,chg);

        double cycletime = sceptar->GetSceptarHit(bhit)->GetTimeStamp() % cycleLength;

        FillHist("betaSinglesCycle_All","#beta vs tape ;charge;cycle time",
            2000,0,5120000,chg,cycleBins,0,cycleLength,cycletime);

        if(chg>beta_threshold) {
          betatimes.push_back(sceptar->GetSceptarHit(bhit)->GetTimeStamp());
          betachans.push_back(sceptar->GetSceptarHit(bhit)->GetChannel()->GetNumber());
          FillHist("BetaSingles_NotClean","BetaSingles above threshold",100,0,100,chan,
              8000,0,5120000,chg);


          FillHist("betaSinglesCycle_NotClean","#beta above threshold vs tape ;charge;Cycle Time [s]",
              2000,0,5120000,chg,cycleBins,0,cycleLength/100000000.,cycletime/100000000.);
 
      if(sceptar->GetMultiplicity()==1) { 
         FillHist("betaSinglesCycle_Clean","#beta above threshold minus doubles vs tape ;charge;Cycle Time [s]",
              2000,0,5120000,chg,cycleBins,0,cycleLength/100000000.,cycletime/100000000.);
 		}


       }


        if(blasttime.count(sceptar->GetSceptarHit(bhit)->GetDetector() )) {
          FillHist("betadead","beta dead time:time between triggers:sceptar det",
              10000,0,10000,sceptar->GetSceptarHit(bhit)->GetTimeStamp() - blasttime[sceptar->GetSceptarHit(bhit)->GetDetector()],
              20,0,20,sceptar->GetSceptarHit(bhit)->GetDetector());
        }
        blasttime[sceptar->GetSceptarHit(bhit)->GetDetector()] = sceptar->GetSceptarHit(bhit)->GetTimeStamp();


      }

      if(betatimes.size()==1) {
        GoodBeta=true;
        betatime=betatimes.at(0) - FindSceptarOffset(betachans.at(0));
        betachan=betachans.at(0);
      }

      if(sceptar->GetMultiplicity()>0)
        FillHist("beta_multi","beta_multi",20,0,20,sceptar->GetMultiplicity());
      if(GoodBeta)     
        FillHist("beta_multi_th","beta_multi_th",20,0,20,betatimes.size());

    }
    FillHist("Trigger","",2000,0,2000,trigger);

    for(int ghit1=0;ghit1<griffin->GetMultiplicity();ghit1++) {
      TGriffinHit *hit1 = griffin->GetGriffinHit(ghit1);




      FillHist("Dirty_gammaSingles","#gamma singles;Energy [keV];xtal id",
          nbins,low,high,hit1->GetEnergy(),70,0,70,hit1->GetArrayNumber());
      if(GoodBeta) {
        FillHist("Dirty_gammaSinglesB","#gamma singles;Energy [keV];xtal id",
            nbins,low,high,hit1->GetEnergy(),70,0,70,hit1->GetArrayNumber());
      }
      FillHist("Dirty_kvalue","",2048,0,2048,hit1->GetKValue());
      if(hit1->GetKValue()==649) {
        FillHist("Dirty_singles_nopileup","no pileup",16000,0,8000,hit1->GetEnergy());
      }else if(hit1->GetKValue()<649) {
        FillHist("Dirty_singles_pileup","pileup",16000,0,8000,hit1->GetEnergy());
      }else {
        FillHist("Dirty_singles_garbage","garbage?",16000,0,8000,hit1->GetEnergy());
      }
    }

    griffin->FixLowGainCrossTalk();


    for(int ghit1=0;ghit1<griffin->GetMultiplicity();ghit1++) {
      TGriffinHit *hit1 = griffin->GetGriffinHit(ghit1);


      FillHist("gamma_summary_wpu","summary;Energy[keV];crystal Id",8000,0,4000,hit1->GetEnergy(),
          70,0,70,hit1->GetArrayNumber());
    } 

    griffin->CleanHits(649,511);
    //std::vector<TGriffinHit> addback = DoAddback(griffin->GetHitVector(TGriffin::kLowGain),prompt,70.);

    for(int ghit1=0;ghit1<griffin->GetMultiplicity();ghit1++) {
      TGriffinHit *hit1 = griffin->GetGriffinHit(ghit1);

      FillHist("gamma_summary","summary;Energy[keV];crystal Id",8000,0,4000,hit1->GetEnergy(),
          70,0,70,hit1->GetArrayNumber());


      FillHist("gamma_singles_charge","charge vs xtal name",70,0,70,hit1->GetChannel()->GetName(),
          16000,0,16000,hit1->GetCharge());


      FillHist("gamma_singles_energy","energy vs xtal name",70,0,70,hit1->GetChannel()->GetName(),
          8000,0,4000,hit1->GetEnergy());

      if((hit1->GetTimeStamp()/1e8)<500) {
        FillHist("gammaSingles_500s","#gamma singles for first 500 seconds;Energy [keV];xtal id",
            nbins,low,high,hit1->GetEnergy(),70,0,70,hit1->GetArrayNumber());
      }

      FillHist("gamma_runtime","#gamma xtal vs runtime:seconds:det number",7200,0,7200,hit1->GetTimeStamp()/1e8,
          70,0,70,hit1->GetArrayNumber());



      double cycletime = hit1->GetTimeStamp() % cycleLength;


      FillHist("gammaSingles_cycle","#gamma singles;Energy [keV];Cycle Time [s]",
          nbins,low,high,hit1->GetEnergy(),cycleBins,0,cycleLength/100000000,cycletime/100000000);

      FillHist(Form("gammaSingles_cyclepart_%i",CheckCycle(hit1)),"#gamma singles;Energy [keV]",
          nbins,low,high,hit1->GetEnergy());



      FillHist("gammaSingles","#gamma singles;Energy [keV];xtal id",
          nbins,low,high,hit1->GetEnergy(),70,0,70,hit1->GetArrayNumber());
      FillHist("kvalue","",2048,0,2048,hit1->GetKValue());
      if(hit1->GetKValue()==649) {
        FillHist("singles_nopileup","no pileup",16000,0,8000,hit1->GetEnergy());
      }else if(hit1->GetKValue()<649) {
        FillHist("singles_pileup","pileup",16000,0,8000,hit1->GetEnergy());
      }else {
        FillHist("singles_garbage","garbage?",16000,0,8000,hit1->GetEnergy());
      }

      if(glasttime.count(hit1->GetArrayNumber())) {
        FillHist("dead_crystal","dead time:time between triggers:xtal id",
            10000,0,10000,hit1->GetTimeStamp() - glasttime[hit1->GetArrayNumber()],
            70,0,70,hit1->GetArrayNumber());
        FillHist("dead_last_energy","dead time:time between triggers:energy",
            10000,0,10000,hit1->GetTimeStamp() - glasttime[hit1->GetArrayNumber()],
            nbins,low,high,glastenergy[hit1->GetArrayNumber()]);
      }
      glasttime[hit1->GetArrayNumber()] = hit1->GetTimeStamp();
      glasttime[hit1->GetArrayNumber()] = hit1->GetEnergy();

      for(int ghit2=ghit1;ghit2<griffin->GetMultiplicity();ghit2++) {
        TGriffinHit *hit2 = griffin->GetGriffinHit(ghit2);
        if(hit1->GetTimeStamp()!=hit2->GetTimeStamp()) {
          FillHist("ggtime_geng","#gamma 1 time - #gamma 2 time;Time Difference [10 ns]",  
              1000,-500,500,(hit1->GetTimeStamp()-hit2->GetTimeStamp()),
              9000,0,9000,hit2->GetEnergy());
        }

	if((fabs(hit1->GetTimeStamp()-hit2->GetTimeStamp())>=ggTlow) && (fabs(hit1->GetTimeStamp()-hit2->GetTimeStamp())<=ggThigh)) {
	  FillHistSym(Form("gg_matrix_cyclepart_%i_prompt",CheckCycle(hit1)),
                Form("gamma-gamma matrix, prompt time for cycle part %i",CheckCycle(hit1)),
                nbins,low,high,hit1->GetEnergy(),
                nbins,low,high,hit2->GetEnergy());

        double ang = (hit1->GetPosition().Angle(hit2->GetPosition()))*TMath::RadToDeg();  // the max ang can be here is 180.
        if(ang>174) {
          FillHistSym("gg_matrix_180","#gamma matrix, prompt time",
              nbins,low,high,hit1->GetEnergy(),
              nbins,low,high,hit2->GetEnergy());
          if(hit1->GetArrayNumber()>16 && hit2->GetArrayNumber()>16) {
            FillHistSym("gg_matrix_180_17_65","#gamma matrix, prompt time",
                nbins,low,high,hit1->GetEnergy(),
                nbins,low,high,hit2->GetEnergy());

            }
           }
	  }


      }
    }


   for(int ahit1=0;ahit1<griffin->GetAddbackMultiplicity();ahit1++) {
      TGriffinHit *hit1 = griffin->GetAddbackHit(ahit1);
//      FillHist("addbackSingles_clover","#addback singles;Energy [keV];Counts per 1 keV",nbins,low,high,hit->GetEnergy(),70,0,70,hit->GetArrayNumber()); 
//}

 //   for(int i=0;i<addback.size();i++) {
 //     TGriffinHit hit1 = addback.at(i);


      double cycletime = hit1->GetTimeStamp() % cycleLength;


      FillHist("addbackSingles_cycle","#gamma singles;Energy [keV];Cycle Time [s]",
          nbins,low,high,hit1->GetEnergy(),cycleBins,0,cycleLength/100000000,cycletime/100000000);

      FillHist(Form("addbackSingles_cyclepart_%i",CheckCycle(hit1)),"#gamma singles;Energy [keV]",
          nbins,low,high,hit1->GetEnergy());

      if((hit1->GetTimeStamp()/1e8)<500) {
        FillHist("addbackSingles_500s","#gamma singles for first 500 seconds;Energy [keV];xtal id",
            nbins,low,high,hit1->GetEnergy(),70,0,70,hit1->GetArrayNumber());
      }


      FillHist("addbackSingles","addback singles;Energy [keV];xtal id",
                nbins,low,high,hit1->GetEnergy(),70,0,70,hit1->GetArrayNumber());

//      for(int j=i+1;j<addback.size();j++) {
//        TGriffinHit hit2 = addback.at(j);

      for(int ahit2=ahit1+1;ahit2<griffin->GetAddbackMultiplicity();ahit2++) {
        TGriffinHit *hit2 = griffin->GetAddbackHit(ahit2);

	if((fabs(hit1->GetTimeStamp()-hit2->GetTimeStamp())>=ggTlow) && (fabs(hit1->GetTimeStamp()-hit2->GetTimeStamp())<=ggThigh)) {

            FillHistSym(Form("aa_matrix_cyclepart_%i_prompt",CheckCycle(hit1)),
                Form("addback-addback matrix, prompt time for cycle part %i",CheckCycle(hit1)),
                nbins,low,high,hit1->GetEnergy(),
                nbins,low,high,hit2->GetEnergy());

        double ang = (hit1->GetPosition().Angle(hit2->GetPosition()))*TMath::RadToDeg();  // the max ang can be here is 180.
        if(ang>174) {
          FillHistSym("aa_matrix_180","addback matrix, prompt time",
              nbins,low,high,hit1->GetEnergy(),
              nbins,low,high,hit2->GetEnergy());
          }
	}

      }
    }



    if(GoodBeta ) {
      //continue;

      ///////////
      //  at this point in the code, we massaged the data.
      //  1) hits are now order in energy from largest to smallest.
      //  2) Pileups have been removed from the hits.
      //  3) sceptar above threshold and size == 1 is required.
      //  4) crosstalk corrections to the hpge has been handled (done in GetEnergy and the calfile) 
      //  5) we are NOT using walk correction -> i.e only GetTimeStamp; GetTime got fired.
      //  6) ....
      ///////////

      for(int ghit1=0;ghit1<griffin->GetMultiplicity();ghit1++) {
        TGriffinHit *hit1 = griffin->GetGriffinHit(ghit1);

        FillHist("gbtime_geng","#gamma time - #beta time;g-b Time Difference [10 ns]",  
            1000,-500,500,(hit1->GetTimeStamp()-betatime),
            9000,0,9000,hit1->GetEnergy());

//        if(!gb_prompt->IsInside(hit1->GetTimeStamp()-betatime,hit1->GetEnergy())) continue;
	if((fabs(hit1->GetTimeStamp()-betatime)>=gbTlow) && (fabs(hit1->GetTimeStamp()-betatime)<=gbThigh)) {

        double cycletime = hit1->GetTimeStamp() % cycleLength;

        FillHist("gammaSinglesB","#gamma singles;Energy [keV];xtal id",
            nbins,low,high,hit1->GetEnergy(),70,0,70,hit1->GetArrayNumber());

        FillHist("gammaSinglesB_cycle","#gamma singles;Energy [keV];Cycle Time [s]",
            nbins,low,high,hit1->GetEnergy(),cycleBins,0,cycleLength/100000000,cycletime/100000000);

        FillHist(Form("gammaSinglesB_cyclepart_%i",CheckCycle(hit1)),
            "#gamma singles;Energy [keV]",
            nbins,low,high,hit1->GetEnergy());


        FillHist("gbtime_xtalid","#gamma time - #beta time;time difference (10ns)",
            1000,-500,500,hit1->GetTimeStamp()-betatime,
            70,0,70,hit1->GetArrayNumber());
        FillHist("gbtime_scepid","#gamma time - #beta time;time difference (10ns)",  
            4000,-2000,2000,hit1->GetTimeStamp()-betatime,
            100,0,100,betachan);
        FillHist("gbtime_geng_gated","#gamma time - #beta time;g-b Time Difference [10 ns]",  
            1000,-500,500,(hit1->GetTimeStamp()-betatime),
            9000,0,9000,hit1->GetEnergy());

        //check if good beta-hit1 time  -> continue;

        for(int ghit2=ghit1+1;ghit2<griffin->GetMultiplicity();ghit2++) {
          TGriffinHit *hit2 = griffin->GetGriffinHit(ghit2);
//          if(!gb_prompt->IsInside(hit2->GetTimeStamp()-betatime,hit2->GetEnergy())) continue;
	  if((fabs(hit2->GetTimeStamp()-betatime)>=gbTlow) && (fabs(hit2->GetTimeStamp()-betatime)<=gbThigh)) {
          //check if good beta-hit2 time  -> continue;

          //check if good hit1-hit2 time  -> continue;
          FillHist("ggbtime_geng","#gamma 1 time - #gamma 2 time with good beta;Time Difference [10 ns]",  
              1000,-500,500,(hit1->GetTimeStamp()-hit2->GetTimeStamp()),
              9000,0,6000,hit2->GetEnergy());
//          if(prompt && prompt->IsInside(hit1->GetTimeStamp()-hit2->GetTimeStamp(),hit2->GetEnergy())) {
	  if((fabs(hit1->GetTimeStamp()-hit2->GetTimeStamp())>=ggTlow) && (fabs(hit1->GetTimeStamp()-hit2->GetTimeStamp())<=ggThigh)) {
            FillHist("ggbtime_geng_gated","#gamma 1 time - #gamma 2 time with good beta and time gate applied;Time Difference [10 ns]",  
                1000,-500,500,(hit1->GetTimeStamp()-hit2->GetTimeStamp()),
                9000,0,6000,hit2->GetEnergy());
            FillHistSym("gg_matrixB","#gamma matrix, prompt time, prompt #beta",
                nbins,low,high,hit1->GetEnergy(),
                nbins,low,high,hit2->GetEnergy());
            FillHistSym(Form("gg_matrixB_cyclepart_%i",CheckCycle(hit1)),
                Form("#gamma-gamma matrix, prompt time, prompt #beta for cycle part %i",CheckCycle(hit1)),
                nbins,low,high,hit1->GetEnergy(),
                nbins,low,high,hit2->GetEnergy());

            double ang = (hit1->GetPosition().Angle(hit2->GetPosition()))*TMath::RadToDeg();  // the max ang can be here is 180.
            if(ang>174) {
              FillHistSym("gg_matrixB_180","#gamma matrix, prompt time, prompt #beta",
                  nbins,low,high,hit1->GetEnergy(),
                  nbins,low,high,hit2->GetEnergy());
            	}
              }
            }
          }
        }
      }

//      for(int i=0;i<addback.size();i++) {
//        TGriffinHit hit1 = addback.at(i);

   for(int ahit1=0;ahit1<griffin->GetAddbackMultiplicity();ahit1++) {
      TGriffinHit *hit1 = griffin->GetAddbackHit(ahit1);

        FillHist("abtime_geng","#addback time - #beta time;b-g Time Difference [10 ns]",  
            1000,-500,500,(hit1->GetTimeStamp()-betatime),
            9000,0,9000,hit1->GetEnergy());
        FillHist("abtime_xtalid","#addack time - #beta time;time difference (10ns)",
            1000,-500,500,hit1->GetTimeStamp()-betatime,
            70,0,70,hit1->GetArrayNumber());

//        if(!gb_prompt->IsInside(hit1->GetTimeStamp()-betatime,hit1->GetEnergy())) continue;
	if((fabs(hit1->GetTimeStamp()-betatime)>=gbTlow) && (fabs(hit1->GetTimeStamp()-betatime)<=gbThigh)) {

        double cycletime = hit1->GetTimeStamp() % cycleLength;
        FillHist("addbackSinglesB","#gamma-addback singles with #beta;Energy [keV];xtal id",
            nbins,low,high,hit1->GetEnergy(),70,0,70,hit1->GetArrayNumber());

        FillHist("addbackSinglesB_cycle","addback singles;Energy [keV];Cycle Time [s]",
            nbins,low,high,hit1->GetEnergy(),cycleBins,0,cycleLength/100000000.,cycletime/100000000.);

        FillHist(Form("addbackSinglesB_cyclepart_%i_prompt",CheckCycle(hit1)),
            "#addback singles;Energy [keV]",
            nbins,low,high,hit1->GetEnergy());

        if(CheckCycleLess(hit1)) {
          FillHist("addbackSinglesB_first3sec","#gamma-addback singles with #beta in the first 3 seconds;Energy [keV];xtal id", 
              nbins,low,high,hit1->GetEnergy(),70,0,70,hit1->GetArrayNumber());
          }
	} 


//        for(int j=i+1;j<addback.size();j++) {
//          TGriffinHit hit2 = addback.at(j);

   for(int ahit2=ahit1+1;ahit2<griffin->GetAddbackMultiplicity();ahit2++) {
        TGriffinHit *hit2 = griffin->GetAddbackHit(ahit2);

          FillHist("aabtime_geng","#gamma 1 time - #gamma 2 time with good beta;Time Difference [10 ns]",  
              1000,-500,500,(hit1->GetTimeStamp()-hit2->GetTimeStamp()),
              9000,0,6000,hit2->GetEnergy());

//          if(!gb_prompt->IsInside(hit2->GetTimeStamp()-betatime,hit2->GetEnergy())) continue;
	  if((fabs(hit2->GetTimeStamp()-betatime)>=gbTlow) && (fabs(hit2->GetTimeStamp()-betatime)<=gbThigh)) {
//          if(prompt && prompt->IsInside(hit1->GetTimeStamp()-hit2->GetTimeStamp(),hit2->GetEnergy())) {
	  if((fabs(hit1->GetTimeStamp()-hit2->GetTimeStamp())>=ggTlow) && (fabs(hit1->GetTimeStamp()-hit2->GetTimeStamp())<=ggThigh)) {
            FillHist("aabtime_geng_gated","#gamma 1 time - #gamma 2 time with good beta;Time Difference [10 ns]",  
                1000,-500,500,(hit1->GetTimeStamp()-hit2->GetTimeStamp()),
                9000,0,6000,hit2->GetEnergy());
            FillHistSym("aa_matrixB","#addback-addback matrix, prompt time, prompt #beta",
                nbins,low,high,hit1->GetEnergy(),
                nbins,low,high,hit2->GetEnergy());
            if(CheckCycleLess(hit1)) {
              FillHistSym("aa_matrixB_first3sec","#addback-addback matrix, prompt time, prompt #beta in the first 3 secs of bbeam on",
                  nbins,low,high,hit1->GetEnergy(),
                  nbins,low,high,hit2->GetEnergy());
            }

            FillHistSym(Form("aa_matrixB_cyclepart_%i_prompt",CheckCycle(hit1)),
                Form("addback-addback matrix, prompt time, prompt #beta for cycle part %i",CheckCycle(hit1)),
                nbins,low,high,hit1->GetEnergy(),
                nbins,low,high,hit2->GetEnergy());

            double ang = (hit1->GetPosition().Angle(hit2->GetPosition()))*TMath::RadToDeg();  // the max ang can be here is 180.
            if(ang>174) {
              FillHistSym("aa_matrixB_180","#addback matrix, prompt time, prompt #beta",
                  nbins,low,high,hit1->GetEnergy(),
                  nbins,low,high,hit2->GetEnergy());
             }
	    }

          } 
        }
      }

    }  //end if good beta

//TIME RANDOM
    for(int ghit1=0;ghit1<griffin->GetMultiplicity();ghit1++) {
      TGriffinHit *hit1 = griffin->GetGriffinHit(ghit1);

	for(int ghit2=ghit1;ghit2<griffin->GetMultiplicity();ghit2++) {
        TGriffinHit *hit2 = griffin->GetGriffinHit(ghit2);

	if(((hit1->GetTimeStamp()-hit2->GetTimeStamp())>=ggBGlow) && ((hit1->GetTimeStamp()-hit2->GetTimeStamp())<=ggBGhigh)) {
            FillHistSym(Form("gg_matrix_cyclepart_%i_rand",CheckCycle(hit1)),
                Form("gamma-gamma matrix, rand time for cycle part %i",CheckCycle(hit1)),
                nbins,low,high,hit1->GetEnergy(),
                nbins,low,high,hit2->GetEnergy());
	  }
       }
    }

   for(int ahit1=0;ahit1<griffin->GetAddbackMultiplicity();ahit1++) {
      TGriffinHit *hit1 = griffin->GetAddbackHit(ahit1);

	for(int ahit2=ahit1+1;ahit2<griffin->GetAddbackMultiplicity();ahit2++) {
        TGriffinHit *hit2 = griffin->GetAddbackHit(ahit2);

//          if(prompt && prompt->IsInside(hit1->GetTimeStamp()-hit2->GetTimeStamp(),hit2->GetEnergy())) {
	  if(((hit1->GetTimeStamp()-hit2->GetTimeStamp())>=ggBGlow) && ((hit1->GetTimeStamp()-hit2->GetTimeStamp())<=ggBGhigh)) {

            FillHistSym(Form("aa_matrix_cyclepart_%i_rand",CheckCycle(hit1)),
                Form("addback-addback matrix, rand time for cycle part %i",CheckCycle(hit1)),
                nbins,low,high,hit1->GetEnergy(),
                nbins,low,high,hit2->GetEnergy());
	   }
          } 
       }

    if(GoodBeta ) {
   for(int ahit1=0;ahit1<griffin->GetAddbackMultiplicity();ahit1++) {
      TGriffinHit *hit1 = griffin->GetAddbackHit(ahit1);

//        if(!gb_prompt->IsInside(hit1->GetTimeStamp()-betatime,hit1->GetEnergy())) continue;
	if((fabs(hit1->GetTimeStamp()-betatime)>=gbBGlow) && (fabs(hit1->GetTimeStamp()-betatime)<=gbBGhigh)) {

        double cycletime = hit1->GetTimeStamp() % cycleLength;

        FillHist(Form("addbackSinglesB_cyclepart_%i_rand",CheckCycle(hit1)),
            "#addback singles;Energy [keV]",
            nbins,low,high,hit1->GetEnergy());
	} 

	for(int ahit2=ahit1+1;ahit2<griffin->GetAddbackMultiplicity();ahit2++) {
        TGriffinHit *hit2 = griffin->GetAddbackHit(ahit2);

//          if(!gb_prompt->IsInside(hit2->GetTimeStamp()-betatime,hit2->GetEnergy())) continue;
	  if((fabs(hit2->GetTimeStamp()-betatime)>=gbTlow) && (fabs(hit2->GetTimeStamp()-betatime)<=gbThigh)) {
//          if(prompt && prompt->IsInside(hit1->GetTimeStamp()-hit2->GetTimeStamp(),hit2->GetEnergy())) {
	  if((fabs(hit1->GetTimeStamp()-hit2->GetTimeStamp())>=ggBGlow) && (fabs(hit1->GetTimeStamp()-hit2->GetTimeStamp())<=ggBGhigh)) {

            FillHistSym(Form("aa_matrixB_cyclepart_%i_rand",CheckCycle(hit1)),
                Form("addback-addback matrix, prompt time, prompt #beta for cycle part %i",CheckCycle(hit1)),
                nbins,low,high,hit1->GetEnergy(),
                nbins,low,high,hit2->GetEnergy());
	   }
          } 
         }
	}
       }


    if((entry%15000)==0) {
      printf("Working on ...   %lu / %lu   entries        \r",entry,maxEntries);
      fflush(stdout);
    }
    //if(entry>1000000)
    //   break;
  }
  printf("Working... on   %lu / %lu   entries        \n",entry,maxEntries);
  fflush(stdout);


  Subtract("gg_matrix_cyclepart_3_prompt","gg_matrix_cyclepart_3_rand","gg_matrix_cyclepart_3_sub","gamma-gamma prompt-rand;Energy [keV];Energy [keV]",ggBGScale);

  Subtract("aa_matrix_cyclepart_3_prompt","aa_matrix_cyclepart_3_rand","aa_matrix_cyclepart_3_sub","addback-addback prompt-rand;Energy [keV];Energy [keV]",ggBGScale);

  Subtract("addbackSinglesB_cyclepart_3_prompt","addbackSinglesB_cyclepart_3_rand","addbackSinglesB_cyclepart_3_sub","beta-addback singles, prompt-rand beta;Energy [keV];Counts per 1 keV",gbBGScale);

  Subtract("aa_matrixB_cyclepart_3_prompt","aa_matrixB_cyclepart_3_rand","aa_matrixB_cyclepart_3_sub","addback-addback-beta prompt-rand;Energy [keV];Energy [keV]",ggBGScale);

  return;
}

std::vector<std::string> RootFiles;
std::vector<std::string> CalFiles;
std::vector<std::string> InfoFiles;

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
//  if(f.Get("prompt")) {
//    prompt = (TCutG*)f.Get("prompt");
//    std::cout << "found prompt banana." << std::endl;
//  }
//  if(f.Get("gb_prompt_gate")) {
//    gb_prompt = (TCutG*)f.Get("gb_prompt_gate");
//    std::cout << "found gb_prompt banana." << std::endl;
//  }
//  if(f.Get("timerandom")) {
//    timerandom = (TCutG*)f.Get("timerandom");
//    std::cout << "found time random banana." << std::endl;
//  }

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



#ifndef __CINT__
int main(int argc,char **argv) {
  for(int i=1;i<argc;i++) 
    AutoFileDetect(argv[i]); 

  if(!gChain || !gChain->GetEntries()) {
    fprintf(stderr,"failed to find stuff to sort, exiting.\n");
    return 0;
  }
  std::string fname = gChain->GetCurrentFile()->GetName();
  int run_number = GetRunNumber(fname.c_str());
  int sub_number = GetSubRunNumber(fname.c_str());

  printf("\nStarting %i with %i files \n\n",run_number,gChain->GetNtrees()); fflush(stdout);

  MakeHistos(); 

  if(outlist && outlist->GetEntries()>0) {
    TFile *file = new TFile(Form("doit_%05i_%03i.root",run_number,sub_number),"recreate");
    outlist->Sort();
    //outlist->Write();
    TIter iter(outlist);
    while(TObject *obj = iter.Next()) {
      printf("obj @ 0x%08x\t",obj); fflush(stdout);
      printf("writing %s\n",obj->GetName()); fflush(stdout);
      obj->Write();
    }
    file->Close();
  }

  return 0;
}

#endif


