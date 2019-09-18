
#include<cstdio>
#include<vector>

#include <TPRegexp.h>

#include<TMidasFile.h>
#include<TMidasEvent.h>
#include<TChannel.h>

#include<TFragment.h>

TFragment *gfrag=0;
TTree     *gtree=0;
TFile     *goutfile=0;









std::string get_run_number(std::string input) {
  int i = input.find("run");
  std::string out = input.substr(input.find("run")+3,9);
  printf(  "%i    %s\nn",i,out.c_str());
  return out;
}


void SetUpTree(const char *infile) {
  std::string rn = get_run_number(infile);
  goutfile = new TFile(Form("fragment%s.root",rn.c_str()),"recreate");
  gfrag    = new TFragment;
  gtree    = new TTree("FragmentTree","FragmentTree");

  gtree->Branch("TFragment",&gfrag);

}


uint32_t lastb=0xffffffff;

struct GRIF16_3 {
  uint32_t data8;
  uint32_t data9;
  uint32_t dataA;
  uint32_t dataB;
  uint32_t dataPH;
  uint32_t dataCFD;
  uint32_t dataE;

  //uint32_t data8;

  uint32_t Address()   { return  (data8&0x000ffff0)>>4;  }
  uint32_t Card()      { return  (data8&0x0e000000)>>25; }
  uint64_t Timestamp() { return  ((( ((uint64_t)dataB)&0x00003fff)<<28) + (((uint64_t)dataA)&0x0fffffff)); }

  uint32_t Counter1()  { return  (data9&0x0fffffff); }
  uint32_t Counter2()  { return  (dataE&0x00003fff); }

  uint32_t Integration() { return (((dataPH&0x7c000000)>>17) + ((dataCFD&0x7fc00000)>>22)); } 
  int32_t Charge()      { return  ((dataPH&0x03ffffff)); } 
  uint32_t Cfd()         { return  ((dataCFD&0x0003ffff)); } 

  void PrintRaw()        {
    //for(int x=0;x<data.size();x++) {
    printf("0x%08x  ",data8);
    printf("0x%08x  ",data9);
    printf("0x%08x  ",dataA);
    printf("0x%08x  ",dataB);
    printf("0x%08x  ",dataPH);
    printf("0x%08x  ",dataCFD);
    printf("0x%08x  ",dataE);
    //}
    TChannel *channel = TChannel::GetChannel(Address());
    if(channel) printf("%s\t%i\n",channel->GetName(),Integration());
    else        printf("\n");
  }
/*
  void Print()  {
    TChannel *channel = TChannel::GetChannel(Address());
    if(channel) printf("%s  ",channel->GetName());
    else        printf("          ");
    printf("0x%08x\t",Address());

    printf("%lu\t",Timestamp());
    printf("%i\t",Integration());
    printf("%i\t",Charge());
    printf("%i\n",Cfd());


  }
*/


  bool IsGood() {
    if(dataB==0xffffffff) {
      dataB=lastb;
    } else {
      lastb=dataB;
    }
    if(((data8&0xf0000000)==0x80000000) &&
        ((data9&0xf0000000)==0x90000000) &&
        ((dataA&0xf0000000)==0xa0000000) &&
        ((dataB&0xf0000000)==0xb0000000) &&
        //((data5&0xf0000000)==0xd0000000) &&
        ((dataE&0xf0000000)==0xe0000000)) return true;
    return false;
  }

};







struct GRIF16_2 {
  uint32_t data8;
  uint32_t data9;
  uint32_t dataA;
  uint32_t dataB;
  uint32_t dataPH;
  uint32_t dataCFD;
  uint32_t dataE;
  //uint32_t data8;

  uint32_t Address()   { return  (data8&0x000ffff0)>>4;  }
  uint32_t Card()      { return  (data8&0x03800000)>>23; }
  uint64_t Timestamp() { return  ((( ((uint64_t)dataB)&0x00003fff)<<28) + (((uint64_t)dataA)&0x0fffffff)); }

  uint32_t Counter1()  { return  (data9&0x0fffffff); }
  uint32_t Counter2()  { return  (dataE&0x00003fff); }

  uint32_t Integration() { return (((dataPH&0x7c000000)>>21) + ((dataCFD&0x7c000000)>>26)); } 
  int32_t Charge()      { return  ((dataPH&0x03ffffff)); } 
  uint32_t Cfd()         { return  ((dataCFD&0x03ffffff)); } 

  void PrintRaw()        {
    //for(int x=0;x<data.size();x++) {
    printf("0x%08x  ",data8);
    printf("0x%08x  ",data9);
    printf("0x%08x  ",dataA);
    printf("0x%08x  ",dataB);
    printf("0x%08x  ",dataPH);
    printf("0x%08x  ",dataCFD);
    printf("0x%08x  ",dataE);
    //}
    TChannel *channel = TChannel::GetChannel(Address());
    if(channel) printf("%s\n",channel->GetName());
    else        printf("\n");
  }

  void Print()  {
    TChannel *channel = TChannel::GetChannel(Address());
    if(channel) printf("%s  ",channel->GetName());
    else        printf("          ");
    printf("0x%08x\t",Address());

    printf("%lu\t",Timestamp());
    printf("%i\t",Integration());
    printf("%i\t",Charge());
    printf("%i\n",Cfd());


  }

  bool IsGood() {
    if(((data8&0xf0000000)==0x80000000) &&
        ((data9&0xf0000000)==0x90000000) &&
        ((dataA&0xf0000000)==0xa0000000) &&
        ((dataB&0xf0000000)==0xb0000000) &&
        //((data5&0xf0000000)==0xd0000000) &&
        ((dataE&0xf0000000)==0xe0000000)) return true;
    return false;
  }

};



uint32_t mserial = 0;
uint32_t fserial = 0;

TFragment MakeFragment(GRIF16_2 *grif,TMidasEvent *event) {
  TFragment fragment;

  fragment.SetAddress(grif->Address());
  fragment.SetTimeStamp(grif->Timestamp());
  fragment.SetKValue(grif->Integration());
  fragment.SetCharge(grif->Charge());
  fragment.SetCfd(grif->Cfd());

  fragment.SetModuleType(grif->Card()); 
  fragment.SetMidasId(event->GetSerialNumber());
  if(event->GetSerialNumber()!=mserial) {
    mserial = event->GetSerialNumber();
    fserial=0;
  }
  fragment.SetFragmentId(fserial++);
  fragment.SetMidasTimeStamp(event->GetTimeStamp());


  return fragment;
}

TFragment MakeFragment(GRIF16_3 *grif,TMidasEvent *event) {
  TFragment fragment;

  fragment.SetAddress(grif->Address());
  fragment.SetTimeStamp(grif->Timestamp());
  fragment.SetKValue(grif->Integration());
  fragment.SetCharge(grif->Charge());
  fragment.SetCfd(grif->Cfd());

  fragment.SetModuleType(grif->Card()); 
  fragment.SetMidasId(event->GetSerialNumber());
  if(event->GetSerialNumber()!=mserial) {
    mserial = event->GetSerialNumber();
    fserial=0;
  }
  fragment.SetFragmentId(fserial++);
  fragment.SetMidasTimeStamp(event->GetTimeStamp());


  return fragment;
}










void PrintVector(std::vector<uint32_t> data) {
  TChannel *channel = 0;//uint32_t address-1; 
  for(int x=0;x<data.size();x++) {
    if((data[x]&0xf000000)==0x8000000) {
      channel = TChannel::GetChannel( (data[x]&0x000ffff0)>>4 ); 
      break;
    }
  }

  for(int x=0;x<data.size();x++) { 
    printf("0x%08x  ",data[x]);
    if(((x+1)%8)==0) printf("\n");
  }
  if(channel) { 
    printf("%s\n",channel->GetName());
  } else {
    printf("\n");
  }
  printf("--------------------------------\n"); 
}


void HandleHit(std::vector<uint32_t> data,TMidasEvent *event) {
    PrintVector(data); 
    return;
  if(data.size()!=7) {
    //PrintVector(data); 
    return;  // 8 appears to be the "normal" size of a grif16 no pileup hit.
  }
  GRIF16_2 *grif = (GRIF16_2*)(data.data());
  if(!grif->IsGood()) {
    //PrintVector(data); 
    return;
  }
  //grif->PrintRaw();
  //grif->Print();
  //TFragment frag = MakeFragment(grif,event);
  //frag.Print();
  TFragment f = MakeFragment(grif,event);
  f.Copy(*gfrag);
  //gfrag->Print();
  gtree->Fill();
  gfrag->Clear();
}

void HandleHits(std::vector<GRIF16_3> grif,TMidasEvent *event) {
  for(size_t x=0;x<grif.size();x++) {
    if(!grif[x].IsGood()) continue;
    TFragment f = MakeFragment(&grif[x],event);
    f.Copy(*gfrag);
    //gfrag->Print();
//    if(f.GetKValue()==700) {
      gtree->Fill();
//    }
    gfrag->Clear();
    //grif[x].PrintRaw();
  }
  return;
}



int main(int argc, char **argv) {

  TMidasFile infile;
  TMidasEvent event;
  infile.Open(argv[1]);

  TChannel::ReadCalFile(argv[2]);  ///"/data1/griffin/nbernier/Cad/132Cd/lon/Cal04496.cal");
  SetUpTree(argv[1]);
  long counter =0;
  long lastmidasid=0;
  long totalbytes=0;
  while(1) {
    int bytes = infile.Read(&event);
    if(bytes<=0) {
      break;
    } else {
      totalbytes+=bytes;
    }

    if((event.GetEventId()&0xf000)==0x8000) continue; 

    int size=0;
    void *ptr=0;

    size = event.LocateBank(0,"GRF3",&ptr);
    if(size) {
      counter++;
      lastmidasid=event.GetSerialNumber();

      //event.SetBankList(); event.Print("all");
      //printf("\n\nstarting breakup\n\n"); 
      bool header=false;
      std::vector<uint32_t> data;
      uint32_t *dptr = (uint32_t*)ptr;
      uint32_t address=0;
      uint32_t card=0;

      std::vector<GRIF16_3> grif;


      //printf("--------------------------\n");
      for(int x=0;x<size;x++) {
        data.push_back((*(dptr+x)));
        if( ((*(dptr+x))&0xf0000000)==0x80000000) {
          if(grif.size() && x>2) {
            if(((*(dptr+x-2))&0x80000000)!=0x80000000) {
              grif.back().dataPH  = (*(dptr+x-2));
            } else {
              grif.back().dataPH  = 0;
            }
            if(((*(dptr+x-1))&0x80000000)!=0x80000000) {
              grif.back().dataCFD = (*(dptr+x-1));
            } else {
              grif.back().dataCFD = 0;
            }
          }

          GRIF16_3 g = { 0xffffffff,0xffffffff,0xffffffff,0xffffffff,
                         0xffffffff,0xffffffff,0xffffffff };
          grif.push_back(g);
          grif.back().data8 = (*(dptr+x));

        } else if( ((*(dptr+x))&0xf0000000)==0x90000000 && grif.size()) { 
          grif.back().data9 = (*(dptr+x));
        } else if( ((*(dptr+x))&0xf0000000)==0xa0000000 && grif.size()) { 
          grif.back().dataA = (*(dptr+x));
        } else if( ((*(dptr+x))&0xf0000000)==0xb0000000 && grif.size()) { 
          grif.back().dataB = (*(dptr+x));
        } else if( ((*(dptr+x))&0xf0000000)==0xe0000000 && grif.size()) {
          grif.back().dataE = (*(dptr+x));
          if(x>2) { 
            if(((*(dptr+x-2))&0x80000000)!=0x80000000) {
              grif.back().dataPH  = (*(dptr+x-2));
            } else {
              grif.back().dataPH  = 0;
            }
            if(((*(dptr+x-1))&0x80000000)!=0x80000000) {
              grif.back().dataCFD = (*(dptr+x-1));
            } else {
              grif.back().dataCFD = 0;
            }
          }      
          for(size_t x=0;x<grif.size()-1;x++) {
            grif[x].dataB = grif.back().dataB;
            grif[x].dataE = grif.back().dataE;
          } 
          //event.SetBankList(); event.Print( Form("all 0x%08x",grif.back().data8));
          HandleHits(grif,&event);



          //printf("-------------------------------\n");
          //PrintVector(data);

          //printf("\n");
          
          data.clear();
          grif.clear();

        }
        else if( ((*(dptr+x))&0xf0000000)==0xc0000000) continue;
        else if( ((*(dptr+x))&0xf0000000)==0xd0000000) continue;
        else if( ((*(dptr+x))&0xfff00000)==0x11100000) continue;
      }
      //printf("--------------------------\n");

    }

    //if(counter>20) break;

    if((counter%5000)==0) {
      printf("    %.02f  bytes     \r",totalbytes/1000000.);
      fflush(stdout);
    }
  }

  printf("\n");
  printf("Counter: %lu       \n",counter);
  printf("MidasId: %lu       \n",lastmidasid);

  TChannel::WriteToRoot(goutfile);
  goutfile->Write();

  return 0;

}

