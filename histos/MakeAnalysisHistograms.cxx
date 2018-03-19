#include "TRuntimeObjects.h"

//#include "TTigress.h"
#include "TGriffin.h"

std::map<int,unsigned long> lasttime;

long largest_current=-1;
long smallest_current=-1;
long largest_last=-1;
long smallest_last=-1;
long counter=0;

extern "C" void MakeAnalysisHistograms(TRuntimeObjects& obj)
{
  auto griffin = obj.GetDetector<TGriffin>();
  if(!griffin) return;
  griffin->FixLowGainCrossTalk();
  counter=0;
  for(int i=0;i<griffin->Size();i++) {
    TGriffinHit *hit = griffin->GetGriffinHit(i); 

    if(counter==0) {
      largest_current=hit->GetTimeStamp();
      smallest_current=hit->GetTimeStamp();
    } else {
      if(hit->GetTimeStamp()>largest_current) {
        largest_current = hit->GetTimeStamp();
      }
      if(hit->GetTimeStamp()<smallest_current) {
        smallest_current = hit->GetTimeStamp();
      }
    }

    obj.FillHistogram("summary",4000,0,2000,hit->GetEnergy(),
                               70,0,70,hit->GetArrayNumber());

    obj.FillHistogram("runtime",3600,0,3600,hit->GetTimeStamp()/1e8,
                               70,0,70,hit->GetArrayNumber());

    TChannel *chan = TChannel::GetChannel(hit->GetAddress());
    if(chan) 
    obj.FillHistogram("xtalmap",70,0,70,hit->GetArrayNumber(),
                                70,0,70,chan->GetNumber());

    //obj.FillHistogram("xtalmap",70,0,70,hit->GetArrayNumber(),
    //                            70,0,70,chan->GetNumber());

    obj.FillHistogram(Form("xtali%02i_drif",hit->GetArrayNumber()),3000,0,3000,hit->GetTimeStamp()/1e8,
                                                                   500,1275,1525,hit->GetEnergy());

    obj.FillHistogram("hitpattern_xtal",70,0,70,hit->GetArrayNumber());
    //obj.FillHistogram("hitpattern_name",70,0,70,hit->GetName());


    if(lasttime.count(hit->GetArrayNumber())) {
      obj.FillHistogram("dead",10000,0,10000,hit->GetTimeStamp()-lasttime[hit->GetArrayNumber()],
                               70,0,70,hit->GetArrayNumber());
    }
    lasttime[hit->GetArrayNumber()] = hit->GetTimeStamp();
  }

  obj.FillHistogram("between_events",10000,0,10000,smallest_current-largest_last,
                                     2,0,2,1);

  largest_last = largest_current;
  smallest_last = smallest_current;

}
