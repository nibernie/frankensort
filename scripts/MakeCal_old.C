{

   cout << "------------------------" << endl;   
   cout << "----- PUNCH! -----" << endl;
   cout << "------------------------" << endl;
   cout << " Num Chans:  "  << TChannel::GetNumberOfChannels() << endl;
   cout << "------------------------" << endl;


   cout << GetRunNumber(_file0->GetName()) << endl;

   GH2D *chg = (GH2D*)_file0->Get("grif_channel_charge");
   //chg->Draw();

   TCalibrator c;
   for(int x=1;x<70;x++) {
      //if(x>1) break;
      GH1D *proj = chg->ProjectionY("_py",x,x);
      TChannel *channel = TChannel::GetChannelByNumber(x-1);


      if(proj->Integral()>10) {
        if(!channel) continue;
         channel->Print();
         cout << proj->Integral() << endl;
         proj->GetXaxis()->SetRangeUser(50,1300);
         c.AddData(proj,"152Eu",2.0,0.05);
         c.Fit();
         //proj->Draw();
         channel->DestroyENGCal();
         channel->AddENGCoefficient(c.GetParameter(0));
         channel->AddENGCoefficient(c.GetParameter(1));
         c.Clear();     
      }

   }
   string filename = Form("Cal%05i_Eu.cal",GetRunNumber(_file0->GetName()));
   TChannel::WriteCalFile(filename.c_str());

}
