void Risoluzione()
{
  //using namespace std;

  string file_input = "RisoluzioneT1.dat";

  int bin_res = 1e6;   // larghezza bin in ns
  int finestra = 1e9; // larghezza della finestrada dare in ns

  int N_bin = finestra/bin_res;

  //istogrammi
  TH1D *h_dt1= new TH1D("h_dt1","Delta t PMT 1 ",N_bin,0,N_bin*bin_res*1e-9);
  TH1D *h_dt2= new TH1D("h_dt2","Delta t PMT 2 ",N_bin,0,N_bin*bin_res*1e-9);
  TH1D *h_dt3= new TH1D("h_dt3","Delta t PMT 3 ",N_bin,0,N_bin*bin_res*1e-9);
  h_dt1->GetXaxis()->SetTitle("Ritardo [s]");
  h_dt2->GetXaxis()->SetTitle("Ritardo [s]");
  h_dt3->GetXaxis()->SetTitle("Ritardo [s]");

  TH1D *h_dt12= new TH1D("h_dt12","Delta t PMT 1 - PMT 2 ",2*N_bin,-N_bin*bin_res*1e-9,N_bin*bin_res*1e-9);
  TH1D *h_dt13= new TH1D("h_dt13","Delta t PMT 1 - PMT 3 ",2*N_bin,-N_bin*bin_res*1e-9,N_bin*bin_res*1e-9);
  TH1D *h_dt23= new TH1D("h_dt23","Delta t PMT 2 - PMT 3 ",2*N_bin,-N_bin*bin_res*1e-9,N_bin*bin_res*1e-9);
  h_dt12->GetXaxis()->SetTitle("Ritardo [s]");
  h_dt13->GetXaxis()->SetTitle("Ritardo [s]");
  h_dt23->GetXaxis()->SetTitle("Ritardo [s]");

  TH1D *h_dt12_abs= new TH1D("h_dt12_abs","Delta t PMT 1 - PMT 2 ",N_bin,0,N_bin*bin_res*1e-9);
  TH1D *h_dt13_abs= new TH1D("h_dt13_abs","Delta t PMT 1 - PMT 3 ",N_bin,0,N_bin*bin_res*1e-9);
  TH1D *h_dt23_abs= new TH1D("h_dt23_abs","Delta t PMT 2 - PMT 3 ",N_bin,0,N_bin*bin_res*1e-9);
  h_dt12_abs->GetXaxis()->SetTitle("Ritardo [s]");
  h_dt13_abs->GetXaxis()->SetTitle("Ritardo [s]");
  h_dt23_abs->GetXaxis()->SetTitle("Ritardo [s]");



  //Variabili
  ifstream file_in;
	TGraph *gr = new TGraph(0);
	int PMT;
	double offset=0;
	double ts,ts_buf=0;
	int i=0;
  double ts_buf_PMT[3]={0.,0.,0.};

	bool salto1=false;
  bool salto2=false;
  bool salto3=false;

  //Lettura file
  file_in.open(file_input);
  if(file_in.is_open())
	{
		while( (!file_in.eof())  )
		{
      file_in >> PMT >> ts;

      ts=ts+offset;

      // sistema il ts che sarebbe a "dente di sega"
      if((ts<ts_buf)&&(ts-ts_buf>-300)&&(ts-ts_buf<-50)){
        offset=offset+85.88;
        ts=ts+85.88;
        salto1=true;
        salto2=true;
        salto3=true;
      }
      if((ts<ts_buf)&&(ts-ts_buf<-300)){
        offset=offset+687.2;
        ts=ts+687.2;
        salto1=true;
        salto2=true;
        salto3=true;
      }


      gr->SetPoint(gr->GetN(), i, ts);
		  i=i+1;

      // differenza di tempo tra i PMT
      //escludo i dt tra triple se c' è il salto nel mezzo
      if ((salto1==false)&&(PMT==1)) {
        if (ts-ts_buf_PMT[PMT-1]>5e-6){
          h_dt1->Fill(ts-ts_buf_PMT[PMT-1]);
          if (salto2==false)  h_dt12->Fill(ts-ts_buf_PMT[2-1]);
          if (salto3==false)  h_dt13->Fill(ts-ts_buf_PMT[3-1]);
        }
      }
      if ((salto2==false)&&(PMT==2)) {
        if (ts-ts_buf_PMT[PMT-1]>5e-6){
          h_dt2->Fill(ts-ts_buf_PMT[PMT-1]);
          if (salto1==false)  h_dt12->Fill(ts_buf_PMT[1-1]-ts);
          if (salto3==false)  h_dt23->Fill(ts-ts_buf_PMT[3-1]);
        }
      }
      if ((salto3==false)&&(PMT==3)) {
        if (ts-ts_buf_PMT[PMT-1]>5e-6){
          h_dt3->Fill(ts-ts_buf_PMT[PMT-1]);
          if (salto1==false)  h_dt13->Fill(ts_buf_PMT[1-1]-ts);
          if (salto2==false)  h_dt23->Fill(ts_buf_PMT[2-1]-ts);
        }
      }


      ts_buf_PMT[PMT-1]=ts;

      ts_buf=ts;

      if ((PMT==1)&&(salto1=true)) salto1=false;
      if ((PMT==2)&&(salto2=true)) salto2=false;
      if ((PMT==3)&&(salto3=true)) salto3=false;



	  }//while endof file
	}//if file is open
}//chiudo tutto
