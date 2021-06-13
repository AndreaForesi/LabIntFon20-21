void Riso()
{
  //using namespace std;

  string file_input = "20210428/5canali.txt";

  int bin_res = 20;   // larghezza bin in ns
  int finestra = 60; // larghezza della finestrada dare in ns

  int N_bin = finestra/bin_res;

  //istogrammi
  TH1D *h_dtD12= new TH1D("h_dtD12","Delta t doppie 12 ",N_bin,0,N_bin*bin_res*1e-9);
  TH1D *h_dtD13= new TH1D("h_dtD13","Delta t doppie 13 ",N_bin,0,N_bin*bin_res*1e-9);
  TH1D *h_dtD23= new TH1D("h_dtD23","Delta t doppie 23 ",N_bin,0,N_bin*bin_res*1e-9);
  h_dtD12->GetXaxis()->SetTitle("Ritardo [s]");
  h_dtD13->GetXaxis()->SetTitle("Ritardo [s]");
  h_dtD23->GetXaxis()->SetTitle("Ritardo [s]");

  TH1D *h_dtD12D13= new TH1D("h_dtD12D13","Delta t doppie 12 - doppie 13 ",2*N_bin,-N_bin*bin_res*1e-9,N_bin*bin_res*1e-9);
  TH1D *h_dtD12D23= new TH1D("h_dtD12D23","Delta t doppie 12 - doppie 23 ",2*N_bin,-N_bin*bin_res*1e-9,N_bin*bin_res*1e-9);
  TH1D *h_dtD13D23= new TH1D("h_dtD13D23","Delta t doppie 13 - doppie 23 ",2*N_bin,-N_bin*bin_res*1e-9,N_bin*bin_res*1e-9);
  h_dtD12D13->GetXaxis()->SetTitle("Ritardo [s]");
  h_dtD12D23->GetXaxis()->SetTitle("Ritardo [s]");
  h_dtD12D23->GetXaxis()->SetTitle("Ritardo [s]");


  TH1F *h_nD12= new TH1F("h_nD12","Conteggi in 5 microsecondi, doppie 12",20,-0.5,19.5);
  TH1F *h_nD13= new TH1F("h_nD13","Conteggi in 5 microsecondi, doppie 13",20,-0.5,19.5);
  TH1F *h_nD23= new TH1F("h_nD23","Conteggi in 5 microsecondi, doppie 23",20,-0.5,19.5);

  TH1D *h_dtMax= new TH1D("h_dtMax","Delta t Max ",2*N_bin,-N_bin*bin_res*1e-9,N_bin*bin_res*1e-9);


  double intervallo= 500e-9; // in secondi

  //Variabili
  ifstream file_in;
	int PMT;
	double offset=0;
	double ts,ts_buf=0;
	int i=0;
  double ts_buf_PMT[8]={0.,0.,0.,0.,0.,0.,0.,0.};
  int conteggi_PMT[8]={0.,0.,0.,0.,0.,0.,0.,0.};
  double ts1,ts3,ts2,ts4,ts5,ts6,ts7,ts8;
  //bool inversione=false;

  double m=0;

  int salto=0;

  //Lettura file
  file_in.open(file_input);
  if(file_in.is_open())
	{
		while( (!file_in.eof())  )
		{
      file_in >> PMT >> ts;

      if (ts<ts_buf-50) salto+=1;

      //Se sono nella finestra di 5 micro
      if ( (ts-ts_buf<intervallo) && (ts-ts_buf>-50) ) {

        //if (ts-ts_buf < 0) inversione=true;

        if (PMT==2){
          if (ts-ts_buf_PMT[PMT-1]>0) h_dtD12->Fill(ts-ts_buf_PMT[PMT-1]);
          if (ts-ts_buf_PMT[PMT-1]<0) h_dtD12->Fill(-ts+ts_buf_PMT[PMT-1]);
          ts2=ts;
        }
        if (PMT==3){
          if (ts-ts_buf_PMT[PMT-1]>0) h_dtD13->Fill(ts-ts_buf_PMT[PMT-1]);
          if (ts-ts_buf_PMT[PMT-1]<0) h_dtD13->Fill(-ts+ts_buf_PMT[PMT-1]);
          ts3=ts;
        }
        if (PMT==4){
          if (ts-ts_buf_PMT[PMT-1]>0) h_dtD23->Fill(ts-ts_buf_PMT[PMT-1]);
          if (ts-ts_buf_PMT[PMT-1]<0) h_dtD23->Fill(-ts+ts_buf_PMT[PMT-1]);
          ts4=ts;
        }

        conteggi_PMT[PMT-1]+=1;

      }// if finestra 5 micro

      else{

        if (conteggi_PMT[1]!=0) h_nD12->Fill(conteggi_PMT[1]);
        if (conteggi_PMT[2]!=0) h_nD13->Fill(conteggi_PMT[2]);
        if (conteggi_PMT[3]!=0) h_nD23->Fill(conteggi_PMT[3]);

        if ((conteggi_PMT[1]==1)&&(conteggi_PMT[2]==1)&&(conteggi_PMT[3]==1)) {   //&&(inversione==false)
          h_dtD12D13->Fill(ts2-ts3);
          h_dtD12D23->Fill(ts2-ts4);
          h_dtD13D23->Fill(ts3-ts4);
          m=ts2-ts3;
          if (abs(ts2-ts4)>abs(m)) m=ts2-ts4;
          if (abs(ts3-ts4)>abs(m)) m=ts3-ts4;
          h_dtMax->Fill(m);
        }

        //inversione=false;
        for(int j=0; j<8; j++)  conteggi_PMT[j]=0;

        ts_buf=ts;
        conteggi_PMT[PMT-1]+=1;
        if(PMT==2) ts2=ts;
        if(PMT==3) ts3=ts;
        if(PMT==4) ts4=ts;


      } // else


      ts_buf_PMT[PMT-1]=ts;

	  }//while endof file
	}//if file is open



  double m4=0;

  double Var=TMath::Power(h_dtMax->GetStdDev(),2);

  for(int i=1;i<=N_bin;i++){
    m4+=(TMath::Power(h_dtMax->GetBinCenter(i)-h_dtMax->GetMean(),4)*(h_dtMax->GetBinContent(i)/h_dtMax->GetEffectiveEntries())    );
  }

  double Var_Var = TMath::Sqrt((m4-((h_dtMax->GetEffectiveEntries()-3)/(h_dtMax->GetEffectiveEntries()-1)*TMath::Power(h_dtMax->GetStdDev(),4)))/h_dtMax->GetEffectiveEntries());


cout << m4 << "\n";

cout << "Max    " << h_dtMax->GetStdDev() << "    Var    " << TMath::Power(h_dtMax->GetStdDev(),2) << "    VarVar    " << Var_Var << "\n";





}//chiudo tutto
