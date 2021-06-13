void Riso()
{
  //using namespace std;

  string file_input = "RisoluzioneT2.dat";

  int bin_res = 20;   // larghezza bin in ns
  int finestra = 1e3; // larghezza della finestrada dare in ns

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


  TH1F *h_nPMT1= new TH1F("h_nPMT1","Conteggi in 5 microsecondi, PMT 1",20,-0.5,19.5);
  TH1F *h_nPMT2= new TH1F("h_nPMT2","Conteggi in 5 microsecondi, PMT 2",20,-0.5,19.5);
  TH1F *h_nPMT3= new TH1F("h_nPMT3","Conteggi in 5 microsecondi, PMT 3",20,-0.5,19.5);


  //Variabili
  ifstream file_in;
	int PMT;
	double offset=0;
	double ts,ts_buf=0;
	int i=0;
  double ts_buf_PMT[3]={0.,0.,0.};
  int conteggi_PMT[3]={0.,0.,0.};
  double ts1,ts3,ts2;
  bool inversione=false;

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
      if ( (ts-ts_buf<5e-6) && (ts-ts_buf>-50) ) {

        if (ts-ts_buf < 0) inversione=true;

        if (PMT==1){
          if (ts-ts_buf_PMT[PMT-1]>0) h_dt1->Fill(ts-ts_buf_PMT[PMT-1]);
          if (ts-ts_buf_PMT[PMT-1]<0) h_dt1->Fill(-ts+ts_buf_PMT[PMT-1]);
          ts1=ts;
        }
        if (PMT==2){
          if (ts-ts_buf_PMT[PMT-1]>0) h_dt2->Fill(ts-ts_buf_PMT[PMT-1]);
          if (ts-ts_buf_PMT[PMT-1]<0) h_dt2->Fill(-ts+ts_buf_PMT[PMT-1]);
          ts2=ts;
        }
        if (PMT==3){
          if (ts-ts_buf_PMT[PMT-1]>0) h_dt3->Fill(ts-ts_buf_PMT[PMT-1]);
          if (ts-ts_buf_PMT[PMT-1]<0) h_dt3->Fill(-ts+ts_buf_PMT[PMT-1]);
          ts3=ts;
        }

        conteggi_PMT[PMT-1]+=1;

      }// if finestra 5 micro

      else{

        if (conteggi_PMT[0]!=0) h_nPMT1->Fill(conteggi_PMT[0]);
        if (conteggi_PMT[1]!=0) h_nPMT2->Fill(conteggi_PMT[1]);
        if (conteggi_PMT[2]!=0) h_nPMT3->Fill(conteggi_PMT[2]);

        if ((conteggi_PMT[0]==1)&&(conteggi_PMT[1]==1)&&(conteggi_PMT[2]==1)&&(inversione==false)) {
          h_dt12->Fill(ts1-ts2);
          h_dt13->Fill(ts1-ts3);
          h_dt23->Fill(ts2-ts3);
        }

        inversione=false;
        for(int j=0; j<3; j++)  conteggi_PMT[j]=0;

        ts_buf=ts;
        conteggi_PMT[PMT-1]+=1;
        if(PMT==1) ts1=ts;
        if(PMT==2) ts2=ts;
        if(PMT==3) ts3=ts;


      } // else


      ts_buf_PMT[PMT-1]=ts;

	  }//while endof file
	}//if file is open



  double m4_12=0;
  double m4_13=0;
  double m4_23=0;

  double Var_12=TMath::Power(h_dt12->GetStdDev(),2);
  double Var_13=TMath::Power(h_dt13->GetStdDev(),2);
  double Var_23=TMath::Power(h_dt23->GetStdDev(),2);

  double v1=0,dv1=0;
  double v2=0,dv2=0;
  double v3=0,dv3=0;


  for(int i=1;i<=N_bin;i++){
    m4_12+=(TMath::Power(h_dt12->GetBinCenter(i)-h_dt12->GetMean(),4)*(h_dt12->GetBinContent(i)/h_dt12->GetEffectiveEntries())    );
    m4_13+=(TMath::Power(h_dt13->GetBinCenter(i)-h_dt13->GetMean(),4)*(h_dt13->GetBinContent(i)/h_dt13->GetEffectiveEntries())    );
    m4_23+=(TMath::Power(h_dt23->GetBinCenter(i)-h_dt23->GetMean(),4)*(h_dt23->GetBinContent(i)/h_dt23->GetEffectiveEntries())    );
  }

  double Var_Var_12 = TMath::Sqrt((m4_12-((h_dt12->GetEffectiveEntries()-3)/(h_dt12->GetEffectiveEntries()-1)*TMath::Power(h_dt12->GetStdDev(),4)))/h_dt12->GetEffectiveEntries());
  double Var_Var_13 = TMath::Sqrt((m4_13-((h_dt13->GetEffectiveEntries()-3)/(h_dt13->GetEffectiveEntries()-1)*TMath::Power(h_dt13->GetStdDev(),4)))/h_dt13->GetEffectiveEntries());
  double Var_Var_23 = TMath::Sqrt((m4_23-((h_dt23->GetEffectiveEntries()-3)/(h_dt23->GetEffectiveEntries()-1)*TMath::Power(h_dt23->GetStdDev(),4)))/h_dt23->GetEffectiveEntries());


cout << m4_12 << "\n";
cout << m4_13 << "\n";
cout << m4_23 << "\n";


cout << "12    " << h_dt12->GetStdDev() << "    " << TMath::Power(h_dt12->GetStdDev(),2) << "    " << Var_Var_12 << "\n";
cout << "13    " << h_dt13->GetStdDev() << "    " << TMath::Power(h_dt13->GetStdDev(),2) << "    " << Var_Var_13 << "\n";
cout << "23    " << h_dt23->GetStdDev() << "    " << TMath::Power(h_dt23->GetStdDev(),2) << "    " << Var_Var_23 << "\n";






v1=(Var_13+Var_12-Var_23)/2;
v2=(Var_12+Var_23-Var_13)/2;
v3=(Var_13+Var_23-Var_12)/2;
dv1=(Var_Var_12+Var_Var_13+Var_Var_23)/2;
dv2=dv1;
dv3=dv1;

cout << "var 1    " << v1 << "    " << dv1 << "\n";
cout << "var 2    " << v2 << "    " << dv2 << "\n";
cout << "var 3    " << v3 << "    " << dv3 << "\n";


}//chiudo tutto
