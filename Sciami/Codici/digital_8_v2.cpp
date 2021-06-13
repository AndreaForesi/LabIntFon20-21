
{
	using namespace std;

	#define start 18.9166667 // è già convertito convertito in ore esempio erano le 18.55
	//#define start 12.08333
	#define days 6. // numero di giorni di durata della presa dati (se sono 4.5->5)
	#define bin 60. // inserire la larghezza del bin in minuti


	//tengo qui i nomi dei file
	//20210430_1855_20210505_911.dat
	//20210505_1115_20210507_1624.dat
	string file_input = "20210430_1855_20210505_911.dat";
	//string file_input = "20210505_1115_20210507_1624.dat";

	ifstream file_in;
	TGraph *gr = new TGraph(0);
	double PMT;
	double offset=0;
	double ts,ts_buf=0;
	int i=0;

	double N_bin=(60./bin)*24.*days;
	double N_bin_day=(60./bin)*24.;

	// flusso su telescopio 1
	TH1D *h_t1= new TH1D("h_t1","t Telescopio 1",N_bin,0,(24*days)); // prendo il tempo in ore
	TH1D *h_dt1= new TH1D("h_dt1","Delta t Telescopio 1 ",1000,0,1000*1000e-9);// qui ci metto la differenza di tempo tra eventi
	h_t1->GetXaxis()->SetTitle("t [ore]");
	h_dt1->GetXaxis()->SetTitle("dt [s]");

	// flusso su telescopio 2
	TH1D *h_t2= new TH1D("h_t2","t Telescopio 2",N_bin,0,(24*days)); // prendo il tempo in ore
	TH1D *h_dt2= new TH1D("h_dt2","Delta t Telescopio 2 ",1000,0,1000*1000e-9);// qui ci metto la differenza di tempo tra eventi
	h_t2->GetXaxis()->SetTitle("t [ore]");
	h_dt2->GetXaxis()->SetTitle("dt [s]");


	//differenza di tempo tra eventi tra i due telescopi 1-2
	TH1D *h_dt12= new TH1D("h_dt12","Delta t Telescopio 1 - Telescopio 2 ",2*10*100,-10*50*20e-9,10*50*20e-9);

	//differenza di tempo tra eventi tra i due telescopi 1-2 (MODULO)
	TH1D *h_dt12_abs= new TH1D("h_dt12_abs","Delta t |Telescopio 1 - Telescopio 2| ",50,0,50*20e-9);

	// flusso degli sciami
	TH1D *h_ts= new TH1D("h_ts","t Sciami",N_bin,0,(24*days));

	//istogrammi delle coincidenze doppie per calcolare le efficienze
	TH1D *h_t12= new TH1D("h_t12","t doppie 12",N_bin,0,(24*days));
	TH1D *h_t13= new TH1D("h_t13","t doppie 13",N_bin,0,(24*days));
	TH1D *h_t23= new TH1D("h_t23","t doppie 23",N_bin,0,(24*days));
	TH1D *h_t45= new TH1D("h_t45","t doppie 45",N_bin,0,(24*days));
	TH1D *h_t46= new TH1D("h_t46","t doppie 46",N_bin,0,(24*days));
	TH1D *h_t56= new TH1D("h_t56","t doppie 56",N_bin,0,(24*days));

	double ts_buf_PMT[8]={0.,0.,0.,0.,0.,0.,0.,0.};
	double dt1,dt2,ts_h,t1,t2,dt12;
	int search=0;
	int salto=0;
	int n_salti=0;
  file_in.open(file_input);
  if(file_in.is_open())
	{
		while( (!file_in.eof())  )
		{
			file_in >> PMT >> ts;
			// sistema il ts che sarebbe a "dente di sega"

			ts=ts+offset;

			if((ts<ts_buf)&&(ts-ts_buf>-300)&&(ts-ts_buf<-50))
			{
				offset=offset+85.88;
				ts=ts+85.88;
				salto=1;
				n_salti+=1;
			}

			if((ts<ts_buf)&&(ts-ts_buf<-300))
			{
				offset=offset+687.2;
				ts=ts+687.2;
				salto=1;
				n_salti+=1;
			}

        // cerco l' 1. trovo PMT =1. c' è il salto-> non calcolo dt(triple).
		if ((search==1)&&(PMT==5)) salto=0;
		if ((search==5)&&(PMT==1)) salto=0;

		ts_buf=ts;

		// differenzaa di tempo tra eveti (triple) sul singolo telescopio.
		if (PMT==1) h_dt1->Fill(ts-ts_buf_PMT[int(PMT)-1]);
		if (PMT==5) h_dt2->Fill(ts-ts_buf_PMT[int(PMT)-1]);

		ts_buf_PMT[int(PMT)-1]=ts;

		ts_h=(ts/(3600))+start; // converto il ts in ore

		// tempi dei telescopi

		if (PMT==1) h_t1->Fill(ts_h);
		if (PMT==5) h_t2->Fill(ts_h);

		// tempi delle coincidenze doppie
		if (PMT==2) h_t12->Fill(ts_h);
		if (PMT==3) h_t13->Fill(ts_h);
		if (PMT==4) h_t23->Fill(ts_h);
		if (PMT==6) h_t45->Fill(ts_h);
		if (PMT==7) h_t46->Fill(ts_h);
		if (PMT==8) h_t56->Fill(ts_h);

		gr->SetPoint(gr->GetN(), i, ts);
		i=i+1;

		// differenza di tempo tra i due telescopi
		//escludo i dt tra triple se c' è il salto nel mezzo
		if ((search==PMT) && (PMT==5) && (salto==0))
		{
			dt12=ts-ts_buf_PMT[1-1];
			h_dt12->Fill(dt12);
			h_dt12_abs->Fill(abs(dt12));
			if((dt12>-100e-9)&&(dt12<20e-9))
			{
				h_ts->Fill(start+((ts+ts_buf_PMT[1-1])/2)/3600);
			}
		}
		if ((search==PMT) && (PMT==1) && (salto==0))
		{
			dt12=-(ts-ts_buf_PMT[5-1]);
			h_dt12->Fill(dt12);
			h_dt12_abs->Fill(abs(dt12));
			if((dt12>-100e-9)&&(dt12<20e-9))
				{
					h_ts->Fill(start+((ts+ts_buf_PMT[5-1])/2)/3600);
				}
		}

		if(PMT==1)
		{
			search=5;
			salto=0; // se trovo PMT 1 o 5 per la volta successiva posso calcolare il dt(triple) quindi salto=0.
		}
		if(PMT==5)
		{
			search=1;
			salto=0; // se trovo PMT 1 o 5 per la volta successiva posso calcolare il dt(triple) quindi salto=0.
		}
	}
	}

	//---------

	// calcolo gli istogrammi delle efficienze

	// efficienze singoli rivelatori
	TH1D *h_e1= new TH1D("h_e1"," efficienza 1",N_bin,0,(24*days));
	TH1D *h_e2= new TH1D("h_e2"," efficienza 2",N_bin,0,(24*days));
	TH1D *h_e3= new TH1D("h_e3"," efficienza 3",N_bin,0,(24*days));
	TH1D *h_e4= new TH1D("h_e4"," efficienza 4",N_bin,0,(24*days));
	TH1D *h_e5= new TH1D("h_e5"," efficienza 5",N_bin,0,(24*days));
	TH1D *h_e6= new TH1D("h_e6"," efficienza 6",N_bin,0,(24*days));


	//efficienze dei due telescopi
	TH1D *h_et1= new TH1D("h_et1"," efficienza telescopio 1",N_bin,0,(24*days));
	TH1D *h_et2= new TH1D("h_et2"," efficienza telescopio 2",N_bin,0,(24*days));

	h_e1->Add(h_t1,h_t23,1.,1.);
	h_e1->Divide(h_t1,h_e1,1.,1.,"B");
	h_e2->Add(h_t1,h_t13,1.,1.);
	h_e2->Divide(h_t1,h_e2,1.,1.,"B");
	h_e3->Add(h_t1,h_t12,1.,1.);
	h_e3->Divide(h_t1,h_e3,1.,1.,"B");

	h_e4->Add(h_t2,h_t56,1.,1.);
	h_e4->Divide(h_t2,h_e4,1.,0.71789,"B");
	h_e5->Add(h_t2,h_t46,1.,1.);
	h_e5->Divide(h_t2,h_e5,1.,1.,"B");
	h_e6->Add(h_t2,h_t45,1.,1.);
	h_e6->Divide(h_t2,h_e6,1.,0.71789,"B");

	h_et1->Multiply(h_e1,h_e2,1.,1.);
	h_et1->Multiply(h_e3);

	h_et2->Multiply(h_e4,h_e5,1.,1.);
	h_et2->Multiply(h_e6);



	//triple senza efficienza
	TH1F *h_t1_noeff = (TH1F*)h_t1->Clone("h_t1_noeff");
	TH1F *h_t2_noeff = (TH1F*)h_t2->Clone("h_t2_noeff");


	// divido il numero di eventi nel tempo per le efficienze.
	h_t1->Divide(h_et1);
	h_t2->Divide(h_et2);

	//istogramma eventi sovrapposti giornalmente

	int i_day=0;
	double err=0;
	double day=0;


	// istogrammi per capire la correlazione tra i flussi dei due telescopi
	//TH2D *h_f1f2= new TH2D("h_f1f2","N eventi in un bin (5 min) telescopio 1 vs Telescopio 2",150,1000,7000,150,1000,9000); // bin 5 minuti
	TH2D *h_f1f2= new TH2D("h_f1f2","N eventi in un bin (5 min) telescopio 1 vs Telescopio 2",150,50e3,60e3,150,60e3,70e3);
	h_f1f2->GetXaxis()->SetTitle("# eventi/5 minuti telescopio 1");
	h_f1f2->GetYaxis()->SetTitle("# eventi/5 minuti telescopio 2");
	for (i=1;i<=N_bin-1;i++)
	{
		if((h_t1->GetBinContent(i)!=0) && (h_t2->GetBinContent(i)!=0))
		{
			h_f1f2->Fill(h_t1->GetBinContent(i),h_t2->GetBinContent(i));
		}


	}


	// vediamo la correlazione tra le efficienze dei due telescopi
	TH2D *h_e1e2= new TH2D("h_e1e2","efficienze (in un Bin) telescopio 1 vs Telescopio 2",1000,0,1,1000,0,1);
	h_e1e2->GetXaxis()->SetTitle("# eff(Bin) telescopio 1");
	h_e1e2->GetYaxis()->SetTitle("# eff(Bin) telescopio 2");
	for (i=1;i<=N_bin-1;i++)
	{
		if((h_et1->GetBinContent(i)!=0) && (h_et2->GetBinContent(i)!=0))
		{
			h_e1e2->Fill(h_et1->GetBinContent(i),h_et2->GetBinContent(i));
		}
	}

	TH2D *h_f1et1= new TH2D("h_f1et1","N eventi in un bin () telescopio 1 vs efficienza",150,50000,70000,1000,0,1);
	h_f1et1->GetXaxis()->SetTitle("# eventi/bin telescopio 1");
	h_f1et1->GetYaxis()->SetTitle("efficienze telescopio1");
	for (i=1;i<=N_bin-1;i++)
	{
		if((h_t1->GetBinContent(i)!=0) && (h_et1->GetBinContent(i)!=0))
		{
			h_f1et1->Fill(h_t1->GetBinContent(i),h_et1->GetBinContent(i));
		}
	}

	TH2D *h_f2et2= new TH2D("h_f2et2","N eventi in un bin () telescopio 2 vs efficienza",150,50000,70000,1000,0,1);
	h_f2et2->GetXaxis()->SetTitle("# eventi/bin telescopio 2");
	h_f2et2->GetYaxis()->SetTitle("efficienze telescopio2");
	for (i=1;i<=N_bin-1;i++)
	{
		if((h_t2->GetBinContent(i)!=0) && (h_et2->GetBinContent(i)!=0))
		{
			h_f2et2->Fill(h_t2->GetBinContent(i),h_et2->GetBinContent(i));
		}
	}

	TH2D *h_f1net1= new TH2D("h_f1net1","N eventi in un bin (no eff) telescopio 1 vs efficienza",100,20000,30000,1000,0,1);
	h_f1et1->GetXaxis()->SetTitle("# eventi/bin telescopio 1");
	h_f1et1->GetYaxis()->SetTitle("efficienze telescopio1");
	for (i=1;i<=N_bin-1;i++)
	{
		if((h_t1->GetBinContent(i)!=0) && (h_et1->GetBinContent(i)!=0))
		{
			h_f1net1->Fill(h_t1_noeff->GetBinContent(i),h_et1->GetBinContent(i));
		}
	}

	TH2D *h_f2net2= new TH2D("h_f2net2","N eventi in un bin (no eff) telescopio 2 vs efficienza",100,27000,37000,1000,0,1);
	h_f2et2->GetXaxis()->SetTitle("# eventi/bin telescopio 2");
	h_f2et2->GetYaxis()->SetTitle("efficienze telescopio2");
	for (i=1;i<=N_bin-1;i++)
	{
		if((h_t2->GetBinContent(i)!=0) && (h_et2->GetBinContent(i)!=0))
		{
			h_f2net2->Fill(h_t2_noeff->GetBinContent(i),h_et2->GetBinContent(i));
		}
	}

	TH2D *h_f1f2ne= new TH2D("h_f1f2ne","N eventi in un bin (no eff) telescopio 1 vs telescopio2",150,20e3,30e3,150,25e3,35e3);
	h_f1f2ne->GetXaxis()->SetTitle("# eventi/bin telescopio 1");
	h_f1f2ne->GetYaxis()->SetTitle("# eventi/bin telescopio 2");
	for (i=1;i<=N_bin-1;i++)
	{
		if((h_t1_noeff->GetBinContent(i)!=0) && (h_t2_noeff->GetBinContent(i)!=0))
		{
			h_f1f2ne->Fill(h_t1_noeff->GetBinContent(i),h_t2_noeff->GetBinContent(i));
		}
	}

	TH2D *h_fD12f2= new TH2D("h_fD12f2","N eventi in un bin Doppie 12 vs telescopio 2",150,64e2,74e2,150,64e3,67e3);
	h_f1f2ne->GetXaxis()->SetTitle(" Doppie 12 / bin");
	h_f1f2ne->GetYaxis()->SetTitle("Triple telescopio 2 / bin");
	for (i=1;i<=N_bin-1;i++)
	{
		if((h_t12->GetBinContent(i)!=0) && (h_t2->GetBinContent(i)!=0))
		{
			h_fD12f2->Fill(h_t12->GetBinContent(i),h_t2->GetBinContent(i));
		}
	}

	cout << "corr telescopio1/telescopio 2 = " << h_f1f2->GetCorrelationFactor() <<"\n";
	cout << "corr  efficienze telescopio1/telescopio 2 = " << h_e1e2->GetCorrelationFactor() <<"\n";

	cout << "corr ev. telescopio1/ eff telescopio 1 = " << h_f1et1->GetCorrelationFactor() <<"\n";
	cout << "corr ev. telescopio2/ eff telescopio 2 = " << h_f2et2->GetCorrelationFactor() <<"\n";

	cout << "corr ev. telescopio1 (noeff)/ eff telescopio 1 = " << h_f1net1->GetCorrelationFactor() <<"\n";
	cout << "corr ev. telescopio2 (noeff)/ eff telescopio 2 = " << h_f2net2->GetCorrelationFactor() <<"\n";

	cout << "corr ev. Doppie 12/ Triple t2 = " << h_fD12f2->GetCorrelationFactor() <<"\n";

}
