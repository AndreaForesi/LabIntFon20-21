double cost1=511/1.28452e5;
double cost2=511/4.32591e4;
double cost3=1274/5.15345e4;

double a1 = -0.030215527720062252;
double b1 =  293.85869580897685;
double o1 =  -12616.48494675589;

double a2 = -0.0027284093341694623;
double b2 =  88.81191865012056;
double o2 =  -1631.6656576431283;


/*
*/
double Enel(double y,double a, double b, double c){
	 double rd = TMath::Sqrt(TMath::Power(b,2)-4*a*(c-y));
	 double x  =   (-b+rd)/(2*a);
	 return x;
 }
TH1F *istC=NULL;
TH1F *istD=NULL;
TH1F *istC1=NULL;
TH1F *istC2=NULL;


void Distri(int bin, double min, double Max, double Delta, bool fit){

istC = new TH1F("istC"," Ritardo tPMT1-tPMT2 ",bin,min,Max);
istD = new TH1F("istD"," Differenza energia PMT-1 - energia PMT-2 ",bin,-Delta-6,Delta-6);
istC1 = new TH1F("istC1"," Energia PMT-1 ",bin,min,Max);
istC2 = new TH1F("istC2"," Misura massa elettrone, PMT-2 ",bin,min,Max);


//TH1F *ist1= new TH1F("ist1"," Triple ",bin,min,Max);
//TH1F *ist2= new TH1F("ist2"," Triple ",bin,min,Max);
//TH1F *ist3= new TH1F("ist3"," Triple ",bin,min,Max);

ifstream inp;
ifstream inp2;
ifstream inp3;
double x,y,z;
double Ts3,Ts1,Ts2;
double Dt12,Dt13,Dt23,c1,c2,c3;

inp.open("1a_misura_me_20210305/Out_CFix_TsRetta2P5Sigma_PMT2_doppie_1725.txt");
inp2.open("1a_misura_me_20210305/Out_CFix_TsRetta2P5Sigma_PMT2_doppie_1740.txt");
inp3.open("1a_misura_me_20210305/Out_CFix_TsRetta2P5Sigma_PMT2_doppie_1759.txt");

while(!inp.eof()){
	inp  >> x >> c1 >> Ts1;
	inp2 >> y >> c2 >> Ts2;
	inp3 >> z >> c3 >> Ts3;
	//istC->Fill(Enel(c1,a1,b1,o1)+Enel(c2,a2,b2,o2));
	//istD->Fill(Enel(c1,a1,b1,o1)-Enel(c2,a2,b2,o2));
	istC1->Fill(Enel(c1,a2,b2,o2));
	istC1->Fill(Enel(c2,a2,b2,o2));
	istC1->Fill(Enel(c3,a2,b2,o2));
	//istC2->Fill(Enel(c2,a2,b2,o2));
}
inp.close();
inp2.close();
inp3.close();
/*
inp.open("1a_misura_me_20210305/Out_CFix_TsRetta2P5Sigma_PMT1_doppie_1740.txt");
inp2.open("1a_misura_me_20210305/Out_CFix_TsRetta2P5Sigma_PMT2_doppie_1740.txt");
//inp3.open("triple180_20210319/TsRetta2p_5SPMT3_triple180_1743.txt");
while(!inp.eof()){
	inp >> x >> c1 >> Ts ;
	inp2 >> y >> c2 >> Ts1;
	//inp3 >> Ts;
	istC->Fill(Enel(c1,a1,b1,o1)+Enel(c2,a2,b2,o2));
	istD->Fill(Enel(c1,a1,b1,o1)-Enel(c2,a2,b2,o2));
	istC1->Fill(Enel(c1,a1,b1,o1));
	istC2->Fill(Enel(c2,a2,b2,o2));
}
inp.close();
inp2.close();


//inp3.close();
inp.open("PuntiLavPMT123_20210224/PMT1_800V_1042.txt");
//inp2.open("1a_misura_me_20210305/Out_Ts5Sigma_PMT2_doppie_1725.txt");
//inp3.open("triple180_20210319/TsRetta2p_5SPMT3_triple180_1837.txt");

while(!inp.eof()){
	inp >> x  >> c1 ;
	//inp2 >> y >> Ts1;
	//inp3 >> Ts;
	istC->Fill(-c1);
	//istC->Fill(Enel(c1,a1,b1,o1)+Enel(c2,a2,b2,o2));
	//istD->Fill(Enel(c1,a1,b1,o1)-Enel(c2,a2,b2,o2));
	//istC1->Fill(Enel(c1,a1,b1,o1));
	//istC2->Fill(Enel(c2,a2,b2,o2));
}
inp.close();
//inp2.close();
//inp3.close();

inp.open("1a_misura_me_20210305/Out_CFix_TsRetta2P5Sigma_PMT2_sodio_1815.txt");

while(!inp.eof()){
	inp >> x >> c1 >> Ts ;
	istC->Fill(c1);
}
inp.close();
*/


istC->GetXaxis()->SetTitle("Energia [keV]");
istD->GetXaxis()->SetTitle("Energia [keV]");
istC1->GetXaxis()->SetTitle("Energia [keV]");
istC2->GetXaxis()->SetTitle("Energia [keV]");
//istC->GetXaxis()->SetTitle("Tempo [ns]");

//istC->Draw();
//istD->Draw();


if (fit==TRUE){
	istC->Fit("gaus","L");
	//istD->Fit("gaus","L");
	//istC1->Fit("gaus","L");
	//istC2->Fit("gaus","L");
}
}

void Dt(int ch){

	int bin=200;
	int min=-400;
	int max=400;

	if (ch==2){
		TH1F *Dt= new TH1F("Dt","Ritardi segnali",bin,min,max);
		ifstream inp;
		double Ev;
		double dt,c1,c2;
		inp.open("Triple_20210310/Dt12_Ts10Sigma_CFix_triple_1125.txt");
		while(!inp.eof()){
			inp >> Ev >> dt >> c1 >> c2;
			Dt->Fill(dt);
		}
		inp.close();
		Dt->GetXaxis()->SetTitle("Dt [ns]");
		Dt->Draw();
	}
	if (ch==3){
		int diffCanali;
		scanf("%d" , &diffCanali);
		if (diffCanali==12){
			TH1F *Dt= new TH1F("Dt","Ritardi segnali 1-2",bin,min,max);
			ifstream inp;
			double Ev;
			double dt12,dt13,dt23,c1,c2,c3;
			inp.open("triple180_20210319/Dt3Ch_Ts5Sigma_CFix_triple180_1707.txt");
			while(!inp.eof()){
				inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
				Dt->Fill(dt12);
			}
			inp.close();
			inp.open("triple180_20210319/Dt3Ch_Ts5Sigma_CFix_triple180_1743.txt");
			while(!inp.eof()){
				inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
				Dt->Fill(dt12);
			}
			inp.close();
			inp.open("triple180_20210319/Dt3Ch_Ts5Sigma_CFix_triple180_1837.txt");
			while(!inp.eof()){
				inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
				Dt->Fill(dt12);
			}
			inp.close();
			/*
			inp.open("tripleor_20210317/Dt12_Ts5Sigma_CFix_tripleOr_1130.txt");
			while(!inp.eof()){
				inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
				Dt->Fill(dt12);
			}
			inp.close();*/
			Dt->GetXaxis()->SetTitle("Dt [ns]");
			Dt->Draw();
		}
		if (diffCanali==13){
			TH1F *Dt= new TH1F("Dt","Ritardi segnali 1-3",bin,min,max);
			ifstream inp;
			double Ev;
			double dt12,dt13,dt23,c1,c2,c3;
			inp.open("triple180_20210319/Dt3Ch_Ts5Sigma_CFix_triple180_1707.txt");
			while(!inp.eof()){
				inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
				//if (dt12>8){
				Dt->Fill(dt13);
				//}
			}
			inp.close();
			inp.open("triple180_20210319/Dt3Ch_Ts5Sigma_CFix_triple180_1743.txt");
			while(!inp.eof()){
				inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
				//if (dt12>8){
				Dt->Fill(dt13);
				//}
			}
			inp.close();
			inp.open("triple180_20210319/Dt3Ch_Ts5Sigma_CFix_triple180_1837.txt");
			while(!inp.eof()){
				inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
				//if ((c1>3000)&&(c1<110000)){
				Dt->Fill(dt13);
				//}
			}
			inp.close();
			/*
			inp.open("tripleor_20210317/Dt3Ch_TsMedia5Sigma_CFix_tripleOr_1130.txt");
			while(!inp.eof()){
				inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
				if ((c1>3000)&&(c1<110000)){
					Dt->Fill(dt13);
				}
			}
			inp.close();
			*/
			Dt->GetXaxis()->SetTitle("Dt [ns]");

			Dt->Draw();
		}
		if (diffCanali==23){

			TH1F *Dt= new TH1F("Dt","Ritardi segnali 2-3",bin,min,max);
			ifstream inp;
			double Ev;
			double dt12,dt13,dt23,c1,c2,c3;
			inp.open("triple180_20210319/Dt3Ch_Ts5Sigma_CFix_triple180_1707.txt");
			while(!inp.eof()){
				inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
				//if ((c2>3000)&&(c2<40000)){
				Dt->Fill(dt23);
				//}
			}
			inp.close();
			inp.open("triple180_20210319/Dt3Ch_Ts5Sigma_CFix_triple180_1743.txt");
			while(!inp.eof()){
				inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
				//if ((c2>3000)&&(c2<40000)){
				Dt->Fill(dt23);
				//}
			}
			inp.close();
			inp.open("triple180_20210319/Dt3Ch_Ts5Sigma_CFix_triple180_1837.txt");
			while(!inp.eof()){
				inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
				//if ((c2>3000)&&(c2<40000)){
				Dt->Fill(dt23);
				//}
			}
			inp.close();
			/*
			inp.open("tripleor_20210317/Dt3Ch_TsMedia5Sigma_CFix_tripleOr_1130.txt");
			while(!inp.eof()){
				inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
				if ((c2>3000)&&(c2<40000)){
					Dt->Fill(dt23);
				}
			}
			inp.close();
			*/
			Dt->GetXaxis()->SetTitle("Dt [ns]");
			Dt->Draw();
		}
	}
}

void Scatter(string name){

	string filePath1 = "triple180_20210319/Dt3Ch_Ts5Sigma_CFix_triple180_1707.txt";
	string filePath2 = "triple180_20210319/Dt3Ch_Ts5Sigma_CFix_triple180_1743.txt";
	string filePath3 = "triple180_20210319/Dt3Ch_Ts5Sigma_CFix_triple180_1837.txt";


	if (name=="C1C2"){

		TH2F *Scatter= new TH2F("Scatter","C1 Vs C2",150,0,450000,100,0,100000);

		ifstream inp;
		double Ev;
		double dt,c1,c2;
		inp.open(filePath1);
		while(!inp.eof()){
			inp >> Ev >> dt >> c1 >> c2;
			Scatter->Fill(c1*cost1,c2*cost2);
		}
		inp.close();
		Scatter->GetXaxis()->SetTitle("Carica sul PMT1");
		Scatter->GetYaxis()->SetTitle("Carica sul PMT2");
		Scatter->Draw("ZColor");
	}
	if (name=="C1C3"){

		TH2F *Scatter= new TH2F("Scatter","C1 Vs C3",100,0,700,100,0,1500);

		ifstream inp;
		double Ev;
		double dt12,dt13,dt23,c1,c2,c3;
		inp.open(filePath1);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(c1*cost1,c3*cost3);
		}
		inp.close();
		inp.open(filePath2);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(c1*cost1,c3*cost3);
		}
		inp.close();
		inp.open(filePath3);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(c1*cost1,c3*cost3);
		}
		inp.close();
		/*
		*/
		Scatter->GetXaxis()->SetTitle("Carica sul PMT1 [keV]");
		Scatter->GetYaxis()->SetTitle("Carica sul PMT3 [keV]");
		Scatter->Draw("ZColor");
	}
	if (name=="C2C3"){

		TH2F *Scatter= new TH2F("Scatter","C2 Vs C3",100,0,700,100,0,1500);

		ifstream inp;
		double Ev;
		double dt12,dt13,dt23,c1,c2,c3;
		inp.open(filePath1);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(c2*cost2,c3*cost3);
		}
		inp.close();
		inp.open(filePath2);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(c2*cost2,c3*cost3);
		}
		inp.close();
		inp.open(filePath3);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(c2*cost2,c3*cost3);
		}
		inp.close();
		/*
		*/
		Scatter->GetXaxis()->SetTitle("Carica sul PMT2 [keV]");
		Scatter->GetYaxis()->SetTitle("Carica sul PMT3 [keV]");
		Scatter->Draw("ZColor");
	}
	if (name=="DtC2"){

		TH2F *Scatter= new TH2F("Scatter","Dt Vs C2",100,-200,200,100,0,200000);

		ifstream inp;
		double Ev;
		double dt,c1,c2;
		inp.open(filePath1);
		while(!inp.eof()){
			inp >> Ev >> dt >> c1 >> c2;
			Scatter->Fill(dt,c2*cost2);
		}
		inp.close();
		Scatter->GetXaxis()->SetTitle("Dt [ns]");
		Scatter->GetYaxis()->SetTitle("Carica sul PMT2 [u.a.]");
		Scatter->Draw("ZColor");
	}
	if (name=="C1-C2"){

		TH1F *TriImp= new TH1F("TriImp","3-Impulso, C1-C2",150,1000,100000);

		ifstream inp;
		double Ev;
		double dt,c1,c2;
		inp.open(filePath1);
		while(!inp.eof()){
			inp >> Ev >> dt >> c1 >> c2;
			TriImp->Fill(c1-c2);
		}
		inp.close();
		TriImp->GetXaxis()->SetTitle("C1-C2, [u.a.]");
		//Scatter->GetYaxis()->SetTitle("");
		TriImp->Draw("");
	}

	if (name=="Dt12C1"){

		TH2F *Scatter= new TH2F("Scatter","Dt12 Vs C1",50,-200,200,100,0,1500);

		ifstream inp;
		double Ev;
		double dt12,dt13,dt23,c1,c2,c3;
		inp.open(filePath1);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(dt12,c1*cost1);
		}
		inp.close();
		inp.open(filePath2);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(dt12,c1*cost1);
		}
		inp.close();
		inp.open(filePath3);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(dt12,c1*cost1);
		}
		inp.close();
		/*
		*/
		Scatter->GetXaxis()->SetTitle("Dt12 [ns]");
		Scatter->GetYaxis()->SetTitle("Carica sul PMT1 [keV]");
		Scatter->Draw("ZColor");
	}
	if (name=="Dt12C2"){

		TH2F *Scatter= new TH2F("Scatter","Dt12 Vs C2",50,-200,200,100,0,1500);

		ifstream inp;
		double Ev;
		double dt12,dt13,dt23,c1,c2,c3;
		inp.open(filePath1);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(dt12,c2*cost2);
		}
		inp.close();
		inp.open(filePath2);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(dt12,c2*cost2);
		}
		inp.close();
		inp.open(filePath3);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(dt12,c2*cost2);
		}
		inp.close();
		/*
		*/
		Scatter->GetXaxis()->SetTitle("Dt12 [ns]");
		Scatter->GetYaxis()->SetTitle("Carica sul PMT2 [keV]");
		Scatter->Draw("ZColor");
	}
	if (name=="Dt12C3"){

		TH2F *Scatter= new TH2F("Scatter","Dt12 Vs C3",50,-200,200,100,0,1500);

		ifstream inp;
		double Ev;
		double dt12,dt13,dt23,c1,c2,c3;
		inp.open(filePath1);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(dt12,c3*cost3);
		}
		inp.close();
		inp.open(filePath2);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(dt12,c3*cost3);
		}
		inp.close();
		inp.open(filePath3);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(dt12,c3*cost3);
		}
		inp.close();
		/*
		*/
		Scatter->GetXaxis()->SetTitle("Dt12 [ns]");
		Scatter->GetYaxis()->SetTitle("Carica sul PMT3 [keV]");
		Scatter->Draw("ZColor");
	}

	if (name=="Dt13C1"){

		TH2F *Scatter= new TH2F("Scatter","Dt13 Vs C1",50,-200,200,100,0,1500);

		ifstream inp;
		double Ev;
		double dt12,dt13,dt23,c1,c2,c3;
		inp.open(filePath1);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(dt13,c1*cost1);
		}
		inp.close();
		inp.open(filePath2);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(dt13,c1*cost1);
		}
		inp.close();
		inp.open(filePath3);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(dt13,c1*cost1);
		}
		inp.close();
		/*
		*/
		Scatter->GetXaxis()->SetTitle("Dt13 [ns]");
		Scatter->GetYaxis()->SetTitle("Carica sul PMT1 [keV]");
		Scatter->Draw("ZColor");
	}
	if (name=="Dt13C3"){

		TH2F *Scatter= new TH2F("Scatter","Dt13 Vs C3",50,-200,200,100,0,2000);

		ifstream inp;
		double Ev;
		double dt12,dt13,dt23,c1,c2,c3;
		inp.open(filePath1);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(dt13,c3*cost3);
		}
		inp.close();
		inp.open(filePath2);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(dt13,c3*cost3);
		}
		inp.close();
		inp.open(filePath3);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(dt13,c3*cost3);
		}
		inp.close();
		Scatter->GetXaxis()->SetTitle("Dt13 [ns]");
		Scatter->GetYaxis()->SetTitle("Carica sul PMT3 [keV]");
		Scatter->Draw("ZColor");
	}
	if (name=="Dt13C2"){

		TH2F *Scatter= new TH2F("Scatter","Dt13 Vs C2",50,-200,200,100,0,2000);

		ifstream inp;
		double Ev;
		double dt12,dt13,dt23,c1,c2,c3;
		inp.open(filePath1);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(dt13,c2*cost2);
		}
		inp.close();
		inp.open(filePath2);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(dt13,c2*cost2);
		}
		inp.close();
		inp.open(filePath3);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(dt13,c2*cost2);
		}
		inp.close();
		Scatter->GetXaxis()->SetTitle("Dt13 [ns]");
		Scatter->GetYaxis()->SetTitle("Carica sul PMT2 [keV]");
		Scatter->Draw("ZColor");
	}

	if (name=="Dt23C2"){

		TH2F *Scatter= new TH2F("Scatter","Dt23 Vs C2",50,-200,200,100,0,2000);

		ifstream inp;
		double Ev;
		double dt12,dt13,dt23,c1,c2,c3;
		inp.open(filePath1);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(dt23,c2*cost2);
		}
		inp.close();
		inp.open(filePath2);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(dt23,c2*cost2);
		}
		inp.close();
		inp.open(filePath3);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(dt23,c2*cost2);
		}
		inp.close();
		Scatter->GetXaxis()->SetTitle("Dt23 [ns]");
		Scatter->GetYaxis()->SetTitle("Carica sul PMT2 [keV]");
		Scatter->Draw("ZColor");
	}
	if (name=="Dt23C3"){

		TH2F *Scatter= new TH2F("Scatter","Dt23 Vs C3",50,-200,200,100,0,2000);

		ifstream inp;
		double Ev;
		double dt12,dt13,dt23,c1,c2,c3;
		inp.open(filePath1);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(dt23,c3*cost3);
		}
		inp.close();
		inp.open(filePath2);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(dt23,c3*cost3);
		}
		inp.close();
		inp.open(filePath3);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(dt23,c3*cost3);
		}
		inp.close();
		Scatter->GetXaxis()->SetTitle("Dt23 [ns]");
		Scatter->GetYaxis()->SetTitle("Carica sul PMT3 [keV]");
		Scatter->Draw("ZColor");
	}
	if (name=="Dt23C1"){

		TH2F *Scatter= new TH2F("Scatter","Dt23 Vs C1",50,-200,200,100,0,2000);

		ifstream inp;
		double Ev;
		double dt12,dt13,dt23,c1,c2,c3;
		inp.open(filePath1);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(dt23,c1*cost1);
		}
		inp.close();
		inp.open(filePath2);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(dt23,c1*cost1);
		}
		inp.close();
		inp.open(filePath3);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(dt23,c1*cost1);
		}
		inp.close();
		Scatter->GetXaxis()->SetTitle("Dt23 [ns]");
		Scatter->GetYaxis()->SetTitle("Carica sul PMT11[keV]");
		Scatter->Draw("ZColor");
	}

	if (name=="Dt23C2C3"){

		TH3F *Scatter= new TH3F("Scatter","C2 Vs C3 Vs Dt23",100,0,700,100,0,1500,5,0,20);

		ifstream inp;
		double Ev;
		double dt12,dt13,dt23,c1,c2,c3;
		inp.open(filePath1);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(c2*cost2,c3*cost3,dt23);
		}
		inp.close();
		inp.open(filePath2);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(c2*cost2,c3*cost3,dt23);
		}
		inp.close();
		inp.open(filePath3);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(c2*cost2,c3*cost3,dt23);
		}
		inp.close();
		/*
		*/
		Scatter->GetXaxis()->SetTitle("Carica sul PMT2 [keV]");
		Scatter->GetYaxis()->SetTitle("Carica sul PMT3 [keV]");
		Scatter->GetYaxis()->SetTitle("Dt23 [ns]");
		Scatter->Draw("ZColor");
	}

	if (name=="Dt13C1C3"){

		TH3F *Scatter= new TH3F("Scatter","C2 Vs C3 Vs Dt23",100,0,700,100,0,1500,5,0,20);

		ifstream inp;
		double Ev;
		double dt12,dt13,dt23,c1,c2,c3;
		inp.open(filePath1);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(c2*cost2,c3*cost3,dt23);
		}
		inp.close();
		inp.open(filePath2);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(c2*cost2,c3*cost3,dt23);
		}
		inp.close();
		inp.open(filePath3);
		while(!inp.eof()){
			inp >> Ev >> dt12 >> dt13 >> dt23 >> c1 >> c2 >> c3;
			Scatter->Fill(c2*cost2,c3*cost3,dt23);
		}
		inp.close();
		/*
		*/
		Scatter->GetXaxis()->SetTitle("Carica sul PMT2 [keV]");
		Scatter->GetYaxis()->SetTitle("Carica sul PMT3 [keV]");
		Scatter->GetYaxis()->SetTitle("Dt23 [ns]");
		Scatter->Draw("ZColor");
	}


}
