
{	
    int N=100000;

	double ly=48;
    double lx=40;
    double h1=20;
    double h2=20;
    double m1=0.6;
    double m2=-0.6;
    double q1=-24;
    double q2=72;
    double q1_2=-24.6;
    double q2_2=72.6;
    double theta;
    double phi;
    double x1,x2,x3,y1,y2,y3;

	gRandom = new TRandom3();
    gRandom->SetSeed (0);
    double pi=TMath::Pi();
    TF1 *f1 = new TF1("f1","cos(x)*cos(x)",0,pi/2);
    double r = f1->GetRandom();
	
    double doppie_12=0;
    double doppie_13=0;
    double triple=0;
    

    //per controllo
    TH2D * hist_xy_1 = new TH2D("hist_1", "xy_1",1000, -100 , 180, 1000,-100,150);
    hist_xy_1->SetTitle("Piano 1");
    hist_xy_1->GetXaxis()->SetTitle("x [cm]");
    hist_xy_1->GetYaxis()->SetTitle("y [cm]");
    TH2D * hist_xy_2 = new TH2D("hist_2", "xy_2",1000, -100 , 180, 1000,-100,150);
    hist_xy_2->SetTitle("Piano 2");
    hist_xy_2->GetXaxis()->SetTitle("x [cm]");
    hist_xy_2->GetYaxis()->SetTitle("y [cm]");
    TH2D * hist_xy_3 = new TH2D("hist_3", "xy_3",1000, -100 , 180, 1000,-100,150);
    hist_xy_3->SetTitle("Piano 3");
    hist_xy_3->GetXaxis()->SetTitle("x [cm]");
    hist_xy_3->GetYaxis()->SetTitle("y [cm]");
    int n=0;
    while(n<N)
    {

       	x1=gRandom->Uniform(lx);
        y1=gRandom->Uniform(ly);
        

        if (((x1>0)&&(y1>0)&&(x1<lx)&&(y1<ly)))
        {
            hist_xy_1->Fill(x1,y1);

            theta=f1->GetRandom();
            phi=gRandom->Uniform(2*pi);
            x2=x1+(h1*TMath::Tan(theta)*TMath::Cos(phi));
            y2=y1+(h1*TMath::Tan(theta)*TMath::Sin(phi));
            x3=x1+((h2+h1)*TMath::Tan(theta)*TMath::Cos(phi));
            y3=y1+((h2+h1)*TMath::Tan(theta)*TMath::Sin(phi));

        
            hist_xy_2->Fill(x2,y2);
            hist_xy_3->Fill(x3,y3);

            int k=0;
            if (((x2<lx+1) && (y2<ly) && (x2>1) && (y2>0)))
            {
                doppie_12+=1;
                k+=1;
            }
            if (((x3<lx) && (y3<ly) && (x3>0) && (y3>0)))
            {
                doppie_13+=1;
                k+=1;
            }
            if (k==2) triple+=1;




            n=n+1;



        }

        

    }

    cout<<"N   "<< N<<"\n";
    cout<<"1 & 2   "<< doppie_12<<"  "<< doppie_12/N<<"\n";
    cout<< "1 & 3   "<< doppie_13<<"  "<< doppie_13/N <<"\n";
    
    cout<< "1 & 2 & 3   "<< triple<<"   (1 & 2 & 3)/N   "<< triple/N << "\n";

    cout << "(1 & 2 & 3)/(1 & 2)   "<<triple/doppie_12<<"\n";
    cout << "(1 & 2 & 3)/(1 & 3)   "<<triple/doppie_13<<"\n";






    



    

   

    

                
        
}