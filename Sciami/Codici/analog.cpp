
{

	using namespace std;

	#include <string>

	//qui i nomi dei file


	//"PMT4_triple1_14521_1719.dat"
	//"test.dat"
	//"triple2_PMT1_145_1840.dat"

	string file_input = "sciame_PMT6_2105_1820.dat";
	//string file_input = "triple2_PMT1_145_1840.dat";
	ifstream file_in;

	double ts,a1,a2;

	N_bin=150;


	TH1D *h_a1= new TH1D("h_a1","ampiezza PMT 1",N_bin,0,1);
	TH1D *h_a1_s= new TH1D("h_a1_s","ampiezza PMT 1",N_bin,0,1);

	vector<double> v1,v2;






	file_in.open(file_input);
    if(file_in.is_open())
	{

	while( (!file_in.eof())  )
	{

		file_in >> ts >> a1 >> a2;
		if (a2<2.5)
		{
			h_a1->Fill(a1);
			v1.push_back (a1);
		}
		else
		{
			h_a1_s->Fill(a1);
			v2.push_back (a2);
		}


	}
	}

	// test ks unbinned 2 sample

	/*
	double * arr1 = &v1[0];
    double * arr2 = &v2[0];

    double pvalue,Dn_1;

    ROOT::Math::GoFTest * goftest = new ROOT::Math::GoFTest(v1.size(), arr1, v2.size(), arr2);

    pvalue = goftest-> KolmogorovSmirnov2SamplesTest();
    Dn_1 = goftest-> KolmogorovSmirnov2SamplesTest("t");

    cout << pvalue <<"    " << Dn_1<<"\n";
	*/



}
