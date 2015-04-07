void plot_bin(const TString & TPol,const int Q2_flag, const int xb_flag, const int t_flag){
	const TString Target = "n";	
	//Q2_Binning:
	const double Q2_Bin[6] = {1.0,1.5,2.0,3.0,4.5,7.0};
	const double Q2_min = Q2_Bin[Q2_flag];
	const double Q2_max = Q2_Bin[Q2_flag+1];
	//xb_Binning:
	const double xb_Bin[6] = {0.1,0.2,0.3,0.4,0.5,0.7};
	const double xb_min = xb_Bin[xb_flag];
	const double xb_max = xb_Bin[xb_flag+1];
	//t_Binning:
	const double t_Bin[7] = {-2.,-0.7,-0.5,-0.4,-0.3,-0.2,-0.1};
	const double t_min = t_Bin[t_flag];
	const double t_max = t_Bin[t_flag+1];
	//phi_Binning:
	//double Phi_Bin[13] ={0,30,60,90,120,150,180,210,240,270,300,330,360};

	TString filename_L = Form("./database_11p8/%s_%s_%d_%d_%d.dat",TPol.Data(), Target.Data(), Q2_flag, xb_flag, t_flag);
	//TString filename_L = Form("./database_11/%s_%s_%d_%d_%d.dat",TPol.Data(), Target.Data(), Q2_flag, xb_flag, t_flag);
	//TString filename_L = Form("./database_8/%s_%s_%d_%d_%d.dat",TPol.Data(), Target.Data(), Q2_flag, xb_flag, t_flag);
	ifstream input_file(filename_L);

	double Q2[12], xb[12], t[12], phi[12], Astat[12], TSA[12], BSA[12], DSA[12], TSA_Err[12], BSA_Err[12], DSA_Err[12];
        int dum;
	for(int i=0;i<12;i++){
		input_file  >> dum >> dum >> Q2[i]  >> xb[i]  >> t[i]  >> phi[i]  >> Astat[i]  >> BSA[i]  >> TSA[i]  >> DSA[i];
		//input_file  >> dum >> dum >> Q2[i]  >> xb[i]  >> t[i]  >> phi[i]  >> Astat[i]  >> BSA[i]  >> TSA[i]  >> DSA[i];

                if(TPol=="Ty"){
                    phi[i] *=180/3.1415926;
                    t[i] /=180./3.1415926;
                }

		if(!(isnan(Astat[i])) && Astat[i]>0){
			phi[i] = phi[i];	  
                        BSA[i] /= 1.0;
                        TSA[i] /= 1.0;
                        DSA[i] /= 1.0;
			BSA_Err[i] = BSA[i] * Astat[i]/1.0; 
			TSA_Err[i] = TSA[i] * Astat[i]/1.0; 
			DSA_Err[i] = DSA[i] * Astat[i]/1.0; 
		}
		else{
			phi[i] = -100;	  
			BSA[i] = -2.;
			TSA[i] = -2.;
			DSA[i] = -2.;
			BSA_Err[i] = 0.0;
			TSA_Err[i] = 0.0;
			DSA_Err[i] = 0.0;
		}
	}
	input_file.close();

	TString titlename;
	TLatex *t1 = new TLatex();
	t1->SetNDC();
        
        gStyle->SetOptStat(0);

        TH2F *h1 = new TH2F("h1","", 100, 0., 360, 100, -0.9, 0.9);
        h1->SetXTitle("#phi");
        h1->SetYTitle(Form("A_{%s}", TPol.Data()));
        h1->Draw();
/*
	Float_t a[2]={0,360};
	Float_t b[2]={-0.9,0.9};
	TGraph *g0 = new TGraph(2,a,b);
	g0->GetYaxis()->SetRangeUser(-1,1);
	g0->GetXaxis()->SetRangeUser(0.1,0.7);
	g0->GetXaxis()->SetTitle("#phi");
	g0->GetYaxis()->SetTitle(Form("A_{%s}", TPol.Data()));
	g0->GetYaxis()->SetTitleOffset(1.3);
	g0->GetXaxis()->SetNdivisions(506);
	g0->SetTitle(0);
	g0->Draw("AP");
*/
	TGraphErrors *g1 = new TGraphErrors(12,phi, BSA,0,BSA_Err);
	g1->Draw("*same");
	g1->SetMarkerColor(2);
	g1->SetMarkerStyle(20);
	g1->SetMarkerSize(1.0);
	g1->SetTitle();

	TGraphErrors *g2 = new TGraphErrors(12,phi,TSA,0,TSA_Err);
	g2->Draw("*same");
	g2->SetMarkerColor(4);
	g2->SetMarkerStyle(21);
	g2->SetMarkerSize(1.0);
	g2->SetTitle();

	TGraphErrors *g3 = new TGraphErrors(12,phi,DSA,0,DSA_Err);
	g3->Draw("*same");
	g3->SetMarkerColor(6);
	g3->SetMarkerStyle(22);
	g3->SetMarkerSize(1.0);
	g3->SetTitle();
	 
  	if(t_flag==0){
	   titlename.Form("#scale[2.0]{%2.1f < Q^{2} < %2.1f}",Q2_min,Q2_max);
	   t1->DrawLatex(0.5,0.90,titlename);
	   titlename.Form("#scale[2.0]{%3.2f < xb < %3.2f}",xb_min,xb_max);
	   t1->DrawLatex(0.5,0.75,titlename);
	   titlename.Form("#scale[2.0]{%3.2f < t < %3.2f}",t_min,t_max);
	   t1->DrawLatex(0.5,0.6,titlename);
	}
	else {
	   titlename.Form("#scale[2.0]{%2.1f < Q^{2} < %2.1f}",Q2_min,Q2_max);
	   t1->DrawLatex(0.15,0.90,titlename);
	   titlename.Form("#scale[2.0]{%3.2f < xb < %3.2f}",xb_min,xb_max);
	   t1->DrawLatex(0.15,0.75,titlename);
	   titlename.Form("#scale[2.0]{%3.2f < t < %3.2f}",t_min,t_max);
	   t1->DrawLatex(0.15,0.6,titlename);
		   }

	TLine *bb = new TLine(0.0,0.0,360,0);
        bb->SetLineStyle(7);
	bb->SetLineColor(3);
	bb->Draw("same");

}
