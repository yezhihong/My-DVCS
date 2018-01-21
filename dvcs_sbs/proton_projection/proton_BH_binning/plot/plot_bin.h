void plot_bin(const TString & energy,const int Q2_flag, const int xb_flag, const int t_flag){
	const TString Target = "p";	
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

	TString filename = Form("../database_%s/%s_%d_%d_%d.dat",energy.Data(),Target.Data(), Q2_flag, xb_flag, t_flag);
	ifstream input_file(filename);

	double Q2[12], xb[12], t[12], phi[12], Astat[12], N[12], XSp[12], XSm[12], XSp_Err[12], XSm_Err[12];
	int dum;
	for(int i=0;i<12;i++){
		input_file  >> dum >> dum >> Q2[i]  >> xb[i]  >> t[i]  >> phi[i]  >> Astat[i] >> N[i]  >> XSp[i]  >> XSm[i];

		if(!(isnan( phi[i])) && !(isnan( Astat[i])) && Astat[i]>0){
			phi[i] = phi[i];	  
            XSp[i] *= 2.0;//now take it as the total
			XSp_Err[i] = XSp[i] * Astat[i]/1.0; 
			XSm_Err[i] = XSm[i] * Astat[i]/1.0; 
		}
		else{
			phi[i] = -100;	  
			XSp[i] = -2.;
			XSm[i] = -2.;
			XSp_Err[i] = 0.0;
			XSm_Err[i] = 0.0;
		}
	}
	input_file.close();

	TString titlename;
	TLatex *t1 = new TLatex();
	t1->SetTextColor(6);
	t1->SetNDC();
	gStyle->SetOptStat(0);
  
	TH2F *h1 = new TH2F("h1","", 100, 0., 360, 100, 1e-5, 1.);
	h1->SetXTitle("#phi");
	h1->SetYTitle(Form("#sigma_{BH} (nb/GeV^{4})"));
	if(t_flag!=0){
		h1->GetYaxis()->SetLabelSize(0.0);
		h1->SetYTitle(""); 
	}

	h1->GetXaxis()->CenterTitle(1);
	h1->GetXaxis()->SetTitleSize(0.12);
	h1->GetXaxis()->SetTitleOffset(0.75);
	h1->GetYaxis()->CenterTitle(1);
	h1->GetYaxis()->SetTitleSize(0.10);
	h1->GetYaxis()->SetTitleOffset(1.00);
	h1->Draw();

	TGraphErrors *g1 = new TGraphErrors(12,phi, XSp,0,XSp_Err);
	g1->Draw("*same");
	g1->SetMarkerColor(4);
	g1->SetMarkerStyle(20);
	g1->SetMarkerSize(1.2);
	g1->SetTitle();
/*
	TGraphErrors *g2 = new TGraphErrors(12,phi,XSm,0,XSm_Err);
	g2->Draw("*same");
	g2->SetMarkerColor(2);
	g2->SetMarkerStyle(21);
	g2->SetMarkerSize(1.2);
	g2->SetTitle();
*/
	if(t_flag==0){
		titlename.Form("#scale[1.8]{%2.1f < Q^{2} < %2.1f}",Q2_min,Q2_max);
		t1->DrawLatex(0.5,0.90,titlename);
		titlename.Form("#scale[1.8]{%3.2f < xb < %3.2f}",xb_min,xb_max);
		t1->DrawLatex(0.5,0.75,titlename);
		titlename.Form("#scale[1.8]{%3.2f < t < %3.2f}",t_min,t_max);
		t1->DrawLatex(0.5,0.6,titlename);
	}
	else {
		titlename.Form("#scale[1.8]{%2.1f < Q^{2} < %2.1f}",Q2_min,Q2_max);
		t1->DrawLatex(0.15,0.90,titlename);
		titlename.Form("#scale[1.8]{%3.2f < xb < %3.2f}",xb_min,xb_max);
		t1->DrawLatex(0.15,0.75,titlename);
		titlename.Form("#scale[1.8]{%3.2f < t < %3.2f}",t_min,t_max);
		t1->DrawLatex(0.15,0.6,titlename);
	}

	TLine *bb = new TLine(0.0,0.0,360,0);
	bb->SetLineStyle(7);
	bb->SetLineColor(3);
	bb->Draw("same");

}
