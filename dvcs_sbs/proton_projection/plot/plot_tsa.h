void plot_tsa(const TString& TPol,const TString& energy,const int Q2_flag, const int xb_flag, const int t_flag){
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
	const double Target_Pol = 0.60;
	const double Beam_Pol = 0.9;
	const double Target_Eff = 0.865;
	const double Dilution_Factor = 0.9;

	TString filename_L = Form("../new_database_%s/%s_%s_%d_%d_%d.dat",energy.Data(),TPol.Data(), Target.Data(), Q2_flag, xb_flag, t_flag);
	ifstream input_file(filename_L);

	double Q2[12], xb[12], t[12], phi[12], Astat[12], TSA[12], BSA[12], DSA[12], TSA_Err[12], BSA_Err[12], DSA_Err[12];
	double N[12], PP[12], PM[12], MP[12], MM[12], AVG[12];
	int dum;
	for(int i=0;i<12;i++){
		input_file  >> dum >> dum >> Q2[i]  >> xb[i]  >> t[i]  >> phi[i]  >> Astat[i] >> N[i]
			>> BSA[i]  >> TSA[i]  >> DSA[i]
			>> PP[i] >> PM[i] >> MP[i] >> MM[i] >> AVG[i];
		//input_file  >> dum >> dum >> Q2[i]  >> xb[i]  >> t[i]  >> phi[i]  >> Astat[i]  >> BSA[i]  >> TSA[i]  >> DSA[i];

		if(!(isnan( phi[i])) && !(isnan( Astat[i])) && Astat[i]>0){
			if(energy!="11p8"){//a mistake I gave in skim.C so it is just a temp solution. will remove later, 04/16/2015 11am
				N[i] *= 2.0;
				Astat[i] /=sqrt(2.0);
			}
			N[i] /= 4.0;
			Astat[i] = 1./sqrt(N[i]);

			phi[i] = phi[i];	  
			TSA[i] /= (Dilution_Factor*Target_Pol*Target_Eff);
			TSA_Err[i] = Astat[i]*sqrt(1-pow(TSA[i]*Dilution_Factor*Target_Pol*Target_Eff,2))/(Dilution_Factor*Target_Pol*Target_Eff); 
		}
		else{
			phi[i] = -100;	  
			TSA[i] = -2.;
			TSA_Err[i] = 0.0;
		}
	}
	input_file.close();

	TString titlename;
	TLatex *t1 = new TLatex();
	t1->SetTextColor(6);
	t1->SetNDC();

	gStyle->SetOptStat(0);

	double Y_Min = -1.01;
	double Y_Max = 1.01;
	if(TPol=="Ty")
	{
	   Y_Min = -0.5;
	   Y_Max = 1.5;
	}

	TH2F *h1 = new TH2F("h1","", 100, 0., 360, 100, Y_Min, Y_Max);

	h1->SetYTitle(""); 
	h1->SetXTitle(""); 
	if(t_flag==2)
		h1->SetXTitle("#phi (degree)" );
	if(xb_flag==2)
		h1->SetYTitle(Form("A_{%s}", TPol.Data()));
	if(t_flag!=0){
		h1->GetYaxis()->SetLabelSize(0.0);
		h1->SetYTitle(""); 
	}
	h1->GetYaxis()->CenterTitle(1);
	h1->GetXaxis()->CenterTitle(1);
	h1->GetXaxis()->SetTitleSize(0.12);
	h1->GetYaxis()->SetTitleSize(0.12);
	h1->GetXaxis()->SetTitleOffset(0.65);
	h1->GetXaxis()->SetTitleOffset(0.65);
	h1->Draw();

	TGraphErrors *g2 = new TGraphErrors(12,phi,TSA,0,TSA_Err);
	g2->Draw("*same");
	g2->SetMarkerColor(4);
	g2->SetMarkerStyle(21);
	g2->SetMarkerSize(1.2);
	g2->SetTitle();

	double y_pos = 0.90;
	if(TPol=="Tx"||TPol=="Ty")
		y_pos = 0.30;
	//	if(Q2_flag<3 && xb_flag==4)
	//		y_pos = 0.55;

	if(xb_flag==4)
		y_pos = 0.5;

	if(t_flag==0){
		titlename.Form("#scale[1.5]{(#%d) %2.1f < Q^{2} < %2.1f}",Q2_flag,Q2_min,Q2_max);
		t1->DrawLatex(0.4,y_pos,titlename);
		titlename.Form("#scale[1.5]{(#%d) %3.2f < xb < %3.2f}",xb_flag,xb_min,xb_max);
		t1->DrawLatex(0.4,y_pos-0.12,titlename);
		titlename.Form("#scale[1.5]{(#%d) %3.2f < t < %3.2f}",t_flag,t_min,t_max);
		t1->DrawLatex(0.4,y_pos-0.24,titlename);
	}
	else {
		titlename.Form("#scale[1.5]{(#%d) %2.1f < Q^{2} < %2.1f}",Q2_flag,Q2_min,Q2_max);
		t1->DrawLatex(0.11,y_pos,titlename);
		titlename.Form("#scale[1.5]{(#%d) %3.2f < xb < %3.2f}",xb_flag,xb_min,xb_max);
		t1->DrawLatex(0.11,y_pos-0.12,titlename);
		titlename.Form("#scale[1.5]{(#%d) %3.2f < t < %3.2f}",t_flag,t_min,t_max);
		t1->DrawLatex(0.11,y_pos-0.24,titlename);
	}

	TLine *bb = new TLine(0.0,0.0,360,0);
	bb->SetLineStyle(7);
	bb->SetLineColor(3);
	bb->Draw("same");

}
