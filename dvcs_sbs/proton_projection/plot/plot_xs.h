void plot_xs(const TString & TPol,const TString & energy,const int Q2_flag, const int xb_flag, const int t_flag){
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
	const double Beam_Pol = 1.0;

	TString filename_L = Form("../new_database_%s/%s_%s_%d_%d_%d.dat",energy.Data(),TPol.Data(), Target.Data(), Q2_flag, xb_flag, t_flag);
	ifstream input_file(filename_L);

	double Q2[12], xb[12], t[12], phi[12], Astat[12], TSA[12], BSA[12], DSA[12], TSA_Err[12], BSA_Err[12], DSA_Err[12];
	double N[12], PP[12], PM[12], MP[12], MM[12], AVG[12];
	double PP_Err[12], PM_Err[12], MP_Err[12], MM_Err[12], AVG_Err[12];
	int dum;
	double min = 1e35, max = -1e35;
	for(int i=0;i<12;i++){
		input_file  >> dum >> dum >> Q2[i]  >> xb[i]  >> t[i]  >> phi[i]  >> Astat[i] >> N[i]
			>> BSA[i]  >> TSA[i]  >> DSA[i]
			>> PP[i] >> PM[i] >> MP[i] >> MM[i] >> AVG[i];

		if(!(isnan( phi[i])) && !(isnan( Astat[i])) && Astat[i]>0){
			phi[i] = phi[i];	  
			BSA[i] /= 1.0;
			TSA[i] /= 1.0;
			DSA[i] /= 1.0;
			BSA_Err[i] = sqrt(1-BSA[i]*BSA[i]) * Astat[i]/1.0 * Target_Pol; 
			TSA_Err[i] = sqrt(1-TSA[i]*TSA[i]) * Astat[i]/1.0 * Beam_Pol; 
			DSA_Err[i] = sqrt(1-DSA[i]*DSA[i]) * Astat[i]/1.0; 
//			DSA[i] /= 1.0;
//			BSA_Err[i] = BSA[i] * Astat[i]/1.0; 
//			TSA_Err[i] = TSA[i] * Astat[i]/1.0; 
//			DSA_Err[i] = DSA[i] * Astat[i]/1.0;
//
            PP_Err[i] = PP[i] * Astat[i];
            PM_Err[i] = PM[i] * Astat[i];
            MP_Err[i] = MP[i] * Astat[i];
            MM_Err[i] = MM[i] * Astat[i];
			AVG[i] *=4.0;//the total cross section
            AVG_Err[i] = AVG[i] * Astat[i];

			if(AVG[i]>max) max=AVG[i];
			if(PP[i]<min) min=PP[i];
			if(PM[i]<min) min=PM[i];
			if(MP[i]<min) min=MP[i];
			if(MM[i]<min) min=MM[i];
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
	t1->SetTextColor(6);
	t1->SetNDC();
        
        gStyle->SetOptStat(0);

        TH2F *h1 = new TH2F("h1","", 100, 0., 360, 100, 5e-5, 1e-0);
		  
		h1->SetYTitle(""); 
		h1->SetXTitle(""); 
        if(t_flag==2)
			h1->SetXTitle("#phi (degree)" );
        if(xb_flag==2)
			h1->SetYTitle(Form("#sigma_{BH+DVCS}^{%s} (nb/GeV^{4})", TPol.Data()));
		if(t_flag!=0){
		  h1->GetYaxis()->SetLabelSize(0.0);
		  h1->SetYTitle(""); 
		}
		h1->GetYaxis()->CenterTitle(1);
		h1->GetXaxis()->CenterTitle(1);
		h1->GetXaxis()->SetTitleSize(0.1);
		h1->GetYaxis()->SetTitleSize(0.1);
		h1->GetXaxis()->SetTitleOffset(0.85);
		h1->GetYaxis()->SetTitleOffset(0.85);

        h1->Draw();
/*
		TGraphErrors *g0 = new TGraphErrors(12,phi, AVG,0,AVG_Err);
		g0->Draw("*same");
		g0->SetMarkerColor(1);
		g0->SetMarkerStyle(20);
		g0->SetMarkerSize(1.2);
		g0->SetTitle();
*/
		TGraphErrors *g1 = new TGraphErrors(12,phi, PP,0,PP_Err);
		g1->Draw("*same");
		g1->SetMarkerColor(2);
		g1->SetMarkerStyle(21);
		g1->SetMarkerSize(1.2);
		g1->SetTitle();

		TGraphErrors *g2 = new TGraphErrors(12,phi, PM,0,PM_Err);
		g2->Draw("*same");
		g2->SetMarkerColor(4);
		g2->SetMarkerStyle(22);
		g2->SetMarkerSize(1.2);
		g2->SetTitle();

		TGraphErrors *g3 = new TGraphErrors(12,phi, MP,0,MP_Err);
		g3->Draw("*same");
		g3->SetMarkerColor(7);
		g3->SetMarkerStyle(23);
		g3->SetMarkerSize(1.2);
		g3->SetTitle();

		TGraphErrors *g4 = new TGraphErrors(12,phi, MM ,0,MM_Err);
		g4->Draw("*same");
		g4->SetMarkerColor(8);
		g4->SetMarkerStyle(24);
		g4->SetMarkerSize(1.2);
		g4->SetTitle();

		double y_pos = 0.90;
		if(Q2_flag<3)
			   y_pos = 0.30;
		if(Q2_flag<3 && xb_flag==4)
			   y_pos = 0.55;

		if(t_flag==0){
			titlename.Form("#scale[1.5]{%2.1f < Q^{2} < %2.1f}",Q2_min,Q2_max);
			t1->DrawLatex(0.5,y_pos,titlename);
			titlename.Form("#scale[1.5]{%3.2f < xb < %3.2f}",xb_min,xb_max);
			t1->DrawLatex(0.5,y_pos-0.12,titlename);
			titlename.Form("#scale[1.5]{%3.2f < t < %3.2f}",t_min,t_max);
			t1->DrawLatex(0.5,y_pos-0.24,titlename);
		}
		else {
			titlename.Form("#scale[1.5]{%2.1f < Q^{2} < %2.1f}",Q2_min,Q2_max);
			t1->DrawLatex(0.15,y_pos,titlename);
			titlename.Form("#scale[1.5]{%3.2f < xb < %3.2f}",xb_min,xb_max);
			t1->DrawLatex(0.15,y_pos-0.12,titlename);
			titlename.Form("#scale[1.5]{%3.2f < t < %3.2f}",t_min,t_max);
			t1->DrawLatex(0.15,y_pos-0.24,titlename);
		}

}
