/*C/C++ Includes{{{*/
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>
#include <cmath>
//#include "Rtypes.h"
//#include "math.h"

/*}}}*/
/*ROOT Includes{{{*/
#include <TSystem.h>
#include <TString.h>
#include <TStyle.h>
#include <Riostream.h>
#include "TObjString.h"
#include <TNamed.h>
#include <TPRegexp.h>
#include <TObjArray.h>
#include <TChain.h>
#include <TMath.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TROOT.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TDatime.h>
#include <TError.h>
#include <TVirtualFitter.h>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TCut.h>
#include <TMultiGraph.h>
#include <TCutG.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TApplication.h>
#include <Rtypes.h>
#include <TTree.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
//#include <TMatrix.h>
/*}}}*/
using namespace std;

void plot_bin(const TString& TPol,const TString& energy,const int Q2_flag, const int xb_flag, const int t_flag){
	/*Define{{{*/
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
	double N_Err[12], PP_Err[12], PM_Err[12], MP_Err[12], MM_Err[12], AVG_Err[12];
	int dum;
	
	double N_Min = 1e35, N_Max = 1e-35;
	/*Define}}}*/
	
	/*Read{{{*/
	for(int i=0;i<12;i++){
		input_file  >> dum >> dum >> Q2[i]  >> xb[i]  >> t[i]  >> phi[i]  >> Astat[i] >> N[i]
			>> BSA[i]  >> TSA[i]  >> DSA[i]
			>> PP[i] >> PM[i] >> MP[i] >> MM[i] >> AVG[i];

		if(!(isnan( phi[i])) && !(isnan( Astat[i])) && Astat[i]>0){
			if(energy!="11p8"){//a mistake I gave in skim.C so it is just a temp solution. will remove later, 04/16/2015 11am
				N[i] *= 2.0;
				Astat[i] /=sqrt(2.0);
			}
			phi[i] = phi[i];	  

			BSA[i] /= Beam_Pol;
			TSA[i] /= (Dilution_Factor*Target_Pol*Target_Eff);
			DSA[i] /= (Dilution_Factor*Target_Pol*Target_Eff*Beam_Pol);
	
			BSA_Err[i] = Astat[i]*sqrt(1-pow(BSA[i]*Beam_Pol,2))/(Beam_Pol); 
			TSA_Err[i] = Astat[i]*sqrt(1-pow(TSA[i]*Dilution_Factor*Target_Pol*Target_Eff,2))/(Dilution_Factor*Target_Pol*Target_Eff); 
			DSA_Err[i] = Astat[i]*sqrt(1-pow(DSA[i]*Dilution_Factor*Target_Pol*Target_Eff*Beam_Pol,2))/(Dilution_Factor*Target_Pol*Target_Eff*Beam_Pol); 
			
			N_Err[i] = 1/sqrt(N[i]);
			if(N[i]>N_Max) N_Max = N[i];
			if(N[i]<N_Min) N_Min = N[i];
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
	/*Read}}}*/

	TString titlename;
	TLatex *t1 = new TLatex();
	t1->SetTextColor(6);
	t1->SetNDC();

	gStyle->SetOptStat(0);

	/*Plot1{{{*/

	TString energy0 = "";
	if(energy=="11p8")
		energy0 = "8.8+11";
	else
		energy0 = energy;

	TLine *bb = new TLine(0.0,0.0,360,0);
	bb->SetLineStyle(7);
	bb->SetLineColor(3);

	double y_pos = 0.30;
	titlename.Form("#scale[1.0]{(#%d) %2.1f < Q^{2} < %2.1f}",Q2_flag,Q2_min,Q2_max);
	t1->DrawLatex(0.20,y_pos,titlename);
	titlename.Form("#scale[1.0]{(#%d) %3.2f < xb < %3.2f}",xb_flag,xb_min,xb_max);
	t1->DrawLatex(0.20,y_pos-0.06,titlename);
	titlename.Form("#scale[1.0]{(#%d) %3.2f < t < %3.2f}",t_flag,t_min,t_max);
	t1->DrawLatex(0.20,y_pos-0.12,titlename);


	TCanvas *c1 = new TCanvas("c1","c1",1400,800);
	c1->Divide(3,1);
	//TH2F *h1 = new TH2F("h1","", 100, 0., 360, 100, min-min*0.2, max+max*0.2);
	TH2F *h1 = new TH2F("h1","", 100, 0., 360, 100, -1.5,1.5);
	h1->SetXTitle("#phi (degree)" );
	h1->SetYTitle(Form("A_{LU}"));
	h1->GetYaxis()->CenterTitle(1);
	h1->GetXaxis()->CenterTitle(1);
	h1->GetXaxis()->SetTitleSize(0.1);
	h1->GetYaxis()->SetTitleSize(0.1);
	h1->GetXaxis()->SetTitleOffset(0.6);
	h1->GetYaxis()->SetTitleOffset(0.6);
    c1->cd(1);
	h1->Draw();
	bb->Draw("same");

	TGraphErrors *g1 = new TGraphErrors(12,phi, BSA,0,BSA_Err);
	g1->Draw("*same");
	g1->SetMarkerColor(2);
	g1->SetMarkerStyle(20);
	g1->SetMarkerSize(1.7);
	g1->SetTitle();
	titlename.Form("#scale[1.0]{(#%d) %2.1f < Q^{2} < %2.1f}",Q2_flag,Q2_min,Q2_max);
	t1->DrawLatex(0.20,y_pos,titlename);
	titlename.Form("#scale[1.0]{(#%d) %3.2f < xb < %3.2f}",xb_flag,xb_min,xb_max);
	t1->DrawLatex(0.20,y_pos-0.06,titlename);
	titlename.Form("#scale[1.0]{(#%d) %3.2f < t < %3.2f}",t_flag,t_min,t_max);
	t1->DrawLatex(0.20,y_pos-0.12,titlename);
 
   	c1->cd(2);
	TH2F *h12 = (TH2F*)h1->Clone();
	h12->SetYTitle(Form("A_{U%s}", TPol.Data()));
	h12->Draw();
	bb->Draw("same");

	TGraphErrors *g2 = new TGraphErrors(12,phi,TSA,0,TSA_Err);
	g2->Draw("*same");
	g2->SetMarkerColor(4);
	g2->SetMarkerStyle(21);
	g2->SetMarkerSize(1.7);
	g2->SetTitle();

	titlename.Form("#scale[1.0]{(#%d) %2.1f < Q^{2} < %2.1f}",Q2_flag,Q2_min,Q2_max);
	t1->DrawLatex(0.20,y_pos,titlename);
	titlename.Form("#scale[1.0]{(#%d) %3.2f < xb < %3.2f}",xb_flag,xb_min,xb_max);
	t1->DrawLatex(0.20,y_pos-0.06,titlename);
	titlename.Form("#scale[1.0]{(#%d) %3.2f < t < %3.2f}",t_flag,t_min,t_max);
	t1->DrawLatex(0.20,y_pos-0.12,titlename);


	c1->cd(3);
	TH2F *h13 = (TH2F*)h1->Clone();
	h13->SetYTitle(Form("A_{L%s}", TPol.Data()));
	h13->Draw();
	bb->Draw("same");

	TGraphErrors *g3 = new TGraphErrors(12,phi,DSA,0,DSA_Err);
	g3->Draw("*same");
	g3->SetMarkerColor(1);
	g3->SetMarkerStyle(22);
	g3->SetMarkerSize(1.7);
	g3->SetTitle();

	titlename.Form("#scale[1.0]{(#%d) %2.1f < Q^{2} < %2.1f}",Q2_flag,Q2_min,Q2_max);
	t1->DrawLatex(0.20,y_pos,titlename);
	titlename.Form("#scale[1.0]{(#%d) %3.2f < xb < %3.2f}",xb_flag,xb_min,xb_max);
	t1->DrawLatex(0.20,y_pos-0.06,titlename);
	titlename.Form("#scale[1.0]{(#%d) %3.2f < t < %3.2f}",t_flag,t_min,t_max);
	t1->DrawLatex(0.20,y_pos-0.12,titlename);

	/*Plot1}}}*/
	
	/*Plot2{{{*/
	TCanvas *c2 = new TCanvas("c2","c2",1200,800);
/*	c2->Divide(1,2);
	c2->cd(1);
    TH2F *h2 = (TH2F*) h1->Clone();
	h2->Draw("");
	g1->Draw("p");
	g2->Draw("p");
	g3->Draw("p");
    bb->Draw("p");
//	l1->Draw();
	c2->cd(2);*/
	TH2F *h3 = new TH2F("h3","", 100, 0., 360, 100, N_Min*0.1, N_Max*2);
	h3->SetXTitle("#phi (degree)" );
	h3->SetYTitle(Form("Counts"));
	h3->GetYaxis()->CenterTitle(1);
	h3->GetXaxis()->CenterTitle(1);
	h3->GetXaxis()->SetTitleSize(0.1);
	h3->GetYaxis()->SetTitleSize(0.1);
	h3->GetXaxis()->SetTitleOffset(0.6);
	h3->GetYaxis()->SetTitleOffset(0.6);
	h3->Draw();

	TGraphErrors *g4 = new TGraphErrors(12,phi, N,0,N_Err);
	g4->Draw("*same");
	g4->SetMarkerColor(2);
	g4->SetMarkerStyle(20);
	g4->SetMarkerSize(1.7);
	g4->SetTitle();

	titlename.Form("#scale[1.2]{(#%d) %2.1f < Q^{2} < %2.1f}",Q2_flag,Q2_min,Q2_max);
	t1->DrawLatex(0.35,y_pos,titlename);
	titlename.Form("#scale[1.2]{(#%d) %3.2f < xb < %3.2f}",xb_flag,xb_min,xb_max);
	t1->DrawLatex(0.35,y_pos-0.06,titlename);
	titlename.Form("#scale[1.2]{(#%d) %3.2f < t < %3.2f}",t_flag,t_min,t_max);
	t1->DrawLatex(0.35,y_pos-0.12,titlename);
	
	/*Plot}}}*/

}
