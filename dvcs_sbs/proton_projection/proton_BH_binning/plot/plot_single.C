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
#include <TLine.h>
//#include <TMatrix.h>
/*}}}*/
using namespace std;

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

	double Q2[12], xb[12], t[12], phi[12], Astat[12], N[12], N_Err[12], XSp[12], XSm[12], XSp_Err[12], XSm_Err[12];
	int dum;
	double min = 1e30, max = -1e30;

	for(int i=0;i<12;i++){
		input_file  >> dum >> dum >> Q2[i]  >> xb[i]  >> t[i]  >> phi[i]  >> Astat[i] >> N[i]  >> XSp[i]  >> XSm[i];

		if(!(isnan( phi[i])) && !(isnan( Astat[i])) && Astat[i]>0){
			phi[i] = phi[i];	  
            XSp[i] *= 2.0;//now take it as the total
			XSp_Err[i] = XSp[i] * Astat[i]/1.0; 
			XSm_Err[i] = XSm[i] * Astat[i]/1.0; 
            N_Err[i] = sqrt(N[i]);
		   	
			if(N[i] > max) max = N[i];
			if(N[i] < min) min = N[i];
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
 
	TCanvas *c1 = new TCanvas("c1","c1", 1200,800);
	c1->Divide(1,2);
    c1->cd(1);
	gPad->SetLogy(1);
	TH2F *h1 = new TH2F("h1","", 100, 0., 360, 100, 1e-3, 1e-1);
	h1->SetXTitle("#phi");
	h1->SetYTitle(Form("#sigma_{BH} (nb/GeV^{4})"));
	h1->GetXaxis()->CenterTitle(1);
	h1->GetXaxis()->SetTitleSize(0.10);
	h1->GetXaxis()->SetTitleOffset(0.50);
	h1->GetYaxis()->CenterTitle(1);
	h1->GetYaxis()->SetTitleSize(0.10);
	h1->GetYaxis()->SetTitleOffset(0.50);
	h1->Draw();

	TGraphErrors *g1 = new TGraphErrors(12,phi, XSp,0,XSp_Err);
	g1->Draw("*same");
	g1->SetMarkerColor(4);
	g1->SetMarkerStyle(20);
	g1->SetMarkerSize(1.2);
	g1->SetTitle();
	if(t_flag==0){
		titlename.Form("#scale[1.8]{%2.1f < Q^{2} < %2.1f}",Q2_min,Q2_max);
		t1->DrawLatex(0.5,0.90,titlename);
		titlename.Form("#scale[1.8]{%3.2f < xb < %3.2f}",xb_min,xb_max);
		t1->DrawLatex(0.5,0.75,titlename);
		titlename.Form("#scale[1.8]{%3.2f < t < %3.2f}",t_min,t_max);
		t1->DrawLatex(0.5,0.6,titlename);
	}
	else {
		titlename.Form("#scale[1.3]{%2.1f < Q^{2} < %2.1f}",Q2_min,Q2_max);
		t1->DrawLatex(0.35,0.80,titlename);
		titlename.Form("#scale[1.3]{%3.2f < xb < %3.2f}",xb_min,xb_max);
		t1->DrawLatex(0.35,0.70,titlename);
		titlename.Form("#scale[1.3]{%3.2f < t < %3.2f}",t_min,t_max);
		t1->DrawLatex(0.35,0.6,titlename);
	}

	TLine *bb = new TLine(0.0,0.0,360,0);
	bb->SetLineStyle(7);
	bb->SetLineColor(3);
	bb->Draw("same");

	c1->cd(2);
	gPad->SetLogy(0);
	TH2F *h2 = new TH2F("h2","", 100, -10.0, 370,100,min-100, max+100);
	h2->SetXTitle("#phi");
	h2->SetYTitle("Counts (48 days)");
	h2->GetXaxis()->CenterTitle(1);
	h2->GetXaxis()->SetTitleSize(0.10);
	h2->GetXaxis()->SetTitleOffset(0.50);
	h2->GetYaxis()->CenterTitle(1);
	h2->GetYaxis()->SetTitleSize(0.10);
	h2->GetYaxis()->SetTitleOffset(0.50);
	h2->Draw();

	TGraphErrors *g2 = new TGraphErrors(12,phi,N,0,N_Err);
	g2->Draw("*same");
	g2->SetMarkerColor(2);
	g2->SetMarkerStyle(21);
	g2->SetMarkerSize(1.2);
	g2->SetTitle();

}
