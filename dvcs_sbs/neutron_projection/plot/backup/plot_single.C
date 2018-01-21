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

void plot_single(const TString & TPol,const int Q2_flag, const int xb_flag, const int t_flag){
	const TString Target = "n";	
	//Q2_Binning:
	const double Q2_Bin[6] = {1.0,1.5,2.0,3.0,4.5,0.7};
	const double Q2_min = Q2_Bin[Q2_flag];
	const double Q2_max = Q2_Bin[Q2_flag+1];
	//xb_Binning:
	const double xb_Bin[6] = {0.1,0.2,0.3,0.4,0.5,0.7};
	const double xb_min = xb_Bin[xb_flag];
	const double xb_max = xb_Bin[xb_flag+1];
	//t_Binning:
	const double t_Bin[6] = {0.1,0.2,0.3,0.4,0.5,0.7};
	const double t_min = t_Bin[t_flag];
	const double t_max = t_Bin[t_flag+1];
	//phi_Binning:
	//double Phi_Bin[13] ={0,30,60,90,120,150,180,210,240,270,300,330,360};

	TString filename_L = Form("./database_11/%s_%s_%d_%d_%d.dat",TPol.Data(), Target.Data(), Q2_flag, xb_flag, t_flag);
	ifstream input_file(filename_L);

	double Q2[12], xb[12], t[12], phi[12], Astat[12], TSA[12], BSA[12], DSA[12], TSA_Err[12], BSA_Err[12], DSA_Err[12];
        int dum;
	for(int i=0;i<12;i++){
		input_file >> dum >> dum  >> Q2[i]  >> xb[i]  >> t[i]  >> phi[i]  >> Astat[i]  >> BSA[i]  >> TSA[i]  >> DSA[i];

                if(TPol=="Ty"){
                    phi[i] *=180/3.1415926;
                    t[i] /=180./3.1415926;
                }

		if(!(isnan(Astat[i])) && Astat[i]>0){
			phi[i] = phi[i];	  
			BSA_Err[i] = BSA[i] * Astat[i]; 
			TSA_Err[i] = TSA[i] * Astat[i]; 
			DSA_Err[i] = DSA[i] * Astat[i]; 

                       cerr<<Form("--- Phi =%f, BSA = %f, TSA = %f, DSA = %f", phi[i], BSA[i], TSA[i], DSA[i])<<endl;
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
	
	TCanvas *c1 = new TCanvas("c1","c1",1600,1200);

        TH2F *h1 = new TH2F("h1","", 100, 0., 360, 100, -0.9, 0.9);
        h1->SetXTitle("#phi");
        h1->SetYTitle(Form("A_{%ss}", TPol.Data()));
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
	g1->SetMarkerSize(0.5);
	g1->SetTitle();

	TGraphErrors *g2 = new TGraphErrors(12,phi,TSA,0,TSA_Err);
	g2->Draw("*same");
	g2->SetMarkerColor(4);
	g2->SetMarkerStyle(21);
	g2->SetMarkerSize(0.5);
	g2->SetTitle();

	TGraphErrors *g3 = new TGraphErrors(12,phi,DSA,0,DSA_Err);
	g3->Draw("*same");
	g3->SetMarkerColor(6);
	g3->SetMarkerStyle(22);
	g3->SetMarkerSize(0.5);
	g3->SetTitle();
	 
  	if(xb_flag==5){
	   titlename.Form("#scale[2.0]{%d < Q^{2} < %d}",(int)Q2_min,(int)Q2_max);
	   t1->DrawLatex(0.15,0.90,titlename);
	   titlename.Form("#scale[2.0]{%3.2f < xb < %3.2f}",xb_min,xb_max);
	   t1->DrawLatex(0.15,0.75,titlename);
	   titlename.Form("#scale[2.0]{%3.2f < -t < %3.2f}",t_min,t_max);
	   t1->DrawLatex(0.15,0.6,titlename);
	}
	else {
	   titlename.Form("#scale[2.0]{%d < Q^{2} < %d}",(int)Q2_min,(int)Q2_max);
	   t1->DrawLatex(0.15,0.90,titlename);
	   titlename.Form("#scale[2.0]{%3.2f < xb < %3.2f}",xb_min,xb_max);
	   t1->DrawLatex(0.15,0.75,titlename);
	   titlename.Form("#scale[2.0]{%3.2f < -t < %3.2f}",t_min,t_max);
	   t1->DrawLatex(0.15,0.6,titlename);
		   }

	TLine *bb = new TLine(0.1,0.0,0.7,0);
	bb->Draw("same");

	c1->Print(Form("%s_%s_%d_%d_%d.png",TPol.Data(), Target.Data(), Q2_flag, xb_flag, t_flag));
}
