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
#include "LHAPDF/LHAPDF.h"
//#include <TMatrix.h>
/*}}}*/

using namespace std;
const double PI = 3.1415926;

double max(double a, double b);
int makebin(TString target, TString energy, int Q2_flag, int x_flag);

/*int main{{{*/
int main(){
	TString target = "x"; cerr<<"--- Target (p or n) = "; cin >> target;
	TString energy = "x"; cerr<<"--- Energy (8 or 11 or 11p8) = "; cin >> energy;
	Int_t x_flag = 0;
	Int_t Q2_flag = 0;
	int err = -1000;
	for(int i=0;i<5;i++){
		for(int j=0;j<5;j++){
			Q2_flag = i;
			x_flag = j;
			err = makebin(target.Data(), energy.Data(), Q2_flag, x_flag);
		}
	}
	return err;
}
/*}}}*/

int makebin(TString target, TString energy, int Q2_flag, int x_flag){
	/*Initial Factors{{{*/
	if(target!="p" && target!="n"){
		cerr<<"I don't know this target flag!"<<endl;
		return -1;
	}
	if(energy!="11" && energy!="8" && energy!="11p8"){
		cerr<<"I don't know this energy flag!"<<endl;
		return -2;
	}

	const double dilute_factor = 0.2;
	const double polarization = 0.6; //60% polarization
	const double det_eff = 0.85; //85% detector efficiency for electrons and hadrons
	const double target_factor = 0.865; //neutron polarization 86.5%
	//	const double Asys = 0.006; 
	//const double Nsys = 2e+6; 

	//Temp, real bin size is determined by the events in each bin
	const double T_MIN  =0.1;
	const double T_MAX = 0.7;

	const double Phi_STEP = 30.;
	const double Phi_MIN = 0.0;
	const double Phi_MAX = 360;
	const int phi_bin = (int)((Phi_MAX-Phi_MIN)/Phi_STEP);
	/*}}}*/

	/*DVCS: Define Rootfile and variables{{{*/
	TString prefix = Form("./rootfiles_%s/", energy.Data());
	TString finalfile = Form("%s/%s_%d_%d_%s.root", prefix.Data(), target.Data(), Q2_flag,x_flag, energy.Data());
	cerr<<Form("--- Reading the rootfile: %s", finalfile.Data())<<endl;

	TFile *file = new TFile(finalfile.Data(),"r");
	TTree *T = (TTree*) file->Get("T");

	Double_t vertexz, E0; 
	Double_t ePx_ini, ePy_ini, ePz_ini;
	//Double_t	hP_ini, hPx_ini, hPy_ini,hPz_ini,hTheta_ini, hPhi_ini;
	Double_t ePx, ePy, ePz, gPx, gPy, gPz, hPx, hPy, hPz;
	Double_t eP_ini,eP, gP,hP;
	Double_t Q2, x, t, phi, XS_P, XS_M, XS_BHp, XS_BHm, PSF;
	Double_t ePhi_ini, ePhi, gPhi,hPhi;
	Double_t eTheta_ini,eTheta, gTheta,hTheta;
	Int_t Ngen,Nacc;

	Double_t eP_res, eTheta_res, ePhi_res, ePx_res, ePy_res, ePz_res;
	Double_t gP_res, gTheta_res, gPhi_res, gPx_res, gPy_res, gPz_res;
	Double_t MM = 0.0,MM_res = 0.0, W=0.0, Wp=0.0;
	Double_t e_acc_f= 0.0, e_acc_l= 0.0, g_acc_f=0, g_acc_l=0;
	Double_t weight, time;

	T->SetBranchAddress("vertexz", &vertexz);
	T->SetBranchAddress("E0", &E0);
	T->SetBranchAddress("eP_ini", &eP_ini);
	T->SetBranchAddress("ePx_ini", &ePx_ini);
	T->SetBranchAddress("ePy_ini", &ePy_ini);
	T->SetBranchAddress("ePz_ini", &ePz_ini);
	T->SetBranchAddress("eTheta_ini",&eTheta_ini);
	T->SetBranchAddress("ePhi_ini",&ePhi_ini);

	T->SetBranchAddress("eP", &eP);
	T->SetBranchAddress("ePx", &ePx);
	T->SetBranchAddress("ePy", &ePy);
	T->SetBranchAddress("ePz", &ePz);
	T->SetBranchAddress("eTheta",&eTheta);
	T->SetBranchAddress("ePhi",&ePhi);

	T->SetBranchAddress("gP", &gP);
	T->SetBranchAddress("gPx", &gPx);
	T->SetBranchAddress("gPy", &gPy);
	T->SetBranchAddress("gPz", &gPz);
	T->SetBranchAddress("gTheta", &gTheta);
	T->SetBranchAddress("gPhi", &gPhi);

	T->SetBranchAddress("hP", &hP);
	T->SetBranchAddress("hPx", &hPx);
	T->SetBranchAddress("hPy", &hPy);
	T->SetBranchAddress("hPz", &hPz);
	T->SetBranchAddress("hTheta",&hTheta);
	T->SetBranchAddress("hPhi",&hPhi);

	T->SetBranchAddress("Q2", &Q2);
	T->SetBranchAddress("x", &x);
	T->SetBranchAddress("t", &t);
	T->SetBranchAddress("phi", &phi);
	T->SetBranchAddress("XS_P", &XS_P);
	T->SetBranchAddress("XS_M", &XS_M);
	T->SetBranchAddress("XS_BHp", &XS_BHp);
	T->SetBranchAddress("XS_BHm", &XS_BHm);
	T->SetBranchAddress("PSF", &PSF);
	T->SetBranchAddress("Ngen", &Ngen);
	T->SetBranchAddress("Nacc", &Nacc);

	T->SetBranchAddress("W", &W);
	T->SetBranchAddress("Wp", &Wp);
	T->SetBranchAddress("MM", &MM);
	T->SetBranchAddress("MM_res", &MM_res);
	T->SetBranchAddress("e_acc_f", &e_acc_f);
	T->SetBranchAddress("e_acc_l", &e_acc_l);
	T->SetBranchAddress("g_acc_f", &g_acc_f);
	T->SetBranchAddress("g_acc_l", &g_acc_l);
	T->SetBranchAddress("eP_res", &eP_res);
	T->SetBranchAddress("ePx_res", &ePx_res);
	T->SetBranchAddress("ePy_res", &ePy_res);
	T->SetBranchAddress("ePz_res", &ePz_res);
	T->SetBranchAddress("ePhi_res", &ePhi_res);
	T->SetBranchAddress("eTheta_res", &eTheta_res);
	T->SetBranchAddress("gP_res", &gP_res);
	T->SetBranchAddress("gPx_res", &gPx_res);
	T->SetBranchAddress("gPy_res", &gPy_res);
	T->SetBranchAddress("gPz_res", &gPz_res);
	T->SetBranchAddress("gPhi_res", &gPhi_res);
	T->SetBranchAddress("gTheta_res", &gTheta_res);

	T->SetBranchAddress("weight", &weight);
	T->SetBranchAddress("time", &time);

	Long64_t N_entries=T->GetEntries();
	T->GetEntry(0);
	Long64_t N_gen = Ngen*2;//11GeV and 8.8 GeV settings should have the same total generated events but please check
	//cout<<"DVCS: total generated events number: "<<N_gen<<"--- and accepted: "<<N_entries<<endl;
	/*}}}*/

	/*Histograms{{{*/
	//t
	TH1F *h1 = new TH1F("h1","h1",1000,T_MIN,T_MAX);
	//Q2
	TH1F *h1Q2  =new TH1F("h1Q2","h1Q2",1000,1.,7.);
	//x 
	TH1F *h1x = new TH1F("h1x","h1x",1000,0.1,0.7);
	// phi
	TH1F *h1phi = new TH1F("h1phi","h1phi",1000.,Phi_MIN,Phi_MAX);
	//t 
	TH1F *h1t = new TH1F("h1t","h1t",1000,-2.0,-0.1);
	//Sigma
	TH1F *h1XSp = new TH1F("h1XSp","h1XSp (log10)",1000,-6.0,3.);
	TH1F *h1XSm = new TH1F("h1XSm","h1XSm (log10)",1000,-6.0,3.);
	/*}}}*/

	/*Binning{{{*/
	TString histoname;
	double dilution = 0.0;
	double N_out = 0.0, Astat=0.0;
	double phimin = 0.0, phimax = 0.0;
	double tmin = 0.0, tmax = 0.0;
	TString cut1,cut2;
	TCut cut;

	const int tbin = 6;
	//const double t1[6] = {0.1, 0.2, 0.3,0.4, 0.5, 0.7};
	const double t1[7] = {-2.0,-0.7, -0.5, -0.4,-0.3, -0.2, -0.1};
	//loop through phi and t bin
	for (Int_t i=0;i<tbin;i++){
		tmin = t1[i];
		tmax = t1[i+1];

		prefix = Form("./database_%s/",energy.Data());
		TString filename = Form("%s_%d_%d_%d.dat",target.Data(),Q2_flag,x_flag, i);
		TString new_filename = prefix + filename;
		ofstream outf(new_filename);

		for (Int_t j=0;j<phi_bin;j++){
			phimin = Phi_MIN + j*Phi_STEP;
			phimax = Phi_MIN + (j+1)*Phi_STEP;

			/*t Binning{{{*/
			h1Q2->Reset();
			h1x->Reset();
			h1phi->Reset();
			h1t->Reset();
			h1XSp->Reset();
			h1XSm->Reset();

			//Sigma_L is the average of ++/+-/-+/--, when calculating the overall statistics, we use the sum not the average
			cut1="(weight*(XS_BHp+XS_BHm))";
			cut2.Form("(phi>=%f&&phi<%f&&t>=%f&&t<%f&&W>2.)",phimin,phimax,tmin,tmax);
			cut2 = cut1 + "*" + cut2;
			cut = (TCut) (cut2);

			T->Project("h1Q2","Q2",cut);
			T->Project("h1x","x",cut);
			T->Project("h1phi","phi",cut);
			T->Project("h1t","t",cut);
			T->Project("h1XSp","log10(XS_BHp)",cut);
			T->Project("h1XSm","log10(XS_BHm)",cut);

			dilution = dilute_factor;//fix the value for now
			N_out = h1phi->GetSum()*(pow(polarization * target_factor * dilution,2) * det_eff);
			if(N_out>10)
				Astat = 1./sqrt(N_out);
			else
				Astat = -1.0;

			outf<<Form("%4d %4d %10.4f %10.4f %10.4f %10.4f %10.4f %12.4e %12.4e %12.4e",
					i,j,
					h1Q2->GetMean(),
					h1x->GetMean(),
					h1t->GetMean(),
					h1phi->GetMean(),
					Astat,
					N_out,
					pow(10., h1XSp->GetMean()),	
					pow(10., h1XSm->GetMean()))	
				<<endl;

			/*}}}*/

		}
		outf.close();
	}
	/*}}}*/

	/*Free{{{*/
	h1->Delete();
	h1Q2->Delete();
	h1x->Delete();
	h1phi->Delete();
	h1t->Delete();
	h1XSp->Delete();
	h1XSm->Delete();
	/*}}}*/
	return 0;
}

double max(double a, double b){
	if(a>b)
		return a;
	else
		return b;
}
