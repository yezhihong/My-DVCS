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
const double DEG = 180./3.1415926;
const double PI = 3.1415926;

double max(double a, double b);
int makebin(int target_flag, int Q2_flag, int t_flag);

/*int main{{{*/
int main(){
	Int_t target_flag = 3; cerr<<"--- Target (1->p, 3->he3) = "; cin >> target_flag;
	Int_t x_flag = 0;
	Int_t Q2_flag = 0;
	int err = -1000;
	for(int i=1;i<=5;i++){
			for(int j=1;j<=8;j++){
				Q2_flag = i;
				x_flag = j;
				err = makebin(target_flag, Q2_flag, x_flag);
			}
		}
		return err;
}
/*}}}*/

int makebin(int target_flag, int Q2_flag, int x_flag){
	TString target = "X";
	if(target_flag==1)
		target ="p";
	else if(target_flag==3)
		target ="n";
	else{
		cerr<<"I don't know this target flag!"<<endl;
		return -1;
	}

	/*DVCS: Define Rootfile and variables{{{*/
	TString prefix = "./skim_rootfiles_11p8/";
	TString finalfile = Form("%s/%s_%d_%d_11p8.root", prefix.Data(), target.Data(), Q2_flag,x_flag);
	cerr<<Form("--- Reading the rootfile: %s", finalfile.Data())<<endl;

	TFile *file = new TFile(finalfile.Data(),"r");
	TTree *T = (TTree*) file->Get("T");

	Double_t vertexz, E0; 
	Double_t ePx_ini, ePy_ini, ePz_ini, hPx_ini, hPy_ini,hPz_ini;
	Double_t ePx, ePy, ePz, gPx, gPy, gPz, hPx, hPy, hPz;
	Double_t eP_ini,hP_ini, eP, gP,hP;
	Double_t Q2, x, t, phi, XS_P, XS_M,PSF;
	Double_t ePhi_ini,hPhi_ini, ePhi, gPhi,hPhi;
	Double_t eTheta_ini,hTheta_ini, eTheta, gTheta,hTheta;
	Int_t Ngen,Nacc;

	Double_t eP_res, eTheta_res, ePhi_res, ePx_res, ePy_res, ePz_res;
	Double_t gP_res, gTheta_res, gPhi_res, gPx_res, gPy_res, gPz_res;
	Double_t MM = 0.0,MM_res = 0.0, W=0.0, Wp=0.0;
	Double_t e_acc_f= 0.0, e_acc_l= 0.0, g_acc_f=0, g_acc_l=0;
	Double_t weight_p, weight_m, time;

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
	
	T->SetBranchAddress("weight_p", &weight_p);
	T->SetBranchAddress("weight_m", &weight_m);
	T->SetBranchAddress("time", &time);

	Long64_t N_entries=T->GetEntries();
	T->GetEntry(0);
	Long64_t N_gen = Ngen*2;//11GeV and 8.8 GeV settings should have the same total generated events but please check
	cout<<"DVCS: total generated events number: "<<N_gen<<"--- and accepted: "<<N_entries<<endl;
	/*}}}*/

	const double dilute_factor = 0.2;
	const double polarization = 0.6; //60% polarization
	const double det_eff = 0.85; //85% detector efficiency for electrons and hadrons
	const double target_factor = 0.865; //neutron polarization 86.5%
//	const double Asys = 0.006; 
	const double Nsys = 2e+6; 

	//Temp, real bin size is determined by the events in each bin
	const double T_MIN  =0.1;
	const double T_MAX = 0.7;
	const double T_STEP = 0.1;
	const int t_bin = (int)((T_MAX-T_MIN)/T_STEP); 

	const double Phi_STEP = 30.;
	const double Phi_MIN = 0.0;
	const double Phi_MAX = 360;
	const int phi_bin = (int)((Phi_MAX-Phi_MIN)/Phi_STEP);

	prefix = "./databases/";
	TString filename = Form("%s_%d_%d.dat",target.Data(),Q2_flag,x_flag);
	TString new_filename = prefix + filename;
	ofstream outf_total(new_filename);

	TString histoname;
	outf_total << phi_bin << endl;
	double phimin = 0.0, phimax = 0.0;
	//loop through phi and t bin
	for (Int_t i=0;i<phi_bin;i++){
		phimin = Phi_MIN + i*Phi_STEP;
		phimax = Phi_MIN + (i+1)*Phi_STEP;
			
		/*Histograms{{{*/
		//phi
		TH1F *h1  =new TH1F("h1","h1",1000,Phi_MIN,Phi_MAX);
		//t
		TH1F *h2 = new TH1F("h2","h2",1000,T_MIN,T_MAX);
		/*}}}*/

		TString cut1,cut2;
		cut1="(weight_p+weight_m)";
		TCut cut(cut1);

		/*X Binning{{{*/
		cut2.Form("(phi>=%f&&phi<%f)",phimin,phimax);
		cut2 = cut1 + "*" + cut2;
		cut = cut2;
		h1->Reset();
		T->Project("h1","phi",cut);
		h2->Reset();
		T->Project("h2","t",cut);

		//double Norm_Factor = (pow(polarization * target_factor*dilute_factor,2) * det_eff);
		//double single_bin_raw = 1.0/Norm_Factor/pow(Asys,2);
		//Double_t total_eve = h2->GetSum(); 
		double single_bin_raw = Nsys/det_eff;

		/*Count How many t bins{{{*/
		double t1[1000];
		int new_bin = 0, nevent_t=0;
		t1[0] = T_MIN;

		int tb_max =0;
		while(tb_max<1000){
			new_bin++;
			nevent_t = 0;

			//	cerr<<Form("--- Start with t=%f, tb=%d",h2->GetBinCenter(tb_max),tb_max)<<endl;
			//while(nevent_t < h2->GetSum()/tbin_fix&&tb_max<1000){
			//while(nevent_t <=max(single_bin_raw,h2->GetSum()/6) &&tb_max<1000){
			while(nevent_t <=h2->GetSum()/6 &&tb_max<1000){
				nevent_t += h2->GetBinContent(tb_max++);
			}
			double tmax = h2->GetBinCenter(tb_max);

			t1[new_bin] = tmax;
			cerr<<Form("    #%d bin: t=%f, tb=%d, N=%e ", new_bin, tmax, tb_max, (double) nevent_t)<<endl;
		}

		int tbin = new_bin;
		/*}}}*/

	//	if (tbin<1) xbin = 1; //Counts as one bin if the statistic is really low

		outf_total << i << "\t" << tbin << endl;
		cerr    << Form("--- Total t-bin in #%d phi bin is %d with N=%e .vs. %e" ,i,tbin,single_bin_raw, h2->GetSum()/100) << endl;
		Double_t tmin = 0.0, tmax = 0.0;
		for (Int_t j=0;j<tbin;j++){
			tmin = t1[j];
			tmax = t1[j+1];

			/*Histograms{{{*/
			//Q2
			TH1F *h1Q2  =new TH1F("h1Q2","h1Q2",1000,1.,7.);
			//x 
			TH1F *h1x = new TH1F("h1x","h1x",1000,0.1,0.7);
			// phi
			TH1F *h1phi = new TH1F("h1phi","h1phi",1000.,Phi_MIN,Phi_MAX);
			//t 
			TH1F *h1t = new TH1F("h1t","h1t",1000,0.1,0.7);
			/*}}}*/

			/*Fill and Save{{{*/
			h1Q2->Reset();
			h1x->Reset();
			h1phi->Reset();
			h1t->Reset();

			cut2.Form("(phi>=%f&&phi<%f&&t>=%f&&t<%f)",phimin,phimax,tmin,tmax);
			cut2 = cut1 + "*" + cut2;
			cut = cut2;

			T->Project("h1Q2","Q2",cut);
			T->Project("h1x","x",cut);
			T->Project("h1phi","phi",cut);
			T->Project("h1t","t",cut);

			double dilution = dilute_factor;//fix the value for now
			double N_out = h1t->GetSum()*(pow(polarization * target_factor * dilution,2) * det_eff);
			outf_total << i << " \t" << j << " \t" 
				<< h1Q2->GetMean() << " \t" 
				<< h1x->GetMean() << " \t" 
				<< h1phi->GetMean() << " \t"
				<< h1t->GetMean() << " \t" 
				<< h1t->GetSum() << " \t"
				<< N_out
				<< endl;

			h1Q2->Delete();
			h1x->Delete();
			h1phi->Delete();
			h1t->Delete();
			/*}}}*/
		}
		/*}}}*/

		h1->Delete();
		h2->Delete();
	}
	outf_total.close();
	return 0;
}

double max(double a, double b){
 if(a>b)
	 return a;
 else
	 return b;
}
