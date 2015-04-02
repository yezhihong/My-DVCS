/*C/C++ Includes{{{*/
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>
#include <cmath>
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
#include <TClass.h>
#include <TPaletteAxis.h>
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

//Here bin on Q2 and x first, and the t and phi binning will be in the next step
Int_t main(Int_t argc, char *argv[]){
	TString particle = "X"; cerr<<"-- What particle (p->proton, n->neutron)? "; cin >> particle;
    TString Target = "";
	if(particle=="n")
		Target= "3he";
	if(particle=="n")
		Target= "NH3";


	/*DVCS: Define Rootfile and variables{{{*/
	TChain *T = new TChain("T");
	if(particle=="p"){
			T->Add("../Proton_8.8GeV.root");
			T->Add("../Proton_11GeV.root");
	}else if(particle=="n"){
			T->Add("../Neutron_11GeV.root");
			T->Add("../Neutron_8.8GeV.root");
	}

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

	Long64_t N_entries=T->GetEntries();
	T->GetEntry(0);
	Long64_t N_gen = Ngen*2;//11GeV and 8.8 GeV settings should have the same total generated events but please check
	cout<<"DVCS: total generated events number: "<<N_gen<<"--- and accepted: "<<N_entries<<endl;
	/*}}}*/
	
	TString prefix = "./skim_rootfiles_11p8/";
	const int Day_11 = 48;//for 11GeV SIDIS run
	const int Day_8 = 21;//for 8.8GeV SIDIS run
	const double Lumi = 1.0e36; // cm-2*s-1, for He3 nuclear not for nucleons
	const double nBcm2 = 1e-33;
//	const Double_t t_MIN=0.1,  t_MAX=0.70;//end at 0.6
//	const Double_t x_MIN=0.1,  x_MAX=0.70;//end at 0.6
//	const Double_t Q2_MIN=1.0, Q2_MAX=7.0;//end at 6.0GeV2
//	const Double_t phi_MIN=0.0,phi_MAX=360.0;//range from 2.0 ~ 358 degree

    const Int_t Q2_bin = 5;
	const Double_t Q2_cut[6] = {1.0, 1.5, 2.0, 3.0, 4.5, 7.0};
    const Int_t x_bin = 8;
	const Double_t x_cut[9] = {0.1, 0.16, 0.22, 0.28, 0.34, 0.40, 0.46, 0.52, 0.70};

	Double_t weight_p = 0.0, weight_m = 0.0, time=0.0, ele_energy = 0.0;
	for(int i=1;i<=Q2_bin;i++){
		for(int j=1;j<=x_bin;j++){
			Double_t Q2min = Q2_cut[i];
			Double_t Q2max = Q2_cut[i+1];
			Double_t xmin = x_cut[j];
			Double_t xmax = x_cut[j+1];

			cerr<<Form("--- Working on Q2_bin#=%d (min=%3.2f,max=%3.2f), x_bin#=%d (min=%3.2f,max=%3.2f) ...",
					i,Q2min,Q2max,j,xmin,xmax)<<endl;

			/*Define new rootfile{{{*/

			TString finalfile = Form("%s/%s_%d_%d_11p8.root", prefix.Data(), particle.Data(), i,j);

			TFile *file = new TFile(finalfile.Data(),"recreate");
			TTree *t1 = new TTree("T","A new tree");

			t1->Branch("vertexz", &vertexz, "vertexz/D");
			t1->Branch("E0", &E0, "E0/D");

			t1->Branch("eP_ini", &eP_ini, "eP_ini/D");
			t1->Branch("ePx_ini", &ePx_ini, "ePx_ini/D");
			t1->Branch("ePy_ini", &ePy_ini, "ePy_ini/D");
			t1->Branch("ePz_ini", &ePz_ini, "ePz_ini/D");
			t1->Branch("eTheta_ini",&eTheta_ini, "eTheta_ini/D");
			t1->Branch("ePhi_ini",&ePhi_ini, "ePhi_ini/D");

			t1->Branch("eP", &eP, "eP/D");
			t1->Branch("ePx", &ePx, "ePx/D");
			t1->Branch("ePy", &ePy, "ePy/D");
			t1->Branch("ePz", &ePz, "ePz/D");
			t1->Branch("eTheta",&eTheta, "eTheta/D");
			t1->Branch("ePhi",&ePhi, "ePhi/D");

			t1->Branch("gP",    &gP,    "gP/D");
			t1->Branch("gPx",   &gPx,   "gPx/D");
			t1->Branch("gPy",   &gPy,   "gPy/D");
			t1->Branch("gPz",   &gPz,   "gPz/D");
			t1->Branch("gTheta",&gTheta,"gTheta/D");
			t1->Branch("gPhi",  &gPhi,  "gPhi/D");

			t1->Branch("hP", &hP, "hP/D");
			t1->Branch("hPx", &hPx, "hPx/D");
			t1->Branch("hPy", &hPy, "hPy/D");
			t1->Branch("hPz", &hPz, "hPz/D");
			t1->Branch("hTheta",&hTheta, "hTheta/D");
			t1->Branch("hPhi",&hPhi, "hPhi/D");

			t1->Branch("Q2", &Q2, "Q2/D");
			t1->Branch("x", &x, "x/D");
			t1->Branch("t", &t, "t/D");
			t1->Branch("phi", &phi, "phi/D");
			t1->Branch("XS_P", &XS_P, "XS_P/D");
			t1->Branch("XS_M", &XS_M, "XS_M/D");
			t1->Branch("PSF", &PSF, "PSF/D");
			t1->Branch("Ngen", &Ngen, "Ngen/I");
			t1->Branch("Nacc", &Nacc, "Nacc/I");

			t1->Branch("MM", &MM, "MM/D");
			t1->Branch("MM_res", &MM_res, "MM_res/D");
			t1->Branch("W", &W, "W/D");
			t1->Branch("Wp", &Wp, "Wp/D");
			t1->Branch("eP_res", &eP_res, "eP_res/D");
			t1->Branch("ePx_res", &ePx_res, "ePx_res/D");
			t1->Branch("ePy_res", &ePy_res, "ePy_res/D");
			t1->Branch("ePz_res", &ePz_res, "ePz_res/D");
			t1->Branch("eTheta_res",&eTheta_res, "eTheta_res/D");
			t1->Branch("ePhi_res",&ePhi_res, "ePhi_res/D");

			t1->Branch("gP_res",    &gP_res,    "gP_res/D");
			t1->Branch("gPx_res",   &gPx_res,   "gPx_res/D");
			t1->Branch("gPy_res",   &gPy_res,   "gPy_res/D");
			t1->Branch("gPz_res",   &gPz_res,   "gPz_res/D");
			t1->Branch("gTheta_res",&gTheta_res,"gTheta_res/D");
			t1->Branch("gPhi_res",  &gPhi_res,  "gPhi_res/D");
			
			t1->Branch("e_acc_f", &e_acc_f, "e_acc_f/D");
			t1->Branch("e_acc_l", &e_acc_l, "e_acc_l/D");
			t1->Branch("g_acc_f", &g_acc_f, "g_acc_f/D");
			t1->Branch("g_acc_l", &g_acc_l, "g_acc_l/D");
			t1->Branch("weight_m", &weight_m, "weight_m/D");
			t1->Branch("weight_p", &weight_p, "weight_p/D");
			t1->Branch("time", &time, "time/D");
			/*}}}*/

			for (Int_t k=0;k<T->GetEntries();k++){
				T->GetEntry(k);

				if (x>=xmin&&x<xmax&&Q2>=Q2min&&Q2<Q2max&&W>2){
					/*Get Beam Time{{{*/
					ele_energy = E0;
					Double_t day = 0; // days
					// nevents = dxs (nbar) * L (nucleons/cm^2/s) * T(days) 
					// 1 day = 24hr*3600s = 86400
					// nbar=10-9 barn = 10^-9*10^-28 m^2=10^-9*10^-28 *10^4 cm^2=10^-33 cm^2
					// so weight  = Lumi * nbar * 86400 *day *acc_ele*acc_had / nsim;
					if(fabs(ele_energy-11.0)<0.1)
						day  = Day_11; // days
					if(fabs(ele_energy-8.80)<0.1)
						day  = Day_8; // days
					time = day * 24 *3600; //sec
					/*}}}*/

				    double total_acceptance = (e_acc_f+e_acc_l) * (g_acc_f+g_acc_l);	
					weight_p=XS_P*PSF/N_gen*nBcm2*Lumi*time*total_acceptance;   //in Count 
					weight_m=XS_M*PSF/N_gen*nBcm2*Lumi*time*total_acceptance;   //in Count 
	//				cerr<<Form("XS=%e,total_acceptance = %f, time=%e, weight_p=%f, weight_m=%f ", (XS_P+XS_M), total_acceptance, time,weight_p, weight_m)<<endl;
		
					t1->Fill();
					if(!(i%10000))
						cerr<<Form("--- Working on Q2_bin#=%d, x_bin#=%d, evt=%d",i,j,k)<<"\r";
				}
			}
			file->Write();
			file->Close();
		}
	}

	delete T;
	return 0;
}
