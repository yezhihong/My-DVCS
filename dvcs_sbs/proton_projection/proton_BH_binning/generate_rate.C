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
  #include <TRandom3.h>
  //#include <TMatrix.h>
  /*}}}*/

#include "/work/halla/solid/yez/dvcs/DVCS_XS_Grid/GetXS/DVCSGrid.h"
using namespace std;

void generate_rato(TString particle,int energy_flag){
	gStyle->SetOptStat(0);

	const double Lumi = 1.0e36; // cm-2*s-1, for He3 nuclear not for nucleons
	const double KHz = 1e-3;
	const double nBcm2 = 1e-33;

	TString Target = "";
	Double_t EBeam = 0;
	/*Define Rootfile and variables{{{*/
	TChain *T=new TChain("T");
	if(particle=="p"){
		 Target = "NH3";	
		if(energy_flag==11){
			EBeam =11.0;
		 T->Add("../DVCS_Proton_11GeV.root");
		}
		else if(energy_flag==8){
			EBeam =8.8;
			T->Add("../DVCS_Proton_8.8GeV.root");
		}
	}else if(particle=="n"){
		 Target = "3he";	
		if(energy_flag==11){
			EBeam =11.0;
			T->Add("../DVCS_Neutron_11GeV.root"); 
		}
		else if(energy_flag==8){
			EBeam =8.8;
			T->Add("../DVCS_Neutron_8.8GeV.root");
		}
	}

	Double_t vertexz, E0, W, Wp; 
	Double_t ePx_ini, ePy_ini, ePz_ini, hPx_ini, hPy_ini,hPz_ini;
	Double_t ePx, ePy, ePz, gPx, gPy, gPz, hPx, hPy, hPz;
	Double_t eP_ini,hP_ini, eP, gP,hP;
	Double_t Q2, x, t, phi, XS_P, XS_M, XS_BHp, XS_BHm, PSF;
	Double_t ePhi_ini,hPhi_ini, ePhi, gPhi,hPhi;
	Double_t eTheta_ini,hTheta_ini, eTheta, gTheta,hTheta;
	Int_t Ngen,Nacc;

	T->SetBranchAddress("vertexz", &vertexz);
	T->SetBranchAddress("E0", &E0);
	T->SetBranchAddress("W", &W);
	T->SetBranchAddress("Wp", &Wp);

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
	/*}}}*/

	Long64_t N_entries=T->GetEntries();
	T->GetEntry(0);
	Long64_t N_gen = Ngen;
	cout<<"total generated events number: "<<N_gen<<"--- and accepted: "<<N_entries<<endl;

    TString Data_Dir = "/work/halla/solid/yez/dvcs/DVCS_XS_Grid/";
	TString TargetPol = ""; //Tx, Ty, L 
	Int_t Debug = 0, err = 0;
	DVCSGrid *grid = new DVCSGrid();

	//deal with acceptance files 
	//read in acceptance histograms
	//electron acceptance	

	/*CLEO Acceptance{{{*/
	TString histoname;
	if(Target=="3he")
		histoname = "acceptance";
	if(Target=="NH3")
		histoname = "acceptance_ThetaP";
	TFile *file_negative=new TFile(Form("../acceptance/acceptance_solid_CLEO_SIDIS_%s_negative_output.root",Target.Data()));
	TH2F *accep_ele_forward=(TH2F*)file_negative->Get(Form("%s_forwardangle",histoname.Data()));       //negative particle forward angle acceptance
	TH2F *accep_ele_large=(TH2F*)file_negative->Get(Form("%s_largeangle",histoname.Data()));           //negative particle large angle acceptance

	TFile* file_neutral=new TFile(Form("../acceptance/acceptance_solid_CLEO_SIDIS_%s_neutral_output.root", Target.Data()));
	TH2F *accep_pho_forward=(TH2F*)file_neutral->Get(Form("%s_forwardangle",histoname.Data()));       //negative particle forward angle acceptance
	TH2F *accep_pho_large=(TH2F*)file_neutral->Get(Form("%s_largeangle",histoname.Data()));       //negative particle forward angle acceptance
  	/*}}}*/

	/*Define Histograms{{{*/
	int Nbinx=100;
	int Nbiny=100;
	int Nbin_theta=100;
	int Nbin_p=100;

	/*Kinematics Variables{{{*/
	//_________x vs t 
	//forward angle
	TH2F *hf_x_t=new TH2F("hf_x_t","x VS t ",Nbinx,0,0.8,Nbiny,-2.10,0.1);
	hf_x_t->GetXaxis()->SetTitle("x_{bj}");
	hf_x_t->GetYaxis()->SetTitle("t (GeV)");
	//large angle
	TH2F *hl_x_t=new TH2F("hl_x_t","x VS t ",Nbinx,0,0.8,Nbiny,-2.10,0.1);
	hl_x_t->GetXaxis()->SetTitle("x_{bj}");
	hl_x_t->GetYaxis()->SetTitle("t (GeV)");
	//total
	TH2F *h_x_t=new TH2F("h_x_t","x VS t ",Nbinx,0,0.8,Nbiny,-2.10,0.1);
	h_x_t->GetXaxis()->SetTitle("x_{bj}");
	h_x_t->GetYaxis()->SetTitle("t (GeV)");

	//_______Q2 vs t  
	//forward angle
	TH2F *hf_Q2_t=new TH2F("hf_Q2_t","Q^{2} VS t ",Nbinx,-2.1,0.1,Nbiny,0,10);
	hf_Q2_t->GetYaxis()->SetTitle("Q^2(GeV^{2})");
	hf_Q2_t->GetXaxis()->SetTitle("t (GeV)");
	//large angle
	TH2F *hl_Q2_t=new TH2F("hl_Q2_t","Q^{2} VS t ",Nbinx,-2.1,0.1,Nbiny,0,10);
	hl_Q2_t->GetYaxis()->SetTitle("Q^2(GeV^{2})");
	hl_Q2_t->GetXaxis()->SetTitle("t (GeV)");
	//total
	TH2F *h_Q2_t=new TH2F("h_Q2_t","Q^{2} VS t ",Nbinx,-2.1,0.1,Nbiny,0,10);
	h_Q2_t->GetYaxis()->SetTitle("Q^2(GeV^{2})");
	h_Q2_t->GetXaxis()->SetTitle("t (GeV)");

	//______Q2 vs x 
	//forward angle
	TH2F *hf_Q2_x=new TH2F("hf_Q2_x","Q^{2} VS x ",Nbinx,0,0.7,Nbiny,0,10);
	hf_Q2_x->GetXaxis()->SetTitle("x_{bj}");
	hf_Q2_x->GetYaxis()->SetTitle("Q^2(GeV^{2})");
	//large angle
	TH2F *hl_Q2_x=new TH2F("hl_Q2_x","Q^{2} VS x ",Nbinx,0,0.7,Nbiny,0,10);
	hl_Q2_x->GetXaxis()->SetTitle("x_{bj}");
	hl_Q2_x->GetYaxis()->SetTitle("Q^2(GeV^{2})");
	//total
	TH2F *h_Q2_x=new TH2F("h_Q2_x","Q^{2} VS x ",Nbinx,0,0.7,Nbiny,0,10);
	h_Q2_x->GetXaxis()->SetTitle("x_{bj}");
	h_Q2_x->GetYaxis()->SetTitle("Q^2(GeV^{2})");

	//______Q2 vs W 
	//forward angle
	TH2F *hf_Q2_W=new TH2F("hf_Q2_W","Q^{2} VS W ",Nbinx,1,4.5,Nbiny,0,10);
	hf_Q2_W->GetXaxis()->SetTitle("W");
	hf_Q2_W->GetYaxis()->SetTitle("Q^2(GeV^{2})");
	//large angle
	TH2F *hl_Q2_W=new TH2F("hl_Q2_W","Q^{2} VS W ",Nbinx,1,4.5,Nbiny,0,10);
	hl_Q2_W->GetXaxis()->SetTitle("W");
	hl_Q2_W->GetYaxis()->SetTitle("Q^2(GeV^{2})");
	//total
	TH2F *h_Q2_W=new TH2F("h_Q2_W","Q^{2} VS W ",Nbinx,1,4.5,Nbiny,0,10);
	h_Q2_W->GetXaxis()->SetTitle("W");
	h_Q2_W->GetYaxis()->SetTitle("Q^2(GeV^{2})");


	//______t vs pe
	//forward angle
	TH2F *hf_t_p=new TH2F("hf_t_p","t VS P (electron) ",Nbinx,0,10.0,Nbiny,-2.1,0.1);
	hf_t_p->GetXaxis()->SetTitle("P (GeV)");
	hf_t_p->GetYaxis()->SetTitle("t(GeV)");
	//large angle
	TH2F *hl_t_p=new TH2F("hl_t_p","t VS P (electron) ",Nbinx,0,10.0,Nbiny,-2.1,0.1);
	hl_t_p->GetXaxis()->SetTitle("P (GeV)");
	hl_t_p->GetYaxis()->SetTitle("t(GeV)");
	//total
	TH2F *h_t_p=new TH2F("h_t_p","t VS P (electron) ",Nbinx,0,10.0,Nbiny,-2.1,0.1);
	h_t_p->GetXaxis()->SetTitle("P (GeV)");
	h_t_p->GetYaxis()->SetTitle("t(GeV)");

	//______x vs pe
	//forward angle
	TH2F *hf_x_p=new TH2F("hf_x_p","x VS P (electron) ",Nbinx,0,10.0,Nbiny,0,0.8);
	hf_x_p->GetYaxis()->SetTitle("x_{bj}");
	hf_x_p->GetXaxis()->SetTitle("P (GeV)");
	//large angle
	TH2F *hl_x_p=new TH2F("hl_x_p","x VS P (electron) ",Nbinx,0,10.0,Nbiny,0,0.8);
	hl_x_p->GetYaxis()->SetTitle("x_{bj}");
	hl_x_p->GetXaxis()->SetTitle("P (GeV)");
	//total
	TH2F *h_x_p=new TH2F("h_x_p","x VS P (electron) ",Nbinx,0,10.0,Nbiny,0,0.8);
	h_x_p->GetYaxis()->SetTitle("x_{bj}");
	h_x_p->GetXaxis()->SetTitle("P (GeV)");

	//______Q2 vs pe 
	//forward angle
	TH2F *hf_Q2_p=new TH2F("hf_Q2_p","Q^{2} VS P (electron) ",Nbinx,0,10.0,Nbiny,0,10);
	hf_Q2_p->GetYaxis()->SetTitle("Q^2(GeV^{2})");
	hf_Q2_p->GetXaxis()->SetTitle("P (GeV)");
	//large angle
	TH2F *hl_Q2_p=new TH2F("hl_Q2_p","Q^{2} VS P (electron) ",Nbinx,0,10.0,Nbiny,0,10);
	hl_Q2_p->GetYaxis()->SetTitle("Q^2(GeV^{2})");
	hl_Q2_p->GetXaxis()->SetTitle("P (GeV)");
	//total
	TH2F *h_Q2_p=new TH2F("h_Q2_p","Q^{2} VS P (electron) ",Nbinx,0,10.0,Nbiny,0,10);
	h_Q2_p->GetYaxis()->SetTitle("Q^2(GeV^{2})");
	h_Q2_p->GetXaxis()->SetTitle("P (GeV)");
    /*}}}*/

	/*electron{{{*/
	/*_____p vs theta for e{{{*/
	TH2F *hf_p_theta=new TH2F("hf_p_theta","P_{e-} VS #theta_{e-}  (forward angle)",Nbin_theta,7.0,26.0,Nbin_p,0,10.0);
	//forward angle
	hf_p_theta->GetYaxis()->SetTitle("P_{e-}");
	hf_p_theta->GetXaxis()->SetTitle("#theta_{e-}");
	//large angle
	TH2F *hl_p_theta=new TH2F("hl_p_theta","P_{e-} VS #theta_{e-} (large angle)",Nbin_theta,7.0,26.0,Nbin_p,0,10.0);
	hl_p_theta->GetYaxis()->SetTitle("P_{e-}");
	hl_p_theta->GetXaxis()->SetTitle("#theta_{e-}");
	//total
	TH2F *h_p_theta=new TH2F("h_p_theta","P_{e-} VS #theta_{e-} (large+forward angle)",Nbin_theta,7.0,26.0,Nbin_p,0,10.0);
	h_p_theta->GetYaxis()->SetTitle("P_{e-}");
	h_p_theta->GetXaxis()->SetTitle("#theta_{e-}");
    /*}}}*/

	/*_____TH1F of theta for e{{{*/
	//forward angle
	TH1F *hf_theta=new TH1F("hf_theta","rate as function of #theta_{e-}  (forward angle)",Nbin_theta,7.0,26.0);
	hf_theta->GetYaxis()->SetTitle("rate (Hz)");
	hf_theta->GetXaxis()->SetTitle("#theta_{e-}");
	//large angle
	TH1F *hl_theta=new TH1F("hl_theta","rate as function of #theta_{e-} (large angle)",Nbin_theta,7.0,26.0);
	hl_theta->GetYaxis()->SetTitle("rate (Hz)");
	hl_theta->GetXaxis()->SetTitle("#theta_{e-}");
	//total
	TH1F *h_theta=new TH1F("h_theta","rate as function of #theta_{e-} (large+forward angle)",Nbin_theta,7.0,26.0);
	h_theta->GetYaxis()->SetTitle("rate (Hz)");
	h_theta->GetXaxis()->SetTitle("#theta_{e-}");
    /*}}}*/

	/*_____TH1F of p e{{{*/
	//forward angle
	TH1F *hf_p=new TH1F("hf_p","rate as function of P_{e-}  (forward angle)",Nbin_p,0.0,10.0);
	hf_p->GetYaxis()->SetTitle("rate (Hz)");
	hf_p->GetXaxis()->SetTitle("P_{e-}");
	//large angle
	TH1F *hl_p=new TH1F("hl_p","rate as function of P_{e-} (large angle)",Nbin_p,0.0,10.0);
	hl_p->GetYaxis()->SetTitle("rate (Hz)");
	hl_p->GetXaxis()->SetTitle("p_{e-}");
	//total
	TH1F *h_p=new TH1F("h_p","rate as function of p_{e-} (large+forward angle)",Nbin_p,0.0,10.0);
	h_p->GetYaxis()->SetTitle("rate (Hz)");
	h_p->GetXaxis()->SetTitle("P_{e-}");
    /*}}}*/
    
	/*}}}*/

	/*neucleon{{{*/
	/*_____p vs theta for n{{{*/
	TH2F *pf_p_theta=new TH2F("pf_p_theta",Form("P_{%s} VS #theta_{%s}  (forward angle)",particle.Data(),particle.Data()),Nbin_theta,7.0,70.0,Nbin_p,0,2.0);
	//forward angle
	pf_p_theta->GetYaxis()->SetTitle("P_{N}");
	pf_p_theta->GetXaxis()->SetTitle("#theta_{N}");
	//large angle
	TH2F *pl_p_theta=new TH2F("pl_p_theta",Form("P_{%s} VS #theta_{%s}  (large angle)",particle.Data(),particle.Data()),Nbin_theta,7.0,70.0,Nbin_p,0,2.0);
	pl_p_theta->GetYaxis()->SetTitle("P_{N}");
	pl_p_theta->GetXaxis()->SetTitle("#theta_{N}");
    //total
	TH2F *p_p_theta=new TH2F("p_p_theta",Form("P_{%s} VS #theta_{%s}  (large+forward angle)",particle.Data(),particle.Data()),Nbin_theta,7.0,70.0,Nbin_p,0,2.0);
	p_p_theta->GetYaxis()->SetTitle("P");
	p_p_theta->GetXaxis()->SetTitle("#theta_{N}");
    /*}}}*/
    
	/*_____TH1F of theta for n {{{*/
	//forward angle
	TH1F *pf_theta=new TH1F("pf_theta",Form("rate as function of #theta_{%s}  (forward angle)",particle.Data()),Nbin_theta,0.0,70.0);
	pf_theta->GetYaxis()->SetTitle("rate (Hz)");
	pf_theta->GetXaxis()->SetTitle("#theta_{N}");
	//large angle
	TH1F *pl_theta=new TH1F("pl_theta",Form("rate as function of #theta_{%s}  (large angle)",particle.Data()),Nbin_theta,0.0,70.0);
	pl_theta->GetYaxis()->SetTitle("rate (Hz)");
	pl_theta->GetXaxis()->SetTitle("#theta_{N}");
	//total
	TH1F *p_theta=new TH1F("p_theta",Form("rate as function of #theta_{%s}  (large+forward angle)",particle.Data()),Nbin_theta,0.0,70.0);
	p_theta->GetYaxis()->SetTitle("rate (Hz)");
	p_theta->GetXaxis()->SetTitle("#theta_{N}");
    /*}}}*/

	/*_____TH1F of p for n{{{*/
	//forward angle
	TH1F *pf_p=new TH1F("pf_p",Form("rate as function of P_{%s}  (forward angle)",particle.Data()),Nbin_p,0.0,2.0);
	pf_p->GetYaxis()->SetTitle("rate (Hz)");
	pf_p->GetXaxis()->SetTitle("P_{N}");
	//large angle
	TH1F *pl_p=new TH1F("pl_p",Form("rate as function of P_{%s}  (forward angle)",particle.Data()),Nbin_p,0.0,2.0);
	pl_p->GetYaxis()->SetTitle("rate (Hz)");
	pl_p->GetXaxis()->SetTitle("P_{N}");
	//total
	TH1F *p_p=new TH1F("p_p",Form("rate as function of P_{%s}  (large+forward angle)",particle.Data()),Nbin_p,0.0,2.0);
	p_p->GetYaxis()->SetTitle("rate (Hz)");
	p_p->GetXaxis()->SetTitle("P_{N}");
	/*}}}*/
    
	/*}}}*/

	/*photon{{{*/
	/*_____p vs theta for g{{{*/
    //forward
	TH2F *gf_p_theta=new TH2F("gf_p_theta","P_{#gamma} VS #theta_#{#gamma} (forward angle)",Nbin_theta,7.0,26.0,Nbin_p,0,10.0);
	gf_p_theta->GetYaxis()->SetTitle("P_{#gamma}");
	gf_p_theta->GetXaxis()->SetTitle("#theta_{#gamma}");
    //large  
	TH2F *gl_p_theta=new TH2F("gl_p_theta","P_{#gamma} VS #theta_#{#gamma} (large angle)",Nbin_theta,7.0,26.0,Nbin_p,0,10.0);
	gl_p_theta->GetYaxis()->SetTitle("P_{#gamma}");
	gl_p_theta->GetXaxis()->SetTitle("#theta_{#gamma}");
    //total
	TH2F *g_p_theta=new TH2F("g_p_theta","P_{#gamma} VS #theta_#{#gamma} (large+forward angle)",Nbin_theta,7.0,26.0,Nbin_p,0,10.0);
	g_p_theta->GetYaxis()->SetTitle("P_{#gamma}");
	g_p_theta->GetXaxis()->SetTitle("#theta_{#gamma}");
    /*}}}*/

	/*_____TH1F of theta g{{{*/
	//total
	TH1F *gf_theta=new TH1F("gf_theta","rate as function of #theta_#{#gamma} (forward angle)",Nbin_theta,7.0,26.0);
	gf_theta->GetYaxis()->SetTitle("rate (Hz)");
	gf_theta->GetXaxis()->SetTitle("#theta_{#gamma}");
	//total
	TH1F *gl_theta=new TH1F("gl_theta","rate as function of #theta_#{#gamma} (large angle)",Nbin_theta,7.0,26.0);
	gl_theta->GetYaxis()->SetTitle("rate (Hz)");
	gl_theta->GetXaxis()->SetTitle("#theta_{#gamma}");
	//total
	TH1F *g_theta=new TH1F("g_theta","rate as function of #theta_#{#gamma} (forward angle+large angle)",Nbin_theta,7.0,26.0);
	g_theta->GetYaxis()->SetTitle("rate (Hz)");
	g_theta->GetXaxis()->SetTitle("#theta_{#gamma}");
    /*}}}*/

	/*_____TH1F of p g{{{*/
	//forard
	TH1F *gf_p=new TH1F("gf_p","rate as function of P_{#gamma} (forward angle)",Nbin_p,0.0,10.0);
	gf_p->GetYaxis()->SetTitle("rate (Hz)");
	gf_p->GetXaxis()->SetTitle("P_{#gamma} (GeV)");
	//forard
	TH1F *gl_p=new TH1F("gl_p","rate as function of P_{#gamma} (large angle)",Nbin_p,0.0,10.0);
	gl_p->GetYaxis()->SetTitle("rate (Hz)");
	gl_p->GetXaxis()->SetTitle("P_{#gamma} (GeV)");
	//total
	TH1F *g_p=new TH1F("g_p","rate as function of P_{#gamma} (forward angle+large angle)",Nbin_p,0.0,10.0);
	g_p->GetYaxis()->SetTitle("rate (Hz)");
	g_p->GetXaxis()->SetTitle("P_{#gamma} (GeV)");
    /*}}}*/
    
	/*}}}*/

	/*_____TH1F of phi e, where phi is the angle between the electron plane and photon plane{{{*/
	//forward angle
	TH1F *hf_phi=new TH1F("hf_phi","rate as function of #Phi (forward angle)",Nbin_theta,0.0,360.0);
	hf_phi->GetYaxis()->SetTitle("rate (Hz)");
	hf_phi->GetXaxis()->SetTitle("#Phi (degree)");
	TH1F *hl_phi=new TH1F("hl_phi","rate as function of #Phi (large angle)",Nbin_theta,0.0,360.0);
	hl_phi->GetYaxis()->SetTitle("rate (Hz)");
	hl_phi->GetXaxis()->SetTitle("#Phi (degree)");
	TH1F *h_phi=new TH1F("h_phi","rate as function of #Phi (large+forward angle)",Nbin_theta,0.0,360.0);
	h_phi->GetYaxis()->SetTitle("rate (Hz)");
	h_phi->GetXaxis()->SetTitle("#Phi (degree)");
	/*}}}*/
	
	/*}}}*/
	
	double ele_rate_forward=0;
	double ele_rate_large=0;
	double pho_rate_forward=0;
	double pho_rate_large=0;

	double rate_forward=0;
	double rate_large=0;

	int ele_theta_bin= 0;
	int ele_p_bin=0;
	double ele_forward_acceptance= 0.0;
	double ele_large_acceptance= 0.0;
			
	for(Long64_t i=0;i<N_entries;i++){
		T->GetEntry(i);
		ele_theta_bin= 0;
		ele_p_bin=0;
		ele_forward_acceptance= 0.0;
		ele_large_acceptance= 0.0;

		//if(1){//any additional cuts should be added in here
		//if(phi>2.0&&phi<358&&W>2.0){//any additional cuts should be added in here
		if(phi>=0.&&phi<360&&W>2.0){//any additional cuts should be added in here

			ele_theta_bin=int(eTheta/0.2)+1;    //0.2 degree per bin
			ele_p_bin=int(eP/0.05)+1;      //0.05 GeV per bin for mom
			ele_forward_acceptance=accep_ele_forward->GetBinContent(ele_theta_bin,ele_p_bin);
			ele_large_acceptance=accep_ele_large->GetBinContent(ele_theta_bin,ele_p_bin);

			if(eP<1.0||eTheta>14.5||eTheta<8.0)//GeV, CLEO
				ele_forward_acceptance=0.0;//Farward-Angle EC Cut at 1 GeV
			if(eP<1.0||eTheta<15.5||eTheta>24)//GeV,CLEO
				ele_large_acceptance=0.0; //Larger-Angle EC Cut at 3 GeV
			if(ele_forward_acceptance>1.) 
				ele_forward_acceptance=1.0; 
			if(ele_large_acceptance>1.) 
				ele_large_acceptance=1.0; 

			double pho_forward_acceptance=0, pho_large_acceptance=0;
			int photon_theta_bin=int(gTheta/0.2)+1;  //0.2 degree per bin
			int photon_p_bin=int(gP/0.05)+1;     //0.05 GeV per bin for mom
			pho_forward_acceptance=accep_pho_forward->GetBinContent(photon_theta_bin,photon_p_bin);
			pho_large_acceptance=accep_pho_large->GetBinContent(photon_theta_bin,photon_p_bin);
			if(gTheta>14.8||gTheta<8.0||gP<0.5||gP>11.)//GeV, CLEO
				pho_forward_acceptance=0.0;
			if(gTheta<16.0||gTheta>24.0||gP<0.5||gP>11.)//GeV, CLEO
				pho_large_acceptance=0.0; 
			if(pho_forward_acceptance>1.) 
				pho_forward_acceptance=1.0; 
			if(pho_large_acceptance>1.) 
				pho_large_acceptance=1.0; 

//			if(gTheta<=14.8&&gTheta>=8.0&&gP>=1.0&&gP<=11.){
//				pho_forward_acceptance=1.0; 
//				pho_large_acceptance=0.0; 
//			}
//			else if(gTheta<=24.0&&gTheta>=16.&&gP>=1.0&&gP<=11.){
//				pho_forward_acceptance=0.0; 
//				pho_large_acceptance=1.0; 
//			}else{
//				pho_forward_acceptance=0.0; 
//				pho_large_acceptance=0.0; 
//			}
//				pho_large_acceptance=0.0; 
				
			double event_weight=(XS_BHp+XS_BHm)*PSF/N_gen*nBcm2*Lumi ;   //in Hz

			double total_ele_forward_acceptance=ele_forward_acceptance*(pho_large_acceptance+pho_forward_acceptance);
			double total_ele_large_acceptance=ele_large_acceptance*(pho_large_acceptance+pho_forward_acceptance);
			double total_pho_forward_acceptance=(ele_forward_acceptance+ele_large_acceptance)*pho_forward_acceptance;
			double total_pho_large_acceptance=(ele_forward_acceptance+ele_large_acceptance)*pho_large_acceptance;
			double total_acceptance=(ele_forward_acceptance+ele_large_acceptance)*(pho_forward_acceptance+pho_large_acceptance);

			/*Rates{{{*/	
			if( eP>=1.0){
				ele_rate_forward+=event_weight*ele_forward_acceptance;
			}
			if( eP>=3.5){
				ele_rate_large+=event_weight*ele_large_acceptance;
			}
			if( gP>=0.0){
				pho_rate_forward+=event_weight*pho_forward_acceptance;
			}
			if( gP>=0.0){
				pho_rate_large+=event_weight*pho_large_acceptance;
			}
		
			if( eP>=1.0 && gP>=0.0){
				rate_forward+=event_weight*total_ele_forward_acceptance;
			}
			if( eP>=3.5 && gP>=0.0){
				rate_large+=event_weight*total_ele_large_acceptance;
			}
			/*}}}*/

			/*fill histogram here{{{*/
            /*Kinematics{{{*/
			//x vs pt
			hf_x_t->Fill(x,t,event_weight*total_ele_forward_acceptance);
			hl_x_t->Fill(x,t,event_weight*total_ele_large_acceptance);
			h_x_t->Fill(x,t,event_weight*(total_acceptance));
			//Q2 vs t
			hf_Q2_t->Fill(t,Q2,event_weight*total_ele_forward_acceptance);
			hl_Q2_t->Fill(t,Q2,event_weight*total_ele_large_acceptance);
			h_Q2_t->Fill(t,Q2,event_weight*(total_acceptance));
			//Q2 vs x
			hf_Q2_x->Fill(x,Q2,event_weight*total_ele_forward_acceptance);
			hl_Q2_x->Fill(x,Q2,event_weight*total_ele_large_acceptance);
			h_Q2_x->Fill(x,Q2,event_weight*(total_acceptance));
			//Q2 vs W
			hf_Q2_W->Fill(W,Q2,event_weight*total_ele_forward_acceptance);
			hl_Q2_W->Fill(W,Q2,event_weight*total_ele_large_acceptance);
			h_Q2_W->Fill(W,Q2,event_weight*(total_acceptance));
			//t vs p
			hf_t_p->Fill(eP,t,event_weight*total_ele_forward_acceptance);
			hl_t_p->Fill(eP,t,event_weight*total_ele_large_acceptance);
			h_t_p->Fill(eP,t,event_weight*(total_acceptance));
			//x vs p
			hf_x_p->Fill(eP,x,event_weight*total_ele_forward_acceptance);
			hl_x_p->Fill(eP,x,event_weight*total_ele_large_acceptance);
			h_x_p->Fill(eP,x,event_weight*(total_acceptance));
			//Q2 vs p
			hf_Q2_p->Fill(eP,Q2,event_weight*total_ele_forward_acceptance);
			hl_Q2_p->Fill(eP,Q2,event_weight*total_ele_large_acceptance);
			h_Q2_p->Fill(eP,Q2,event_weight*(total_acceptance));			//x vs z
            /*}}}*/

            /*Electrons{{{*/
			//p vs theta e
			hf_p_theta->Fill(eTheta,eP,event_weight*total_ele_forward_acceptance);
			hl_p_theta->Fill(eTheta,eP,event_weight*total_ele_large_acceptance);
			h_p_theta->Fill( eTheta,eP, event_weight*(total_acceptance));
			//TH1F theta
			hf_theta->Fill(eTheta,event_weight*total_ele_forward_acceptance);
			hl_theta->Fill(eTheta,event_weight*total_ele_large_acceptance);
			h_theta->Fill( eTheta,event_weight*(total_acceptance));
			//TH1F p
			hf_p->Fill(eP,event_weight*total_ele_forward_acceptance);
			hl_p->Fill(eP,event_weight*total_ele_large_acceptance);
			h_p->Fill( eP,event_weight*(total_acceptance));
            /*}}}*/

            /*Photons{{{*/
			//p vs theta g
			gf_p_theta->Fill( gTheta,gP, event_weight*(total_pho_forward_acceptance));
			gl_p_theta->Fill( gTheta,gP, event_weight*(total_pho_large_acceptance));
			g_p_theta->Fill( gTheta,gP, event_weight*(total_acceptance));
			//TH1F theta
			gf_theta->Fill( gTheta,event_weight*(total_pho_forward_acceptance));
			gl_theta->Fill( gTheta,event_weight*(total_pho_large_acceptance));
			g_theta->Fill( gTheta,event_weight*(total_acceptance));
			//TH1F p
			gf_p->Fill( gP,event_weight*(total_pho_forward_acceptance));
			gl_p->Fill( gP,event_weight*(total_pho_large_acceptance));
			g_p->Fill( gP,event_weight*(total_acceptance));
            /*}}}*/

            /*Nucleons{{{*/
			//p vs theta p/n
			pf_p_theta->Fill(hTheta,hP,event_weight);
			pl_p_theta->Fill(hTheta,hP,event_weight);
			p_p_theta->Fill( hTheta,hP, event_weight);
			//TH1F theta p/n
			pf_theta->Fill(hTheta,event_weight);
			pl_theta->Fill(hTheta,event_weight);
			p_theta->Fill( hTheta,event_weight);
			//TH1F p/n
			pf_p->Fill(hP,event_weight);
			pl_p->Fill(hP,event_weight);
			p_p->Fill( hP,event_weight);
            /*}}}*/
	
			hf_phi->Fill(phi,event_weight*total_ele_forward_acceptance);
			hl_phi->Fill(phi,event_weight*total_ele_large_acceptance);
			h_phi->Fill( phi,event_weight*(total_acceptance));
			/*}}}*/
		}

	}// events loop ends here

	/*Print&Save{{{*/
	ofstream outf(Form("%s_%d_he3_rate_8.dat", particle.Data(), energy_flag));

	outf<<"______Single integral rate____________________"<<endl;
	outf<<"eletron rate_forward: "<<ele_rate_forward*KHz<<endl;
	outf<<"eletron rate_large: "<<ele_rate_large*KHz<<endl;
	outf<<"photon rate_forward: "<<pho_rate_forward*KHz<<endl;
	outf<<"photon rate_large: "<<pho_rate_large*KHz<<endl;

	outf<<"______forward and large angle integral rate____________________"<<endl;
	outf<<"rate_forward: "<<rate_forward*KHz<<endl;
	outf<<"rate_large: "<<rate_large*KHz<<endl;

	cout<<"______Single integral rate____________________"<<endl;
	cout<<"eletron rate_forward: "<<ele_rate_forward*KHz<<endl;
	cout<<"eletron rate_large: "<<ele_rate_large*KHz<<endl;
	cout<<"photon rate_forward: "<<pho_rate_forward*KHz<<endl;
	cout<<"photon rate_large: "<<pho_rate_large*KHz<<endl;

	cout<<"______forward and large angle integral rate____________________"<<endl;
	cout<<"rate_forward: "<<rate_forward*KHz<<endl;
	cout<<"rate_large: "<<rate_large*KHz<<endl;
	/*}}}*/

	/*Plot-New{{{*/

	/*Kinematics{{{*/
	
	TCanvas *can_forward=new TCanvas("can_forward","forward angle", 800,600);
	can_forward->Divide(2,2);
	can_forward->cd(1);
	hf_x_t->Draw("colz");
	can_forward->cd(2);
	hf_Q2_x->Draw("colz");
	can_forward->cd(3);
	hf_Q2_t->Draw("colz");
	can_forward->cd(4);
	hf_Q2_W->Draw("colz");
//	can_forward->cd(2);
//	hf_x_p->Draw("colz");
//	can_forward->cd(4);
//	hf_t_p->Draw("colz");
//	can_forward->cd(6);
//	hf_Q2_p->Draw("colz");

	can_forward->Print(Form("CLEO_DVCS_%s_%d_Forward.png", particle.Data(),energy_flag));
	can_forward->Print(Form("CLEO_DVCS_%s_%d_Forward.pdf", particle.Data(),energy_flag));

	TCanvas *can_large=new TCanvas("can_large","large angle", 800,600);
	can_large->Divide(2,2);
	can_large->cd(1);
	hl_x_t->Draw("colz");
	can_large->cd(2);
	hl_Q2_x->Draw("colz");
	can_large->cd(3);
	hl_Q2_t->Draw("colz");
	can_large->cd(4);
	hl_Q2_W->Draw("colz");
//	can_large->cd(2);
//	hl_x_p->Draw("colz");
//	can_large->cd(4);
//	hl_t_p->Draw("colz");
//	can_large->cd(6);
//	hl_Q2_p->Draw("colz");
	can_large->Print(Form("CLEO_DVCS_%s_%d_Large.png", particle.Data(),energy_flag));
	can_large->Print(Form("CLEO_DVCS_%s_%d_Large.pdf", particle.Data(),energy_flag));

	TCanvas *can_total=new TCanvas("can_total","total", 800,600);
	can_total->Divide(2,2);
	can_total->cd(1);
	hf_x_t->Draw();
	hl_x_t->SetMarkerColor(8); hl_x_t->Draw("same");
	can_total->cd(2);
	hf_Q2_x->Draw();
	hl_Q2_x->SetMarkerColor(8); hl_Q2_x->Draw("same");
	can_total->cd(3);
	hf_Q2_t->Draw();
	hl_Q2_t->SetMarkerColor(8); hl_Q2_t->Draw("same");
	can_total->cd(4);
	hf_Q2_W->Draw();
	hl_Q2_W->SetMarkerColor(8); hl_Q2_W->Draw("same");

//	can_total->cd(2);
//	hf_x_p->Draw();
//	hl_x_p->SetMarkerColor(8); hl_x_p->Draw("same");
//	can_total->cd(4);
//	hf_t_p->Draw();
//	hl_t_p->SetMarkerColor(8); hl_t_p->Draw("same");
//	can_total->cd(6);
//	hf_Q2_p->Draw();
//	hl_Q2_p->SetMarkerColor(8); hl_Q2_p->Draw("same");

	can_total->Print(Form("CLEO_DVCS_%s_%d_Total.png", particle.Data(), energy_flag));
	can_total->Print(Form("CLEO_DVCS_%s_%d_Total.pdf", particle.Data(), energy_flag));
   
	TCanvas *can_total2=new TCanvas("can_total2","total2", 800,600);
	can_total2->Divide(2,2);
	can_total2->cd(1);
	h_x_t->Draw("colz");
	can_total2->cd(2);
	h_Q2_x->Draw("colz");
	can_total2->cd(3);
	h_Q2_t->Draw("colz");
	can_total2->cd(4);
	h_Q2_W->Draw("colz");

//	can_total2->cd(2);
//	h_x_p->Draw("colz");
//	can_total2->cd(4);
//	h_t_p->Draw("colz");
//	can_total2->cd(6);
//	h_Q2_p->Draw("colz");

	can_total2->Print(Form("CLEO_DVCS_%s_%d_Total2.png", particle.Data(), energy_flag));
	can_total2->Print(Form("CLEO_DVCS_%s_%d_Total2.pdf", particle.Data(), energy_flag));
    
	/*}}}*/

	TCanvas *c1=new TCanvas("c1","#Phi plots", 800,800);
	c1->Divide(1,3);
	c1->cd(1);
    hf_phi->Draw();
	c1->cd(2);
    hl_phi->Draw();
	c1->cd(3);
    h_phi->Draw();
	c1->Print(Form("CLEO_DVCS_%s_%d_Phi.png", particle.Data(), energy_flag));
	c1->Print(Form("CLEO_DVCS_%s_%d_Phi.pdf", particle.Data(), energy_flag));

	TCanvas *c2=new TCanvas("c2","theta and P_{e-} plots", 800,800);
	c2->Divide(1,3);
	c2->cd(1);
	h_p_theta->Draw("COLZ");
	c2->cd(2);
	hf_theta->SetLineColor(2); hf_theta->Draw();
	hl_theta->SetLineColor(4); hl_theta->Draw("same");
	c2->cd(3);
	hf_p->SetLineColor(2); hf_p->Draw();
	hl_p->SetLineColor(4); hl_p->Draw("same");

	c2->Print(Form("CLEO_DVCS_%s_%d_Electron.png", particle.Data(), energy_flag));
	c2->Print(Form("CLEO_DVCS_%s_%d_Electron.pdf", particle.Data(), energy_flag));

	TCanvas *c3=new TCanvas("c3","Hadron plots (No acceptance cut)", 800,800);
	c3->Divide(1,3);
	c3->cd(1);
    p_p_theta->Draw("colz");
	c3->cd(2);
    pf_theta->SetLineColor(2); pf_theta->Draw("");
    pl_theta->SetLineColor(4); pl_theta->Draw("same");
	c3->cd(3);
	pf_p->SetLineColor(2); pl_p->Draw("");
	pl_p->SetLineColor(4); pl_p->Draw("same");

	c3->Print(Form("CLEO_DVCS_%s_%d_Hadron.png", particle.Data(), energy_flag));
	c3->Print(Form("CLEO_DVCS_%s_%d_Hadron.pdf", particle.Data(), energy_flag));
	
	TCanvas *c4=new TCanvas("c4","Photon plots (forward angle only)", 800,800);
	c4->Divide(1,3);
	c4->cd(1);
    g_p_theta->Draw("colz");
	c4->cd(2);
    gf_theta->SetLineColor(2); gf_theta->Draw("");
    gl_theta->SetLineColor(4); gl_theta->Draw("same");
	c4->cd(3);
	gf_p->SetLineColor(2); gf_p->Draw("");
	gl_p->SetLineColor(4); gl_p->Draw("same");

	c4->Print(Form("CLEO_DVCS_%s_%d_Photon.png", particle.Data(), energy_flag));
	c4->Print(Form("CLEO_DVCS_%s_%d_Photon.pdf", particle.Data(), energy_flag));
	/*}}}*/

	/*Save Histo{{{*/
	TString outputfilename="test.root";
	if(particle=="p"){
		if(energy_flag==11){
			outputfilename = "proton_11_histo.root";
		}
		if(energy_flag==8){
			outputfilename = "proton_8_histo.root";
		}
	}

	if(particle=="n"){
		if(energy_flag==11){
			outputfilename = "neutron_11_histo.root";
		}
		if(energy_flag==8){
			outputfilename = "neutron_8_histo.root";
		}
	}

	TFile *fout = new TFile(outputfilename.Data(),"recreate");
	fout->cd();
	h_x_t->Write();
	h_Q2_t->Write();
	h_Q2_x->Write();
	h_Q2_W->Write();
	h_t_p->Write();
	h_x_p->Write();
	h_Q2_p->Write();

	hl_x_t->Write();
	hl_Q2_t->Write();
	hl_Q2_x->Write();
	hl_Q2_W->Write();
	hl_t_p->Write();
	hl_x_p->Write();
	hl_Q2_p->Write();
	hf_x_t->Write();
	hf_Q2_t->Write();
	hf_Q2_x->Write();
	hf_Q2_W->Write();
	hf_t_p->Write();
	hf_x_p->Write();
	hf_Q2_p->Write();

	hf_p_theta->Write();
	hf_theta->Write();
	hf_p->Write();
	hl_p_theta->Write();
	hl_theta->Write();
	hl_p->Write();
	h_p_theta->Write();
	h_theta->Write();
	h_p->Write();

	gf_p_theta->Write();
	gl_p_theta->Write();
	g_p_theta->Write();
	gf_theta->Write();
	gl_theta->Write();
	g_theta->Write();
	gf_p->Write();
	gl_p->Write();
	g_p->Write();
    
	pf_p->Write();
	pf_theta->Write();
	pf_p_theta->Write();
	pl_p->Write();
	pl_theta->Write();
	pl_p_theta->Write();
	p_p->Write();
	p_theta->Write();
	p_p_theta->Write();
	
	hf_phi->Write();
	hl_phi->Write();
	h_phi->Write();

	fout->Close();
	/*}}}*/
//	file_negative->Close();
//	file_neutral->Close();
	
}

