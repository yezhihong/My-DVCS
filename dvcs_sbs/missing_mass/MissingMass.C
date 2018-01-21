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
  #include <TLatex.h>
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
using namespace std;

const double Deg2Rad = TMath::DegToRad(); 
const double Sigma_EC_E = 0.02; //2%, energy resolution for electron from GEM tracking
const double Sigma_EC_G = 0.05; //10%, energy resolution for photon
const double Sigma_Theta_E = 0.6/1000.; //0.6mrad, Angular resolution for electron, determined by GEM tracking
const double Sigma_Phi_E =  5.0/1000.; //5mrad, Angular resolution for electron, determined by GEM tracking
const double Sigma_X_G = 1.0; //cm, Position resolution for photon, determined by EC cluster reconstruction 
const double Sigma_Y_G = 1.0; //cm, Position resolution for photon, determined by EC cluster reconstruction
const double Sigma_VZ = 0.5; //cm, VertexZ resolution for photon, determined by electron GEM tracking
const double Length = 790.0; //cm, FAEC to target distance

const double eMass = 0.511/1000;//electron mass 
const double piMass = 0.13957018;
const double pMass = 0.938272;
const double nMass = 0.939565;

Int_t CheckLaws(TLorentzVector* P_t,TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_g, TLorentzVector* P_h);
Int_t CheckLawsPi0(TLorentzVector* P_t,TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_g1, TLorentzVector* P_g2, TLorentzVector* P_h);
Double_t GetMM( TLorentzVector* P_t,TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_g);
Double_t GetPi0MM(TLorentzVector* P_t,TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_g);

int MissingMass(TString particle,int energy_flag){
	gStyle->SetOptStat(1);
	const double Lumi = 1.0e36; // cm-2*s-1, for He3 nuclear not for nucleons
	const double KHz = 1e-3;
	const double nBcm2 = 1e-33;
	double hMass = 0.0;
	TString Target = "";
	TString Energy = "";
	gRandom->SetSeed(0);
	/* Target and Energy variables{{{*/
	if(particle=="p"){
		 Target = "NH3"; hMass = 0.938272;	
		if(energy_flag==11){
		 Energy = "11GeV";
		}
		else if(energy_flag==8){
		 Energy = "8.8GeV";
		}
	}else if(particle=="n"){
		 Target = "3he"; hMass = 0.939565;	
		if(energy_flag==11){
		 Energy = "11GeV";
		}
		else if(energy_flag==8){
		 Energy = "8.8GeV";
		}
	}/*}}}*/
	
	/*CLEO Acceptance{{{*/
	TFile *file_negative=new TFile(Form("../acceptance/acceptance_solid_CLEO_SIDIS_%s_negative_output.root",Target.Data()));
	TH2F *accep_ele_forward=(TH2F*)file_negative->Get("acceptance_forwardangle");       //negative particle forward angle acceptance
	TH2F *accep_ele_large=(TH2F*)file_negative->Get("acceptance_largeangle");           //negative particle large angle acceptance

	TFile* file_neutral=new TFile(Form("../acceptance/acceptance_solid_CLEO_SIDIS_%s_neutral_output.root", Target.Data()));
	TH2F *accep_pho_forward=(TH2F*)file_neutral->Get("acceptance_forwardangle");       //negative particle forward angle acceptance
	TH2F *accep_pho_large=(TH2F*)file_neutral->Get("acceptance_largeangle");       //negative particle forward angle acceptance
  	/*}}}*/

	/*Vectors & Histograms{{{*/
	TLorentzVector *P_E0 = new TLorentzVector();//incoming electron
	TLorentzVector *P_e_i = new TLorentzVector();//scattered electron
	TLorentzVector *P_e_f = new TLorentzVector();//scattered electron with Eloss
	TLorentzVector *P_g = new TLorentzVector();//photon
	TLorentzVector *P_h = new TLorentzVector();//proton or neutron
	TLorentzVector *P_t = new TLorentzVector();//target, either proton or neutron
	TLorentzVector *P_e_res = new TLorentzVector();//scattered electron with resolution
	TLorentzVector *P_g_res = new TLorentzVector();//photon with resolution
	
	TH1F *hMM = new TH1F("hMM",Form("Missing Mass of %s at %s", particle.Data(),Energy.Data()), 200,0.0, 2.8);
	hMM->SetXTitle("Hadron Missing Mass (GeV)");
	hMM->GetXaxis()->CenterTitle(1);
	hMM->SetYTitle("Rate (Hz)");
	hMM->GetYaxis()->CenterTitle(1);
	TH1F *hMM_res = new TH1F("hMM_res",Form("Missing Mass of %s at %s", particle.Data(),Energy.Data()), 200,0.0, 2.8);
	hMM_res->SetXTitle("Hadron Missing Mass (GeV)");
	hMM_res->GetXaxis()->CenterTitle(1);
	hMM_res->SetYTitle("Rate (Hz)");
	hMM_res->GetYaxis()->CenterTitle(1);
  
	TH1F *piMM = new TH1F("piMM",Form("Missing Mass of n+#gamma at %s",Energy.Data()), 200,0.0, 2.8);
	piMM->SetXTitle("Hadron+#gamma Missing Mass (GeV)");
	piMM->GetXaxis()->CenterTitle(1);
	piMM->SetYTitle("Rate (Hz)");
	piMM->GetYaxis()->CenterTitle(1);
	TH1F *piMM_res = new TH1F("piMM_res",Form("Missing Mass of n+#gamma at %s",Energy.Data()), 200,0.0, 2.8);
	piMM_res->SetXTitle("Hadron+#gamma Missing Mass (GeV)");
	piMM_res->GetXaxis()->CenterTitle(1);
	piMM_res->SetYTitle("Rate (Hz)");
	piMM_res->GetYaxis()->CenterTitle(1);

	TH1F *piMM_g12 = new TH1F("piMM_g12",Form("Missing Mass of n at %s",Energy.Data()), 200,0.0, 2.8);
	piMM_g12->SetXTitle("Hadron Missing Mass (GeV)");
	piMM_g12->GetXaxis()->CenterTitle(1);
	piMM_g12->SetYTitle("Rate (Hz)");
	piMM_g12->GetYaxis()->CenterTitle(1);

	Double_t vertexz, E0; 
	Double_t ePx_ini, ePy_ini, ePz_ini, hPx_ini, hPy_ini,hPz_ini;
	Double_t ePx, ePy, ePz, gPx, gPy, gPz, hPx, hPy, hPz;
	Double_t eP_ini,hP_ini, eP, gP,hP;
	Double_t Q2, x, t, phi, Sigma_L, XS_P, XS_M,PSF;
	Double_t ePhi_ini,hPhi_ini, ePhi, gPhi,hPhi;
	Double_t eTheta_ini,hTheta_ini, eTheta, gTheta,hTheta;
	Int_t Ngen,Nacc;
    /*}}}*/

	//Define new measured quantities that include detector resolutions
	double eE_i, eE_f, gE, hE;
	double eE_i_res, eE_f_res,eP_res, eTheta_res, ePhi_res, ePx_res, ePy_res, ePz_res;
	double gP_res, gTheta_res, gPhi_res, gPx_res, gPy_res, gPz_res;
	double MM = 0.0,MM_res = 0.0;

	/*DVCS: Define Rootfile and variables{{{*/
	TChain *T=new TChain("T");
	if(particle=="p"){
		if(energy_flag==11)
			T->Add("../proton_projection/DVCS_Proton_11GeV_New0.root");
		else if(energy_flag==8)
			T->Add("../proton_projection/DVCS_Proton_8.8GeV_New0.root");
	}else if(particle=="n"){
		if(energy_flag==11)
			T->Add("../neutron_projection/DVCS_Neutron_11GeV_New0.root"); 
		else if(energy_flag==8)
			T->Add("../neutron_projection/DVCS_Neutron_8.8GeV_New0.root");
	}

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
	T->SetBranchAddress("Sigma_L", &Sigma_L);
	T->SetBranchAddress("XS_P", &XS_P);
	T->SetBranchAddress("XS_M", &XS_M);
	T->SetBranchAddress("PSF", &PSF);
	T->SetBranchAddress("Ngen", &Ngen);
	T->SetBranchAddress("Nacc", &Nacc);

	Long64_t N_entries=T->GetEntries();
	T->GetEntry(0);
	Long64_t N_gen = Ngen;
	cout<<"DVCS: total generated events number: "<<N_gen<<"--- and accepted: "<<N_entries<<endl;

	/*}}}*/
   
	/*Get DVCS Missing Mass{{{*/  
	int ele_theta_bin= 0;
	int ele_p_bin=0;
	double ele_forward_acceptance= 0.0;
	double ele_large_acceptance= 0.0;
	double g_forward_acceptance=0, g_large_acceptance=0;
	int g_theta_bin=0, g_p_bin=0;

	for(int i=0;i<N_entries; i++){
		T->GetEntry(i);
		ele_theta_bin= 0;
		ele_p_bin=0;
		ele_forward_acceptance= 0.0;
		ele_large_acceptance= 0.0;

		if(phi>2.0&&phi<358&&(Sigma_L)>0){//any additional cuts should be added in here
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

			g_theta_bin=int(gTheta/0.2)+1;  //0.2 degree per bin
			g_p_bin=int(gP/0.05)+1;     //0.05 GeV per bin for mom
			g_forward_acceptance=accep_pho_forward->GetBinContent(g_theta_bin,g_p_bin);
			g_large_acceptance=accep_pho_large->GetBinContent(g_theta_bin,g_p_bin);
			/*Additional cuts{{{*/
			/*
			   if(gTheta<=14.8&&gTheta>=8.0&&gP>=1.0&&gP<=11.){
			   g_forward_acceptance=1.0; 
			   g_large_acceptance=0.0; 
			   }
			   else if(gTheta<=24.0&&gTheta>=16.&&gP>=1.0&&gP<=11.){
			   g_forward_acceptance=0.0; 
			   g_large_acceptance=1.0; 
			   }else{
			   g_forward_acceptance=0.0; 
			   g_large_acceptance=0.0; 
			   }
			   */
			/*}}}*/

			double ele_acceptance = ele_forward_acceptance+ele_large_acceptance;
			double g_acceptance = g_forward_acceptance+g_large_acceptance;
			double eg_acceptance = ele_acceptance*g_acceptance;

			double event_weight=(Sigma_L)*PSF/N_gen*Lumi*nBcm2;   //put into Hz, Note, XS in nb
			//double event_weight=(XS_P+XS_M)*PSF/N_gen*Lumi*nBcm2;   //put into Hz, Note, XS in nb

			eE_i = sqrt(eP_ini*eP_ini + eMass*eMass);
			eE_f = sqrt(eP*eP + eMass*eMass);
			hE = sqrt(hP*hP + hMass*hMass);
			gE = gP;

			P_E0->SetPxPyPzE(0.,0.,E0, E0);
			P_t->SetPxPyPzE(0.,0.,0., hMass);
			P_e_i->SetPxPyPzE(ePx_ini, ePy_ini, ePz_ini, eE_i);
			P_e_f->SetPxPyPzE(ePx, ePy, ePz, eE_f);	
			P_g->SetPxPyPzE(gPx, gPy, gPz, gE);	
			P_h->SetPxPyPzE(hPx, hPy, hPz, hE);

			int err = CheckLaws(P_t, P_E0, P_e_i, P_g, P_h);//Check whether momentum and energy conserve first
			if (err < 1e-33){
				cerr<<"---- Momentum and Energy Conservation Laws are broken!! Something is wrong!!!"<<endl;
				return -222;
			}

			MM = GetMM(P_t, P_E0, P_e_i, P_g);	
			hMM->Fill(MM, event_weight*eg_acceptance);

			//////////////////////////////////////////////////////////////////////////
			//Now consider the detector resolution here
			//////////////////////////////////////////////////////////////////////////
			//
			/////////////////////////////////////////
			//Electron w/o E-loss
			eP_res = gRandom->Gaus(eP_ini, Sigma_EC_E*eP_ini);//GeV, for electron, E ~= P
			eTheta_res = gRandom->Gaus(eTheta_ini*Deg2Rad, Sigma_Theta_E);//rad
			ePhi_res = gRandom->Gaus(ePhi_ini*Deg2Rad, Sigma_Phi_E);//rad

			ePx_res = eP_res * sin(eTheta_res)*cos(ePhi_res); 
			ePy_res = eP_res * sin(eTheta_res)*sin(ePhi_res); 
			ePz_res = eP_res * cos(eTheta_res);

			eE_i_res = sqrt(eP_res*eP_res + eMass*eMass);
			P_e_res->SetPxPyPzE(ePx_res, ePy_res, ePz_res, eE_i_res);	

			/////////////////////////////////////////
			//Electron with E-loss
			/*
			   eP_res = gRandom->Gaus(eP, Sigma_EC_E*sqrt(eP));//GeV, for electron, E ~= P
			   eTheta_res = gRandom->Gaus(eTheta*Deg2Rad, Sigma_Theta_E);//rad
			   ePhi_res = gRandom->Gaus(ePhi*Deg2Rad, Sigma_Phi_E);//rad

			   ePx_res = eP_res * sin(eTheta_res)*cos(ePhi_res); 
			   ePy_res = eP_res * sin(eTheta_res)*sin(ePhi_res); 
			   ePz_res = eP_res * cos(eTheta_res);

			   eE_f_res = sqrt(eP_res*eP_res + eMass*eMass);
			//P_e_res->SetPxPyPzE(ePx_res, ePy_res, ePz_res, eE_i_res);	
			*/
			/////////////////////////////////////////
			//Photon
			double gZ = Length-vertexz; 
			double gX = gZ * tan(gTheta*Deg2Rad)*cos(gPhi*Deg2Rad); 
			double gY = gZ * tan(gTheta*Deg2Rad)*sin(gPhi*Deg2Rad); 
			double gX_res =gRandom->Gaus(gX, Sigma_X_G);//cm
			double gY_res =gRandom->Gaus(gY, Sigma_Y_G);//cm
			double gZ_res =gRandom->Gaus(gZ, Sigma_VZ);//cm 

			gP_res = gRandom->Gaus(gP, Sigma_EC_G*sqrt(gP));//GeV, for photon, E = P
			gPhi_res = atan2(gY_res, gX_res);
			gTheta_res = atan2(sqrt(gX_res*gX_res+gY_res*gY_res), gZ_res);

			gPx_res = gP_res * sin(gTheta_res)*cos(gPhi_res); 
			gPy_res = gP_res * sin(gTheta_res)*sin(gPhi_res); 
			gPz_res = gP_res * cos(gTheta_res);
			P_g_res->SetPxPyPzE(gPx_res, gPy_res, gPz_res, gP_res);	

			/////////////////////////////////////////
			MM_res = GetMM(P_t, P_E0, P_e_res, P_g_res);	
			hMM_res->Fill(MM_res, event_weight*eg_acceptance);	
			/////////////////////////////////////////
		}
	}
	delete T;
	/*}}}*/

	/*Pi0: Define Rootfile and variables{{{*/
	Double_t g1Px, g1Py, g1Pz, g1P, g1Theta, g1Phi,g1E;
	Double_t g2Px, g2Py, g2Pz, g2P, g2Theta, g2Phi,g2E;
	double g1P_res, g1Theta_res, g1Phi_res, g1Px_res, g1Py_res, g1Pz_res;
	double g2P_res, g2Theta_res, g2Phi_res, g2Px_res, g2Py_res, g2Pz_res;
	TLorentzVector *P_g1 = new TLorentzVector();//photon#1
	TLorentzVector *P_g2 = new TLorentzVector();//photon#1
	TLorentzVector *P_g1_res = new TLorentzVector();//photon with resolution
	TLorentzVector *P_g2_res = new TLorentzVector();//photon with resolution
	TChain *T1=new TChain("T");
	if(particle=="p"){
		if(energy_flag==11)
			T1->Add("../Pi0_Proton_11GeV.root");
		else if(energy_flag==8)
			T1->Add("../Pi0_Proton_8.8GeV.root");
	}else if(particle=="n"){
		if(energy_flag==11)
			T1->Add("../Pi0_Neutron_11GeV.root"); 
		else if(energy_flag==8)
			T1->Add("../Pi0_Neutron_8.8GeV.root");
	}

	T1->SetBranchAddress("vertexz", &vertexz);
	T1->SetBranchAddress("E0", &E0);

	T1->SetBranchAddress("eP_ini", &eP_ini);
	T1->SetBranchAddress("ePx_ini", &ePx_ini);
	T1->SetBranchAddress("ePy_ini", &ePy_ini);
	T1->SetBranchAddress("ePz_ini", &ePz_ini);
	T1->SetBranchAddress("eTheta_ini",&eTheta_ini);
	T1->SetBranchAddress("ePhi_ini",&ePhi_ini);

	T1->SetBranchAddress("eP", &eP);
	T1->SetBranchAddress("ePx", &ePx);
	T1->SetBranchAddress("ePy", &ePy);
	T1->SetBranchAddress("ePz", &ePz);
	T1->SetBranchAddress("eTheta",&eTheta);
	T1->SetBranchAddress("ePhi",&ePhi);

	T1->SetBranchAddress("g1P",   &g1P);
	T1->SetBranchAddress("g1Px",  &g1Px);
	T1->SetBranchAddress("g1Py",  &g1Py);
	T1->SetBranchAddress("g1Pz",  &g1Pz);
	T1->SetBranchAddress("g1Theta", &g1Theta);
	T1->SetBranchAddress("g1Phi", &g1Phi);

	T1->SetBranchAddress("g2P",   &g2P);
	T1->SetBranchAddress("g2Px",  &g2Px);
	T1->SetBranchAddress("g2Py",  &g2Py);
	T1->SetBranchAddress("g2Pz",  &g2Pz);
	T1->SetBranchAddress("g2Theta", &g2Theta);
	T1->SetBranchAddress("g2Phi", &g2Phi);

	T1->SetBranchAddress("hP", &hP);
	T1->SetBranchAddress("hPx", &hPx);
	T1->SetBranchAddress("hPy", &hPy);
	T1->SetBranchAddress("hPz", &hPz);
	T1->SetBranchAddress("hTheta",&hTheta);
	T1->SetBranchAddress("hPhi",&hPhi);

	T1->SetBranchAddress("Q2", &Q2);
	T1->SetBranchAddress("x", &x);
	T1->SetBranchAddress("t", &t);
	T1->SetBranchAddress("phi", &phi);
	T1->SetBranchAddress("XS_P", &XS_P);
	T1->SetBranchAddress("XS_M", &XS_M);
	T1->SetBranchAddress("PSF", &PSF);
	T1->SetBranchAddress("Ngen", &Ngen);
	T1->SetBranchAddress("Nacc", &Nacc);

	N_entries=T1->GetEntries();
	T1->GetEntry(0);
	N_gen = Ngen;
	cout<<"Pi0: total generated events number: "<<N_gen<<"--- and accepted: "<<N_entries<<endl;

	/*}}}*/

	/*Get Pi0 Decay Missing Mass{{{*/
	//Define new measured quantities that include detector resolutions
	double g1_forward_acceptance=0, g1_large_acceptance=0;
	int g1_theta_bin=0, g1_p_bin=0;
	double g2_forward_acceptance=0, g2_large_acceptance=0;
	int g2_theta_bin=0, g2_p_bin=0;

	for(int i=0;i<N_entries; i++){
		T1->GetEntry(i);
		ele_theta_bin= 0;ele_p_bin=0;
		ele_forward_acceptance= 0.0;ele_large_acceptance= 0.0;
		g1_forward_acceptance = 0; g2_forward_acceptance = 0; 
		g1_large_acceptance = 0; g2_large_acceptance = 0;
		g1_theta_bin=0; g2_theta_bin = 0; g1_p_bin =0; g2_p_bin=0;	

		if(phi/Deg2Rad>2.0&&phi/Deg2Rad<358){//any additional cuts should be added in here
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

			g1_theta_bin=int(g1Theta/0.2)+1;  //0.2 degree per bin
			g1_p_bin=int(g1P/0.05)+1;     //0.05 GeV per bin for mom
			g1_forward_acceptance=accep_pho_forward->GetBinContent(g1_theta_bin,g1_p_bin);
			g1_large_acceptance=accep_pho_large->GetBinContent(g1_theta_bin,g1_p_bin);

			g2_theta_bin=int(g2Theta/0.2)+1;  //0.2 degree per bin
			g2_p_bin=int(g2P/0.05)+1;     //0.05 GeV per bin for mom
			g2_forward_acceptance=accep_pho_forward->GetBinContent(g2_theta_bin,g2_p_bin);
			g2_large_acceptance=accep_pho_large->GetBinContent(g2_theta_bin,g2_p_bin);
			/*Additional cuts{{{*/
			/*
			   if(g1Theta<=14.8&&g1Theta>=8.0&&g1P>=1.0&&g1P<=11.){
			   g1_forward_acceptance=1.0; 
			   g1_large_acceptance=0.0; 
			   }
			   else if(g1Theta<=24.0&&g1Theta>=16.&&g1P>=1.0&&g1P<=11.){
			   g1_forward_acceptance=0.0; 
			   g1_large_acceptance=1.0; 
			   }else{
			   g1_forward_acceptance=0.0; 
			   g1_large_acceptance=0.0; 
			   }
			   if(g2Theta<=14.8&&g2Theta>=8.0&&g2P>=1.0&&g2P<=11.){
			   g2_forward_acceptance=1.0; 
			   g2_large_acceptance=0.0; 
			   }
			   else if(g2Theta<=24.0&&g2Theta>=16.&&g2P>=1.0&&g2P<=11.){
			   g2_forward_acceptance=0.0; 
			   g2_large_acceptance=1.0; 
			   }else{
			   g2_forward_acceptance=0.0; 
			   g2_large_acceptance=0.0; 
			   }
			   */
			/*}}}*/

			double ele_acceptance = ele_forward_acceptance+ele_large_acceptance;
			double g1_acceptance = g1_forward_acceptance+g1_large_acceptance;
			double g2_acceptance = g2_forward_acceptance+g2_large_acceptance;
			double eg1_acceptance = ele_acceptance*g1_acceptance;
			double eg12_acceptance = ele_acceptance*g1_acceptance*g2_acceptance;

			double Norm_Fact = 0.01;
			double event_weight=Norm_Fact*PSF/N_gen*Lumi*nBcm2;   //put into Hz, Note, XS in nb

			eE_i = sqrt(eP_ini*eP_ini + eMass*eMass);
			eE_f = sqrt(eP*eP + eMass*eMass);
			hE = sqrt(hP*hP + hMass*hMass);

			g1E = g1P;	g2E = g2P;

			P_E0->SetPxPyPzE(0.,0.,E0, E0);
			P_t->SetPxPyPzE(0.,0.,0., hMass);
			P_e_i->SetPxPyPzE(ePx_ini, ePy_ini, ePz_ini, eE_i);
			P_e_f->SetPxPyPzE(ePx, ePy, ePz, eE_f);	
			P_h->SetPxPyPzE(hPx, hPy, hPz, hE);
			P_g1->SetPxPyPzE(g1Px, g1Py, g1Pz, g1E);	
			P_g2->SetPxPyPzE(g2Px, g2Py, g2Pz, g2E);	

			int err = CheckLawsPi0(P_t, P_E0, P_e_i, P_g1, P_g2, P_h);//Check whether momentum and energy conserve first
			// int err = 1;
			if (err < 1e-33){
				cerr<<"---- Momentum and Energy Conservation Laws are broken!! Something is wrong!!!"<<endl;
				return -222;
			}

			MM = GetPi0MM(P_t, P_E0, P_e_i, P_g1);	
			piMM->Fill(MM, event_weight*eg1_acceptance);

			//////////////////////////////////////////////////////////////////////////
			//Now consider the detector resolution here
			//////////////////////////////////////////////////////////////////////////
			//
			/////////////////////////////////////////
			//Electron w/o E-loss
			eP_res = gRandom->Gaus(eP_ini, Sigma_EC_E*eP_ini);//GeV, for electron, E ~= P
			eTheta_res = gRandom->Gaus(eTheta_ini*Deg2Rad, Sigma_Theta_E);//rad
			ePhi_res = gRandom->Gaus(ePhi_ini*Deg2Rad, Sigma_Phi_E);//rad

			ePx_res = eP_res * sin(eTheta_res)*cos(ePhi_res); 
			ePy_res = eP_res * sin(eTheta_res)*sin(ePhi_res); 
			ePz_res = eP_res * cos(eTheta_res);

			eE_i_res = sqrt(eP_res*eP_res + eMass*eMass);
			P_e_res->SetPxPyPzE(ePx_res, ePy_res, ePz_res, eE_i_res);	

			/////////////////////////////////////////
			//Electron with E-loss
			/*
			   eP_res = gRandom->Gaus(eP, Sigma_EC_E*sqrt(eP));//GeV, for electron, E ~= P
			   eTheta_res = gRandom->Gaus(eTheta*Deg2Rad, Sigma_Theta_E);//rad
			   ePhi_res = gRandom->Gaus(ePhi*Deg2Rad, Sigma_Phi_E);//rad

			   ePx_res = eP_res * sin(eTheta_res)*cos(ePhi_res); 
			   ePy_res = eP_res * sin(eTheta_res)*sin(ePhi_res); 
			   ePz_res = eP_res * cos(eTheta_res);

			   eE_f_res = sqrt(eP_res*eP_res + eMass*eMass);
			//P_e_res->SetPxPyPzE(ePx_res, ePy_res, ePz_res, eE_i_res);	
			*/
			/////////////////////////////////////////
			//Photon
			double g1Z = Length-vertexz; 
			double g1X = g1Z * tan(g1Theta*Deg2Rad)*cos(g1Phi*Deg2Rad); 
			double g1Y = g1Z * tan(g1Theta*Deg2Rad)*sin(g1Phi*Deg2Rad); 
			double g1X_res =gRandom->Gaus(g1X, Sigma_X_G);//cm
			double g1Y_res =gRandom->Gaus(g1Y, Sigma_Y_G);//cm
			double g1Z_res =gRandom->Gaus(g1Z, Sigma_VZ);//cm 

			g1P_res = gRandom->Gaus(g1P, Sigma_EC_G*sqrt(g1P));//GeV, for photon, E = P
			g1Phi_res = atan2(g1Y_res, g1X_res);
			g1Theta_res = atan2(sqrt(g1X_res*g1X_res+g1Y_res*g1Y_res), g1Z_res);

			g1Px_res = g1P_res * sin(g1Theta_res)*cos(g1Phi_res); 
			g1Py_res = g1P_res * sin(g1Theta_res)*sin(g1Phi_res); 
			g1Pz_res = g1P_res * cos(g1Theta_res);
			P_g1_res->SetPxPyPzE(g1Px_res, g1Py_res, g1Pz_res, g1P_res);	

			/////////////////////////////////////////
			MM_res = GetPi0MM(P_t, P_E0, P_e_res, P_g1_res);	
			piMM_res->Fill(MM_res, event_weight*eg1_acceptance);	
			piMM_g12->Fill(MM_res, event_weight*(eg1_acceptance-eg12_acceptance));	
			/////////////////////////////////////////
		}
	}
	delete T1;
	/* }}}*/

	delete P_g, P_t, P_e_i, P_e_f, P_h, P_e_res, P_g_res;
	delete P_g1, P_g2, P_g1_res, P_g2_res;   

	/*Plot{{{*/
	TCanvas *c1 = new TCanvas("c1","c1", 1000,600);
	hMM_res->SetLineColor(4); hMM_res->Draw();// hMM_res->Fit("gaus");
	hMM->SetLineColor(8); hMM->Draw("same");
	piMM->SetLineColor(7); piMM->Draw("same");
	piMM_res->SetLineColor(2); piMM_res->Draw("same");;
	piMM_g12->SetLineColor(6); piMM_g12->Draw("same");;
	TLatex *t1=new TLatex();
	t1->SetNDC();
	t1->SetTextColor(6);
	t1->SetTextSize(0.040);
	t1->DrawLatex(0.18,0.90, Form("#deltaE_{e}=%2.1f%, #deltaE_{g}=%2.1f%/#sqrt{E}", Sigma_EC_E*100, Sigma_EC_G*100.));
	t1->DrawLatex(0.18,0.85, Form("#delta#theta_{e}=%2.1fmrad, #delta#phi_{e}=%2.1fmrad", Sigma_Theta_E*1000, Sigma_Phi_E*1000.));
	t1->DrawLatex(0.18,0.80, Form("#deltaX_{g}=%2.1fcm, #deltaY_{g}=%2.1fcm", Sigma_X_G, Sigma_Y_G));
	t1->DrawLatex(0.18,0.75, Form("#deltaVZ=%2.1fcm",Sigma_VZ));
	/*}}}*/
}

/*Double_t GetMM(TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_g){{{*/
Double_t GetMM(TLorentzVector* P_t,TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_g){
	Double_t kMM = 0.0;

	double E_miss = (P_t->E() + P_E0->E()) - (P_e->E()+P_g->E());
	double px_miss = (P_t->Px() + P_E0->Px()) - (P_e->Px()+P_g->Px()); 
	double py_miss = (P_t->Py() + P_E0->Py()) - (P_e->Py()+P_g->Py()); 
	double pz_miss = (P_t->Pz() + P_E0->Pz()) - (P_e->Pz()+P_g->Pz()); 

	kMM = sqrt( E_miss*E_miss - (pow(px_miss,2)+pow(py_miss,2)+pow(pz_miss,2)) );

	return kMM;
}
/*}}}*/

/*Double_t CheckLaws(TLorentzVector* P_t,TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_g, TLorentzVector* P_h){{{*/
Int_t CheckLaws(TLorentzVector* P_t,TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_g, TLorentzVector* P_h){
	Int_t err = -1;

	double energy_check = (P_t->E() + P_E0->E()) - (P_e->E()+P_g->E()+P_h->E());
	double px_check =(P_t->Px() + P_E0->Px()) - (P_e->Px()+P_g->Px()+P_h->Px()); 
	double py_check =(P_t->Py() + P_E0->Py()) - (P_e->Py()+P_g->Py()+P_h->Py()); 
	double pz_check =(P_t->Pz() + P_E0->Pz()) - (P_e->Pz()+P_g->Pz()+P_h->Pz()); 

	if(fabs(energy_check)<0.01 && fabs(px_check)<0.01 && fabs(py_check)<0.01 && fabs(pz_check)<0.01)
		err = 1;
	else{

		cerr<<Form("*** dE = %f,  dPx = %f, dPy = %f, dPz = %f", energy_check, px_check, py_check, pz_check)<<endl;

		err = -1;
	}
	return err;

}
/*}}}*/

/*Double_t GetPi0MM(TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_g){{{*/
Double_t GetPi0MM(TLorentzVector* P_t,TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_g){
	Double_t kMM = 0.0;

	double E_miss = (P_t->E() + P_E0->E()) - (P_e->E()+P_g->E());
	double px_miss = (P_t->Px() + P_E0->Px()) - (P_e->Px()+P_g->Px()); 
	double py_miss = (P_t->Py() + P_E0->Py()) - (P_e->Py()+P_g->Py()); 
	double pz_miss = (P_t->Pz() + P_E0->Pz()) - (P_e->Pz()+P_g->Pz()); 

	kMM = sqrt( E_miss*E_miss - (pow(px_miss,2)+pow(py_miss,2)+pow(pz_miss,2)) ) + (piMass+pMass-nMass);

	return kMM;
}
/*}}}*/

/*Double_t CheckLawsPi0(TLorentzVector* P_t,TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_g1, TLorentzVector* P_g2, TLorentzVector* P_h){{{*/
Int_t CheckLawsPi0(TLorentzVector* P_t,TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_g1, TLorentzVector* P_g2, TLorentzVector* P_h){
	Int_t err = -1;

	double energy_check = (P_t->E() + P_E0->E()) - (P_e->E()+P_g1->E()+P_g2->E()+P_h->E());
	double px_check =(P_t->Px() + P_E0->Px()) - (P_e->Px()+P_g1->Px()+P_g2->Px()+P_h->Px()); 
	double py_check =(P_t->Py() + P_E0->Py()) - (P_e->Py()+P_g1->Py()+P_g2->Py()+P_h->Py()); 
	double pz_check =(P_t->Pz() + P_E0->Pz()) - (P_e->Pz()+P_g1->Pz()+P_g2->Pz()+P_h->Pz()); 

	if(fabs(energy_check)<0.01 && fabs(px_check)<0.01 && fabs(py_check)<0.01 && fabs(pz_check)<0.01)
		err = 1;
	else{

		cerr<<Form("*** dE = %f,  dPx = %f, dPy = %f, dPz = %f", energy_check, px_check, py_check, pz_check)<<endl;

		err = -1;
	}
	return err;

}
/*}}}*/
