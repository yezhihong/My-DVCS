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
#include "/work/halla/solid/yez/dvcs/DVCS_XS_Grid/GetXS/DVCSGrid.h"
using namespace std;

const double Deg2Rad = TMath::DegToRad(); 
const double Sigma_EC_E = 0.02; //2%, energy resolution for electron from GEM tracking
const double Sigma_EC_G = 0.05; //10%, energy resolution for photon
const double Sigma_Theta_E = 0.6/1000.; //0.6mrad, Angular resolution for electron, determined by GEM tracking
const double Sigma_Phi_E =  5.0/1000.; //5mrad, Angular resolution for electron, determined by GEM tracking
const double Sigma_X_G = 1.0; //cm, Position resolution for photon, determined by EC cluster reconstruction 
const double Sigma_Y_G = 1.0; //cm, Position resolution for photon, determined by EC cluster reconstruction
const double Sigma_VZ = 0.5; //cm, VertexZ resolution for photon, determined by electron GEM tracking
const double Length_FAEC = 790.0; //cm, FAEC to target distance
const double Length_LAEC = 310.0; //cm, LAEC to target distance

const double eMass = 0.511/1000;//electron mass 
const double piMass = 0.13957018;
const double pMass = 0.938272;
const double nMass = 0.939565;

Int_t CheckLaws(TLorentzVector* P_t,TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_g, TLorentzVector* P_h);
Int_t CheckLawsPi0(TLorentzVector* P_t,TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_g1, TLorentzVector* P_g2, TLorentzVector* P_h);
Double_t GetMM( TLorentzVector* P_t,TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_g);
Double_t GetPi0MM(TLorentzVector* P_t,TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_g);

//int add_branch(TString particle,int energy_flag){
int main(){
	TString particle; cerr<<"--- Particle (n or p)"; cin>>particle;
	Int_t energy_flag; cerr<<"--- Energy Flag (11 or 8)"; cin>>energy_flag;

	gStyle->SetOptStat(1);
	//const double Lumi = 1.0e36; // cm-2*s-1, for He3 nuclear not for nucleons
	//const double KHz = 1e-3;
	//const double nBcm2 = 1e-33;
	double hMass = 0.0;
	double EBeam = 0.0;
	TString Target = "";
	TString TargetName = "";
	TString Energy = "";
	gRandom->SetSeed(0);
	/* Target and Energy variables{{{*/
	if(particle=="p"){
		TargetName = "Proton";
		Target = "NH3"; hMass = 0.938272;	
		if(energy_flag==11){
			Energy = "11GeV";
			EBeam = 11.0;
		}
		else if(energy_flag==8){
			Energy = "8.8GeV";
			EBeam = 8.80;
		}
	}else if(particle=="n"){
		TargetName = "Neutron";
		Target = "3he"; hMass = 0.939565;	
		if(energy_flag==11){
			Energy = "11GeV";
			EBeam = 11.0;
		}
		else if(energy_flag==8){
			Energy = "8.8GeV";
			EBeam = 8.80;
		}
	}/*}}}*/

	/*CLEO Acceptance{{{*/
	TFile *file_negative=new TFile(Form("../acceptance/acceptance_solid_CLEO_SIDIS_%s_negative_output.root",Target.Data()));
	TH2F *accep_e_forward=(TH2F*)file_negative->Get("acceptance_forwardangle");       //negative particle forward angle acceptance
	TH2F *accep_e_large=(TH2F*)file_negative->Get("acceptance_largeangle");           //negative particle large angle acceptance

	TFile* file_neutral=new TFile(Form("../acceptance/acceptance_solid_CLEO_SIDIS_%s_neutral_output.root", Target.Data()));
	TH2F *accep_pho_forward=(TH2F*)file_neutral->Get("acceptance_forwardangle");       //negative particle forward angle acceptance
	TH2F *accep_pho_large=(TH2F*)file_neutral->Get("acceptance_largeangle");       //negative particle forward angle acceptance
	/*}}}*/

	/*Vectors & Variables{{{*/
	TLorentzVector *P_E0 = new TLorentzVector();//incoming electron
	TLorentzVector *P_e_i = new TLorentzVector();//scattered electron
	TLorentzVector *P_e_f = new TLorentzVector();//scattered electron with Eloss
	TLorentzVector *P_g = new TLorentzVector();//photon
	TLorentzVector *P_h = new TLorentzVector();//proton or neutron
	TLorentzVector *P_t = new TLorentzVector();//target, either proton or neutron
	TLorentzVector *P_e_res = new TLorentzVector();//scattered electron with resolution
	TLorentzVector *P_g_res = new TLorentzVector();//photon with resolution

	Double_t vertexz, E0; 
	Double_t ePx_ini, ePy_ini, ePz_ini; 
	//Double_t hP_ini, hPx_ini, hPy_ini,hPz_ini,hTheta_ini,hPhi_ini;
	Double_t ePx, ePy, ePz, gPx, gPy, gPz, hPx, hPy, hPz;
	Double_t eP_ini,eP, gP,hP;
	Double_t Q2, x, t, phi, XS_P, XS_M, XS_BHp, XS_BHm, PSF;
	Double_t ePhi_ini, ePhi, gPhi,hPhi;
	Double_t eTheta_ini, eTheta, gTheta,hTheta;
	Int_t Ngen,Nacc;

	//Define new measured quantities that include detector resolutions
	double eE_i, eE_f, gE, hE;
	double eE_i_res, eE_f_res,eP_res, eTheta_res, ePhi_res, ePx_res, ePy_res, ePz_res;
	double gP_res, gTheta_res, gPhi_res, gPx_res, gPy_res, gPz_res;
	double Length=0.0,MM = 0.0,MM_res = 0.0, W=0.0, Wp=0.0;

	int e_theta_bin= 0;
	int e_p_bin=0;
	double e_acc_f= 0.0;
	double e_acc_l= 0.0;
	double g_acc_f=0, g_acc_l=0;
	int g_theta_bin=0, g_p_bin=0;
	/*}}}*/

	/*DVCS: Define Rootfile and variables{{{*/
	TString filename = "";
	if(particle=="p"){
		if(energy_flag==11)
			filename = "../DVCS_Proton_11GeV.root";
		else if(energy_flag==8)
			filename = "../DVCS_Proton_8.8GeV.root";
	}else if(particle=="n"){
		if(energy_flag==11)
			filename = "../DVCS_Neutron_11GeV.root";
		else if(energy_flag==8)
			filename = "../DVCS_Neutron_8.8GeV.root";
	}

	TFile *file = new TFile(filename.Data(), "r");
	TTree *T = (TTree*) file->Get("T");

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

	Long64_t N_entries=T->GetEntries();
	T->GetEntry(0);
	Long64_t N_gen = Ngen;
	cout<<"DVCS: total generated events number: "<<N_gen<<"--- and accepted: "<<N_entries<<endl;
	/*}}}*/

	/*Add new branches{{{*/
	filename.ReplaceAll(".root","_New0.root");
	TFile *newfile = new TFile(filename.Data(), "recreate");
	TTree *NewT = T->CloneTree(0);;

	NewT->Branch("W", &W, "W/D");
	NewT->Branch("Wp", &Wp, "Wp/D");
	NewT->Branch("MM", &MM, "MM/D");
	NewT->Branch("MM_res", &MM_res, "MM_res/D");
	NewT->Branch("e_acc_f", &e_acc_f, "e_acc_f/D");
	NewT->Branch("e_acc_l", &e_acc_l, "e_acc_l/D");
	NewT->Branch("g_acc_f", &g_acc_f, "g_acc_f/D");
	NewT->Branch("g_acc_l", &g_acc_l, "g_acc_l/D");

	NewT->Branch("eP_res", &eP_res, "eP_res/D");
	NewT->Branch("ePx_res", &ePx_res, "ePx_res/D");
	NewT->Branch("ePy_res", &ePy_res, "ePy_res/D");
	NewT->Branch("ePz_res", &ePz_res, "ePz_res/D");
	NewT->Branch("eTheta_res", &eTheta_res, "eTheta_res/D");
	NewT->Branch("ePhi_res", &ePhi_res, "ePhi_res/D");

	NewT->Branch("gP_res", &gP_res, "gP_res/D");
	NewT->Branch("gPx_res", &gPx_res, "gPx_res/D");
	NewT->Branch("gPy_res", &gPy_res, "gPy_res/D");
	NewT->Branch("gPz_res", &gPz_res, "gPz_res/D");
	NewT->Branch("gTheta_res", &gTheta_res, "gTheta_res/D");
	NewT->Branch("gPhi_res", &gPhi_res, "gPhi_res/D");

	Double_t tmin, BSA_L, TSA_L, DSA_L;
	Double_t Sigma_L, SigmaPP_L,SigmaPM_L, SigmaMP_L, SigmaMM_L;
	Double_t BSA_Tx, TSA_Tx, DSA_Tx;
	Double_t Sigma_Tx, SigmaPP_Tx,SigmaPM_Tx, SigmaMP_Tx, SigmaMM_Tx;
	Double_t BSA_Ty, TSA_Ty, DSA_Ty;
	Double_t Sigma_Ty, SigmaPP_Ty,SigmaPM_Ty, SigmaMP_Ty, SigmaMM_Ty;

	NewT->Branch("tmin",      &tmin,     "tmin/D");
	NewT->Branch("Sigma_L",   &Sigma_L,   "Sigma_L/D");
	NewT->Branch("SigmaPP_L", &SigmaPP_L, "SigmaPP_L/D");
	NewT->Branch("SigmaPM_L", &SigmaPM_L, "SigmaPM_L/D");
	NewT->Branch("SigmaMP_L", &SigmaMP_L, "SigmaMP_L/D");
	NewT->Branch("SigmaMM_L", &SigmaMM_L, "SigmaMM_L/D");
	NewT->Branch("BSA_L",     &BSA_L,     "BSA_L/D");
	NewT->Branch("TSA_L",     &TSA_L,     "TSA_L/D");
	NewT->Branch("DSA_L",     &DSA_L,     "DSA_L/D");

	NewT->Branch("Sigma_Tx",   &Sigma_Tx,   "Sigma_Tx/D");
	NewT->Branch("SigmaPP_Tx", &SigmaPP_Tx, "SigmaPP_Tx/D");
	NewT->Branch("SigmaPM_Tx", &SigmaPM_Tx, "SigmaPM_Tx/D");
	NewT->Branch("SigmaMP_Tx", &SigmaMP_Tx, "SigmaMP_Tx/D");
	NewT->Branch("SigmaMM_Tx", &SigmaMM_Tx, "SigmaMM_Tx/D");
	NewT->Branch("BSA_Tx",     &BSA_Tx,     "BSA_Tx/D");
	NewT->Branch("TSA_Tx",     &TSA_Tx,     "TSA_Tx/D");
	NewT->Branch("DSA_Tx",     &DSA_Tx,     "DSA_Tx/D");

	NewT->Branch("Sigma_Ty",   &Sigma_Ty,   "Sigma_Ty/D");
	NewT->Branch("SigmaPP_Ty", &SigmaPP_Ty, "SigmaPP_Ty/D");
	NewT->Branch("SigmaPM_Ty", &SigmaPM_Ty, "SigmaPM_Ty/D");
	NewT->Branch("SigmaMP_Ty", &SigmaMP_Ty, "SigmaMP_Ty/D");
	NewT->Branch("SigmaMM_Ty", &SigmaMM_Ty, "SigmaMM_Ty/D");
	NewT->Branch("BSA_Ty",     &BSA_Ty,     "BSA_Ty/D");
	NewT->Branch("TSA_Ty",     &TSA_Ty,     "TSA_Ty/D");
	NewT->Branch("DSA_Ty",     &DSA_Ty,     "DSA_Ty/D");
	/*}}}*/ 

	TString Data_Dir = "/work/halla/solid/yez/dvcs/DVCS_XS_Grid/";
	TString TargetPol = ""; //Tx, Ty, L 
	Int_t Debug = 0, err = 0;

	DVCSGrid *grid = new DVCSGrid();
	grid->Init(EBeam,TargetName.Data(), Data_Dir.Data(), Debug);

	/*Loop DVCS Events{{{*/  

	for(int i=0;i<N_entries; i++){
		T->GetEntry(i);

		grid->FindBin(Q2, x, -t, phi);// here -t was used instead of t
		//grid->PrintRanges();
		tmin = grid->GetTMin();	

		/*Reading the Longitudinal Target Spin{{{*/
		SigmaPP_L = 1e-36; SigmaPM_L = 1e-36;  SigmaMP_L = 1e-36;  SigmaMM_L = 1e-36; 
		Sigma_L = 1e-36;   BSA_L = 1e-36; TSA_L = 1e-36; DSA_L = 1e-36;
		TargetPol = "L";
		err = grid->LoadXS(TargetPol);
		if(err>=0){
			//grid->PrintXS();
			SigmaPP_L = grid->GetSigmaPP();//++
			SigmaPM_L = grid->GetSigmaPM();//+-
			SigmaMP_L = grid->GetSigmaMP();//-+
			SigmaMM_L = grid->GetSigmaMM();//--
			Sigma_L = grid->GetSigma();//average of four XSs above
			BSA_L = grid->GetBSA();//beam spin asym
			TSA_L = grid->GetTSA();//target spin asym
			DSA_L = grid->GetDSA();//double spin asym
		}
		/*Reading the Longitudinal Target Spin}}}*/	

		/*Reading the Transverse Target Spin on x{{{*/	
		SigmaPP_Tx = 1e-36;SigmaPM_Tx = 1e-36; SigmaMP_Tx = 1e-36; SigmaMM_Tx = 1e-36; 
		Sigma_Tx = 1e-36;  BSA_Tx = 1e-36;TSA_Tx = 1e-36;DSA_Tx = 1e-36;
		TargetPol = "Tx";
		err = grid->LoadXS(TargetPol);
		if(err>=0){
			//grid->PrintXS();
			SigmaPP_Tx = grid->GetSigmaPP();//++
			SigmaPM_Tx = grid->GetSigmaPM();//+-
			SigmaMP_Tx = grid->GetSigmaMP();//-+
			SigmaMM_Tx = grid->GetSigmaMM();//--
			Sigma_Tx = grid->GetSigma();//average of four XSs above
			BSA_Tx = grid->GetBSA();//beam spin asym
			TSA_Tx = grid->GetTSA();//target spin asym
			DSA_Tx = grid->GetDSA();//double spin asym
		//	cerr<<Form("-Tx-- BS = %e (%e)", BSA_Tx, (SigmaPP_Tx+SigmaPM_Tx - SigmaMP_Tx-SigmaMM_Tx)/4. )<<endl;
	//		cerr<<Form("-Tx-- TS = %e (%e)", TSA_Tx, (SigmaPP_Tx+SigmaMP_Tx - SigmaPM_Tx-SigmaMM_Tx)/4. )<<endl;
	//		cerr<<Form("-Tx-- DS = %e (%e)", DSA_Tx, (SigmaPP_Tx+SigmaMM_Tx - SigmaPM_Tx-SigmaMP_Tx)/4. )<<endl;
	//		cerr<<Form("-Tx-- AVG = %e ", Sigma_Tx)<<endl;
		}
		/*Reading the Transverse Target Spin on x}}}*/	

		/*Reading the Transverse Target Spin on y{{{*/	
		SigmaPP_Ty = 1e-36;SigmaPM_Ty = 1e-36; SigmaMP_Ty = 1e-36; SigmaMM_Ty = 1e-36; 
		Sigma_Ty = 1e-36;  BSA_Ty = 1e-36;TSA_Ty = 1e-36;DSA_Ty = 1e-36;

		TargetPol = "Ty";
		err = grid->LoadXS(TargetPol);
		if(err>=0){
			//grid->PrintXS();
			SigmaPP_Ty = grid->GetSigmaPP();//++
			SigmaPM_Ty = grid->GetSigmaPM();//+-
			SigmaMP_Ty = grid->GetSigmaMP();//-+
			SigmaMM_Ty = grid->GetSigmaMM();//--
			Sigma_Ty = grid->GetSigma();//average of four XSs above
			BSA_Ty = grid->GetBSA();//beam spin asym
			TSA_Ty = grid->GetTSA();//target spin asym
			DSA_Ty = grid->GetDSA();//double spin asym
		}
		/*Reading the Transverse Target Spin on y}}}*/	

		/*Acceptance{{{*/
		e_theta_bin= 0;	e_p_bin=0;
		e_acc_f= -1e6; e_acc_l= -1e6;	g_acc_f= -1e6; g_acc_l= -1e6;
		MM = -1e6; MM_res = -1e6; W = -1e6;; Wp = -1e6;
		eP_res = -1e6; ePx_res = -1e6; ePy_res = -1e6; ePz_res = -1e6; eTheta_res = -1e6; ePhi_res = -1e6;
		gP_res = -1e6; gPx_res = -1e6; gPy_res = -1e6; gPz_res = -1e6; gTheta_res = -1e6; gPhi_res = -1e6;

		e_theta_bin=int(eTheta/0.2)+1;    //0.2 degree per bin
		e_p_bin=int(eP/0.05)+1;      //0.05 GeV per bin for mom
		e_acc_f=accep_e_forward->GetBinContent(e_theta_bin,e_p_bin);
		e_acc_l=accep_e_large->GetBinContent(e_theta_bin,e_p_bin);

		if(eP<1.0||eTheta>14.5||eTheta<8.0)//GeV, CLEO
			e_acc_f=0.0;//Farward-Angle EC Cut at 1 GeV
		if(eP<1.0||eTheta<15.5||eTheta>24)//GeV,CLEO
			e_acc_l=0.0; //Larger-Angle EC Cut at 3 GeV
		if(e_acc_f>1.) 
			e_acc_f=1.0; 
		if(e_acc_l>1.) 
			e_acc_l=1.0; 

		g_theta_bin=int(gTheta/0.2)+1;  //0.2 degree per bin
		g_p_bin=int(gP/0.05)+1;     //0.05 GeV per bin for mom
		g_acc_f=accep_pho_forward->GetBinContent(g_theta_bin,g_p_bin);
		g_acc_l=accep_pho_large->GetBinContent(g_theta_bin,g_p_bin);
		/*Additional cuts{{{*/
		   if(gTheta>14.8||gTheta<8.0||gP<1.0||gP>11.){
               g_acc_f=0.0; 
		   }
		   else if(gTheta>24.0||gTheta<16.||gP<1.0||gP>11.){
               g_acc_l=0.0; 
           }
           if(g_acc_f>1.0) g_acc_f = 1.0;
           if(g_acc_l>1.0) g_acc_l = 1.0;
           /*}}}*/

		//	double e_acceptance = e_acc_f+e_acc_l;
		//	double g_acceptance = g_acc_f+g_acc_l;
		//	double eg_acceptance = e_acceptance*g_acceptance;
		//double event_weight=(XS_P+XS_M)*PSF/N_gen*Lumi*nBcm2;   //put into Hz, Note, XS in nb
		/*Acceptance}}}*/

		/*Missing Mass w/o resolutions{{{*/
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
		/*}}}*/

		/*Missing Mass w resolutions{{{*/
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
		if(g_acc_f>g_acc_l)
           Length = Length_FAEC;
		else
           Length = Length_LAEC;

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
		/////////////////////////////////////////
		/*}}}*/

		W = (*P_E0 + *P_t - *P_e_i)*(*P_E0 + *P_t - *P_e_i);
		Wp = (*P_E0 + *P_t - *P_e_i - *P_h)*(*P_E0 + *P_t - *P_e_i - *P_h);
		W = sqrt(W);
		Wp = sqrt(Wp);


		if(!(i%1000)) cerr<<Form("--- Processing #event = %d/%d", i,(int)(N_entries) ) <<"\r";

		NewT->Fill();
	}
	/*}}}*/

	/*Save and Free{{{*/
	NewT->Write("",TObject::kOverwrite);
	newfile->Close();
	file->Close();
	delete grid;
	delete P_g; delete P_t; delete P_e_i; delete P_e_f; delete P_h; delete P_e_res; delete P_g_res;
	/*Save and Free}}}*/
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
