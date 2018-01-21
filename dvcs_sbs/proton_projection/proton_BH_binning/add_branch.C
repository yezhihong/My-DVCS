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

//int add_branch(TString particle,int energy_flag){
int main(){
	TString particle; cerr<<"--- Particle (n or p)"; cin>>particle;
	Int_t energy_flag; cerr<<"--- Energy Flag (11 or 8)"; cin>>energy_flag;

    gStyle->SetOptStat(1);
	//const double Lumi = 1.0e36; // cm-2*s-1, for He3 nuclear not for nucleons
	//const double KHz = 1e-3;
	//const double nBcm2 = 1e-33;
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
	TString histoname;
	if(Target=="3he")
		histoname = "acceptance";
	if(Target=="NH3")
		histoname = "acceptance_ThetaP";
	TFile *file_negative=new TFile(Form("../acceptance/acceptance_solid_CLEO_SIDIS_%s_negative_output.root",Target.Data()));
	TH2F *accep_e_forward=(TH2F*)file_negative->Get(Form("%s_forwardangle",histoname.Data()));       //negative particle forward angle acceptance
	TH2F *accep_e_large=(TH2F*)file_negative->Get(Form("%s_largeangle",histoname.Data()));           //negative particle large angle acceptance

	TFile* file_neutral=new TFile(Form("../acceptance/acceptance_solid_CLEO_SIDIS_%s_neutral_output.root", Target.Data()));
	TH2F *accep_pho_forward=(TH2F*)file_neutral->Get(Form("%s_forwardangle",histoname.Data()));       //negative particle forward angle acceptance
	TH2F *accep_pho_large=(TH2F*)file_neutral->Get(Form("%s_largeangle",histoname.Data()));       //negative particle forward angle acceptance
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
	Double_t Q2, x, t, phi, XS_P, XS_M,XS_BHp, XS_BHm,PSF;
	Double_t ePhi_ini, ePhi, gPhi,hPhi;
	Double_t eTheta_ini, eTheta, gTheta,hTheta;
	Int_t Ngen,Nacc;

	//Define new measured quantities that include detector resolutions
	double eE_i, eE_f, gE, hE;
	double eE_i_res, eE_f_res,eP_res, eTheta_res, ePhi_res, ePx_res, ePy_res, ePz_res;
	double gP_res, gTheta_res, gPhi_res, gPx_res, gPy_res, gPz_res;
	double MM = 0.0,MM_res = 0.0, W=0.0, Wp=0.0;
    
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
	
	TFile *file = new TFile(filename.Data(), "update");
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
    TBranch *br_W = T->Branch("W", &W, "W/D");
    TBranch *br_Wp = T->Branch("Wp", &Wp, "Wp/D");
	TBranch *br_MM = T->Branch("MM", &MM, "MM/D");
	TBranch *br_MM_res = T->Branch("MM_res", &MM_res, "MM_res/D");
    TBranch *br_e_acc_f = T->Branch("e_acc_f", &e_acc_f, "e_acc_f/D");
    TBranch *br_e_acc_l = T->Branch("e_acc_l", &e_acc_l, "e_acc_l/D");
    TBranch *br_g_acc_f = T->Branch("g_acc_f", &g_acc_f, "g_acc_f/D");
    TBranch *br_g_acc_l = T->Branch("g_acc_l", &g_acc_l, "g_acc_l/D");
    
	TBranch *br_eP_res = T->Branch("eP_res", &eP_res, "eP_res/D");
    TBranch *br_ePx_res = T->Branch("ePx_res", &ePx_res, "ePx_res/D");
    TBranch *br_ePy_res = T->Branch("ePy_res", &ePy_res, "ePy_res/D");
    TBranch *br_ePz_res = T->Branch("ePz_res", &ePz_res, "ePz_res/D");
    TBranch *br_eTheta_res = T->Branch("eTheta_res", &eTheta_res, "eTheta_res/D");
    TBranch *br_ePhi_res = T->Branch("ePhi_res", &ePhi_res, "ePhi_res/D");
    
	TBranch *br_gP_res = T->Branch("gP_res", &gP_res, "gP_res/D");
    TBranch *br_gPx_res = T->Branch("gPx_res", &gPx_res, "gPx_res/D");
    TBranch *br_gPy_res = T->Branch("gPy_res", &gPy_res, "gPy_res/D");
    TBranch *br_gPz_res = T->Branch("gPz_res", &gPz_res, "gPz_res/D");
    TBranch *br_gTheta_res = T->Branch("gTheta_res", &gTheta_res, "gTheta_res/D");
    TBranch *br_gPhi_res = T->Branch("gPhi_res", &gPhi_res, "gPhi_res/D");

	/*}}}*/

	/*Get DVCS Missing Mass{{{*/  

	for(int i=0;i<N_entries; i++){
		T->GetEntry(i);
		e_theta_bin= 0;	e_p_bin=0;
		e_acc_f= -1e6; e_acc_l= -1e6;	g_acc_f= -1e6; g_acc_l= -1e6;
		MM = -1e6; MM_res = -1e6; W = -1e6;; Wp = -1e6;
        eP_res = -1e6; ePx_res = -1e6; ePy_res = -1e6; ePz_res = -1e6; eTheta_res = -1e6; ePhi_res = -1e6;
        gP_res = -1e6; gPx_res = -1e6; gPy_res = -1e6; gPz_res = -1e6; gTheta_res = -1e6; gPhi_res = -1e6;

		if(phi>=0.0&&phi<360){//any additional cuts should be added in here

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
			/*
			   if(gTheta<=14.8&&gTheta>=8.0&&gP>=1.0&&gP<=11.){
			   g_acc_f=1.0; 
			   g_acc_l=0.0; 
			   }
			   else if(gTheta<=24.0&&gTheta>=16.&&gP>=1.0&&gP<=11.){
			   g_acc_f=0.0; 
			   g_acc_l=1.0; 
			   }else{
			   g_acc_f=0.0; 
			   g_acc_l=0.0; 
			   }
			   */
			/*}}}*/

		//	double e_acceptance = e_acc_f+e_acc_l;
		//	double g_acceptance = g_acc_f+g_acc_l;
		//	double eg_acceptance = e_acceptance*g_acceptance;
			//double event_weight=(XS_P+XS_M)*PSF/N_gen*Lumi*nBcm2;   //put into Hz, Note, XS in nb

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
		}
		/*Fill Branches{{{*/
	
		br_e_acc_f->Fill();
		br_e_acc_l->Fill();
		br_g_acc_f->Fill();
		br_g_acc_l->Fill();

		br_W->Fill();
		br_Wp->Fill();
		br_MM->Fill();
		br_MM_res->Fill();
	
		br_eP_res->Fill();
		br_ePx_res->Fill();
		br_ePy_res->Fill();
		br_ePz_res->Fill();
		br_eTheta_res->Fill();
		br_ePhi_res->Fill();

		br_gP_res->Fill();
		br_gPx_res->Fill();
		br_gPy_res->Fill();
		br_gPz_res->Fill();
		br_gTheta_res->Fill();
		br_gPhi_res->Fill();
		/*}}}*/
	}
	/*}}}*/

	T->Write("",TObject::kOverwrite);
	file->Close();
	delete P_g; delete P_t; delete P_e_i; delete P_e_f; delete P_h; delete P_e_res; delete P_g_res;
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
