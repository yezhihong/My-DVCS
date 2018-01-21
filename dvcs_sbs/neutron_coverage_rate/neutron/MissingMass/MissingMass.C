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

const double DEG = TMath::DegToRad(); 
const double Sigma_EC_E = 0.00; //10%, energy resolution for electron
const double Sigma_EC_G = 0.00; //10%, energy resolution for photon
const double Sigma_Theta_E = 0.6/1000.; //0.6mrad, Angular resolution for electron, determined by GEM tracking
const double Sigma_Phi_E =  5.0/1000.; //5mrad, Angular resolution for electron, determined by GEM tracking
const double Sigma_Theta_G = 0.0; //Rad, Angular resolution for photon, determined by EC cluster reconstruction 
const double Sigma_Phi_G = 0.0; //Rad, Angular resolution for electron, determined by GEM tracking

Int_t CheckLaws(TLorentzVector* P_t,TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_g, TLorentzVector* P_h);
Double_t GetMM( TLorentzVector* P_t,TLorentzVector* P_E0, TLorentzVector* P_e, TLorentzVector* P_g);

int MissingMass(TString particle,int energy_flag){
	gStyle->SetOptStat(1);

	const double Lumi = 1.0e36; // cm-2*s-1, for He3 nuclear not for nucleons
	const double KHz = 1e-3;
	const double nBcm2 = 1e-33;

	const double eMass = 0.511/1000.;//electron mass 
    double hMass = 0.0;

	TString Target = "";
	/*Define Rootfile and variables{{{*/
	TChain *T=new TChain("T");
	if(particle=="p"){
		 Target = "NH3"; hMass = 0.938272;	
		if(energy_flag==11){
		 T->Add("./Proton_11GeV.root");
		}
		else if(energy_flag==8){
			T->Add("./Proton_8.8GeV.root");
		}
	}else if(particle=="n"){
		 Target = "3he"; hMass = 0.939565;	
		if(energy_flag==11){
			T->Add("./Neutron_11GeV.root"); 
		}
		else if(energy_flag==8){
			T->Add("./Neutron_8.8GeV.root");
		}
	}

	Double_t vertexz, E0; 
	Double_t ePx_ini, ePy_ini, ePz_ini, hPx_ini, hPy_ini,hPz_ini;
	Double_t ePx, ePy, ePz, gPx, gPy, gPz, hPx, hPy, hPz;
	Double_t eP_ini,hP_ini, eP, gP,hP;
	Double_t Q2, x, t, phi, XS_P, XS_M,PSF;
	Double_t ePhi_ini,hPhi_ini, ePhi, gPhi,hPhi;
	Double_t eTheta_ini,hTheta_ini, eTheta, gTheta,hTheta;
	Int_t Ngen,Nacc;

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
	/*}}}*/
	Long64_t N_entries=T->GetEntries();
	T->GetEntry(0);
	Long64_t N_gen = Ngen;
	cout<<"total generated events number: "<<N_gen<<"--- and accepted: "<<N_entries<<endl;

	TLorentzVector *P_E0 = new TLorentzVector();//incoming electron
	TLorentzVector *P_e_i = new TLorentzVector();//scattered electron
	TLorentzVector *P_e_f = new TLorentzVector();//scattered electron with Eloss
	TLorentzVector *P_g = new TLorentzVector();//photon
	TLorentzVector *P_h = new TLorentzVector();//proton or neutron
	TLorentzVector *P_t = new TLorentzVector();//target, either proton or neutron
	TLorentzVector *P_e_res = new TLorentzVector();//scattered electron with resolution
	TLorentzVector *P_g_res = new TLorentzVector();//photon with resolution

	TH1F *hMM = new TH1F("hMM",Form("Missing Mass of %s", particle.Data()), 200, 0.0, 2.0);
	hMM->SetXTitle("Nucleon Missing Mass (GeV)");
	hMM->SetYTitle("Rate (Hz)");
	TH1F *hMM_res = new TH1F("hMM_res",Form("Missing Mass of %s (with detector resolutions)", particle.Data()), 200, 0.0, 2.0);
	hMM_res->SetXTitle("Nucleon Missing Mass (GeV)");
	hMM_res->SetYTitle("Rate (Hz)");
   
	gRandom->SetSeed(0);
	//Define new measured quantities that include detector resolutions
	double eE_i, eE_f, gE, hE;
	double eE_i_res, eE_f_res,eP_res, eTheta_res, ePhi_res, ePx_res, ePy_res, ePz_res;
	double gP_res, gTheta_res, gPhi_res, gPx_res, gPy_res, gPz_res;
	double MM = 0.0,MM_res = 0.0;
	for(int i=0;i<N_entries; i++){
        T->GetEntry(i);
		if(phi/DEG>2.0&&phi/DEG<358){//any additional cuts should be added in here
		double event_weight=(XS_P+XS_M)*PSF/N_gen*Lumi*nBcm2;   //put into Hz, Note, XS in nb
      
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
	    hMM->Fill(MM, event_weight);

		//////////////////////////////////////////////////////////////////////////
		//Now consider the detector resolution here
		//////////////////////////////////////////////////////////////////////////
		//
		/////////////////////////////////////////
		//Electron w/o E-loss
		eP_res = gRandom->Gaus(eP_ini, Sigma_EC_E*sqrt(eP_ini));//GeV, for electron, E ~= P
	    eTheta_res = gRandom->Gaus(eTheta_ini*DEG, Sigma_Theta_E);//rad
	    ePhi_res = gRandom->Gaus(ePhi_ini*DEG, Sigma_Phi_E);//rad

        ePx_res = eP_res * sin(eTheta_res)*cos(ePhi_res); 
        ePy_res = eP_res * sin(eTheta_res)*sin(ePhi_res); 
        ePz_res = eP_res * cos(eTheta_res);
		
		eE_i_res = sqrt(eP_res*eP_res + eMass*eMass);
	    P_e_res->SetPxPyPzE(ePx_res, ePy_res, ePz_res, eE_i_res);	

		/////////////////////////////////////////
		//Electron with E-loss
		/*
		eP_res = gRandom->Gaus(eP, Sigma_EC_E*sqrt(eP));//GeV, for electron, E ~= P
	    eTheta_res = gRandom->Gaus(eTheta*DEG, Sigma_Theta_E);//rad
	    ePhi_res = gRandom->Gaus(ePhi*DEG, Sigma_Phi_E);//rad

        ePx_res = eP_res * sin(eTheta_res)*cos(ePhi_res); 
        ePy_res = eP_res * sin(eTheta_res)*sin(ePhi_res); 
        ePz_res = eP_res * cos(eTheta_res);
		
		eE_f_res = sqrt(eP_res*eP_res + eMass*eMass);
	    //P_e_res->SetPxPyPzE(ePx_res, ePy_res, ePz_res, eE_i_res);	
        */
		/////////////////////////////////////////
		//Photon
		gP_res = gRandom->Gaus(gP, Sigma_EC_G*sqrt(gP));//GeV, for photon, E = P
	    gTheta_res = gRandom->Gaus(gTheta*DEG, Sigma_Theta_G);//rad
	    gPhi_res = gRandom->Gaus(gPhi*DEG, Sigma_Phi_G);//rad

        gPx_res = gP_res * sin(gTheta_res)*cos(gPhi_res); 
        gPy_res = gP_res * sin(gTheta_res)*sin(gPhi_res); 
        gPz_res = gP_res * cos(gTheta_res);
	    P_g_res->SetPxPyPzE(gPx_res, gPy_res, gPz_res, gP_res);	

		/////////////////////////////////////////
        MM_res = GetMM(P_t, P_E0, P_e_res, P_g_res);	
	    hMM_res->Fill(MM_res, event_weight);	
		/////////////////////////////////////////
		}
	}
   delete P_g, P_t, P_e_i, P_e_f, P_h, P_e_res, P_g_res;	

   TCanvas *c1 = new TCanvas("c1","c1", 800,600);
   hMM_res->Draw(); hMM_res->Fit("gaus");
   hMM->SetLineColor(8); hMM->Draw("same");
   TLatex *t1=new TLatex();
   t1->SetNDC();
   t1->SetTextColor(4);
   t1->DrawLatex(0.25,0.40, Form("#deltaE_{e}=%2.1f%/#sqrt{E}, #deltaE_{g}=%2.1f%/#sqrt{E}", Sigma_EC_E*100, Sigma_EC_G*100.));
   t1->DrawLatex(0.25,0.30, Form("#delta#theta_{e}=%2.1fmrad, #delta#phi_{e}=%2.1fmrad", Sigma_Theta_E*1000, Sigma_Phi_E*1000.));
   t1->DrawLatex(0.25,0.20, Form("#delta#theta_{g}=%2.1fmrad, #delta#phi_{g}=%2.1fmrad", Sigma_Theta_G*1000, Sigma_Phi_G*1000.));

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
