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

void get_rate(TString particle,int energy_flag){
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
		 T->Add("../proton_projection/DVCS_Proton_11GeV.root");
		}
		else if(energy_flag==8){
			EBeam =8.8;
			T->Add("../proton_projection/DVCS_Proton_8.8GeV.root");
		}
	}else if(particle=="n"){
		 Target = "3he";	
		if(energy_flag==11){
			EBeam =11.0;
			T->Add("../neutron_projection/DVCS_Neutron_11GeV_New0.root"); 
		}
		else if(energy_flag==8){
			EBeam =8.8;
			T->Add("../neutron_projection/DVCS_Neutron_8.8GeV_New0.root"); 
		}
	}

	Double_t vertexz, E0, W, Wp; 
	Double_t ePx_ini, ePy_ini, ePz_ini, hPx_ini, hPy_ini,hPz_ini;
	Double_t ePx, ePy, ePz, gPx, gPy, gPz, hPx, hPy, hPz;
	Double_t eP_ini,hP_ini, eP, gP,hP;
	Double_t Q2, x, t, phi, XS_P, XS_M,PSF;
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
	T->SetBranchAddress("PSF", &PSF);
	T->SetBranchAddress("Ngen", &Ngen);
	T->SetBranchAddress("Nacc", &Nacc);

	Double_t SigmaPP_L, SigmaPM_L, SigmaMP_L, SigmaMM_L, Sigma_L, BSA_L, TSA_L, DSA_L;
	Double_t SigmaPP_Tx, SigmaPM_Tx, SigmaMP_Tx, SigmaMM_Tx, Sigma_Tx, BSA_Tx, TSA_Tx, DSA_Tx;
	Double_t SigmaPP_Ty, SigmaPM_Ty, SigmaMP_Ty, SigmaMM_Ty, Sigma_Ty, BSA_Ty, TSA_Ty, DSA_Ty;

	T->SetBranchAddress("SigmaPP_L", &SigmaPP_L);
	T->SetBranchAddress("SigmaPM_L", &SigmaPM_L);
	T->SetBranchAddress("SigmaMP_L", &SigmaMP_L);
	T->SetBranchAddress("SigmaMM_L", &SigmaMM_L);
	T->SetBranchAddress("Sigma_L", &Sigma_L);
	T->SetBranchAddress("BSA_L", &BSA_L);
	T->SetBranchAddress("TSA_L", &TSA_L);
	T->SetBranchAddress("DSA_L", &DSA_L);

	T->SetBranchAddress("SigmaPP_Tx", &SigmaPP_Tx);
	T->SetBranchAddress("SigmaPM_Tx", &SigmaPM_Tx);
	T->SetBranchAddress("SigmaMP_Tx", &SigmaMP_Tx);
	T->SetBranchAddress("SigmaMM_Tx", &SigmaMM_Tx);
	T->SetBranchAddress("Sigma_Tx", &Sigma_Tx);
	T->SetBranchAddress("BSA_Tx", &BSA_Tx);
	T->SetBranchAddress("TSA_Tx", &TSA_Tx);
	T->SetBranchAddress("DSA_Tx", &DSA_Tx);

	T->SetBranchAddress("SigmaPP_Ty", &SigmaPP_Ty);
	T->SetBranchAddress("SigmaPM_Ty", &SigmaPM_Ty);
	T->SetBranchAddress("SigmaMP_Ty", &SigmaMP_Ty);
	T->SetBranchAddress("SigmaMM_Ty", &SigmaMM_Ty);
	T->SetBranchAddress("Sigma_Ty", &Sigma_Ty);
	T->SetBranchAddress("BSA_Ty", &BSA_Ty);
	T->SetBranchAddress("TSA_Ty", &TSA_Ty);
	T->SetBranchAddress("DSA_Ty", &DSA_Ty);

    Double_t tmin;
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

    //CLEO
	TFile *file_negative=new TFile(Form("../acceptance/acceptance_solid_CLEO_SIDIS_%s_negative_output.root",Target.Data()));
	TH2F *accep_ele_forward=(TH2F*)file_negative->Get("acceptance_forwardangle");       //negative particle forward angle acceptance
	TH2F *accep_ele_large=(TH2F*)file_negative->Get("acceptance_largeangle");           //negative particle large angle acceptance

	TFile* file_neutral=new TFile(Form("../acceptance/acceptance_solid_CLEO_SIDIS_%s_neutral_output.root", Target.Data()));
	TH2F *accep_pho_forward=(TH2F*)file_neutral->Get("acceptance_forwardangle");       //negative particle forward angle acceptance
	TH2F *accep_pho_large=(TH2F*)file_neutral->Get("acceptance_largeangle");       //negative particle forward angle acceptance
	
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
		//if(phi>2.0&&phi<358&&W>2.0 && Sigma_L >1e-25 && gTheta>8.5){//any additional cuts should be added in here
		if(phi>=0.0&&phi<360&&W>2.0 && Sigma_L >1e-25 && gTheta>8.5){//any additional cuts should be added in here

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
				
			double event_weight=(XS_P+XS_M)/2.0*PSF/N_gen*nBcm2*Lumi ;   //in Hz
	//		double event_weight=(Sigma_Tx*1.0)*PSF/N_gen*nBcm2*Lumi ;   //in Hz

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
//	file_negative->Close();
//	file_neutral->Close();
	
}

