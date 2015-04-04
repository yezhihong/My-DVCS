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
#include "/work/halla/solid/yez/dvcs/VGG/GetXS/DVCSGrid.h"
using namespace std;
const double DEG = 180./3.1415926;
const double PI = 3.1415926;

//Here bin on Q2 and x first, and the t and phi binning will be in the next step
Int_t main(Int_t argc, char *argv[]){
	TString particle = "X"; cerr<<"-- What particle (p->proton, n->neutron)? "; cin >> particle;
	Double_t EBeam= 0.0;   cerr<<"-- What energy (11, 8.8 (GeV))? "; cin >> EBeam;  
    TString Target = "";
	if(particle=="n")
		Target= "Neutron";
	if(particle=="n")
		Target= "Proton";

	/*DVCS: Define Rootfile and variables{{{*/
	TString filename = "";
	if(particle=="p"){
		if(fabs(EBeam-11)<1.0)
			filename = "./Proton_11GeV.root";
		else if(fabs(EBeam-8.8)<1.)
			filename = "./Proton_8.8GeV.root";
	}else if(particle=="n"){
		if(fabs(EBeam-11)<1.0)
			filename = "./Neutron_11GeV.root";
		else if(fabs(EBeam-8.8)<1.)
			filename = "./Neutron_8.8GeV.root";
	}
	TFile *file = new TFile(filename.Data(), "update");
    TTree *T = (TTree*) file->Get("T");

	Double_t vertexz, E0; 
	Double_t ePx_ini, ePy_ini, ePz_ini;
    //Double_t	hP_ini, hPx_ini, hPy_ini,hPz_ini,hTheta_ini,hPhi_ini;
	Double_t ePx, ePy, ePz, gPx, gPy, gPz, hPx, hPy, hPz;
	Double_t eP_ini,eP, gP,hP;
	Double_t Q2, x, t, phi, XS_P, XS_M,PSF;
	Double_t ePhi_ini, ePhi, gPhi,hPhi;
	Double_t eTheta_ini, eTheta, gTheta,hTheta;
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

	/*Add new branches{{{*/
	Double_t BSA_L, TSA_L, DSA_L;
	Double_t Sigma_L, SigmaPP_L,SigmaPM_L, SigmaMP_L, SigmaMM_L;
	Double_t BSA_Tx, TSA_Tx, DSA_Tx;
	Double_t Sigma_Tx, SigmaPP_Tx,SigmaPM_Tx, SigmaMP_Tx, SigmaMM_Tx;
	Double_t BSA_Ty, TSA_Ty, DSA_Ty;
	Double_t Sigma_Ty, SigmaPP_Ty,SigmaPM_Ty, SigmaMP_Ty, SigmaMM_Ty;
			
	TBranch *br_Sigma_L   = T->Branch("Sigma_L",   &Sigma_L,   "Sigma_L/D");
	TBranch *br_SigmaPP_L = T->Branch("SigmaPP_L", &SigmaPP_L, "SigmaPP_L/D");
	TBranch *br_SigmaPM_L = T->Branch("SigmaPM_L", &SigmaPM_L, "SigmaPM_L/D");
	TBranch *br_SigmaMP_L = T->Branch("SigmaMP_L", &SigmaMP_L, "SigmaMP_L/D");
	TBranch *br_SigmaMM_L = T->Branch("SigmaMM_L", &SigmaMP_L, "SigmaMM_L/D");
	TBranch *br_BSA_L     = T->Branch("BSA_L",     &BSA_L,     "BSA_L/D");
	TBranch *br_TSA_L     = T->Branch("TSA_L",     &TSA_L,     "TSA_L/D");
	TBranch *br_DSA_L     = T->Branch("DSA_L",     &DSA_L,     "DSA_L/D");

	TBranch *br_Sigma_Tx   = T->Branch("Sigma_Tx",   &Sigma_Tx,   "Sigma_Tx/D");
	TBranch *br_SigmaPP_Tx = T->Branch("SigmaPP_Tx", &SigmaPP_Tx, "SigmaPP_Tx/D");
	TBranch *br_SigmaPM_Tx = T->Branch("SigmaPM_Tx", &SigmaPM_Tx, "SigmaPM_Tx/D");
	TBranch *br_SigmaMP_Tx = T->Branch("SigmaMP_Tx", &SigmaMP_Tx, "SigmaMP_Tx/D");
	TBranch *br_SigmaMM_Tx = T->Branch("SigmaMM_Tx", &SigmaMP_Tx, "SigmaMM_Tx/D");
	TBranch *br_BSA_Tx     = T->Branch("BSA_Tx",     &BSA_Tx,     "BSA_Tx/D");
	TBranch *br_TSA_Tx     = T->Branch("TSA_Tx",     &TSA_Tx,     "TSA_Tx/D");
	TBranch *br_DSA_Tx     = T->Branch("DSA_Tx",     &DSA_Tx,     "DSA_Tx/D");

	TBranch *br_Sigma_Ty   = T->Branch("Sigma_Ty",   &Sigma_Ty,   "Sigma_Ty/D");
	TBranch *br_SigmaPP_Ty = T->Branch("SigmaPP_Ty", &SigmaPP_Ty, "SigmaPP_Ty/D");
	TBranch *br_SigmaPM_Ty = T->Branch("SigmaPM_Ty", &SigmaPM_Ty, "SigmaPM_Ty/D");
	TBranch *br_SigmaMP_Ty = T->Branch("SigmaMP_Ty", &SigmaMP_Ty, "SigmaMP_Ty/D");
	TBranch *br_SigmaMM_Ty = T->Branch("SigmaMM_Ty", &SigmaMP_Ty, "SigmaMM_Ty/D");
	TBranch *br_BSA_Ty     = T->Branch("BSA_Ty",     &BSA_Ty,     "BSA_Ty/D");
	TBranch *br_TSA_Ty     = T->Branch("TSA_Ty",     &TSA_Ty,     "TSA_Ty/D");
	TBranch *br_DSA_Ty     = T->Branch("DSA_Ty",     &DSA_Ty,     "DSA_Ty/D");
	/*}}}*/ 
    TString Data_Dir = "/work/halla/solid/yez/dvcs/VGG/";
	TString TargetPol = ""; //Tx, Ty, L 
	Int_t Debug = 0, err = 0;

	DVCSGrid *grid = new DVCSGrid();
    grid->Init(EBeam, Target.Data(), Data_Dir.Data(), Debug);

	for(int i=0;i<N_entries;i++){
		T->GetEntry(i);

		grid->FindBin(Q2, x, -t, phi);// here -t was used instead of t
		//grid->PrintRanges();

	    /*Reading the Longitudinal Target Spin{{{*/	
   		SigmaPP_L = -1000; SigmaPM_L = -1000;  SigmaMP_L = -1000;  SigmaMM_L = -1000; 
		Sigma_L = -1000;   BSA_L = -1000; TSA_L = -1000; DSA_L = -1000;
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
		br_SigmaPP_L->Fill();  br_SigmaPM_L->Fill();  br_SigmaMP_L->Fill();  br_SigmaMM_L->Fill();
        br_Sigma_L->Fill();    br_BSA_L->Fill();  br_TSA_L->Fill();  br_DSA_L->Fill();
	    /*Reading the Longitudinal Target Spin}}}*/	
	
	    /*Reading the Transverse Target Spin on x{{{*/	
   		SigmaPP_Tx = -1000;SigmaPM_Tx = -1000; SigmaMP_Tx = -1000; SigmaMM_Tx = -1000; 
		Sigma_Tx = -1000;  BSA_Tx = -1000;TSA_Tx = -1000;DSA_Tx = -1000;
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
		}
		br_SigmaPP_Tx->Fill(); br_SigmaPM_Tx->Fill(); br_SigmaMP_Tx->Fill(); br_SigmaMM_Tx->Fill();
		br_Sigma_Tx->Fill();   br_BSA_Tx->Fill(); br_TSA_Tx->Fill(); br_DSA_Tx->Fill();
	    /*Reading the Transverse Target Spin on x}}}*/	
	
	    /*Reading the Transverse Target Spin on y{{{*/	
		SigmaPP_Ty = -1000;SigmaPM_Ty = -1000; SigmaMP_Ty = -1000; SigmaMM_Ty = -1000; 
		Sigma_Ty = -1000;  BSA_Ty = -1000;TSA_Ty = -1000;DSA_Ty = -1000;

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
		br_SigmaPP_Ty->Fill(); br_SigmaPM_Ty->Fill(); br_SigmaMP_Ty->Fill(); br_SigmaMM_Ty->Fill();
		br_Sigma_Ty->Fill();   br_BSA_Ty->Fill(); br_TSA_Ty->Fill(); br_DSA_Ty->Fill();
	    /*Reading the Transverse Target Spin on y}}}*/	
	}
	
	T->Write("",TObject::kOverwrite);
    file->Close();

	delete grid;
	return 0;
}
