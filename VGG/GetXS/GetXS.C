#include <iostream>
#include <stdlib.h>
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
#include "DVCSGrid.h"
using namespace std;

int main(){

	const Double_t E0 = 11.0;//GeV
	const TString Target = "Neutron";
	const TString Data_Dir = "/work/halla/solid/yez/dvcs/VGG/";
	Int_t Debug = 3;

	Double_t Q2, xb, t, phi;
	TString TargetPol;//Tx, Ty, L

	DVCSGrid *grid = new DVCSGrid();
    grid->Init(E0, Target.Data(), Data_Dir.Data(), Debug);

	//reset the kinematic ranges and sizes if they have been changed from the default values
	//           Min, Max, Step
	//SetQ2Range(1.0, 9.0, 1.0); //GeV2
	//SetXbRange(0.05, 0.75, 0.05);
	//SetTRange(0.1, 2.5, 0.1); //GeV2
	//SetPhiRange(0.0, 360.0, 15.0); //degree
	
    Int_t err = -100;

	cerr<<"--- Q2 (GeV2) = "; cin >> Q2;
	cerr<<"---        xb = "; cin >> xb;
	cerr<<"---         t = "; cin >> t;
	cerr<<"--- phi (Deg) = "; cin >> phi;
	cerr<<"--- Target Polarization (L, Tx, Ty) = "; cin >> TargetPol;

	grid->FindBin(Q2, xb, t, phi);//last input, >0 to print the finding result
	grid->PrintRanges();
	err = grid->LoadXS(TargetPol);
	if(err>=0){
		grid->PrintXS();

		Double_t XS_PP = grid->GetSigmaPP();//++
		Double_t XS_PM = grid->GetSigmaPM();//+-
		Double_t XS_MP = grid->GetSigmaMP();//-+
		Double_t XS_MM = grid->GetSigmaMM();//--
		Double_t XS = grid->GetSigma();//average of four XSs above
		Double_t XS_BS = grid->GetSigmaBS();//beam spin asym
		Double_t XS_TS = grid->GetSigmaTS();//target spin asym
		Double_t XS_DS = grid->GetSigmaDS();//double spin asym
	}

	delete grid;
}
