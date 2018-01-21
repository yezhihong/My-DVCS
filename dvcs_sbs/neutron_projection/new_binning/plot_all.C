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
#include <TLatex.h>
#include <TLine.h>
//#include <TMatrix.h>
/*}}}*/
#include "plot_bin.h"
using namespace std;

void plot_all(){
	gStyle->SetLabelSize(0.1,"Y");
	gStyle->SetTitleSize(0.1,"Y");
	gStyle->SetTitleOffset(0.8,"Y");
	gStyle->SetLabelSize(0.1,"X");
	gStyle->SetTitleSize(0.1,"X"); 
	gStyle->SetPadBorderMode(0);
	
	TString energy="x";
	cerr<<"--- Energy Setting (8, 11 or 11p8) = "; cin >> energy;

	TCanvas *c1 = new TCanvas("c1","c1",1600,1200);
	c1->SetFillColor(10);
	c1->Divide(6,5,0,0);
	for (Int_t k=0;k!=5;k++){
		for (Int_t i=0;i!=5;i++){
			for (Int_t j=0;j!=6;j++){
				c1->cd(i*6+j+1);
				if (i!=4){
					gPad->SetBottomMargin(0);
				}else{
					gPad->SetBottomMargin(0.2);
				}
				gPad->SetTopMargin(0);

				if (j!=0){
					gPad->SetLeftMargin(0);
				}else{
					gPad->SetLeftMargin(0.3);
				}
				if (j!=5){
					gPad->SetRightMargin(0);
				}else{
					gPad->SetRightMargin(0.3);
				}

				gPad->SetFillColor(10);
				plot_bin("L",energy.Data(),k,i,j);
			}
		}
		c1->Print(Form("n_L_E%s_Q2_%d.png", energy.Data(), k));
		c1->Print(Form("n_L_E%s_Q2_%d.pdf", energy.Data(), k));
	}

	TCanvas *c2 = new TCanvas("c2","c2",1600,1200);
	c2->SetFillColor(10);
	c2->Divide(6,5,0,0);
	for (Int_t k=0;k!=5;k++){
		for (Int_t i=0;i!=5;i++){
			for (Int_t j=0;j!=6;j++){
				c2->cd(i*6+j+1);
				if (i!=4){
					gPad->SetBottomMargin(0);
				}else{
					gPad->SetBottomMargin(0.2);
				}
				gPad->SetTopMargin(0);

				if (j!=0){
					gPad->SetLeftMargin(0);
				}else{
					gPad->SetLeftMargin(0.3);
				}
				if (j!=5){
					gPad->SetRightMargin(0);
				}else{
					gPad->SetRightMargin(0.3);
				}

				gPad->SetFillColor(10);
				plot_bin("Tx",energy.Data(),k,i,j);
			}
		}
		c2->Print(Form("n_Tx_E%s_Q2_%d.png", energy.Data(), k));
		c2->Print(Form("n_Tx_E%s_Q2_%d.pdf", energy.Data(), k));
		}

	TCanvas *c3 = new TCanvas("c3","c3",1600,1200);
	c3->SetFillColor(10);
	c3->Divide(6,5,0,0);
	for (Int_t k=0;k!=5;k++){
		for (Int_t i=0;i!=5;i++){
			for (Int_t j=0;j!=6;j++){
				c3->cd(i*6+j+1);
				if (i!=4){
					gPad->SetBottomMargin(0);
				}else{
					gPad->SetBottomMargin(0.2);
				}
				gPad->SetTopMargin(0);

				if (j!=0){
					gPad->SetLeftMargin(0);
				}else{
					gPad->SetLeftMargin(0.3);
				}
				if (j!=5){
					gPad->SetRightMargin(0);
				}else{
					gPad->SetRightMargin(0.3);
				}

				gPad->SetFillColor(10);
				plot_bin("Ty",energy.Data(),k,i,j);
			}
		}
		c3->Print(Form("n_Ty_E%s_Q2_%d.png", energy.Data(), k));
		c3->Print(Form("n_Ty_E%s_Q2_%d.pdf", energy.Data(), k));
		}


}
