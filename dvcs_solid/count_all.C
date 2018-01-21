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
using namespace std;

double count_single(const TString & TPol,const TString & energy,const int Q2_flag, const int xb_flag, const int t_flag){
	const TString Target = "n";	
	const double Target_Pol = 0.65;
	const double Beam_Pol = 1.0;

	TString filename_L = Form("../database_%s/%s_%s_%d_%d_%d.dat",energy.Data(),TPol.Data(), Target.Data(), Q2_flag, xb_flag, t_flag);
	ifstream input_file(filename_L);
	cerr<<"--- Reading file from"<<filename_L.Data()<<endl;

	double Q2[12], xb[12], t[12], phi[12], Astat[12], TSA[12], BSA[12], DSA[12], TSA_Err[12], BSA_Err[12], DSA_Err[12];
	int dum;
	double count_bsa=0, count_tsa=0.0, count_dsa=0.0;
	for(int i=0;i<12;i++){
		input_file  >> dum >> dum >> Q2[i]  >> xb[i]  >> t[i]  >> phi[i]  >> Astat[i]  >> BSA[i]  >> TSA[i]  >> DSA[i];

		if(!(isnan( phi[i])) && !(isnan( Astat[i])) && Astat[i]>0){
			count_bsa += pow(1./(Astat[i]*Target_Pol), 2);
			count_tsa += pow(1./(Astat[i]*Beam_Pol), 2);
			count_dsa += pow(1./(Astat[i]), 2);

		}
	}
	input_file.close();

	return count_dsa;
}

void count_all(){
   int Q2_bin = 5;
   int x_bin = 5;
   int t_bin = 6;

   double total_count = 0.0;

   for(int i=0;i<Q2_bin;i++){
	   for(int j=0;j<Q2_bin;j++){
		   for(int k=0;k<Q2_bin;k++){

			   double count_temp = count_single("L", "11", i,j,k);

			   total_count += count_temp;
		   }
	   }
   }

   cerr<<"-- for 11 GeV at L, Total Count = "<< total_count<<endl;
}
