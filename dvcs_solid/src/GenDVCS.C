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

#include "TGenSoLIDDVCS.h"
using namespace std;

int main(int argc, char *argv[])
{
	cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
	cout << "&&&& SoLID-DVCS Generator &&&&" << endl;
	cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
	if(argc!=4) {
		cout << "Usage: GenSoLID <Ebeam(GeV)> <Target(0=P, 1=N)> <Total_Event>" << endl;
		return -1;
	}
	const Double_t Ebeam = atof(argv[1]);// or 8.8 GeV
	const Int_t Target_Type = atoi(argv[2]); //target:0->P,1->N,2->D
	const Int_t nev=atoi(argv[3]); // Good for Proton/Neutron
	cout << Form("--> Beam = %4.2f (GeV), and Target = %d (0=P, 1=N, 2=D, Event_Gen = %d)",Ebeam,Target_Type, nev)<<endl;

	TGenSoLIDDVCS *gEv=new TGenSoLIDDVCS(Ebeam, Target_Type, 0, 0);//last two parameters unused

	Double_t Target_Length = 0.0; //cm
	Double_t Target_Density = 0.0;//g/cm3
	if(Target_Type==0){ // Proton target, DNP-NH3,817kg/m3 for solid state
	  Target_Length = 3.0;//cm
	  Target_Density = 0.82;//g/cm3
	}else if(Target_Type==1){ // Neutron target, He3, note that N will be 1/3 from the luminosity
	  Target_Length = 40.;//cm
	  Target_Density = 0.001345;//g/cm3
	}
	else{
	  cerr<<"---- I don't know the target!"<<endl;
	  return -100;
	}

	gEv->SetTargetParam(Target_Length, Target_Density);

	////////////////////////////////////////////////////////
	// Geometry of spectrometer for electron and the calorimeter for photon detection
    // Use both FAEC and LAEC to detect electrons
	////////////////////////////////////////////////////////

	//Note that for SoLID only it is (Theta-Min, Theta-Max, Phi-Min, Phi-Max, Mom-Min, Mom-Max), in rad and GeV
	//     for other configurations, it is alway "Max" first then "Min"..
	//
	//Detect photons at FAEC from  6 to 16 where the actual range is 8 to 14.8,
	//           and at LAEC from 16 to 26 where the actual range is 16 to 24.0
	gEv->SetSoLIDAcceptance(6.*TMath::DegToRad(), 26*TMath::DegToRad(),0*TMath::DegToRad(), 360*TMath::DegToRad(),0.5,11.);

	//Only detect photons at FAEC from 6 to 16 where the actuall values is 8 to 14.8
	//Note that it is (Theta-Max, Theta-Min, Phi-Max, Phi-Min), in rad
	gEv->SetSoLIDCaloAcceptance(6.*TMath::DegToRad(), 26*TMath::DegToRad(),0*TMath::DegToRad(), 360*TMath::DegToRad());
    //////////////////////////////////////////////////////

	////////////////////////
	//gEv->SetFermi(kTRUE); // For Fermi only
	////////////////////////

	Int_t j=0;
	cout<<endl<<"--> Generating "<<nev<<" events..."<<endl<<endl;

	for(Int_t i=0;i<nev;i++){
		if ((i%1000)==0) cout<<i<<"/"<<nev<<"\r";
		gEv->Settmin(-2.); //Use this method to change tmin (default -2 GeV)
		gEv->Settmax(0.);  //Use this method to change tmax (default 0 GeV)
		gEv->GenerateVertex();//Generates the vertex
		gEv->ExtBrem();//Make external pre-vertex radiative corrections
		gEv->GenKinSoLID();//Computes scattered electron kinematics, 
		//-------------------for SoLID(wider), 6<Theta<26, 0.5<Ee<11., 0.05<xb<0.8, 0.1<Q2<13.0

		gEv->IntRCBef();//Internal radiative corrections (before vertex)
		/////////////////////
	//	gEv->GenFermiIni();//To use if fFermi is set
		/////////////////////
		
		if(gEv->ComputeElectron()){//Compute the scattered electron 4-vector
			gEv->IntRCAft();//Internal radiative corrections (after vertex)
			if(gEv->HitsSoLID(gEv->GetScatteredElectron())){//Check hor spectro accep.
				gEv->ComputeDVCS();//Compute the gamma*p->gamma p collision
				gEv->ApplySoLIDAcc();//Rotates all final vectors around the beam axis
				if(gEv->HitsSoLIDCalo(gEv->GetFinalPhoton())){//Checks calo acceptance
					gEv->XSecSum(0);//(1)-->Get only BH term, othetwise BH+VCS+I, unit is in nb/....
					gEv->Write2File();//This writes the event in a temporary file
					j++;
				}
			}
		}
	}
	gEv->CloseTmpFile();
	cout<<endl<<"--> Accepted Events = "<<j<<endl<<endl;

	//gEv->DumpFinalFile("outdvcs.txt",nev,j);
	
	TString Target_Name = "Test";
	TString Energy = "11";
	if(Target_Type==0)
	   Target_Name = "Proton"; 
	if(Target_Type==1)
		Target_Name = "Neutron"; 
	if(Target_Type==2)
		Target_Name = "Deutron"; 
	if(abs(Ebeam-11)<0.001)
		Energy="11";
	if(abs(Ebeam-8.8)<0.001)
		Energy="8.8";
    cout<<Form("Saving data into file: DVCS_%s_%sGeV.root",Target_Name.Data(), Energy.Data())<<endl;
	gEv->DumpRootFile(Form("DVCS_%s_%sGeV.root",Target_Name.Data(), Energy.Data()),nev,j);

	return 0;
}
