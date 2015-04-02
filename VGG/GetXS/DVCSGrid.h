//DVCSGrid.h, v1.0, Wed Apr 1 15:26
// Author: Z. Ye
//
//-----------------------------------------
//   I generated the cross sections with VGG,
// by looping the Q2, xb, t and phi within
// certain ranges and step size. Each table
// stores a grid point of Q2, xb and t, and 
// it contains XS values with different phi bins.
// 
//   This code is designed to obtain the XSs 
// from the tables with known Q2, xb, t and phi
// values. Linear correlation is assumed between
// two grid point. e.g., for fixed xb, t and phi,
// Q2 falls between Q2_low and Q2_high, and  its
// XS value is the linear conbination of the XSs
// at Q2_low and Q2_high.  
//   

#ifndef __DVCSGrid__
#define __DVCSGrid__
class DVCSGrid{
	public:

		/*Constructor and Destructor{{{*/
		DVCSGrid(){
			//constructor
			fQ2_Min = 1.0;//GeV
			fQ2_Max = 9.0;//GeV
			fQ2_Step = 1.0;//GeV
			fXb_Min = 0.05;//
			fXb_Max = 0.75;//
			fXb_Step = 0.05;//
			fT_Min = 0.1;//, -t, GeV
			fT_Max = 2.5;//, -t, GeV
			fT_Step = 0.1;//, -t, GeV
			fPhi_Min = 0.0;//, degree
			fPhi_Max = 360.;//, degree
			fPhi_Step = 15.;//, degreee

			fE0 = 11.0; //GeV
			fTarget = "Neutron";
			fTablesDir = "/work/halla/solid/yez/dvcs/VGG/";

			fMass_Neutron = 0.939565378;
		}

		DVCSGrid(const DVCSGrid&){
		}

		~DVCSGrid(){
			//destructor
		}
		/*Constructor and Destructor}}}*/

		/*Init & Set Range{{{*/
		void Init( const Double_t kE0, const TString& kTarget, const TString& kTablesDir){
			fE0 = kE0;
			fTarget = kTarget;
			fTablesDir = kTablesDir;
		}

		void SetEnergy(const Double_t kE0){
			fE0 = kE0;
		}

		void SetTarget(const TString& kTarget){
			fTarget = kTarget;
		}

		void SetTablesDir(const TString& kTablesDir){
			fTablesDir = kTablesDir;
		}

		void SetQ2Range(const Double_t kQ2_Min, const Double_t kQ2_Max, const Double_t kQ2_Step){
			fQ2_Max = kQ2_Max;
			fQ2_Min = kQ2_Min;
			fQ2_Step = kQ2_Step;
		}

		void SetXbRange(const Double_t kXb_Min, const Double_t kXb_Max, const Double_t kXb_Step){
			fXb_Max = kXb_Max;
			fXb_Min = kXb_Min;
			fXb_Step = kXb_Step;
		}

		void SetTRange(const Double_t kT_Min, const Double_t kT_Max, const Double_t kT_Step){
			fT_Max = kT_Max;
			fT_Min = kT_Min;
			fT_Step = kT_Step;
		}

		void SetPhiRange(const Double_t kPhi_Min, const Double_t kPhi_Max, const Double_t kPhi_Step){
			fPhi_Max = kPhi_Max;
			fPhi_Min = kPhi_Min;
			fPhi_Step = kPhi_Step;
		}
		/*Set Range}}}*/

		/*Load and Fine XS{{{*/
		void FindBin(const Double_t kQ2, const Double_t kXb, const Double_t kT, const Double_t kPhi){

			fQ2 = kQ2;
			fXb = kXb;
			fT = kT;
			fPhi = kPhi;

			fQ2_Bin = (int)( (fQ2 - fQ2_Min) / fQ2_Step );
			fXb_Bin = (int)( (fXb - fXb_Min) / fXb_Step );

			Double_t fEp = fE0 -(fQ2/(2.*fMass_Neutron*fXb));
			Double_t fNu = fE0-fEp;
			Double_t tmin = -fQ2-(fQ2*(1.-fXb)*(fNu-sqrt(fQ2+fNu*fNu))/(fXb*(fMass_Neutron+fNu-sqrt(fQ2+fNu*fNu))));
			tmin = tmin+0.02*tmin;
			tmin = int(1000*tmin)/1000;
			tmin=-1.*tmin+0.005;

			fT_Min = tmin;
			fT_Bin = (int)( (fT - fT_Min)/fT_Step );

			fPhi_Bin = (int)( (fPhi - fPhi_Min) / fPhi_Step);
		}

		void LoadXS(const TString& kTargetPol){
			TString filename = "";
			double kPhi0,kXS_PP0,kXS_PM0, kXS_MP0, kXS_MM0, kXS0;
			double kPhi1,kXS_PP1,kXS_PM1, kXS_MP1, kXS_MM1, kXS1;
			double kPhi2,kXS_PP2,kXS_PM2, kXS_MP2, kXS_MM2, kXS2;
			double kPhi3,kXS_PP3,kXS_PM3, kXS_MP3, kXS_MM3, kXS3;
			Double_t phi, xs_pp, xs_pm, xs_mp, xs_mm,xs;
			Double_t phi1, xs_pp1, xs_pm1, xs_mp1, xs_mm1,xs1;

			Double_t kQ2_low =  fQ2_Min + fQ2_Bin*fQ2_Step;
			Double_t kQ2_high = kQ2_low + fQ2_Step;
			Double_t kXb_low =  fXb_Min + fXb_Bin*fXb_Step;
			Double_t kXb_high = kXb_low + fXb_Step;
			Double_t kT_low =  fT_Min + fT_Bin*fT_Step;
			Double_t kT_high = kT_low + fT_Step;
			Double_t kPhi_low =  fPhi_Min + fPhi_Bin*fPhi_Step;
			//Double_t kPhi_high = kPhi_low + fPhi_Step;

			/*Closest{{{*/
			filename = Form("%s/%s_%d/E%d-%s-Q2-%d-xb-%d-t-%d.dat", 
					fTablesDir.Data(),
					fTarget.Data(),
					(int)(fE0),
					(int)(fE0),
					kTargetPol.Data(),
					fQ2_Bin,
					fXb_Bin,
					fT_Bin);
			ifstream file_temp0(filename.Data());

			while(!file_temp0.eof()){
				file_temp0 >> phi >> xs_pp >> xs_pm >> xs_mp >> xs_mm >> xs; 
				if(fabs(phi-kPhi_low)<fPhi_Step*0.1){
					file_temp0 >> phi1 >> xs_pp1 >> xs_pm1 >> xs_mp1 >> xs_mm1 >> xs1; 
				
					kPhi0 = phi;
					kXS_PP0 = xs_pp * GetCorrLinear(fPhi, phi, phi1, xs_pp, xs_pp1); 
					kXS_PM0 = xs_pm * GetCorrLinear(fPhi, phi, phi1, xs_pm, xs_pm1); 
					kXS_MP0 = xs_mp * GetCorrLinear(fPhi, phi, phi1, xs_mp, xs_mp1); 
					kXS_MM0 = xs_mm * GetCorrLinear(fPhi, phi, phi1, xs_mm, xs_mm1); 
					kXS0 = xs * GetCorrLinear(fPhi, phi, phi1, xs, xs1); 
					break;
				}
			}
			file_temp0.close();
			/*}}}*/

			/*Q2-up{{{*/
			filename = Form("%s/%s_%d/E%d-%s-Q2-%d-xb-%d-t-%d.dat", 
					fTablesDir.Data(),
					fTarget.Data(),
					(int)(fE0),
					(int)(fE0),
					kTargetPol.Data(),
					fQ2_Bin+1,
					fXb_Bin,
					fT_Bin);
			ifstream file_temp1(filename.Data());

			while(!file_temp1.eof()){
				file_temp1 >> phi >> xs_pp >> xs_pm >> xs_mp >> xs_mm >> xs; 
				if(fabs(phi-kPhi_low)<fPhi_Step*0.1){
					file_temp1 >> phi1 >> xs_pp1 >> xs_pm1 >> xs_mp1 >> xs_mm1 >> xs1; 
				
					kPhi1 = phi;
					kXS_PP1 = xs_pp * GetCorrLinear(fPhi, phi, phi1, xs_pp, xs_pp1); 
					kXS_PM1 = xs_pm * GetCorrLinear(fPhi, phi, phi1, xs_pm, xs_pm1); 
					kXS_MP1 = xs_mp * GetCorrLinear(fPhi, phi, phi1, xs_mp, xs_mp1); 
					kXS_MM1 = xs_mm * GetCorrLinear(fPhi, phi, phi1, xs_mm, xs_mm1); 
					kXS1 = xs * GetCorrLinear(fPhi, phi, phi1, xs, xs1); 
					break;
				}
			}
			file_temp1.close();
			/*}}}*/

			/*Xb-up{{{*/
			filename = Form("%s/%s_%d/E%d-%s-Q2-%d-xb-%d-t-%d.dat", 
					fTablesDir.Data(),
					fTarget.Data(),
					(int)(fE0),
					(int)(fE0),
					kTargetPol.Data(),
					fQ2_Bin,
					fXb_Bin+1,
					fT_Bin);
			ifstream file_temp2(filename.Data());

			while(!file_temp2.eof()){
				file_temp2 >> phi >> xs_pp >> xs_pm >> xs_mp >> xs_mm >> xs; 
				if(fabs(phi-kPhi_low)<fPhi_Step*0.1){
					file_temp2 >> phi1 >> xs_pp1 >> xs_pm1 >> xs_mp1 >> xs_mm1 >> xs1; 
				
					kPhi2 = phi;
					kXS_PP2 = xs_pp * GetCorrLinear(fPhi, phi, phi1, xs_pp, xs_pp1); 
					kXS_PM2 = xs_pm * GetCorrLinear(fPhi, phi, phi1, xs_pm, xs_pm1); 
					kXS_MP2 = xs_mp * GetCorrLinear(fPhi, phi, phi1, xs_mp, xs_mp1); 
					kXS_MM2 = xs_mm * GetCorrLinear(fPhi, phi, phi1, xs_mm, xs_mm1); 
					kXS2 = xs * GetCorrLinear(fPhi, phi, phi1, xs, xs1); 
					break;
				}
			}
			file_temp2.close();
			/*}}}*/

			/*T-up{{{*/
			filename = Form("%s/%s_%d/E%d-%s-Q2-%d-xb-%d-t-%d.dat", 
					fTablesDir.Data(),
					fTarget.Data(),
					(int)(fE0),
					(int)(fE0),
					kTargetPol.Data(),
					fQ2_Bin,
					fXb_Bin+1,
					fT_Bin);
			ifstream file_temp3(filename.Data());

			while(!file_temp3.eof()){
				file_temp3 >> phi >> xs_pp >> xs_pm >> xs_mp >> xs_mm >> xs; 
				if(fabs(phi-kPhi_low)<fPhi_Step*0.1){
					file_temp3 >> phi1 >> xs_pp1 >> xs_pm1 >> xs_mp1 >> xs_mm1 >> xs1; 
				
					kPhi3 = phi;
					kXS_PP3 = xs_pp * GetCorrLinear(fPhi, phi, phi1, xs_pp, xs_pp1); 
					kXS_PM3 = xs_pm * GetCorrLinear(fPhi, phi, phi1, xs_pm, xs_pm1); 
					kXS_MP3 = xs_mp * GetCorrLinear(fPhi, phi, phi1, xs_mp, xs_mp1); 
					kXS_MM3 = xs_mm * GetCorrLinear(fPhi, phi, phi1, xs_mm, xs_mm1); 
					kXS3 = xs * GetCorrLinear(fPhi, phi, phi1, xs, xs1); 
					break;
				}
			}
			file_temp3.close();
			/*}}}*/

			Double_t cor1 = 0.0, cor2 = 0.0, cor3 = 0.0;

			cor1 = GetCorrLinear(fQ2, kQ2_low, kQ2_high, kXS_PP0, kXS_PP1);
			cor2 = GetCorrLinear(fXb, kXb_low, kXb_high, kXS_PP0, kXS_PP2);
			cor3 = GetCorrLinear(fT, kT_low, kT_high, kXS_PP0, kXS_PP3);
			fSigmaPP = kXS_PP0 * cor1 * cor2 * cor3;

			cor1 = GetCorrLinear(fQ2, kQ2_low, kQ2_high, kXS_PM0, kXS_PM1);
			cor2 = GetCorrLinear(fXb, kXb_low, kXb_high, kXS_PM0, kXS_PM2);
			cor3 = GetCorrLinear(fT, kT_low, kT_high, kXS_PM0, kXS_PM3);
			fSigmaPM = kXS_PM0 * cor1 * cor2 * cor3;

			cor1 = GetCorrLinear(fQ2, kQ2_low, kQ2_high, kXS_MP0, kXS_MP1);
			cor2 = GetCorrLinear(fXb, kXb_low, kXb_high, kXS_MP0, kXS_MP2);
			cor3 = GetCorrLinear(fT, kT_low, kT_high, kXS_MP0, kXS_MP3);
			fSigmaMP = kXS_MP0 * cor1 * cor2 * cor3;

			cor1 = GetCorrLinear(fQ2, kQ2_low, kQ2_high, kXS_MM0, kXS_MM1);
			cor2 = GetCorrLinear(fXb, kXb_low, kXb_high, kXS_MM0, kXS_MM2);
			cor3 = GetCorrLinear(fT, kT_low, kT_high, kXS_MM0, kXS_MM3);
			fSigmaMM = kXS_MM0 * cor1 * cor2 * cor3;

			cor1 = GetCorrLinear(fQ2, kQ2_low, kQ2_high, kXS0, kXS1);
			cor2 = GetCorrLinear(fXb, kXb_low, kXb_high, kXS0, kXS2);
			cor3 = GetCorrLinear(fT, kT_low, kT_high, kXS0, kXS3);
			fSigma = kXS0 * cor1 * cor2 * cor3;
		}

		Double_t GetCorrLinear(const Double_t x, const Double_t x1, const Double_t x2, const Double_t y1, const Double_t y2){
			Double_t A = y2 - y1;
			Double_t B = (x - x1) / (x2-x1);
			Double_t y = A*B + y1;
			return y/y1;//return the correction for y1, not y itself
		}
		/*Load and Fine XS}}}*/

		/*return XS{{{*/
		Double_t GetSigmaBS(){
			Double_t xs = 1./4. * (fSigmaPP + fSigmaPM - fSigmaMP - fSigmaMM);
			return xs;
		}

		Double_t GetSigmaTS(){
			Double_t xs = 1./4. * (fSigmaPP - fSigmaPM + fSigmaMP - fSigmaMM);
			return xs;
		}

		Double_t GetSigmaDS(){
			Double_t xs = 1./4. * (fSigmaPP - fSigmaPM - fSigmaMP + fSigmaMM);
			return xs;
		}

		Double_t GetSigma(){
			return fSigma;
		}

		Double_t GetSigmaPP(){
			return fSigmaPM;
		}

		Double_t GetSigmaPM(){
			return fSigmaPM;
		}

		Double_t GetSigmaMP(){
			return fSigmaPM;
		}

		Double_t GetSigmaMM(){
			return fSigmaPM;
		}
		/*Return XS}}}*/
		
		/*Print{{{*/
		void PrintRanges(){
	        cerr<<"==============================================="<<endl;
		    cerr<<Form("%8s %8s %8s %8s","VAR","Min","Max","STEP")<<endl;	
		    cerr<<Form("%8s %8.4f %8.4f %8.4f","Q2", fQ2_Min, fQ2_Max, fQ2_Step)<<endl;	
		    cerr<<Form("%8s %8.4f %8.4f %8.4f","xb", fXb_Min, fXb_Max, fXb_Step)<<endl;	
		    cerr<<Form("%8s %8.4f %8.4f %8.4f","t",  fT_Min,  fT_Max,  fT_Step)<<endl;	
		    cerr<<Form("%8s %8.4f %8.4f %8.4f","phi",fPhi_Min,fPhi_Max,fPhi_Step)<<endl;	
	        cerr<<"==============================================="<<endl;
		}

		void PrintXS(){
	        cerr<<"==============================================="<<endl;
	        cerr<<Form("    Q2 = %8.4f, xb = %8.4f, t = %8.4f, phi = %8.4f", fQ2, fXb, fT, fPhi)<<endl;    
		    cerr<<Form("@@@ XS_PP = %10.6e, XS_PM = %10.6e, XS_MP = %10.6e, XS_MM = %10.6e, XS_All = %10.6e",
					        fSigmaPP, fSigmaPM, fSigmaMP, fSigmaMM, fSigma)<<endl;	
		    cerr<<Form("&&& XS_BS = %10.6e, XS_TS = %10.6e, XS_DS = %10.6e",
					        GetSigmaBS(), GetSigmaTS(), GetSigmaDS())<<endl;	
				cerr<<"==============================================="<<endl;
		}

		/*Print}}}*/

	protected:
		/*variables{{{*/
		Double_t fE0; //incoming beam energy, 11GeV or 8.8 GeV, GeV
		Double_t fMass_Neutron; //should be a constant but I still initialize it
		Double_t fQ2; //GeV2
		Double_t fXb;
		Double_t fT; //GeV2
		Double_t fPhi; //Degree

		Double_t fQ2_Min; //Q2 low range, GeV2
		Double_t fQ2_Max; //Q2 low range, GeV2
		Double_t fQ2_Step; //Q2 step in the generator, GeV2
		Double_t fXb_Min; //xbj low range,
		Double_t fXb_Max; //xbj low range,
		Double_t fXb_Step; //xbj step in the generator, 
		Double_t fT_Min; //-t low range, GeV2
		Double_t fT_Max; //-t low range, GeV2
		Double_t fT_Step; //-t step in the generator, GeV2
		Double_t fPhi_Min; //phitlow range, GeV
		Double_t fPhi_Max; //phi low range, GeV
		Double_t fPhi_Step; //phi step in the generator, GeV
        Int_t fQ2_Bin;
        Int_t fXb_Bin;
        Int_t fT_Bin;
        Int_t fPhi_Bin;

		Double_t fSigmaPP;//nb/GeV, beam+, target+
		Double_t fSigmaPM;//nb/GeV, beam+, target-
		Double_t fSigmaMP;//nb/GeV, beam-, target+
		Double_t fSigmaMM;//nb/GeV, beam-, target-
		Double_t fSigma;//nb/GeV, average of the four cross section values

		Double_t fSigmaBS;//beam spin-dependent xs
		Double_t fSigmaTS;//target spin-dependent xs
		Double_t fSigmaDS;//beam+target double spin-dependent xs

		TString fTablesDir; //The directory to store the subdirectories of tables
		TString fTarget; // target name, Neutron or Proton
		/*variables}}}*/

};
#endif
