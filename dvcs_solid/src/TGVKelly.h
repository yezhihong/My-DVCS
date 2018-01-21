/*---------------------------------------- TGVKelly.h ---------------------------------------*
 |                                                                                      |
 | Author : H. Moutarde (CEA-Saclay, IRFU/SPhN).                                        |
 | v1.3, July, 25th 2008.                                                               |
 |                                                                                      |
 | Computation of the differential cross sections of the Bethe-Heitler and Deeply       |
 | Virtual Compton Scattering processes (and also of the interference of the two        |
 | processes).                                                                          |
 |                                                                                      |
 | The formulas are derived from a Mathematica package written by P. Guichon            |
 | (CEA-Saclay, IRFU/SPhN) and M. Vanderhaegen (W&M College), see arXiv:0808.????       |
 |                                                                                      | 
 | This library is not standalone, since some Root librairies are necessary.            |
 |                                                                                      |
 | All expressions are written in a class named TGVKelly.                                    | 
 |                                                                                      |
 *--------------------------------------------------------------------------------------*/ 
 
 
 
/*------------------------------ Compilation of the class ------------------------------*/ 

#ifndef __TGVKelly__ // Checks that the class is compiled only once
#define __TGVKelly__



/*----------------------------- Necessary Root Libraries -------------------------------*/

#include "TObject.h"
#include "TLorentzVector.h"
#include "TMath.h"



/*------------------------------ Definition of the class TGVKelly ---------------------------*/

class TGVKelly
{
	private :


/*---------------------------------- Tests variables -----------------------------------*/

		// Printouts

		Bool_t Validation;	
		Bool_t NoPrint;	
	
	
		// Initialisation flags
	
		Bool_t InitBeamEnergy;	
		Bool_t InitKinematics;
		Bool_t InitMvcs;
		Bool_t InitExactBHCrossSections;
		Bool_t InitExactVCSAndInterfCrossSections;
		Bool_t InitLeadingBHCrossSections;
		Bool_t InitLeadingVCSAndInterfCrossSections;
			
	

/*---------------------------------------- Kinematics ----------------------------------*/

		// Frame dependent scalars	
	
    	Double_t ELab; // Beam energy
    	Double_t thetag; // Angle between real and virtual photons
    	Double_t qpPerp; // Component (here x-axis) of the real photon 3-momentum orthogonal 
    					// to the virtual photon trajectory (here z-axis) in the 
    					// hadronic plane (here xz-plane)
		Double_t Ur[100]; // Coefficients of the expansion of the interference cross section
						// wrt (combinations of) helicity amplitudes
		Double_t Omega; // (Fonction of) the linear polarization rate
    
        
    	// Invariant scalars
    	
		Double_t xB; // Bjorken variable
		Double_t Q2; // Virtuality of the photon
		Double_t t; // Mandelstam variable (square of the 4-momentum transfer)
		Double_t s; // Mandelstam variable (square of the total incoming 4-momentum)
		Double_t Q;  // Photon virtual mass i.e. square root of Q2   
    	Double_t PhaseSpace; // Differential element of cross section
    	Double_t M; // Proton mass
//    	const static Double_t M = 0.93827231; // Proton mass
    
       
    	// 4-vectors defined in the CM frame :
    
    	TLorentzVector qCM; // Virtual photon (propagates along z-axis)
    	TLorentzVector pCM; // Incoming proton (propagates along z-axis)
    	TLorentzVector qpCM; // Real photon (defines hadronic plane xz)
    	TLorentzVector ppCM; // Outgoing proton



/*------------------------ (Combinations of) helicity amplitudes -----------------------*/

		Double_t Jem[4][3]; // Helicity amplitudes of the interference process assuming the 
								// real photon has helicity +1.
		Double_t RMvcs[4][3]; //  Real part of the helicity amplitudes of the VCS process 
								//assuming the real photon has helicity +1.
		Double_t IMvcs[4][3]; //  Imaginary part of the helicity amplitudes of the VCS process 
								// assuming the real photon has helicity +1.
	
	
	
/*------------------ Expansion of cross sections for harmonic analysis -----------------*/

		// Bethe Heitler process
	
		Double_t SigmaBHPol0[4]; // coefficients for the unpolarized cross section
		Double_t SigmaBHPolX[2]; // coefficients for the x-polarized cross section 
		Double_t SigmaBHPolY; // coefficient for the y-polarized cross section
		Double_t SigmaBHPolZ[2]; // coefficients for the z-polarized cross section
	
	
		// Virtual Compton Scattering process
	
		Double_t SigmaVCSPol0[5]; // coefficients for the unpolarized cross section
		Double_t SigmaVCSPolX[4]; // coefficients for the x-polarized cross section	
		Double_t SigmaVCSPolY[5]; // coefficients for the y-polarized cross section	
		Double_t SigmaVCSPolZ[4]; // coefficients for the z-polarized cross section
	
	
		// Interference
	
		Double_t SigmaIPol0[8]; // coefficients for the unpolarized cross section	
		Double_t SigmaIPolX[8]; // coefficients for the x-polarized cross section	
		Double_t SigmaIPolY[8]; // coefficients for the x-polarized cross section
		Double_t SigmaIPolZ[8]; // coefficients for the x-polarized cross section
	
	
	
/*-------------------------------------- Methods ---------------------------------------*/	

		// Initialisations
	
		void SetBeamEnergy( Double_t EBeam );
				// Sets beam energy
	
		void MakeKinematics( Double_t Q2Input, Double_t xBInput, Double_t tInput );	
				// Sets all the kinematics except the Ur's
				
		void MakeVCSHelicityAmplitudes( Double_t ReH, Double_t ImH, Double_t ReE, Double_t ImE, Double_t ReHT, Double_t ImHT, Double_t ReET, Double_t ImET );
				// Fills the RMvcs and IMvcs arrays
	  
		void MakeExactBHCrossSections();
				// Fills the Jem and SigmaBHPol. arrays
	  
		void MakeExactVCSAndInterfCrossSections();
				// Fills the Ur, SigmaVCSPol. and SigmaIPol. arrays
	  
		void MakeLeadingBHCrossSections();
				// Fills the Jem and SigmaBHPol. arrays at leading order
				// in the 1/Q expansion.

	  
		void MakeLeadingVCSAndInterfCrossSections();
				// Fills the Ur, SigmaVCSPol. and SigmaIPol. arrays at leading order
				// in the 1/Q expansion.
				
	
		// Misc.
	
		Double_t DdirectDcrossed( Double_t phi );
				// Denominator of the Bethe Heitler cross section

		// Form Factors

		Double_t KellyE(Double_t q2);
		Double_t KellyM(Double_t q2);


	
	public :
	
		// Cross sections
	
		Double_t SqrAmplBH( Double_t Q2Input, Double_t xBInput, Double_t tInput, Double_t phi, Int_t BeamHeli, Int_t TargetPolar, Bool_t Exact ); 
			// Bethe Heitler cross section divided by phase space
		
		Double_t SqrAmplVCS( Double_t Q2Input, Double_t xBInput, Double_t tInput, Double_t phi, Int_t BeamHeli, Int_t TargetPolar, Double_t ReH, Double_t ImH, Double_t ReE, Double_t ImE, Double_t ReHT, Double_t ImHT, Double_t ReET, Double_t ImET, Bool_t Exact );
			// Virtual Compton Scattering cross section divided by phase space	
			
		Double_t SqrAmplInterf( Double_t Q2Input, Double_t xBInput, Double_t tInput, Double_t phi, Int_t BeamHeli, Int_t TargetPolar, Int_t BeamCharge, Double_t ReH, Double_t ImH, Double_t ReE, Double_t ImE, Double_t ReHT, Double_t ImHT, Double_t ReET, Double_t ImET, Bool_t Exact );
			// Interference cross section divided by phase space
	
		Double_t CrossSectionBH( Double_t Q2Input, Double_t xBInput, Double_t tInput, Double_t phi, Int_t BeamHeli, Int_t TargetPolar, Bool_t Exact ); 
			// Bethe Heitler cross section
		
		Double_t CrossSectionVCS( Double_t Q2Input, Double_t xBInput, Double_t tInput, Double_t phi, Int_t BeamHeli, Int_t TargetPolar, Double_t ReH, Double_t ImH, Double_t ReE, Double_t ImE, Double_t ReHT, Double_t ImHT, Double_t ReET, Double_t ImET, Bool_t Exact );
			// Virtual Compton Scattering cross section	
			
		Double_t CrossSectionInterf( Double_t Q2Input, Double_t xBInput, Double_t tInput, Double_t phi, Int_t BeamHeli, Int_t TargetPolar, Int_t BeamCharge, Double_t ReH, Double_t ImH, Double_t ReE, Double_t ImE, Double_t ReHT, Double_t ImHT, Double_t ReET, Double_t ImET, Bool_t Exact );
			// Interference cross section
	
	
		// Constructor, destructor
	
		TGVKelly(); // default constructor is mandatory for Root
		TGVKelly( Double_t EBeam, Bool_t Valid, Bool_t NoPrint );
		~TGVKelly();

}; // end of the TGVKelly class definition 
	
#endif
