#!/usr/bin/perl
$WORK_DIR = "/work/halla/solid/yez/dvcs/VGG";
$VGG_CMD = "$WORK_DIR/env_ld.sh";
$run_number = $ARGV[0];
$split = $ARGV[1];
$Mn =0.939565378;
$Energy = 8.8;
$E0 = int($Energy);
$Q2_Min = 1;
$Q2_Max = 9;
$Q2_Step = 1.0;
$Q2_Bin = ($Q2_Max - $Q2_Min)/$Q2_Step;
$xb_Min = 0.05;
$xb_Max = 0.75;
$xb_Step = 0.05;
$xb_Bin = ($xb_Max - $xb_Min)/$xb_Step;
#$t_Min = 0.1;#no in use here, tmin will be recalculated
$t_Max = 2.5;
$t_Step = 0.1;
#$t_Bin = ($t_Max - $t_Min)/$t_Step;
$phi_Min = 0.0;
$phi_Max = 360.;
$phi_Step = 15.;
#$phi_Bin = ($phi_Max - $phi_Min)/$phi_Step;

$TargetPol = 2; #1->Tx, 2->Ty, 3->Long
$TP = "Ty";
$DATA_DIR = "$WORK_DIR/Neutron_$E0";
#########################################################################3
system("echo $WORK_DIR");

for ( $iQ2=0;$iQ2<$Q2_Bin;$iQ2++){#from 1.0 upto 10.0
	$Q2 = $Q2_Min+$iQ2 * $Q2_Step;
	for ( $ixb=0;$ixb<$xb_Bin;$ixb++){#from 0.1 upto 0.9
		$xb = $xb_Min+$ixb * $xb_Step;

		$Kp = $Energy -($Q2/(2.*$Mn*$xb));
		$nu = $Energy-$Kp;
        
		$tmax = ($Q2*$Mn+2.*$Mn*$nu*($nu+sqrt($nu*$nu+$Q2))) / (-sqrt($nu*$nu+$Q2)-$nu-$Mn);
		$tmax = -1*$tmax;
		if ($tmax < $t_Max) {
			$t_Max = $tmax;
		}

		##This definition is the same as the other one but I prefer to use that one
		#$tmin = -$Q2-($Q2*(1.-$xb)*($nu-sqrt($Q2+$nu*$nu))/($xb*($Mn+$nu-sqrt($Q2+$nu*$nu))));
	
		$tmin = ($Q2*$Mn+2.*$Mn*$nu*($nu-sqrt($nu*$nu+$Q2))) / (sqrt($nu*$nu+$Q2)-$nu-$Mn);
		$tmin = $tmin+0.02*$tmin;
		$tmin = int(1000*$tmin)/1000;

#		system("echo max=$t_Max");
#		system("echo min=$tmin");

		$t=-1.*$tmin;
#		for ($it=0;$it<20;$it++){#upto <62 before
        $it=0;
	    while ($t<=$t_Max){

			#$filename = "$WORK_DIR/inputs_np/E$E0-$TP-Q2-$iQ2-xb-$ixb-t-$it";
			#unless (-e $filename) {
			$datafilename = "$DATA_DIR/E$E0-$TP-Q2-$iQ2-xb-$ixb-t-$it.dat";
			unless (-e $datafilename) {
				#unless (0) {
				open(ENV_FILE,"> $filename");#_$ARGV[0]_$ARGV[1]_$ARGV[2]");

#Choose one of the following menu items :
# 1 : 5-fold DIFFERENTIAL CROSS SECTION for ep -> ep gamma as function of t
# 2 : e+e- asymmetry for DVCS
# 3 : electron single spin asymmetry (SSA)
# 4 : 2-body DOUBLY POLARIZED cross sections for (D)DVCS
#     polarized electron, polarized target
# 5 : DIFFERENTIAL CROSS SECTION for DVCS as function of Q2
# 6 : DIFFERENTIAL CROSS SECTION for DVCS as function of t
# 11 : parton distributions
# 12 : GPDs
# 13 : sum rules and positivity tests
# 14 : form factors from GPDs and WACS
# 15 : accuracy test of integrals for twist-2 amplitude
# 16 : accuracy test of integrals for twist-3 amplitude
# 17 : integration of BH over bin size
# 18 : fit of kinematics
# 20 : 3-body DOUBLY POLARIZED cross sections for DDVCS
#     polarized electron, polarized target
# 21 : 2-body DIFFERENTIAL CROSS SECTION for Rho as function of Q2 or theta_lab
# 22 : Timelike Compton scattering
# 30 : WACS
# 40 : Fit of DVCS data
# 50 : Transverse densities
# 60 : Resolutions
				print ENV_FILE "4\n";#DDVCS

#Make a choice for the MECHANISM :
# 1 : Bethe-Heitler contribution
# 2 : DVCS contribution
# 3 : Bethe-Heitler + DVCS contribution
# 4 : rho contribution
# 5 : DVCS + rho contribution
# 6 : Bethe-Heitler + DVCS + rho contribution
# 7 : Bethe-Heitler + rho contribution
				print ENV_FILE "3\n";

#Choice to calculate on proton or on neutron (default = 1)
# 1 = proton
				print ENV_FILE "2\n";#neutron

#Choose model for GPD
# 35 : xi dependent parametrization with MRST02 NNLO distribution at mu^2 = 1 GeV^2
				print ENV_FILE "35\n";

#Give the value for the power b in the profile function for the valence contribution to H (e.g. 1.)
				print ENV_FILE "1\n";

#Give the value for the power b in the profile function for the sea contribution to H (e.g. 1.)
				print ENV_FILE "1\n";# b input 6

#Choose the model for the t-dependence of the GPD H
#1 = Factorized model for the t-dependence
#2 = Regge inspired ansatz for the t-dependence
#3 = Experimental exponential input by user
#4 = Experimental exponential fit (fit of all data)
#5 = Experimental exponential fit (fit of subrange -by Mick-)
#6 = Hybrid model (FF vor val and Regge for sea)
#7 = Diehl et al. model (t-dep in DDs)
#8 = R2 Regge ansatz model (t-dep in DDs)
#9 = Diehl et al. model (t-dep out of DDs)
#10 = R2 Regge ansatz model (t-dep out of DDs)
				print ENV_FILE "8\n";

#Enter slope alphap (GeV-2)
				print ENV_FILE "1.098\n";

#Do you want to evaluate the D-term contribution to the GPD H?
#1 = Yes
#2 = No
				print ENV_FILE "2\n";

#Do you want to evaluate the GPD E?
#1 = only D-term contribution
#2 = double distribution contribution + D-term contribution
#3 = No, E = 0
				print ENV_FILE "2\n";

#Give the model for the double distribution part of the GPD E
# 1 = only valence quark contribution (factorized t-dependence)
# 2 = valence quark + VM contribution (factorized t-dependence)
# 3 = only VM contribution (factorized t-dependence)
# 4 = Reggeized Double Distribution based on e(x) = (1 - x)^eta q_v(x)
# (only valence contribution)
# 6 = Reggeized Double Distribution based on e(x) = (1 - x)^eta q_v(x)
# (valence + VM contribution)
				print ENV_FILE "2\n";

#Give the value of Ju (e.g. 0.3)
				print ENV_FILE "0.3\n";#Ju

#Give the value of Jd (e.g. 0.1)
				print ENV_FILE "0.0\n";

#Do you want to evaluate the pi0 pole contribution (i.e. SPD Etilde)?
#1 = Yes
#2 = No
				print ENV_FILE "1\n";

#Do you want to include twist-3 corrections ?
#1: Do not include twist-3 corrections for L photon,
#   but include the minimal corrections to restore gauge invariance for T photon
#2: Include twist-3 corrections for L photon in Wandzura-Wilczek approximation
				print ENV_FILE "2\n"; #twist 3

# With (1) or without (2) Htilde ?
				print ENV_FILE "1\n";# Htilde

#Give the polarization of the target proton
# 1 : proton polarized along x-axis
# 2 : proton polarized along y-axis (perpendicular to lepton plane)
# 3 : proton polarized along z-axis (along the virtual photon direction)
				print ENV_FILE "$TargetPol\n"; # zpol target

#Calculation for what LEPTON charge ?
# 1 : negatively charged lepton (JLab)
# 2 : positively charged lepton (HERMES, COMPASS)
				print ENV_FILE "1\n"; # e- beam, jlab

#Give the value of beam energy in GeV (e.g. 27.)
				print ENV_FILE "$Energy\n"; # Energy

# 1 : As a function of t
# 2 : As a function of xB
# 3 : As a function of Phi
# 4 : As a function of Q2
				print ENV_FILE "3\n"; # thetat dist

#Give the value of Q^2 in GeV^2 (e.g. 5.0)
# WARNING !!! If >20, thetae, not Q2
				print ENV_FILE "$Q2\n"; #Q2

#Give the value of x_B (e.g. 0.3)
				print ENV_FILE "$xb\n"; #xb

#Give the value of Qprime^2 in GeV^2 (e.g. 2.0)
				print ENV_FILE "0.0\n"; #Q'

#Give the value of -t (in GeV^2)
				print ENV_FILE "$t\n"; #t
#cos(theta) = 0.948207, ga_en = 2.496229 GeV, ga_mom = 2.496229 GeV, theta_lab = 18.521140 deg

#Give the first value for the angle phi (in deg) to calculate
				print ENV_FILE "$phi_Min\n";

#Give the step in the angle phi (in deg) (e.g. 10.)
				print ENV_FILE "$phi_Step\n";# delta phi
#Give the last value in the angle phi (in deg) (e.g. 180.)
				print ENV_FILE "$phi_Max\n";#final phi / inv

				close(ENV_FILE);
				system("echo created file: $filename");
#########################################################################3

				$script_name = "$WORK_DIR/inputs_np/jsub_E$E0-$TP-Q2-$iQ2-xb-$ixb-t-$it";#_$ARGV[0]_$ARGV[1]_$ARGV[2]";

#			$mail_to     = $ENV{'USER'};
				print "Openning file $script_name\n";
#			print "Running Job as user : $mail_to\n";

				open(OUT_FILE,">$script_name");
				print OUT_FILE "JOBNAME: VGG_E$E0-$TP-Q2-$iQ2-xb-$ixb-t-$it\n";
				print OUT_FILE "PROJECT: solid\n";
#			print OUT_FILE "MAIL: $mail_to\@jlab.org\n";
				print OUT_FILE "COMMAND: cat $WORK_DIR/inputs_np/E$E0-$TP-Q2-$iQ2-xb-$ixb-t-$it | $VGG_CMD\n";
				print OUT_FILE "OUTPUT_DATA: dvcs.out\n";#_$ARGV[0]_$ARGV[1]_$ARGV[2].dat\n";
				print OUT_FILE "OUTPUT_TEMPLATE:$DATA_DIR/E$E0-$TP-Q2-$iQ2-xb-$ixb-t-$it.dat\n";#_$ARGV[0]_$ARGV[1]_$ARGV[2].dat\n";
				print OUT_FILE "MEMORY: 256 MB\n";
				print OUT_FILE "TRACK: simulation\n";
				print OUT_FILE "OS: centos65\n";
				print OUT_FILE "TIME: 4320\n";

				close(OUT_FILE);

				system("echo jsub $script_name");
				system("/site/bin/jsub $script_name");
			}

			$t = $t+$t_Step;
			$t =int(1000*$t)/1000;
			$it = $it + 1;
		}
	}
}

