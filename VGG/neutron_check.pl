#!/usr/bin/perl
$WORK_DIR = "/work/halla/solid/yez/dvcs/VGG";
$VGG_CMD = "$WORK_DIR/env_ld.sh";
$Energy = $ARGV[0];
$TP = $ARGV[1];
$Mn =0.939565378;
#$Energy = 11.;
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

$TargetPol = 1; #1->Tx, 2->Ty, 3->Long
#$TP = "Tx";
$DATA_DIR = "$WORK_DIR/Neutron_$E0";
#########################################################################3
system("echo $WORK_DIR");

$filename = "$WORK_DIR/check/E$E0-$TP.txt";
open(ENV_FILE,"> $filename");#_$ARGV[0]_$ARGV[1]_$ARGV[2]");

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
		$tmin = ($Q2*$Mn+2.*$Mn*$nu*($nu-sqrt($nu*$nu+$Q2))) / (sqrt($nu*$nu+$Q2)-$nu-$Mn);
		$tmin = $tmin+0.02*$tmin;
		$tmin = int(1000*$tmin)/1000;
		$t=-1.*$tmin;
		$it=0;
		while ($t<=$t_Max){
			$datafilename = "$DATA_DIR/E$E0-$TP-Q2-$iQ2-xb-$ixb-t-$it.dat";
			unless (-e $datafilename) {
					$datafilename = "E$E0-$TP-Q2-$iQ2-xb-$ixb-t-$it.dat";
					print ENV_FILE "$datafilename\n";#DDVCS
				}
			$t = $t+$t_Step;
			$t =int(1000*$t)/1000;
			$it = $it + 1;
		}
	}
#}
}
close(ENV_FILE);
system("echo created file: $filename");
#################

