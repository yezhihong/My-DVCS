#!/usr/bin/perl
#
$pgm      = 'aao_rad';
$infile   = 'aao-pi0-1.6.inp';
$outfile  = '1.6-aao-5000evts';
$OS       = $ENV{"OS_NAME"};
$TOP_DIR  = $ENV{"TOP_DIR"};
$USER     = $ENV{"USER"};
$workdir = "/work/clas/disk5/$USER";
chdir("$workdir");
$cmd = "$TOP_DIR/bin/$OS/$pgm < $infile";
print  "cmd: $cmd\n";
system $cmd;
rename("$pgm.evt","$outfile.evt") ||
    warn "cannot rename aaoradgen.evt to $outfile.evt: $!";
rename("$pgm.sum","$outfile.sum") ||
    warn "cannot rename aaoradgen.sum to $outfile.sum: $!";
rename("$pgm.rz","$outfile.rz") ||
    warn "cannot rename aaoradgen.rz to $outfile.rz: $!";
rename("$pgm.out","$outfile.out") ||
    warn "cannot rename aaoradgen.out to $outfile.out: $!";
