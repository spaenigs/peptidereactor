#!/usr/bin/perl -w

#get perl executive path
$perlPath=`which perl`;
chomp($perlPath);
die "No perl found in searching folders.\n" if ($perlPath=~/no perl/);
#die "./setup.pl [ -cpu [CPU number = 1 (default)] ] " if @ARGV<1;
#%params=@ARGV;
$cnfhome=`pwd`;chomp($cnfhome);


# genPSP.sh -> BLAST
open fhRUN,"<".$cnfhome."/util/genPSP.sh";
open fhRUNNEW,">".$cnfhome."/util/genPSP.sh.new";
while(<fhRUN>){
    if(/^RaptorX_HOME=/){
	print fhRUNNEW "RaptorX_HOME=".$cnfhome."\n";
	next;
    }
    print fhRUNNEW "$_";
}
close fhRUN;
close fhRUNNEW;
$cmd="mv ".$cnfhome."/util/genPSP.sh.new ".$cnfhome."/util/genPSP.sh";
`$cmd`;
chmod(0755, $cnfhome."/util/genPSP.sh");


# genSS3.sh -> PSIPRED
open fhRUN,"<".$cnfhome."/util/genSS3.sh";
open fhRUNNEW,">".$cnfhome."/util/genSS3.sh.new";
while(<fhRUN>){
    if(/^RaptorX_HOME=/){
	print fhRUNNEW "RaptorX_HOME=".$cnfhome."\n";
        next;
    }
    print fhRUNNEW "$_";
}
close fhRUN;
close fhRUNNEW;
$cmd="mv ".$cnfhome."/util/genSS3.sh.new ".$cnfhome."/util/genSS3.sh";
`$cmd`;
chmod(0755, $cnfhome."/util/genSS3.sh");


# genDIS.sh -> DISOPRED
open fhRUN,"<".$cnfhome."/util/genDIS.sh";
open fhRUNNEW,">".$cnfhome."/util/genDIS.sh.new";
while(<fhRUN>){
	if(/^RaptorX_HOME=/){
		print fhRUNNEW "RaptorX_HOME=".$cnfhome."\n";
		next;
	}
	print fhRUNNEW "$_";
}
close fhRUN;
close fhRUNNEW;
$cmd="mv ".$cnfhome."/util/genDIS.sh.new ".$cnfhome."/util/genDIS.sh";
`$cmd`;
chmod(0755, $cnfhome."/util/genDIS.sh");


# genMTX.sh
open fhRUN,"<".$cnfhome."/util/genMTX.sh";
open fhRUNNEW,">".$cnfhome."/util/genMTX.sh.new";
while(<fhRUN>){
    if(/^RaptorX_HOME=/){
	print fhRUNNEW "RaptorX_HOME=".$cnfhome."\n";
        next;
    }
    print fhRUNNEW "$_";
}
close fhRUN;
close fhRUNNEW;
$cmd="mv ".$cnfhome."/util/genMTX.sh.new ".$cnfhome."/util/genMTX.sh";
`$cmd`;
chmod(0755, $cnfhome."/util/genMTX.sh");


# genHMM.sh
open fhRUN,"<".$cnfhome."/util/genHMM.sh";
open fhRUNNEW,">".$cnfhome."/util/genHMM.sh.new";
while(<fhRUN>){
    if(/^RaptorX_HOME=/){
	print fhRUNNEW "RaptorX_HOME=".$cnfhome."\n";
        next;
    }
    print fhRUNNEW "$_";
}
close fhRUN;
close fhRUNNEW;
$cmd="mv ".$cnfhome."/util/genHMM.sh.new ".$cnfhome."/util/genHMM.sh";
`$cmd`;
chmod(0755, $cnfhome."/util/genHMM.sh");


# genSS8.sh
open fhRUN,"<".$cnfhome."/util/genSS8.sh";
open fhRUNNEW,">".$cnfhome."/util/genSS8.sh.new";
while(<fhRUN>){
    if(/^RaptorX_HOME=/){
        print fhRUNNEW "RaptorX_HOME=".$cnfhome."\n";
	next;
    }
    print fhRUNNEW "$_";
}
close fhRUN;
close fhRUNNEW;
$cmd="mv ".$cnfhome."/util/genSS8.sh.new ".$cnfhome."/util/genSS8.sh";
`$cmd`;
chmod(0755, $cnfhome."/util/genSS8.sh");


# genACC.sh
open fhRUN,"<".$cnfhome."/util/genACC.sh";
open fhRUNNEW,">".$cnfhome."/util/genACC.sh.new";
while(<fhRUN>){
    if(/^RaptorX_HOME=/){
        print fhRUNNEW "RaptorX_HOME=".$cnfhome."\n";
        next;
    }
    print fhRUNNEW "$_";
}
close fhRUN;
close fhRUNNEW;
$cmd="mv ".$cnfhome."/util/genACC.sh.new ".$cnfhome."/util/genACC.sh";
`$cmd`;
chmod(0755, $cnfhome."/util/genACC.sh");



# ======== now we run into directory =========
# HMMpred -> buildali2.pl
open fhRUN,"<".$cnfhome."/util/HHpred/buildali2.pl";
open fhRUNNEW,">".$cnfhome."/util/HHpred/buildali2.pl.new";
while(<fhRUN>){
    if(/^my \$raptorx=/){
	print fhRUNNEW "my \$raptorx=\"".$cnfhome."\";\n";
        next;
    }
    print fhRUNNEW "$_";
}
close fhRUN;
close fhRUNNEW;
$cmd="mv ".$cnfhome."/util/HHpred/buildali2.pl.new ".$cnfhome."/util/HHpred/buildali2.pl";
`$cmd`;
chmod(0755, $cnfhome."/util/HHpred/buildali2.pl");


# HMMpred -> alignhits.pl
open fhRUN,"<".$cnfhome."/util/HHpred/alignhits.pl";
open fhRUNNEW,">".$cnfhome."/util/HHpred/alignhits.pl.new";
while(<fhRUN>){
    if(/^my \$raptorx=/){
        print fhRUNNEW "my \$raptorx=\"".$cnfhome."\";\n";
        next;
    }
    print fhRUNNEW "$_";
}
close fhRUN;
close fhRUNNEW;
$cmd="mv ".$cnfhome."/util/HHpred/alignhits.pl.new ".$cnfhome."/util/HHpred/alignhits.pl";
`$cmd`;
chmod(0755, $cnfhome."/util/HHpred/alignhits.pl");



#SS8_Predict -> run_raptorx-ss8.pl
open fhRUN,"<".$cnfhome."/util/SS8_Predict/bin/run_raptorx-ss8.pl";
open fhRUNNEW,">".$cnfhome."/util/SS8_Predict/bin/run_raptorx-ss8.pl.new";
while(<fhRUN>){
    if(/^\$raptorx=/){
        print fhRUNNEW "\$raptorx=\"".$cnfhome."\";\n";
        next;
    }
    print fhRUNNEW "$_";
}
close fhRUN;
close fhRUNNEW;
$cmd="mv ".$cnfhome."/util/SS8_Predict/bin/run_raptorx-ss8.pl.new ".$cnfhome."/util/SS8_Predict/bin/run_raptorx-ss8.pl";
`$cmd`;
chmod(0755, $cnfhome."/util/SS8_Predict/bin/run_raptorx-ss8.pl");

