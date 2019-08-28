#!/usr/bin/perl -w
# Nov. 2008, by Eshel Faraggi
# 
# Send error messages to the user, not system log
   open(STDERR,'<&STDOUT'); $| = 1;   

   $outdir = "spXout";
   if (($#ARGV >= 1)&&($ARGV[0] ne "-h"))
     {
     $flist = $ARGV[0];
     $prfdir=$ARGV[1];
     if ($#ARGV >= 2)
       {
       $outdir = $ARGV[2];
       }
     }
   else
     {
     print "
     1) Use \"$0 protein_list_file input_directory [output_dir.]\" to run program on list_file 
     2) protein_list_file should contain protein id's
     3) input_directory should contain either the psiblast profiles (*.mat ending) or fasta sequences (no ending)
     4) You must set the environment variable spineXcodir pointing to the location of the spineX code\ne.g., export spineXcodir=/path/\n\nDepending on your local configuration you may also need to set spineXblast\nto point to the blast root directory\n
     5) output_dir is optional dir. name for output (default is spXout)\n\n";
     exit;
     }

   $t1 =  "~/servers/spineX/";
   if (-d $t1)
     {
     $spcodir = $t1;
     }
   else
     {
     $t1 = `readlink -f $0`;
     chomp($t1 = `dirname "$t1"`);
     if (-e "$t1/bin/phipsi_phipsi0.e")
       {
       $spcodir = "$t1/";
       }
     else
       {
       $spcodir = $ENV{'spineXcodir'};
       }
     }
   if ($spcodir)
     {
     do
       {
       $rndm = rand();
       $irnd = int($rndm*1000000);
       $pwdir = `pwd`;
       chomp($pwdir);
       if (-d "_subspXbiglist.")
         {
         $workdir = $pwdir . "/_subspXbiglist./spxtemp$irnd/";
         }
       else
         {
         $workdir = $pwdir . "/spxtemp$irnd/";
         }
       }
     while (-d $workdir);
     mkdir $workdir or die "Can't make temp. dir. ($workdir) aborting\n";
     $codir = "$spcodir/bin/"
     }
   else
     {
     die "ABORTING:\nYou must set the environment variable spineXcodir\npointing to the location of the spineX perl code\ne.g., export spineXcodir=/path/\n";
     }

   print "Program index (in temp directory): $irnd\n";

main: {

   open(FL,$flist) || die "Could not open $flist, aborting";
   @pdat = <FL>;
   close(FL);
   for (@pdat)
     {
     chomp;
     $infl = "$prfdir/$_.mat";
     if (! -e $infl)
       {
       $infl = "$prfdir/$_";
       if (-e $infl)
         {

   $blastdir = $ENV{'spineXblast'};
   if ($blastdir)
     {
     print "You have set spineXblast\n";
     if (!(-e "$blastdir/bin/blastpgp"))
       {
       die "Aborting: can't find blastpgp in $blastdir/bin\nset spineXblast environment variable to blast directory\n";
       }
     }
   else
     {
     $blastdir="~/blast/";
     if (! -e "$blastdir/bin/blastpgp")
       {
       die "Aborting: can't find bin/blastpgp in $blastdir \nset spineXblast environment variable to its location\n";
       }
     }
         $tmprf = "$prfdir/$_.mat";
         print "Note: Using spineXblast=$blastdir for blastpgp\n";
         system("$blastdir/bin/blastpgp -d ~/Downloads/uniref50.db -j 3 -i $infl -Q $tmprf -a 4 > $workdir/_tmp2.$irnd")==0 or die "Aborting(spineX): Can't do psiblast\n";
         print "Psiblast profile for $_ saved in $prfdir\n";
         }
       else
         {
         die "Can't find profile ($_.mat) or fasta ($_) file for $_ in $prfdir, aborting";
         }
       }
     }

   print "Please wait, calculating prediction ...\n";
   `mkdir -p $outdir`;
   die "Can't find $outdir, aborting" if (! -d "$outdir");
   system("$codir/phipsi_ss0.e $irnd $spcodir $workdir $flist $prfdir")==0 or die "Aborting: Can't run first SS predictor\n";
   system("$codir/phipsi_rsa.e $irnd $spcodir $workdir $flist")==0 or die "Aborting: Can't run rsa predictor\n";
   system("$codir/phipsi_phipsi0.e $irnd $spcodir $workdir $flist")==0 or die "Aborting: Can't run first phipsi predictor\n";
   system("$codir/phipsi_ss1.e $irnd $spcodir $workdir $flist")==0 or die "Aborting: Can't run second SS predictor\n";
   system("$codir/phipsi_phipsi1.e $irnd $spcodir $workdir $flist $outdir")==0 or die "Aborting: Can't run second phipsi predictor\n";
   sleep(1);
   system("rm $workdir/_nres.out  $workdir/out_*")==0 or print "Couldn't remove content of temp. dir. $workdir\n";
   system("rmdir $workdir")==0 or print "Couldn't remove temp. dir. $workdir\n";
   print "Predictions in $outdir dir.\n\n";
  }
__END__
