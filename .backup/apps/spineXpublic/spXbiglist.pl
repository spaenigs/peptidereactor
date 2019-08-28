#!/usr/bin/perl -w
# Nov. 2008, by Eshel Faraggi
# 
# Send error messages to the user, not system log
# nohup time ./spXbiglist.pl list.23678 profiles/ 8 2.5 &
   open(STDERR,'<&STDOUT'); $| = 1;   

   $outdir = "spXout";
   if ($#ARGV >= 3)
     {
     $flist = $ARGV[0];
     $prfdir=$ARGV[1];
     $tncpu=$ARGV[2];
     $runtime=$ARGV[3];
     if ($#ARGV >= 4)
       {
       $outdir = $ARGV[4];
       }
     }
   else
     {
     print "
     1) Use \"$0 protein_list_file input_directory number_cpu run_time [output_dir.]\" to run program on list_file 
     2) protein_list_file should contain protein id's
     3) input_directory should contain either the psiblast profiles (*.mat ending) or fasta sequences (*.fasta ending)
     4) You must set the environment variable spineXcodir pointing to the location of the spineX code\ne.g., export spineXcodir=/path/\n\nDepending on your local configuration you may also need to set spineXblast\nto point to the blast root directory\n
     5) output_dir is optional dir. name for output (default is spXout)
     6) number_cpu : number of separate processes to start
     7) run_time : amount of run time (hours, aspen estimate). May not finish all protein list if too small. \n\n";
     exit;
     }

main: {

  $ncpu = 0;
  open(FL,"$flist") || die "Could not open $flist, aborting";
  chomp($bflist = `basename $flist`);
  @pdat2 = <FL>;
  close(FL);
  $m1 = 3600 * $runtime / 0.008;
  if ($m1 < 1999999)
    {
    $maxj = $m1;
    }
  else
    {
    $maxj = 1999999;
    }
  $jt = 0;
  for (@pdat2)
    {
    chomp ($prftmp = "$prfdir/$pdat2[0]");
    $jtt = `cat $prftmp.mat | wc -l`;
    $jt += $jtt - 9;
    }
  $m1 = $jt * 1.0 / $tncpu;
  if ($m1 < $maxj)
    {
    $maxj = $m1 + 1;
    }
  `mkdir -p _subspXbiglist.`;
  chomp ($xbig = `readlink -f _subspXbiglist.`);
  do 
   {
   $i=0; $jt=0; @pdat = ();
   while (($i<7999)&&($jt<$maxj)&&($pdat2[0]))
     {
     chomp ($prftmp = "$prfdir/$pdat2[0]");
     $jtt = `cat $prftmp.mat | wc -l`;
     $jt += $jtt - 9;
     if (($jt<$maxj)&&($pdat2[0]))
       {
       $pdat[$i++] = $pdat2[0];
       shift(@pdat2);
       $j = $jt;
       }
     }
   $premain = $#pdat2 + 1;
   $tmplist = "$xbig/$bflist.spXbiglist.$premain";
   open(FL,">$tmplist") || die "Could not write to $tmplist, aborting";
   print "Remaining $premain proteins $i residues $j file $tmplist\n";
   for (@pdat)
     {
     if ($_)
       {
       print FL "$_";
       }
     }
   close(FL);
   system("nohup time /data3/efaraggi/server/spX/spineX/spX.pl $tmplist $prfdir $outdir &");
   $ncpu++;
   } while (($premain > 0)&&($ncpu<$tncpu));
   
   if ($premain > 0)
     {
     $tmplist = "$xbig/$bflist.spXbiglist.remain";
     open(FL,">$tmplist") || die "Could not write to $tmplist, aborting";
     for (@pdat2)
       {
       if ($_)
         {
         print FL "$_";
         }
       }
     close(FL);
     print "Did not finish list, remaining in $xbig/$bflist.spXbiglist_remain\n";
     }
     
   print "Predictions in $outdir dir.\n\n";
  }
__END__
