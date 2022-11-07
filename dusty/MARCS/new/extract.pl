#!/usr/bin/perl
#
#$ng   = 11;
#$g[0] = '00'; 
#$g[1] = '05'; 
#$g[2] = '10'; 
#$g[3] = '15'; 
#$g[4] = '20'; 
#$g[5] = '25'; 
#$g[6] = '30'; 
#$g[7] = '35'; 
#$g[8] = '40'; 
#$g[9] = '45'; 
#$g[10]= '50'; 
$nfe = 5;
$fe[0] = 'p02';
$fe[1] = 'p05';
$fe[2] = 'm05';
$fe[3] = 'm10';
$fe[4] = 'p00';

# see http://www.oact.inaf.it/castelli/castelli/grids.html for list of models
open(GETTEMP,"<temp.dat");
while(<GETTEMP>) {
  $temp = $_;
  $temp =~ s/^\s*//;
  $temp =~ s/\s*\n//;
  for ($i=0; $i<$nfe; $i++) {
    $file = 'http://www.oact.inaf.it/castelli/castelli/grids/grid' . $fe[$i] . 'k2odfnew/f' . $fe[$i] . 't' . $temp . 'g40k2odfnew.dat';
    #$file = 'http://wwwuser.oat.ts.astro.it/castelli/grids/gridp00k2odfnew/fp00t' . $temp . 'g' . $g[$i] . 'k2odfnew.dat';
    print "$file\n";
    open(GETFILE,"wget $file |");
    while(<GETFILE>) {
      $line = $_;
      printf "$line";
      }  
    close(GETFILE);
    }
  }
