#!/usr/bin/perl


# typical filter width is about 15% in wavelength
# want to average  lam f_lam dlam

# kurucz is fnu
$n = 0;
open(INFILE,"grep -v \# kurucz_Ostar.dat|");
while(<INFILE>) { 
  $line = $_;
  $line =~ s/\n//;
  @entries = split(' ',$line);
  $wave[$n] = $entries[0];
  $flux[$n] = $entries[1];
  # convert to flam
  $flux[$n] /= ($wave[$n]**2);
  $n++;
  }
close(INFILE);
printf "read $n entries from MARCS \n";

$output2 = $output . '2';

open(OUTPUT2,">kurucz.dat2");
open(INPUT,"<../dustymc3600/lambda_grid.dat");
while(<INPUT>) {
  $lwant = $_;
  $lwant =~ s/\n//;
  $bar0 = 0.0;
  $bar1 = 0.0;
  $bar2 = 0.0;
  $nbin = 0;
  $lmin = 0.92*$lwant;
  $lmax = 1.08*$lwant;
  for ($i=0; $i<$n; $i++) {
    if (($wave[$i] >= $lmin) && ($wave[$i] <= $lmax)) {
      $nbin++;
      $bar0 += $wave[$i];
      $bar1 += $wave[$i]*$flux[$i];
      $bar2 += $wave[$i]*$wave[$i]*$flux[$i];
      }
    }
  if ($nbin > 0) { 
    $value = $bar1/$bar0;
    $bar2 /= $bar1;
    printf OUTPUT2 "%13.6f %13.6f %4d \n",$lwant,$bar2,$nbin;
    }
  }

