#!/usr/local/bin/perl


# typical filter width is about 15% in wavelength
# want to average  lam f_lam dlam
$input  = $ARGV[0];
$output = $ARGV[1];

open(INPUT1,"<$input");
open(INPUT2,"<flx_wavelengths.vac");
$n = 0;
while(<INPUT1>) {
  $flux[$n] = $_;
  $wave[$n] = <INPUT2>;
  $flux[$n] =~ s/\n//;
  $wave[$n] =~ s/\n//;
  $wave[$n] /= 10000.0;
  $n++;
  }
close(INPUT1);
close(INPUT2);
printf "read $n entries from MARCS \n";

$output2 = $output . '2';

open(OUTPUT1,">$output");
open(OUTPUT2,">$output2");
open(INPUT,"<../dustymc3600/lambda_grid.dat");
while(<INPUT>) {
  $lwant = $_;
  $lwant =~ s/\n//;
  $bar0 = 0.0;
  $bar1 = 0.0;
  $bar2 = 0.0;
  $nbin = 0;
  # average width of the broadband filters (FWHM/cetner) is 20%
  $lmin = 0.90*$lwant;
  $lmax = 1.10*$lwant;
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
    printf OUTPUT1 "%13.6f %13.6f\n",$lwant,$value;
    printf OUTPUT2 "%13.6f %13.6f\n",$lwant,$bar2;
    }
  }
