#!/usr/local/bin/perl

# 2015-09-25: Changed output to scientific notation
#       flux wasn't smooth with decimal output
#       at low flux values
#  Also added interpolation for lambda > 10 microns

# typical filter width is about 15% in wavelength
# want to average  lam f_lam dlam
$input  = $ARGV[0];
$output = $ARGV[1];
$finput = $ARGV[2];
$foutput = $ARGV[3];

# number of steps to interpolate between each spectral point beyond 10-microns:
$intsteps = 100;

# the marcs models have a separate file flx_wavelengths with the
# wavelengths, so this reads in the spectrum and the wavelengths
# the MARCS models are angstroms flambda
#       convert angstroms --> microns
if ($input =~ /s\d+_g/) {
  printf "I BELIEVE THIS IS A MARCS MODEL\n";
  open(INPUT1,"<$input");
  open(INPUT2,"<flx_wavelengths.vac");
  $n = 0;
  while(<INPUT1>) {
    $flux[$n] = $_;
    $wave[$n] = <INPUT2>;
    $flux[$n] =~ s/\n//;
    $wave[$n] =~ s/\n//;
    $wave[$n] /= 10000.0;  # convert angstroms to microns
    # normalize by the value close to 0.5 microns
    if ($n > 0) {
      if (($wave[$n] >= 0.5) && ($wave[$n-1] < 0.5)) {
        $flux5 = $flux[$n-1] + ($flux[$n]-$flux[$n-1])*(0.5-$wave[$n-1])/($wave[$n]-$wave[$n-1]);
        }
      }
    $n++;
    }
  close(INPUT1);
  close(INPUT2);
  printf "read $n entries from MARCS \n";
  }

# the kurucz files in /home/bhima/ckochanek/blackhole/dustform/kurucz/new
#         are in F_nu -- column 3 is wavelength in nm
#                        column 5 is F_nu
#         so nm need to be converted to microns
#            F_nu has to be converted to F_lambda = F_nu/lambda^2
#            dusty does not care about the absolute normalization
# things to look out for -- two leading lines
#                           final trailing line
#                           at end, columns 2/3 run together for 4 entries
#                           flux can be zero
if ($input =~ /f*\d+t/) {
  printf "I BELIEVE THIS IS A KURUCZ MODEL\n";
  open(INPUT1,"<$input");
  # there are two leading lines
  $junk = <INPUT1>;
  $junk = <INPUT1>;
  $k    = 0;
  while(<INPUT1>) {
    $line     = $_;
    if ($line =~ /FLUX\s*\d/) { 
      # in the last four lines, the columns get run together -- fix
      # which would break the split command
      $line =~ s/FLUX 1218/FLUX 1218 /;
      $line =~ s/FLUX 1219/FLUX 1219 /;
      $line =~ s/FLUX 1220/FLUX 1220 /;
      $line =~ s/FLUX 1221/FLUX 1221 /;
      $line     =~ s/\n//;
      @entries  = split(' ',$line);
      #can have entries like 1.3559-139 (no E) which confuses perl
      if ($entries[4] =~ /\d\-\d/) { $entries[4] = 0.0; }
      #useds 1.e-50 for zero even though then tries to fill in other numbers set to zero
      if ($entries[4] < 1.e-50) { $entries[4] = 0.0;}
      $origwave[$k] = $entries[2]/1000;          # convert nanometers to microns
      $origflux[$k] = $entries[4]/$origwave[$k]**2;  # convert F_nu to F_\lambda
      # normalize by the value close to 0.5 microns
      if ($k > 0) {
        if (($origwave[$k] >= 0.5) && ($origwave[$k-1] < 0.5)) {  
          $flux5 = $origflux[$k-1] + ($origflux[$k]-$origflux[$k-1])*(0.5-$origwave[$k-1])/($origwave[$k]-$origwave[$k-1]);
          }
        }
      $k++;
      }
    }
  close(INPUT1);
  printf "read $k entries from KURUCZ flux5 is $flux5\n";
  }

  # (linear) interpolate spectra (on log scale) for lambda>1-micron since the Kurucz resolution is coarser than 
  #   lambda_grid.dat for these wavelengths
  printf "interpolating spectra for lambda > 10-micron\n";
  $n = 0;
  for ($i=0; $i<$k-1; $i++) {
    if ($origwave[$i] > 1.0) {
      for ($int=0; $int<$intsteps; $int++) {
        $wave[$n] = 10**($int/$intsteps*(log10($origwave[$i+1])-log10($origwave[$i]))+log10($origwave[$i]));
        $flux[$n] = 10**($int/$intsteps*(log10($origflux[$i+1])-log10($origflux[$i]))+log10($origflux[$i]));
        #printf "$int $origwave[$i] $origflux[$i] $wave[$n] $flux[$n]\n";
        $n++;
        }
    } else {
      $wave[$n] = $origwave[$i];
      $flux[$n] = $origflux[$i];
      #printf "$int $origwave[$i] $origflux[$i] $wave[$n] $flux[$n] no interp\n";
      $n++;
    }
  }
  # add last point without interpolation
  $wave[$n] = $origwave[$i];
  $flux[$n] = $origflux[$i];
  $n++;
  #printf "$int $origwave[$i] $origflux[$i] $wave[$n] $flux[$n] last point\n";
  #for ($i=0; $i<$n; $i++) {
  #  printf "$wave[$i] $flux[$i]\n";
  #}


# while DUSTY does not care about normalization, if we are going
# to interpolate between models they should be on some common 
# amplitude scale -- 
for ($i=0; $i<$n; $i++) {
  $flux[$i] /= $flux5;
  }

# this reads in a filter -- "temp" if you are just smoothing
# an input spectrum for dusty rather than generating a filter model
open(INFILT,"<$finput");
$nf = 0;
while(<INFILT>) {
  $line = $_;
  $line =~ s/\n//;
  @entries = split(' ',$line);
  $flam[$nf] = $entries[0]/10000.0;
  $tran[$nf] = $entries[1];
  $nf++;
  }
printf "read $nf filter entries\n";

$output2 = $output . '2';

# this reads in the wavelength grid that will be used by dusty
open(INPUT,"<../../dustymc4000/lambda_grid.dat");
$line = <INPUT>;
$nw   = 0;
while(<INPUT>) { 
 $line = $_;
 $line =~ s/\n//;
 $dwave[$nw] = $line; 
 $nw++;
 }
printf "got $nw wavelengths from dusty\n";

# loop over the dusty wavelength grid
open(OUTPUT1,">$output");
open(OUTPUT2,">$output2");
open(FOUTPUT,">$foutput");
for ($k=0; $k<$nw; $k++) {
  # set the max/min wavelengths for this bin
  if ($k==0) {
    $dl   = $dwave[1]-$dwave[0];
    $lmin = $dwave[0]-$dl/2.0;
    $lmax = $dwave[0]+$dl/2.0;
  } else {
    if (($k+1)==$nw) {
       $dl   = $dwave[$nw-1]-$dwave[$nw-2];
       $lmin = $dwave[$nw-1]-$dl/2.0;
       $lmax = $dwave[$nw-1]+$dl/2.0;
    } else {
      $lmin = 0.5*($dwave[$k]+$dwave[$k-1]);
      $lmax = 0.5*($dwave[$k]+$dwave[$k+1]);
      }
    }
  # now loop over the high resolution spectrum and 
  # add up its contributions to this bin
  # bar0 = integrating lambda dlambda
  # bar1 = integrating lambda F_lambda dlambda
  # bar2 = integrating lambda^2 F_lambda dlambda
  $bar0 = 0.0;
  $bar1 = 0.0;
  $bar2 = 0.0;
  $nbin = 0;
  # the MARCS models are flambda
  for ($i=0; $i<$n-1; $i++) {
    if (($wave[$i] < $lmax) && ($wave[$i+1] >$lmin)) {
      $nbin++;
      $low   = ($wave[$i]   > $lmin) ? $wave[$i]   : $lmin;
      $hig   = ($wave[$i+1] < $lmax) ? $wave[$i+1] : $lmax;
      $dlam  = $hig-$low;
      if ($dlam < 0) {
        printf "error in getting bin width $dlam\n";
        exit(-1);
        }
      $bar0 += 0.5*$dlam*($wave[$i]+$wave[$i+1]);
      $bar1 += 0.5*$dlam*($wave[$i]*$flux[$i]+$wave[$i+1]*$flux[$i+1]);
      $bar2 += 0.5*$dlam*($wave[$i]*$wave[$i]*$flux[$i]+$wave[$i+1]*$wave[$i+1]*$flux[$i+1]);
      #printf "   %3d %3d %8.4f %8.4f %8.4f %13.6g %13.6g %13.6g %13.6g %13.6g\n",$i,$n,$dlam,$wave[$i],$wave[$i+1],$flux[$i],$flux[$i+1],$bar0,$bar1,$bar2;
      }
    }
  # done integrating, output to
  # OUTPUT1 = integrate lambda F_lambda dlambda / integrate dlambda
  # OTUPUT2 = integrate lambda^2 F_lambda dlambda / integrate F_lambda lambda dlambda
  if (($nbin > 0) && ($bar0 > 0) && ($bar1 > 0)) { 
    $value = $bar1/$bar0;
    $bar2 /= $bar1;
  } else {
    $value = 0.0;
    $bar2  = 0.0;
    }
  #printf "FINAL %13.6f %13.6f $nbin $bar0 $bar1 $bar2 \n",$dwave[$k],$value;
  printf OUTPUT1 "%13.6g %13.6g $nbin $bar0 $bar1 $bar2 \n",$dwave[$k],$value;
  printf OUTPUT2 "%13.6g %13.6g\n",$dwave[$k],$bar2;
  # now average the filter band pass over the dust wavelength grid 
  # bar0 = integrate lambda*dlambda
  # bar1 = integrate lambda Transmission_lambda dlambda
  # output bar1/bar0
  $bar0 = 0.0;
  $bar1 = 0.0;
  $nbin = 0;
  for ($i=0; $i<$nf-1; $i++) {
    if (($flam[$i] < $lmax) && ($flam[$i+1] >$lmin)) {
      $low   = ($wave[$i]   > $lmin) ? $wave[$i]   : $lmin;
      $hig   = ($wave[$i+1] < $lmax) ? $wave[$i+1] : $lmax;
      $dlam  = $hig-$low;
      $nbin++;
      $bar0 += 0.5*$dlam*($flam[$i]+$flam[$i+1]);
      $bar1 += 0.5*$dlam*($flam[$i]*$tran[$i]+$flam[$i+1]*$tran[$i+1]);
      }
    }
  if ($nbin > 0) { $bar1 /= $bar0; }
  printf FOUTPUT "%8.4g %8.4g \n",$dwave[$k],$bar1;
  }

sub log10 {
  my $n = shift;
  return log($n)/log(10);
}
