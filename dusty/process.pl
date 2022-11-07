#!/usr/bin/perl
# hst filters added by Scott - May 2014 - but only with "mag" support (not "mJy")
#   added filter name as the last column of the output file
$clight = 2.99e10;
$mp     = 1.66e-24;
$sigmat = 6.65e-25;
$au     = 1.496e13;
$msun   = 1.989e33;
$day    = 3600*24;
$kms    = 1.0e5;
$year   = 365*$day;
$lsun   = 3.83e33;
$sigmabb= 5.66e-5;
$kboltz = 1.38e-16;
$pc     = 3.09e18;
$planck = 6.626e-27;
$mpc    = 1.e6*$pc;
$jansky = 1.e-23;
$microjansky = $jansky/1.e6;
$micron = 1.0e-4;

$zerou = 1884;
$zerob = 4620;
$zerov = 3590;
$zerosdssuVega = 1568.54;
$zerosdssu = 3631;
$zerosdssg = 3631;
$zerosdssr = 3631;
$zerosdssi = 3631;
$zerosdssz = 3631;
$zeroDECamgAB = 3631;
$zeroDECamrAB = 3631;
$zeroDECamiAB = 3631;
$zeroDECamzAB = 3631;
$zeroDECamYAB = 3631;
$zeror = 3009.9;
$zeroi = 2408.8;
$zeroj = 1594.0;
$zeroh = 1024;
$zerok =  667;
$zeroz = 2261.4; #UKIDSS Z
$zeroy = 2057.2; #UKIDSS Y
$zerol =  278;
$zerolp=  233;
$zeromp=  167;
$zeron =   39.8;
$zero36 = 280.9;
$zero45 = 179.7;
$zero58 = 115.0;
$zero80 =  64.9;
$zero24=   7.14;
$zero70=   0.775;
$zero160=  0.159;
# WISE http://wise2.ipac.caltech.edu/docs/release/prelim/expsup/figures/sec4_3gt4.gif Fnu*iso
$zero34w = 306.682;
$zero46w = 170.663;
$zero12w =  29.045;
$zero22w =   8.284;
# Stromgren mags -- on Vega system -- Fabregat and Reig 1996 
$zerostrmu = 4680.0;
$zerostrmv = 4890.0;
$zerostrmb = 4250.0;
$zerostrmy = 3700.0;
# Swift on AB system
$zeroSV   = 3631.0;
$zeroSB   = 3631.0;
$zeroSU   = 3631.0;
$zeroSW1  = 3631.0;
$zeroSWM  = 3631.0;
$zeroSW2  = 3631.0;
$zeroPS1r = 3151.44; # Vega
$zeroPS1i = 2584.56; # Vega

# HST filters added by Scott 30 Apr 2014 from http://www.stsci.edu/hst/wfc3/phot_zp_lbn
#  see also: /home/poseidon/sadams/impostors/dusty/NOTES  
$zeroWFC3uvF275W = 913.495;
$zeroWFC3uvF336W = 1219.51;
$zeroWFC3uvF438W = 4181.25;
$zeroWFC3uvF475W = 3941.67;
$zeroWFC3uvF555W = 3716.95;
$zeroWFC3uvF606W = 3365.32;
$zeroWFC3uvF775W = 2541.32;
$zeroWFC3uvF814W = 2470.28;
$zeroWFC3irF110W = 1803.96;
$zeroWFC3irF160W = 1146.74;
$zeroWFPC2F606W = 3274.2;
$zeroWFPC2F439W = 4171.25;
$zeroWFPC2F450W = 3929.99;
$zeroWFPC2F675W = 2879.46;
$zeroNIC3F160W = 1071.23;
# JWST filters
$zeroNIRCamF115W = 1753.8;
$zeroNIRCamF150W = 1183.1;
$zeroNIRCamF200W = 765.7;

$lamNIC3F160W = 1.578588; # microns
$lamPS1i = 0.74846; # microns
$lamPS1r = 0.612953; # microns
$lamWFC3uvF275W = 0.27102; # microns
$lamWFC3uvF336W = 0.33548;
$lamWFC3uvF438W = 0.43265;
$lamWFC3uvF475W = 0.46961;
$lamWFC3uvF555W = 0.53081;
$lamWFC3uvF606W = 0.58874;
$lamWFC3uvF775W = 0.75867;
$lamWFC3uvF814W = 0.80295;
$lamWFC3irF110W = 1.15340;
$lamWFC3irF160W = 1.53690;
$lamWFPC2F606W = 0.57344;
$lamWFPC2F439W = 0.429134; 
$lamWFPC2F450W = 0.444413;
$lamWFPC2F675W = 0.670297;

# More HST filters added by Scott 30 Oct 2014 from http://svo2.cab.inta-csic.es/theory/fps/index.php?mode=browse&gname=HST&gname2=WFC3_UVIS1
$zeroACSWFCF435W = 4039.8; # Jy
$zeroACSWFCF555W = 3621.9;
$zeroACSWFCF814W = 2414.9;
$zeroWFPC2F555W = 3661.0;
$zeroWFPC2F814W = 2372.1;

$lamACSWFCF435W = 0.43319; # microns
$lamACSWFCF555W = 0.53309;
$lamACSWFCF814W = 0.79854;
$lamWFPC2F555W = 0.53301;
$lamWFPC2F814W = 0.83695;

# Even more HST filters added 20 May 2015 from http://svo2.cab.inta-csic.es/theory/fps/index.php?mode=browse&gname=HST&gname2=ACS_WFC
$zeroACSWFCF475W = 3931.7; # Jy
$zeroACSWFCF606W = 3233.0;
$zeroACSWFCF658N = 3079.9;

$lamACSWFCF475W = 0.47082; # microns
$lamACSWFCF606W = 0.58082;
$lamACSWFCF658N = 0.65841;
# JWST
$lamNIRCamF115W = 1.14;
$lamNIRCamF150W = 1.50;
$lamNIRCamF200W = 1.96;
#
$zeroks = 666.8;
$lamks = 2.159; # 2MASS Ks

$lamu  =  0.36;
$lamb  =  0.44;
$lamv  =  0.55;
$lamr  =  0.658;
$lami  =  0.806;
$lamsdssuVega = 0.359493;
$lamsdssu = 0.354;
$lamsdssg = 0.477;
$lamsdssr = 0.623;
$lamsdssi = 0.763;
$lamsdssz = 0.913;
$lamDECamgAB = 0.4734;
$lamDECamrAB = 0.63452;
$lamDECamiAB = 0.77496;
$lamDECamzAB = 0.91382;
$lamDECamYAB = 0.98758;
$lamz  =  0.881;
$lamy  =  1.03;
$lamj  =  1.24;
$lamh  =  1.62;
$lamk  =  2.19;
$laml  =  3.4;
$lamlp =  3.8;
$lammp =  4.7;
$lamn  = 10.2;
$lam36 =  3.6;
$lam45 =  4.5;
$lam58 =  5.8;
$lam80 =  8.0;
$lam24 = 23.68;
$lam70 = 71.42;
$lam160=155.90;
$lam11 =  11.2;
$lamnuv=  0.232;
$lamfuv=  0.154;
#akari
$lam24a= 2.4;
$lam32a= 3.2;
$lam41a= 4.1;
#Wise
$lam34w = 3.35;
$lam46w = 4.60;
$lam12w = 11.56;
$lam22w = 22.08;
# fake MIPS
$lam7mips = 7.0;
$lam9mips = 9.0;
$lam10mips = 10.0;
$lam13mips = 13.0;
$lam14mips = 14.0;
# SPITZER IRS
$lam15irs = 15.6;
$lam22irs = 22.0;
# SWIFT UVOT zero points
# these are from the WWW site, but Berger gives 0.252, 0.188
# http://www.mssl.ucl.ac.uk/www_astro/uvot/uvot_instrument/filterwheel/filterwheel.html
$lamSV  = 0.543;
$lamSB  = 0.434;
$lamSU  = 0.344;
$lamSW1 = 0.2910;
$lamSW2 = 0.2120;
$lamSWM = 0.2310;
# these are from some other swift www site
#$lamSW1 = 0.2634;  0.2910;
#$lamSW2 = 0.2030;  0.2120;
#$lamSWM = 0.2231;  0.2310;
# Stromgren
$lamstrmu = 0.345;
$lamstrmv = 0.411;
$lamstrmb = 0.467;
$lamstrmy = 0.548;
# Hershchel
$lamh70   =  70.0;
$lamh100  = 100.0;
$lamh160  = 160.0;

$ln10 = log(10.0);
$pi   = 3.141592654;

$open = 0;
open(INPUT,"grep -v \# $ARGV[0] |");
#open(INPUT,"grep -v \# DATA |");
while(<INPUT>) { 
  $line    = $_;
  $line   =~ s/\n//;
  @entries = split(' ',$line);
  $done    = 0;
  $flux    = 0;
  printf "doing $line\n";
  if ( ($line =~ /dist/) && ($done == 0)) {
    $dist = $entries[1]*$mpc;
    printf "  setting distance to $entries[1]\n";
    $done = 1;
    }
  if ( ($line =~ /egal/) && ($done == 0)) {
    $ebv = $entries[1];
    printf "  setting extinction  to E(B-V) = $ebv\n";
    $done = 1;
    }
  if ( ($line =~ /dmod/ && ($done == 0)) ) {
    $dist = 10.0*$pc*10.0**($entries[1]/5.0);
    printf "  setting distance to $dist for distance modulus $entries[0]\n";
    $done = 1;
    }
  if ( ($line =~ /file/) && ($done == 0) ) {
    if ($open == 1) { close(OUTPUT); }
    printf "  opening $entries[1] for output \n";
    open(OUTPUT,">$entries[1]");
    $done = 1;
    }
  if ( ($line =~ /mag/i) && ($done == 0)) {
    printf "  processing as magnitude\n";
    $limit = 0;
    if ($entries[1] < 0) { $limit = 1;}

    if ( !($line =~ /strm/) && !($line =~ /sdss/) && !($line =~ /S[UBV]/i) ) {  
      if ($line =~ /U/)   {  $jy = $zerou*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamu;  $done = 1;}
      if ($line =~ /B/)   {  $jy = $zerob*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamb;  $done = 1;}
      if ($line =~ /V/)   {  $jy = $zerov*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamv;  $done = 1;}
      if ($line =~ /R/)   {  $jy = $zeror*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamr;  $done = 1;}
      if ($line =~ /I/)   {  $jy = $zeroi*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lami;  $done = 1;}
      }
    if ($line =~ /sdssu/){ $jy = $zerosdssu*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamsdssu;  $done = 1;}
    if ($line =~ /sdssg/){ $jy = $zerosdssg*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamsdssg;  $done = 1;}
    if ($line =~ /sdssr/){ $jy = $zerosdssr*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamsdssr;  $done = 1;}
    if ($line =~ /sdssi/){ $jy = $zerosdssi*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamsdssi;  $done = 1;}
    if ($line =~ /sdssz/){ $jy = $zerosdssz*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamsdssz;  $done = 1;}
    if ($line =~ /Z/)   {  $jy = $zeroz*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamz;  $done = 1;}
    if ($line =~ /Y/)   {  $jy = $zeroy*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamy;  $done = 1;}
    if ($line =~ /J/)   {  $jy = $zeroj*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamj;  $done = 1;}
    if ($line =~ /H/)   {  $jy = $zeroh*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamh;  $done = 1;}
    if ($line =~ /K/)   {  $jy = $zerok*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamk;  $done = 1;}
    if ($line =~ /Ks/)   {  $jy = $zeroks*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamks;  $done = 1;}
    if (($line =~ /L/) && ($line =~ /Lp/))   {  $jy = $zerol*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $laml; $filer = $line;  $done = 1;}
    if ($line =~ /Lp/)   {  $jy = $zerolp*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamlp;  $done = 1;}
    if ($line =~ /Mp/)   {  $jy = $zeromp*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lammp;  $done = 1;}
    if (($line =~ /N/ && !($line =~ /UV/)) )   {  $jy = $zeron*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamn;  $done = 1;}
    if ($line =~ /F36/) {  $jy = $zero36*10**(-0.4*$entries[0]); $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lam36; $done = 1;}
    if ($line =~ /F45/) {  $jy = $zero45*10**(-0.4*$entries[0]); $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lam45; $done = 1;}
    if ($line =~ /F58/) {  $jy = $zero58*10**(-0.4*$entries[0]); $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lam58; $done = 1;}
    if ($line =~ /F80/) {  $jy = $zero80*10**(-0.4*$entries[0]); $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lam80; $done = 1;}
    if ( ($line =~ /F24/) && !($line =~ /F24A/) ) {  
                           $jy = $zero24*10**(-0.4*$entries[0]); $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lam24; $done = 1;}
    if ($line =~ /F70/) {  $jy = $zero70*10**(-0.4*$entries[0]); $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lam70; $done = 1;}
    if ($line =~ /F160/) {  $jy = $zero160*10**(-0.4*$entries[0]); $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lam160; $done = 1;}
    if ($line =~ /F24A/) {  printf "NO AKARI MAGNITUDES \n"; exit(-1); }
    if ($line =~ /F32A/) {  printf "NO AKARI MAGNITUDES \n"; exit(-1); }
    if ($line =~ /F41A/) {  printf "NO AKARI MAGNITUDES \n"; exit(-1); }
    #if ($line =~ /F11/)  {  printf "NO Michelle MAGNITUDES \n"; exit(-1); }
    if ($line == /F11/)  {  printf "NO Michelle MAGNITUDES \n"; exit(-1); }
    if ($line =~ /F34W/) {  $jy = $zero34w*10**(-0.4*$entries[0]); $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lam34w; $done = 1;}
    if ($line =~ /F46W/) {  $jy = $zero46w*10**(-0.4*$entries[0]); $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lam46w; $done = 1;}
    if ($line =~ /F12W/) {  $jy = $zero12w*10**(-0.4*$entries[0]); $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lam12w; $done = 1;}
    if ($line =~ /F22W/) {  $jy = $zero22w*10**(-0.4*$entries[0]); $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lam22w; $done = 1;}
    if ($line =~ /MIPS/) {  printf "NO MIPS MAGNITUDES \n"; exit(-1); }
    if ($line =~ /FH/)   {  printf "NO Heschel MAGNITUDES \n"; exit(-1); }
    if ($line =~ /SV/)   {  $jy = $zeroSV*10**(-0.4*$entries[0]); $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamSV;  $done = 1;}
    if ($line =~ /SB/)   {  $jy = $zeroSB*10**(-0.4*$entries[0]); $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamSB;  $done = 1;}
    if ($line =~ /SU/)   {  $jy = $zeroSU*10**(-0.4*$entries[0]); $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamSU;  $done = 1;}
    if ($line =~ /SW1/)  {  $jy = $zeroSW1*10**(-0.4*$entries[0]); $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamSW1;  $done = 1;}
    if ($line =~ /SW2/)  {  $jy = $zeroSW2*10**(-0.4*$entries[0]); $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamSW2;  $done = 1;}
    if ($line =~ /SWM/)  {  $jy = $zeroSWM*10**(-0.4*$entries[0]); $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamSWM;  $done = 1;}
    if ($line =~ /IRS/i) {  printf "NO IRS MAGNITUDES \n"; exit(-1); }
    # Stromgren
    if ($line =~ /strmu/) {  $jy = $zerostrmu*10**(-0.4*$entries[0]); $ejy = 0.4*$ln10*$jy*$entries[1];$lam = $lamstrmu;  $done = 1;}
    if ($line =~ /strmv/) {  $jy = $zerostrmv*10**(-0.4*$entries[0]); $ejy = 0.4*$ln10*$jy*$entries[1];$lam = $lamstrmv;  $done = 1;}
    if ($line =~ /strmb/) {  $jy = $zerostrmb*10**(-0.4*$entries[0]); $ejy = 0.4*$ln10*$jy*$entries[1];$lam = $lamstrmb;  $done = 1;}
    if ($line =~ /strmy/) {  $jy = $zerostrmy*10**(-0.4*$entries[0]); $ejy = 0.4*$ln10*$jy*$entries[1];$lam = $lamstrmy;  $done = 1;}

    # Addition for HST Filters by Scott Adams - 1 May 2014 
    #   Note: I had issues before when placing these HST filters earlier because the lam & jy were being rewritten by I & U band interpretations
    #     since "=~" just means string contains, not string identical to.
    if ($line =~ /NIC3F160W/)   {  $jy = $zeroNIC3F160W*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamNIC3F160W;  $done = 1;}
    if ($line =~ /WFC3uvF275W/)   {  $jy = $zeroWFC3uvF275W*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamWFC3uvF275W;  $done = 1;}
    if ($line =~ /WFC3uvF336W/)   {  $jy = $zeroWFC3uvF336W*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamWFC3uvF336W;  $done = 1;}
    if ($line =~ /WFC3uvF438W/)   {  $jy = $zeroWFC3uvF438W*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamWFC3uvF438W;  $done = 1;}
    if ($line =~ /WFC3uvF475W/)   {  $jy = $zeroWFC3uvF475W*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamWFC3uvF475W;  $done = 1;}
    if ($line =~ /WFC3uvF555W/)   {  $jy = $zeroWFC3uvF555W*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamWFC3uvF555W;  $done = 1;}
    if ($line =~ /WFC3uvF606W/)   {  $jy = $zeroWFC3uvF606W*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamWFC3uvF606W;  $done = 1;}
    if ($line =~ /WFC3uvF775W/)   {  $jy = $zeroWFC3uvF775W*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamWFC3uvF775W;  $done = 1;}
    if ($line =~ /WFC3uvF814W/)   {  $jy = $zeroWFC3uvF814W*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamWFC3uvF814W;  $done = 1;}
    if ($line =~ /WFC3irF110W/)   {  $jy = $zeroWFC3irF110W*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamWFC3irF110W;  $done = 1;}
    if ($line =~ /WFC3irF160W/)   {  $jy = $zeroWFC3irF160W*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamWFC3irF160W;  $done = 1;}
    if ($line =~ /WFPC2F606W/)   {  $jy = $zeroWFPC2F606W*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamWFPC2F606W;  $done = 1;}

    if ($line =~ /ACSWFCF435W/)   {  $jy = $zeroACSWFCF435W*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamACSWFCF435W;  $done = 1;}
    if ($line =~ /ACSWFCF475W/)   {  $jy = $zeroACSWFCF475W*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamACSWFCF475W;  $done = 1;}
    if ($line =~ /ACSWFCF555W/)   {  $jy = $zeroACSWFCF555W*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamACSWFCF555W;  $done = 1;}
    if ($line =~ /ACSWFCF606W/)   {  $jy = $zeroACSWFCF606W*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamACSWFCF606W;  $done = 1;}
    if ($line =~ /ACSWFCF814W/)   {  $jy = $zeroACSWFCF814W*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamACSWFCF814W;  $done = 1;}
    if ($line =~ /ACSWFCF658N/)   {  $jy = $zeroACSWFCF658N*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamACSWFCF658N;  $done = 1;}
    if ($line =~ /WFPC2F555W/)   {  $jy = $zeroWFPC2F555W*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamWFPC2F555W;  $done = 1;}
    if ($line =~ /WFPC2F814W/)   {  $jy = $zeroWFPC2F814W*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamWFPC2F814W;  $done = 1;}
    if ($line =~ /WFPC2F439W/)   {  $jy = $zeroWFPC2F439W*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamWFPC2F439W;  $done = 1;}
    if ($line =~ /WFPC2F450W/)   {  $jy = $zeroWFPC2F450W*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamWFPC2F450W;  $done = 1;}
    if ($line =~ /WFPC2F675W/)   {  $jy = $zeroWFPC2F675W*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamWFPC2F675W;  $done = 1;}
    if ($line =~ /NIRCamF115W/)   {  $jy = $zeroNIRCamF115W*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamNIRCamF115W;  $done = 1;}
    if ($line =~ /NIRCamF150W/)   {  $jy = $zeroNIRCamF150W*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamNIRCamF150W;  $done = 1;}
    if ($line =~ /NIRCamF200W/)   {  $jy = $zeroNIRCamF200W*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamNIRCamF200W;  $done = 1;}

    if ($line =~ /lbcbBesselU/)   {  $jy = $zerou*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamu;  $done = 1;}
    if ($line =~ /lbcbBesselB/)   {  $jy = $zerob*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamb;  $done = 1;}
    if ($line =~ /lbcbBesselV/)   {  $jy = $zerov*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamv;  $done = 1;}
    if ($line =~ /lbcbBesselR/)   {  $jy = $zeror*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamr;  $done = 1;}
    if ($line =~ /PS1r/)   {  $jy = $zeroPS1r*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamPS1r;  $done = 1;}
    if ($line =~ /PS1i/)   {  $jy = $zeroPS1i*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamPS1i;  $done = 1;}

    if ($line =~ /DECamgAB/)   {  $jy = $zeroDECamgAB*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamDECamgAB;  $done = 1;}
    if ($line =~ /DECamrAB/)   {  $jy = $zeroDECamrAB*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamDECamrAB;  $done = 1;}
    if ($line =~ /DECamiAB/)   {  $jy = $zeroDECamiAB*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamDECamiAB;  $done = 1;}
    if ($line =~ /DECamzAB/)   {  $jy = $zeroDECamzAB*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamDECamzAB;  $done = 1;}
    if ($line =~ /DECamYAB/)   {  $jy = $zeroDECamYAB*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamDECamYAB;  $done = 1;}
    if ($line =~ /sdssuVega/)   {  $jy = $zerosdssuVega*10**(-0.4*$entries[0]);  $ejy = 0.4*$ln10*$jy*$entries[1]; $lam = $lamsdssuVega;  $done = 1;}

    if ($done == 1) { $flux = 1; }
    }   

 
  # generic mjy
  if ( ($line =~ /default/i) && ($done == 0)) {
    printf "  processing as default mJy \n";
    $limit = 0;
    if ($entries[1] < 0) { $limit = 1;}
    $jy    = $entries[0]/1000;
    $ejy   = $entries[1]/1000;
    $rval  = 0;
    $lam   = $entries[3];
    if (!($entries[2] =~ /default/i)) {
      printf "default entry must by Jy e(Jy) default Lambda \n";
      exit(-1);
      }
    $flux  = 1;
    $done  = 1;
    }
  if ( ($line =~ /mJy/i) && ($done == 0)) {
    printf "  processing as mJy \n";
    $limit = 0;
    if ($entries[1] < 0) { $limit = 1;}
    $jy    = $entries[0]/1000;
    $ejy   = $entries[1]/1000;
    if ( !($line =~ /strm/) && !($line =~ /sdss/) && !($line =~ /S[UBV]/i) ) {  
      if ($line =~ /U/)   { $lam = $lamu;  $done = 1;}
      if ($line =~ /B/)   { $lam = $lamb;  $done = 1;}
      if ($line =~ /V/)   { $lam = $lamv;  $done = 1;}
      if ($line =~ /R/)   { $lam = $lamr;  $done = 1;}
      }
    if ($line =~ /sdssu/){$lam = $lamsdssu;  $done = 1;}
    if ($line =~ /sdssg/){$lam = $lamsdssg;  $done = 1;}
    if ($line =~ /sdssr/){$lam = $lamsdssr;  $done = 1;}
    if ($line =~ /sdssi/){$lam = $lamsdssi;  $done = 1;}
    if ($line =~ /sdssz/){$lam = $lamsdssz;  $done = 1;}
    if ($line =~ /I/)   { $lam = $lami; $done = 1;}
    if ($line =~ /Z/)   { $lam = $lamz; $done = 1;}
    if ($line =~ /Y/)   { $lam = $lamy; $done = 1;}
    if ($line =~ /J/)   { $lam = $lamj; $done = 1;}
    if ($line =~ /H/)   { $lam = $lamh; $done = 1;}
    if ($line =~ /K/)   { $lam = $lamk; $done = 1;}
    if ($line =~ /L/)   { $lam = $laml; $done = 1;}
    if (($line =~ /N/ && !($line =~ /UV/)) )   { $lam = $laml;  $done = 1;}
    if ($line =~ /F36/) { $lam = $lam36;  $done = 1;}
    if ($line =~ /F45/) { $lam = $lam45;  $done = 1;}
    if ($line =~ /F58/) { $lam = $lam58;  $done = 1;}
    if ($line =~ /F80/) { $lam = $lam80;  $done = 1;}
    if (($line =~ /F24/) && !($line =~ /F 24A/)) 
                        { $lam = $lam24;  $done = 1;}
    if ($line =~ /F70/) { $lam = $lam70;  $done = 1;}
    if ($line =~ /F24A/){ $lam = $lam24a;  $done = 1;}
    if ($line =~ /F32A/){ $lam = $lam32a;  $done = 1;}
    if ($line =~ /F41A/){ $lam = $lam41a;  $done = 1;}
    if ($line =~ /F11/) { $lam = $lam11;  $done = 1;}
    if ($line =~/F160/) { $lam =$lam160;  $done = 1;}
    if ($line =~ /NUV/) { $lam =$lamnuv;  $done = 1;}
    if ($line =~ /FUV/) { $lam =$lamfuv;  $done = 1;}
    #Herschel
    if ($line =~ /FH70/)  { $lam = $lamh70;  $done = 1;}
    if ($line =~ /FH100/) { $lam = $lamh100; $done = 1;}
    if ($line =~ /FH160/) { $lam = $lamh160; $done = 1;}
    # WISE
    if ($line =~ /F34W/) { $lam =$lam34w;  $done = 1;}
    if ($line =~ /F46W/) { $lam =$lam46w;  $done = 1;}
    if ($line =~ /F12W/) { $lam =$lam12w;  $done = 1;}
    if ($line =~ /F22W/) { $lam =$lam22w;  $done = 1;}
    # MIPS
    if ($line =~ /F7MIPS/)  { $lam =$lam7mips;  $done = 1;}
    if ($line =~ /F9MIPS/)  { $lam =$lam9mips;  $done = 1;}
    if ($line =~ /F10MIPS/) { $lam =$lam10mips; $done = 1;}
    if ($line =~ /F13MIPS/) { $lam =$lam13mips; $done = 1;}
    if ($line =~ /F14MIPS/) { $lam =$lam14mips; $done = 1;}
    # IRS
    if ($line =~ /F15IRS/) { $lam =$lam15irs; $done = 1;}
    if ($line =~ /F22IRS/) { $lam =$lam22irs; $done = 1;}
    # SWIFT UVOT
    if ($line =~ /SV/)      {  $lam = $lamSV;  $done = 1;}
    if ($line =~ /SB/)      {  $lam = $lamSB;  $done = 1;}
    if ($line =~ /SU/)      {  $lam = $lamSU;  $done = 1;}
    if ($line =~ /SW1/)     {  $lam = $lamSW1; $done = 1;}
    if ($line =~ /SW2/)     {  $lam = $lamSW2; $done = 1;}
    if ($line =~ /SWM/)     {  $lam = $lamSWM; $done = 1;}
    # Stromgren
    if ($line =~ /strmu/)     {  $lam = $lamstrmu; $done = 1;}
    if ($line =~ /strmv/)     {  $lam = $lamstrmv; $done = 1;}
    if ($line =~ /strmb/)     {  $lam = $lamstrmb; $done = 1;}
    if ($line =~ /strmy/)     {  $lam = $lamstrmy; $done = 1;}

    # HST -- added by Scott 2015-05-15
    if ($line =~ /WFC3uvF275W/)   { $lam = $lamWFC3uvF275W;  $done = 1;}
    if ($line =~ /WFC3uvF336W/)   { $lam = $lamWFC3uvF336W;  $done = 1;}
    if ($line =~ /WFC3uvF438W/)   { $lam = $lamWFC3uvF438W;  $done = 1;}
    if ($line =~ /WFC3uvF475W/)   { $lam = $lamWFC3uvF475W;  $done = 1;}
    if ($line =~ /WFC3uvF555W/)   { $lam = $lamWFC3uvF555W;  $done = 1;}
    if ($line =~ /WFC3uvF775W/)   { $lam = $lamWFC3uvF775W;  $done = 1;}
    if ($line =~ /WFC3uvF814W/)   { $lam = $lamWFC3uvF814W;  $done = 1;}
    if ($line =~ /WFC3irF110W/)   { $lam = $lamWFC3irF110W;  $done = 1;}
    if ($line =~ /WFC3irF160W/)   { $lam = $lamWFC3irF160W;  $done = 1;}
    if ($line =~ /WFPC2F606W/)   { $lam = $lamWFPC2F606W;  $done = 1;}


    if ($line =~ /ACSWFCF435W/)   { $lam = $lamACSWFCF435W;  $done = 1;}
    if ($line =~ /ACSWFCF555W/)   { $lam = $lamACSWFCF555W;  $done = 1;}
    if ($line =~ /ACSWFCF814W/)   { $lam = $lamACSWFCF814W;  $done = 1;}
    if ($line =~ /WFPC2F555W/)   { $lam = $lamWFPC2F555W;  $done = 1;}
    if ($line =~ /WFPC2F814W/)   { $lam = $lamWFPC2F814W;  $done = 1;}
    if ($line =~ /ACSWFCF475W/)   { $lam = $lamACSWFCF475W;  $done = 1;}
    if ($line =~ /ACSWFCF606W/)   { $lam = $lamACSWFCF606W;  $done = 1;}
    if ($line =~ /ACSWFCF658N/)   { $lam = $lamACSWFCF658N;  $done = 1;}

    if ($line =~ /lbcbBesselR/)   { $lam = $lamr;  $done = 1;}
    if ($line =~ /lbcbBesselV/)   { $lam = $lamv;  $done = 1;}
    if ($line =~ /lbcbBesselB/)   { $lam = $lamb;  $done = 1;}
    if ($line =~ /lbcbBesselU/)   { $lam = $lamu;  $done = 1;}

    if ($done == 1) { $flux = 1; }
    }
 
 
  if ($flux == 1) { 
    if ($entries[2] =~/mjy/) {$filter = $entries[3];}
    if ($entries[3] =~/mag/) {$filter = $entries[2];}
    printf "  flux flag is $flux, writing luminosities flux for wavelength $lam is $jy janskys \n";
    $lum = ($jy*$jansky/$lsun)*(4*$pi*$dist*$dist)*($clight/($lam*$micron));
    printf "  luminosity $lum\n";
    $elum= $lum*$ejy/$jy;
    # correct for extinction
    $rval = &rl($lam);
    $ecor= 10.0**(-0.4*$rval*$ebv);
    $lumold = $lum;
    $lum /= $ecor;
    $elum/= $ecor;
    printf "  using R $rval for wavelength $lam ebv is $ebv luminosity correction factor is $ecor old lum $lumold new lum $lum\n";
    if ($limit == 0) { 
      printf "  Writing luminosity $lam $lum $elum ecor $ecor \n";
      $minerr = 0.01;
      #$minerr = 0.2;	# This was Kochanek's original setting
      #$minerr = 0.1;	# Lets go for 10% uncertainty
      if ($elum < $minerr*$lum) {
        printf "resetting error bar to minimum fractinal error of $minerr\n";
        $elum = $minerr*$lum;
        }
      printf OUTPUT "%6.3f %9.3f %10.4f %s\n", $lam,$lum,$elum,$filter;
    } else {
      printf "  Writing limit $lam $lum \n";
      printf OUTPUT "%6.3f -1  %12.6f %s\n", $lam,$lum,$filter;
      } 
    }
  if ($done == 0) { 
    printf "  failed to parse $line\n";
    }
  }
close(OUTPUT);

#$lam = 10;
#temp = &rl($lam);
#rintf "for wavelength $lam rl = $temp\n";

sub rl {
  $x = 1.0/@_[0];
  printf "x is $x\n";
  $rv = 3.1;

  if ($x < 1.1) {
    $a =  0.574*$x**1.61;
    $b = -0.527*$x**1.61;
  } else {
    if ($x < 3.3) {
      $y = $x-1.82;
      $a = 1.0+0.17699*$y-0.50447*$y*$y-0.02427*$y**3+0.72085*$y**4 +
               0.01979*$y**5-0.77530*$y**6+0.32999*$y**7;
      $b = 1.41338*$y+2.28305*$y**2+1.07233*$y**3-5.38434*$y**4-
                      0.62251*$y**5+5.30260*$y**6-2.09002*$y**7;
    } else {
      if ($x < 5.9) {
        $fa = 0.0;
        $fb = 0.0;
      } else {
        $fa = -0.04473*($x-5.9)**2 - 0.009779*($x-5.9)**3;
        $fb =  0.21300*($x-5.9)**2 + 0.120700*($x-5.9)**3;
        }
      $a =  1.752 - 0.316*$x - 0.104/(($x-4.67)**2+0.341) + $fa;
      $b = -3.090 + 1.825*$x + 1.206/(($x-4.67)**2+0.263) + $fb;
      }
    }

    $rl = $rv*($a+$b/$rv);


  }
