import corner
import matplotlib.pyplot as plt
from astropy.table import Table, Column, vstack
from astropy.io import ascii
import numpy as np
from os import system
import os
import sys

def geninput(x,blackbody=True,tstarmin=3500,tstarmax=48999):
    #x = {}
    #for n in range(len(varxnames)):
    #    x[varxnames[n]] = varx[n]
    #for n in range(len(constxnames)):
    #    x[constxnames[n]] = constx[n]
    tstar = x['tstar']
    td = x['td']
    tau = x['tau']
    thick = x['thick']
    ebv = x['ebv']
    idtype = x['idtype']
    fileuse = x['fileuse']
    output = open('foo1.inp','w')
    if (tstar < tstarmin or tstar > tstarmax) or (blackbody):
        output.write('Spectrum = 1\n')
        output.write('Number of BB = 1\n')
        output.write('Temperature = %0.2f\n' % (tstar))
    else:
        output.write('Spectrum = 5   \n')
        output.write('%s\n' % (fileuse))
    output.write('   optical properties index = 1 \n')
    output.write('   #   Sil-Ow  Sil-Oc  Sil-DL  grf-DL  amC-Hn  SiC-Pg \n')
    if idtype == 0:
        output.write('    x = 0.00    0.00   0.00    1.00    0.00    0.00 \n')
    else:
        output.write('    x = 0.00    0.00   1.00    0.00    0.00    0.00 \n')
    if custom_grain_distribution:
        output.write('- size distribution = 2  % custom       \n')
        output.write('  q = 3.5, a(min) = %s micron, a(max) = %s micron\n' % (aminnew,amaxnew))
    else:
        output.write('- size distribution = 1  % standard MRN    \n')
    output.write('- temperature = %s K \n' % (td))
    output.write('- density type = 1                   \n')
    output.write('- number of powers = 1              \n')
    output.write('- shells relative thickness = %s\n' % (thick))
    output.write('- power = 2 \n')
    output.write('- grid type = 1                  % linear grid \n')
    output.write('- lambda0 = 0.55 micron          % optical depth specified  \n')
    output.write('- tau(min) = '+str(tau)+' ; tau(max) = 1000.0   % for the visual wavelength \n')
    output.write('- number of models = 1           \n')
    output.write('- accuracy for flux conservation = 0.05             \n')
    output.write('- verbosity flag;                              verbose = 1  \n')
    output.write('- properties of emerging spectra;            fname.spp = 1  \n')
    output.write('- detailed spectra for each model;          fname.s### = 1  \n')
    output.write('- images at specified wavelengths;          fname.i### = 1  \n')
    output.write('     number of wavelengths = 5  \n')
    output.write('     wavelengths = 3.5, 4.5, 6.0, 8.0, 24.0 micron  \n')
    output.write('- radial profiles for each model;           fname.r### = 1  \n')
    output.write('- detailed run-time messages;               fname.m### = 1  \n')
    output.write('- visibility function at spec. wavelengths; fname.v### = 0  \n')
    output.close()


def plotdata(filename, color='black',ferrs = False):
    olam, olum, oerr = np.loadtxt(filename, usecols=(0,1,2), unpack=True)
    measured = (olum>0)
    llamuse = np.log10(olam[measured])
    llumuse = np.log10(olum[measured])
    lerruse = oerr[measured]/olum[measured]/np.log(10) # double check that this should be a natural log...
    if not ferrs:
        plt.errorbar(10**llamuse, 10**llumuse, yerr=oerr[measured], fmt='.', color=color, capsize=1.5, elinewidth=0.8, markeredgewidth=0.8, markersize=12)
    else :
         plt.errorbar(10**llamuse, 10**llumuse, yerr=10**(llumuse-0.5), fmt='.', color=color, capsize=1.5, elinewidth=0.8, markeredgewidth=0.8, markersize=12)
    limit = (olum<0)
    llamuse = np.log10(olam[limit])
    llumuse = np.log10(oerr[limit])


def run_dusty(pd,verbose=True):
    if verbose: print('calling dusty with tstar = %0.1f; tau = %0.2f; td = %0.1f; thick = %0.2f; E(B-V) = %0.4f' % (pd['tstar'],pd['tau'],pd['td'],pd['thick'],pd['ebv']))
    system('./dusty')

    ierror = 0
      #try:
    data = ascii.read('foo1.stb')
    lam = data['col1']
    flx = data['col2']
    npt = len(lam)

    i = 0
    with open('foo1.out','r') as f:
        for line in f:
            if i == 42:
                # print 'line read in from foo1.out:', line
                line_s = line.split()
                id = int(line_s[0])
                tau0 = float(line_s[1])
                f1 = float(line_s[2])
                r1 = float(line_s[3])
                r1torstar = float(line_s[4])
                theta1 = float(line_s[5])
                tdout = float(line_s[6])
                break
            i += 1   
    #except:
    #  ierror = 1
    #  lam = np.nan
    #  flx = np.nan
    #  npt = np.nan
    #  r1 = np.nan
    return lam, flx, npt, r1, ierror

def mol_abs_fractions(T,filters,abs_path='/scr2/viraj/Gattini/dusty_modeling/mol_abs'):
    print('Interpolating to calculate molecular absorption fraction')
    abs_fracs = []
    for f in filters:
        if os.path.exists('%s/%s_molecular_absorption_fractions.dat'%(abs_path,f)):
            fracs = ascii.read('%s/%s_molecular_absorption_fractions.dat'%(abs_path,f))
            abs_fracs.append(np.interp(T,fracs['col1'],fracs['col2']))
        else:
            print('Could not find file for filter %s. Assuming no molecular absorption.'%(f))
            abs_fracs.append(1)
    return np.array(abs_fracs)

object_photometry_file = sys.argv[1]#'lums_AOHer_ebv0.017_nodisterr_thickvary_1.dat'
chains_filename = object_photometry_file + '_results.dat'
colnames = ['tstar','tdust','tau','lstar','r1','lnl']

t1 = ascii.read('proc_0/%s'%(chains_filename),data_start=4,names=colnames)
t2= ascii.read('proc_1/%s'%(chains_filename),data_start=4,names=colnames)
t3 = ascii.read('proc_2/%s'%(chains_filename),data_start=4,names=colnames)
t4 = ascii.read('proc_3/%s'%(chains_filename),data_start=4,names=colnames)
print(len(t1),len(t2),len(t3),len(t4))
t =vstack([t1,t2,t3,t4])
print('Writing full results file of length %s'%(len(t)))
#t.write('%s_full_results.dat'%(chains_filename.split('.dat')[0]))
t['r1'] = np.log10(t['r1']) + 0.5*(t['lstar']-4.0)

labels = [r'T$_{*}$ [K]',r'T$_{d}$ [K]',r'$\tau$']#,r'thick'],'log($L_{*}/L_{\odot}$)','log(r1/cm)']
s = np.array([t[x] for x in t.columns[:len(labels)]])
s = s.transpose()

fig = corner.corner(s,labels=labels,show_titles=True, title_kwargs={"fontsize": 12})
plt.savefig('%s.pdf'%(chains_filename.split('.dat')[0]))

#Best-fit
x = t[np.argmax(t['lnl'])]
pdict = {}
x = t[np.argmax(t['lnl'])]
pdict['tstar'] = x['tstar']
pdict['td'] = x['tdust']
pdict['tau'] = x['tau']
pdict['ebv'] = 0.0
pdict['thick'] = 10#x['thick']#1000
pdict['fileuse'] = ''
pdict['idtype'] = 0
custom_grain_distribution = False
geninput(pdict)
lam, flx, npt, r1, ierror = run_dusty(pdict)
filters_names = np.array(['B','g','r','i','z','y','J','H','K','W1','W2','W3','W4'])
abs_fracs = mol_abs_fractions(pdict['tstar'],filters_names)
filt_lams = np.array([0.43,0.48,0.63,0.76,0.87,0.96,1.23,1.66,2.16,3.4,4.6,12,22])
abs_fracs_interp = np.interp(lam,filt_lams,abs_fracs)

bfit = Table()
bfit.add_column(Column(name='lam',data=lam))
bfit.add_column(Column(name='flx',data=flx*abs_fracs_interp))
bfit.add_column(Column(name='scaled',data=10**x['lstar']*flx*abs_fracs_interp))
bfit.write('bestfit_%s.csv'%(object_photometry_file))
plt.figure()
plotdata(object_photometry_file,color='black')
plt.plot(lam,(10**x['lstar'])*flx*abs_fracs_interp,c='blue',linestyle='-',alpha=0.7)
plt.xlim(0.2,25)
plt.ylim(0.001,1e3)
plt.xscale('log')
plt.yscale('log')
plt.ylabel(r'$\lambda$L$_{\lambda}$ [L$_{\odot}$]')
plt.xlabel(r'$\lambda$')
plt.savefig('%s_sed_bestfit.pdf'%(chains_filename.split('.dat')[0]))
