#!/usr/bin/env python

# Scott Adams (OSU) November 2011

# weighted least-squares
def wlsq(x,y,sigma=''):
  """weighted least-squares 
       inputs: x, y, (sigma)
       outputs: m, sigma_m, b, sigma_b"""
  import numpy as np
  # if no y-errors given assume uniform uncertainties
  if sigma == '': sigma = np.ones(len(x))
  # make sure data is formatted as arrays
  x = np.array(x)
  y = np.array(y)
  sigma = np.array(sigma)
  # calculate sums
  sumx = np.sum(x/sigma**2.)
  sumy = np.sum(y/sigma**2.)
  sumxx = np.sum(x*x/sigma**2.)
  sumxy = np.sum(x*y/sigma**2.)
  delta = np.sum(1./sigma**2.)*np.sum(sumxx) - sumx**2.
  # find slope and intercept
  m = ( np.sum(1./sigma**2) * sumxy - sumx*sumy ) / delta
  b = (sumxx*sumy - sumx*sumxy) / delta
  # find uncertainties in fit
  sigma_m = np.sqrt( np.sum(1./sigma**2) / delta )
  sigma_b = np.sqrt( sumxx / delta )
  return(m,sigma_m,b,sigma_b)

# Calculate linear-least squares fit
def lsq(x,y):
  """Inputs: x, y\nOutputs: m, b, r2, chi^2"""
  # initialize sums
  sumx = 0.
  sumy = 0.
  sumxx = 0.
  sumxy = 0.
  # calculate sums
  for i in range(len(x)):
    sumx = sumx + x[i]
    sumy = sumy + y[i]
    sumxx = sumxx + x[i]*x[i]
    sumxy = sumxy + x[i]*y[i]
  n = float(len(x))
  xave = sumx / n
  yave = sumy / n
  # find slope and intercept
  m = ( n*sumxy - sumx*sumy ) / ( n*sumxx - sumx*sumx )
  b = (sumxx*sumy - sumx*sumxy)/ ( n*sumxx - sumx*sumx )
  # find R^2
  sstot = 0
  ssreg = 0
  chisq = 0
  for i in range(len(x)):
    sstot = sstot + ( y[i] - yave )**2.
    ssreg = ssreg + ( y[i] - (m*x[i]+b))**2.
    chisq = chisq + ( y[i] - b - m*x[i] )**2.
  if sstot == 0:
    r2 = 1
  else:
    r2 = ssreg / sstot
  return(m,b,r2,chisq)

def lsq_array(x,y):
  print 'mylib.py lsq_array: this might not be working correctly'
  exit(1)
  from numpy import sum, average
  # initialize sums
  sumx = 0.
  sumy = 0.
  sumxx = 0.
  sumxy = 0.
  # calculate sums
  sumx = sum(x)
  sumy = sum(y)
  sumxx = sum(x*x)
  xumxy = sum(x*y)
  xave = average(x)
  yave = average(y)
  n = len(x)
  # find slope and intercept
  m = ( n*sumxy - sumx*sumy ) / ( n*sumxx - sumx*sumx )
  b = (sumxx*sumy - sumx*sumxy)/ ( n*sumxx - sumx*sumx )
  # find R^2
  sstot = 0
  ssreg = 0
  chisq = 0
  sstot = sum((y-yave)**2.)
  ssreg = sum((y-(m*x+b))**2.)
  chisq = sum((y-b-m*x)**2.)
  r2 = ssreg / sstot
  return(m,b,r2,chisq)

#from numpy import loadtxt
#x,y = loadtxt('test.dat', usecols=(0,1), unpack=True)
#m,b,chisq = lsq(x,y)
#print 'y = '+str(m)+'x + '+str(b)
#print 'chi^2 = ', chisq

# Routine for deleting old files
def del_old(filename):
  import os
  if os.path.isfile(filename) == 1:
    #print 'deleting old ' + filename
    os.remove(filename)

# Module for sigma-clipping
def sigma_clip(x,sigma,maxiter,verbose=0):
  import numpy as np
  x = np.array(x)
  ncut = 10.
  niter = 0
  while ncut > 0 and niter < maxiter:
    std1 = np.std(x)
    med1 = np.median(x)
    length1 = len(x)
    mask = (x>med1-sigma*std1)
    mask = mask*(x<med1+sigma*std1)
    x1 = x[mask]
    xcut = x[np.invert(mask)]
    length2 = len(x1)
    ncut = length1-length2
    niter += 1
    x = x1
    if verbose:
      std2 = np.std(x)
      med2 = np.std(x)
      print med1, std1, med2, std2, length1, length2, niter
  return x

# Module for sigma-clipping array mask
def sigma_clip_array(x,sigma,maxiter,verbose):
  import numpy as np
  import math
  x = np.array(x)
  ncut = 10.
  niter = 0
  mask = (x==x)
  while ncut > 0 and niter < maxiter:
    std1 = np.std(x[mask])
    med1 = np.median(x[mask])
    length1 = len(x[mask])
    mask = mask*(x>med1-sigma*std1)
    mask = mask*(x<med1+sigma*std1)
    length2 = len(x[mask])
    ncut = length1-length2
    niter += 1
    if verbose:
      std2 = np.std(x[mask])
      med2 = np.std(x[mask])
      print med1, std1, med2, std2, length1, length2, niter
  return mask

def sigma_clip_lsq(x,y,sigma,maxiter,verbose):
  import numpy as np
  x = np.array(x)
  y = np.array(y)
  ncut = 10.
  niter = 0
  while ncut > 0 and niter < maxiter:
    m,b,r2,chisq = lsq(x,y)
    if verbose: print 'y =', m, '* x +', b, '; r2 =', r2, '; chisq =', chisq
    yfit = m*x+b
    diff = y-yfit
    std1 = np.std(diff)
    med1 = np.median(diff)
    length1 = len(diff)
    mask = (abs(diff)<=sigma*std1)
    if verbose:
      cutmask = np.invert(mask)
      xcut = x[cutmask]
      ycut = y[cutmask]
    x,y,diff = x[mask],y[mask],diff[mask]
    length2 = len(diff)
    ncut = length1-length2
    niter += 1
    if verbose:
      med2 = np.median(diff)
      std2 = np.std(diff)
      print med1, std1, med2, std2, length1, length2, niter
      print xcut, ycut
  return x, y

# Routine for calculating time since Big Bang for a given redshift
#  for Concordance Model
def cosmotime(z):
  from math import sqrt, log
  omega_m0 = 0.3
  omega_lambda0 = 0.7
  #H0 = 2.2977E-18  # [s^-1] -> equivalent to 71km/s/Mpc
  H0 = 0.072185 # [Gyr^-1] -> equivalent to 71km/s/Mpc
  a = 1. / (1.+z)
  a_ml = (omega_m0/omega_lambda0)**(1./3.)
  t = 2./(3.*H0*sqrt(1.-omega_m0)) * log( (a/a_ml)**1.5 + sqrt( 1+(a/a_ml)**3. ) )
  t0 = 2./(3.*H0*sqrt(1.-omega_m0)) * log( (1./a_ml)**1.5 + sqrt( 1+(1./a_ml)**3. ) )
  lookback = t0-t
  return(round(t,3),round(lookback,3))

# module to add a keword and value to a file
#   updating the value if the keyword exists
def update_file(filename, keyword, value):
  data = []
  keyword_added = 0
  for line in file(filename):
    line = line.split()
    if line[0] == keyword:
      line[1] = value
      keyword_added = 1
    sline = ''
    for i in xrange(len(line)):
      sline += str(line[i])+'\t'
    data.append(sline)
  if keyword_added==0:
    data.append(keyword+'\t'+str(value))
  outfile = open(filename,'w')
  for item in data:
    outfile.write(str(item)+'\n')
  outfile.close()

def update_optfile(filename, keyword, value):
  data = []
  keyword_added = 0
  for line in file(filename):
    line = line.split('=')
    if line[0] == keyword:
      line[1] = value+'\n'
      keyword_added = 1
    sline = line[0]+'='+line[1]
    #for i in xrange(len(line)):
    #  sline += str(line[i])+'\t'
    data.append(sline)
  #if keyword_added==0:
  #  data.append(keyword+'\t'+str(value))
  outfile = open(filename,'w')
  for item in data:
    outfile.write(str(item))
  outfile.close()

def read_keyword(filename, keyword):
  from sys import exit
  keyword_found = 0
  for line in file(filename):
    line = line.split()
    if len(line)==0: continue
    if line[0] == keyword:
      keyword_found = line[1]
  if keyword_found == 0:
    print 'ERROR: '+keyword+' not found'
    exit(1)
  return keyword_found  

def return_ascii_starting_header(filename, comment_marker):
  header = []
  for line in file(filename):
    if line.startswith('#'):
      header.append(line)
    else:
      return header
  return header

def read_header(header,field):
  """for HEADER with content:
  FIELD1 = VALUE1 % DESCRIPTION1
  FIELD2 = VALUE2 % DESCRIPTION2
  read_header(HEADER, FIELD1) -> VALUE1"""
  from sys import exit
  for i in xrange(len(header)):
    line = header[i].strip('# ')
    sline = line.split(' = ')
    if sline[0] == field:
      trimmed_line = sline[1].strip('\n')
      return trimmed_line.split('%')[0]
  print "ERROR: '%s' not found in header" % (field)
  exit(1)


def progress_bar(n,p,nl):
  from sys import stdout
  n += 1
  p += 1
  if float(p)/float(nl) > 0.01:
    p = 0
    stdout.write('\r'+str(int(float(n)/float(nl)*100))+'% complete')
    stdout.flush()
  return n, p

def trim_dictionary(dict,mask):
  new_dict = {}
  for item in dict:
    new_dict[item] = dict[item][mask]
  return new_dict

def append_dictionary_orig(dict,addition):
  for item in addition:
    for i in xrange(len(addition[item])):
      if len(dict[item])==0:
        dict[item] = addition[item][i]
      else:
        dict[item].append(addition[item][i])
  return dict

def append_dictionary(dict,addition):
  for item in addition:
      if len(dict[item])==0:
        dict[item] = addition[item]
      else:
        dict[item].append(addition[item])
  return dict

def dictionary_line(dict,i):
  new_dict = {}
  for item in dict:
    new_dict[item] = dict[item][i]
  return new_dict

def append_line(dict,addition,i):
  for item in addition:
    try:
      dict[item].append(addition[item][i])
    except:
      dict[item] = []
      dict[item].append(addition[item][i])
  return dict

# verify that the field is listed in the dictionary
def check_dictionary(dict,field):
  names = dict.dtype.names
  for name in names:
    if name == field:
      return 1
  return 0

__all__ = ['lsq','lsq_array','del_old','cosmotime','sigma_clip','sigma_clip_array','update_file','read_keyword','progress_bar','trim_dictionary','check_dictionary','return_ascii_starting_header']
