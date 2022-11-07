#!/usr/bin/env python

"""Convert ACS/WFC Throughput fits files from http://www.stsci.edu/hst/acs/analysis/reference_files/synphot_tables.html into the format expected by smooth2.pl"""

import argparse
#
from astropy.io import fits

parser = argparse.ArgumentParser()
parser.add_argument("input", help="Name on input file")
parser.add_argument("output", help="Name of output file")
args = parser.parse_args()

data = fits.open(args.input)
table = data[1].data

output = open(args.output, 'w')
for i in xrange(len(table['WAVELENGTH'])):
  output.write('%s %s 0\n' % (table['WAVELENGTH'][i], table['THROUGHPUT'][i]))
output.close()
