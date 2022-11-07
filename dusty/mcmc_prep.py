#!/usr/bin/env python

import os
import argparse
from sys import exit

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Create and prepare directory to run dusty_wrapper.py")
parser.add_argument("directory", help="directory to create")
parser.add_argument("magfile", help="file containing luminosity constraints")
parser.add_argument("--path", help="relative path to dusty_wrapper.py", default='/scr/sma/software/dusty/')
args = parser.parse_args()

#if args.path:
#  path = args.path
#else:
#  path = ''
path = args.path
if os.path.isdir(args.directory):
  print 'ERROR: directory already exists'
  exit(-1)
os.mkdir(args.directory)
os.system('cp '+path+'dusty '+args.directory)
os.system('cp '+path+'dusty.inp '+str(args.directory))
os.system('cp '+path+'lambda_grid.dat '+str(args.directory))
os.system('cp '+str(args.magfile)+' '+str(args.directory))
os.system('cp '+path+'dusty_wrapper.py '+str(args.directory))
os.system('cp '+path+'plot_sed.py '+str(args.directory))
