import os
import logging
from pydusty.parameters import DustyParameters
from astropy.io import ascii

logger = logging.getLogger(__name__)


class BaseDusty:

    def __init__(self,
                 parameters: DustyParameters,
                 dusty_working_directory: str,
                 file_basename='foo1',
                 ):
        self.parameters = parameters
        self.file_basename = file_basename
        self.dusty_working_directory = dusty_working_directory
        self.setup_dusty_dir()

    def setup_dusty_dir(self):
        if not os.path.exists(self.dusty_working_directory):
            os.makedirs(self.dusty_working_directory)
        os.system(f'cp dusty_files/dusty {self.dusty_working_directory}/')
        os.system(f'cp dusty.inp {self.dusty_working_directory}/')
        os.system(f'cp dusty.f {self.dusty_working_directory}/')
        os.system(f'cp lambda_grid.dat {self.dusty_working_directory}/')
        logger.info(f'Copied dusty files to {self.dusty_working_directory}')

    def generate_input(self):
        pass

    def run(self):
        # if verbose:
        #     pd = {}
        #     for n in range(len(varxnames)):
        #         pd[varxnames[n]] = varx[n]
        #     for n in range(len(constxnames)):
        #         pd[constxnames[n]] = constx[n]
        #     print('calling dusty with tstar = %0.1f; tau = %0.2f; td = %0.1f; thick = %0.2f; E(B-V) = %0.4f' % (
        #         pd['tstar'], pd['tau'], pd['td'], pd['thick'], pd['ebv']))
        os.system('./dusty')

    def get_results(self):
        ierror = 0
        # try:
        data = ascii.read(f'{self.file_basename}.stb')
        lam = data['col1']
        flx = data['col2']
        npt = len(lam)

        i = 0
        with open(f'{self.file_basename}.out', 'r') as f:
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
                # except:
        #  ierror = 1
        #  lam = np.nan
        #  flx = np.nan
        #  npt = np.nan
        #  r1 = np.nan
        return lam, flx, npt, r1, ierror


class Dusty(BaseDusty):

    def generate_input(self):

        output = open(f'{self.file_basename}.inp', 'w')

        if (self.parameters.tstar.value < self.parameters.tstarmin.value or self.parameters.tstar.value > self.parameters.tstarmax.value) or (
        self.parameters.blackbody.value):
            output.write('Spectrum = 1\n')
            output.write('Number of BB = 1\n')
            output.write(f'Temperature = {round(self.parameters.tstar.value, 2)}\n')
        else:
            output.write('Spectrum = 5   \n')
            output.write('\n')
        output.write('   optical properties index = 1 \n')
        output.write('   #   Sil-Ow  Sil-Oc  Sil-DL  grf-DL  amC-Hn  SiC-Pg \n')
        if self.parameters.dust_type.value == 0:
            output.write('    x = 0.00    0.00   0.00    1.00    0.00    0.00 \n')
        else:
            output.write('    x = 0.00    0.00   1.00    0.00    0.00    0.00 \n')
        if self.parameters.custom_grain_distribution.value:
            output.write('- size distribution = 2  % custom       \n')
            output.write(f'  q = 3.5, a(min) = {self.parameters.min_grain_size.value} micron, a(max) = {self.parameters.max_grain_size.value} micron\n')
        else:
            output.write('- size distribution = 1  % standard MRN    \n')
        output.write(f'- temperature = {self.parameters.tdust.value} K \n')
        output.write('- density type = 1                   \n')
        output.write('- number of powers = 1              \n')
        output.write(f'- shells relative thickness = {self.parameters.shell_thickness.value}\n')
        output.write('- power = 2 \n')
        output.write('- grid type = 1                  % linear grid \n')
        output.write('- lambda0 = 0.55 micron          % optical depth specified  \n')
        output.write('- tau(min) = ' + str(
            self.parameters.tau.value) + ' ; tau(max) = 1000.0   % for the visual wavelength \n')
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
        logger.info(f"Writing dusty file with {self.parameters.get_printable_string()}")
