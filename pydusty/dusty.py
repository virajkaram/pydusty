import os
import logging
import numpy as np
from pydusty.parameters import DustyParameters
from astropy.io import ascii
from glob import glob
import subprocess

logger = logging.getLogger(__name__)


class BaseDusty:

    def __init__(self,
                 parameters: DustyParameters,
                 dusty_working_directory: str,
                 file_basename='foo1',
                 dusty_file_directory='data/dusty_files',
                 ):
        self.parameters = parameters
        self.file_basename = file_basename
        self.dusty_working_directory = dusty_working_directory
        self.dusty_file_directory = dusty_file_directory

        self.setup_dusty_dir()

    def setup_dusty_dir(self):
        if not os.path.exists(self.dusty_working_directory):
            os.makedirs(self.dusty_working_directory)

        dusty_required_files = ['dusty', 'dusty.inp', 'dusty.f', 'lambda_grid.dat']
        for req_file in dusty_required_files:
            filepath = os.path.join(self.dusty_file_directory, req_file)
            if not os.path.exists(filepath):
                err = f'Required dusty file : {filepath} not found in dusty_file_directory {self.dusty_file_directory}'
                raise ValueError(err)

        dusty_files = glob(f'{self.dusty_file_directory}/*')
        logger.info(f'Copying following files from {self.dusty_file_directory} '
                    f'to {self.dusty_working_directory}')
        logger.info(f'{dusty_files}')
        for filepath in dusty_files:
            os.system(f'cp {filepath} {self.dusty_working_directory}/')

        logger.info(f'Copied dusty files to {self.dusty_working_directory}')

    def generate_input(self):
        raise NotImplementedError

    def run(self):
        # if verbose:
        #     pd = {}
        #     for n in range(len(varxnames)):
        #         pd[varxnames[n]] = varx[n]
        #     for n in range(len(constxnames)):
        #         pd[constxnames[n]] = constx[n]
        #     print('calling dusty with tstar = %0.1f; tau = %0.2f; td = %0.1f; thick = %0.2f; E(B-V) = %0.4f' % (
        #         pd['tstar'], pd['tau'], pd['td'], pd['thick'], pd['ebv']))
        curdir = os.getcwd()
        try:
            os.chdir(self.dusty_working_directory)
            logger.info(f"Executing command {self.dusty_working_directory}/dusty")
            status = subprocess.run(f"{self.dusty_working_directory}/dusty", check=True, capture_output=True).stdout.decode()
        except subprocess.CalledProcessError as e:
            logger.error(e)
            raise e
        os.chdir(curdir)

    def get_results(self):
        ierror = 0
        # try:
        data = ascii.read(f'{self.file_basename}.stb')
        lam = data['col1']
        flx = data['col2']
        npt = len(lam)

        i = 0
        with open(f'{self.file_basename}.out', 'r') as f:
            lines = f.readlines()
        lines = np.array(lines)
        mask = np.array(['RESULTS:' in x for x in lines])
        result_ind = np.where(mask)[0][0]
        result_line = lines[result_ind+5]

        line_s = result_line.split()
        id = int(line_s[0])
        tau0 = float(line_s[1])
        f1 = float(line_s[2])
        r1 = float(line_s[3])
        r1torstar = float(line_s[4])
        theta1 = float(line_s[5])
        tdout = float(line_s[6])

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
        if self.parameters.dust_type.value == 'graphite':
            output.write('    x = 0.00    0.00   0.00    1.00    0.00    0.00 \n')
        elif self.parameters.dust_type.value == 'silicate':
            output.write('    x = 0.00    0.00   1.00    0.00    0.00    0.00 \n')
        elif self.parameters.dust_type.value == 'amorphous_carbon':
            output.write('    x = 0.00    0.00   0.00    0.00    1.00    0.00 \n')
        elif self.parameters.dust_type.value == 'silicate_carbide':
            output.write('    x = 0.00    0.00   0.00    0.00    0.00    1.00 \n')
        else:
            raise ValueError(f'Unknown dust type {self.parameters.dust_type.value}')
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
        output.write(f'- lambda0 = {self.parameters.tau_wavelength_microns.value} micron          % optical depth specified  \n')
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


class Dusty_Alumina_SilOW(BaseDusty):

    def generate_input(self):
        assert self.parameters.al_com_abundance is not None, "Alumina abundance must be specified"
        output = open(f'{self.file_basename}.inp', 'w')

        if (self.parameters.tstar.value < self.parameters.tstarmin.value or self.parameters.tstar.value > self.parameters.tstarmax.value) or (
        self.parameters.blackbody.value):
            output.write('Spectrum = 1\n')
            output.write('Number of BB = 1\n')
            output.write(f'Temperature = {round(self.parameters.tstar.value, 2)}\n')
        else:
            output.write('Spectrum = 5   \n')
            output.write('\n')
        output.write('   optical properties index = 2 \n')
        output.write('Abundances for supported grain types:\n')
        output.write('   #   Sil-Ow  Sil-Oc  Sil-DL  grf-DL  amC-Hn  SiC-Pg \n')
        output.write(f'    x = {1 - self.parameters.al_com_abundance.value}    0.00   0.00    0.00    0.00    0.00 \n')
        output.write('Number of additional components = 1, properties listed in:\n')
        output.write('Al2O3-comp.nk\n')
        output.write(f'Abundances for these components = {self.parameters.al_com_abundance.value} \n')
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
        output.write(f'- lambda0 = {self.parameters.tau_wavelength_microns.value} micron          % optical depth specified  \n')
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


class Dusty_Alumina_SilDL(BaseDusty):

    def generate_input(self):
        assert self.parameters.al_com_abundance is not None, "Alumina abundance must be specified"
        output = open(f'{self.file_basename}.inp', 'w')

        if (self.parameters.tstar.value < self.parameters.tstarmin.value or self.parameters.tstar.value > self.parameters.tstarmax.value) or (
        self.parameters.blackbody.value):
            output.write('Spectrum = 1\n')
            output.write('Number of BB = 1\n')
            output.write(f'Temperature = {round(self.parameters.tstar.value, 2)}\n')
        else:
            output.write('Spectrum = 5   \n')
            output.write('\n')
        output.write('   optical properties index = 2 \n')
        output.write('Abundances for supported grain types:\n')
        output.write('   #   Sil-Ow  Sil-Oc  Sil-DL  grf-DL  amC-Hn  SiC-Pg \n')
        output.write(f'    x = 0.00    0.00   {1 - self.parameters.al_com_abundance.value}    0.00    0.00    0.00 \n')
        output.write('Number of additional components = 1, properties listed in:\n')
        output.write('Al2O3-comp.nk\n')
        output.write(f'Abundances for these components = {self.parameters.al_com_abundance.value} \n')
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
        output.write(f'- lambda0 = {self.parameters.tau_wavelength_microns.value} micron          % optical depth specified  \n')
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


class Dusty_Two_Component_Silicate(BaseDusty):

    def generate_input(self):
        assert self.parameters.al_com_abundance is not None, "Alumina abundance must be specified"
        output = open(f'{self.file_basename}.inp', 'w')

        if self.parameters.al_com_abundance.name == 'alumina':
            dust_comp_file = 'Al2O3-comp.nk'
        elif self.parameters.al_com_abundance.name == 'crystaline_silicate':
            dust_comp_file = 'crMgFeSil.nk'
        elif self.parameters.al_com_abundance.name == 'glassy_olivine':
            dust_comp_file = 'gloliMg50.nk'
        elif self.parameters.al_com_abundance.name == 'glassy_pyroxene':
            dust_comp_file = 'glpyrMg50.nk'
        elif self.parameters.al_com_abundance.name == 'feo':
            dust_comp_file = 'FeO.nk'
        else:
            raise ValueError(f'Unknown dust component {self.parameters.al_com_abundance.name}')
        if (self.parameters.tstar.value < self.parameters.tstarmin.value or self.parameters.tstar.value > self.parameters.tstarmax.value) or (
        self.parameters.blackbody.value):
            output.write('Spectrum = 1\n')
            output.write('Number of BB = 1\n')
            output.write(f'Temperature = {round(self.parameters.tstar.value, 2)}\n')
        else:
            output.write('Spectrum = 5   \n')
            output.write('\n')
        output.write('   optical properties index = 2 \n')
        output.write('Abundances for supported grain types:\n')
        output.write('   #   Sil-Ow  Sil-Oc  Sil-DL  grf-DL  amC-Hn  SiC-Pg \n')
        output.write(f'    x = 0.00    0.00   {1 - self.parameters.al_com_abundance.value}    0.00    0.00    0.00 \n')
        output.write('Number of additional components = 1, properties listed in:\n')
        output.write(f'{dust_comp_file}\n')
        output.write(f'Abundances for these components = {self.parameters.al_com_abundance.value} \n')
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
        output.write(f'- lambda0 = {self.parameters.tau_wavelength_microns.value} micron          % optical depth specified  \n')
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


external_elements_filenames_dict = {'alumina': 'Al2O3-comp.nk',
                                    'crystalline_silicate': 'crMgFeSil.nk',
                                    'glassy_olivine': 'gloliMg50.nk',
                                    'glassy_pyroxene': 'glpyrMg50.nk',
                                    'feo': 'FeO.nk',
                                    'glassy_silicate': 'glSil.nk'
                                    }


class Dusty_Multi_Composition(BaseDusty):

    def generate_input(self):
        inbuilt_abundances = {}
        for elem in ['silow', 'siloc', 'sildl', 'grfdl', 'amchn', 'sicpg']:
            inbuilt_abundances[elem] = 0
            if elem in self.parameters.abundances_dict.keys():
                inbuilt_abundances[elem] = self.parameters.abundances_dict[elem]

        external_abundances = {}
        for elem in self.parameters.abundances_dict.keys():
            if elem not in inbuilt_abundances.keys():
                external_abundances[elem] = self.parameters.abundances_dict[elem]

        output = open(f'{self.file_basename}.inp', 'w')

        if (self.parameters.tstar.value < self.parameters.tstarmin.value or self.parameters.tstar.value > self.parameters.tstarmax.value) or (
        self.parameters.blackbody.value):
            output.write('Spectrum = 1\n')
            output.write('Number of BB = 1\n')
            output.write(f'Temperature = {round(self.parameters.tstar.value, 2)}\n')
        else:
            output.write('Spectrum = 5   \n')
            output.write('\n')
        output.write('   optical properties index = 2 \n')
        output.write('Abundances for supported grain types:\n')
        output.write('   #   Sil-Ow  Sil-Oc  Sil-DL  grf-DL  amC-Hn  SiC-Pg \n')
        output.write(f"    x = {inbuilt_abundances['silow']}    {inbuilt_abundances['siloc']}   {inbuilt_abundances['sildl']}    {inbuilt_abundances['grfdl']}    {inbuilt_abundances['amchn']}    {inbuilt_abundances['sicpg']} \n")
        if len(external_abundances) > 0:
            output.write(f'Number of additional components = {len(external_abundances)}, properties listed in:\n')
            for external_elem in external_abundances.keys():
                dust_comp_file = external_elements_filenames_dict[external_elem]
                output.write(f'{dust_comp_file}\n')
            output.write(f'''Abundances for these components = {", ".join([str(x) for x in external_abundances.values()])} \n''')
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
        output.write(f'- lambda0 = {self.parameters.tau_wavelength_microns.value} micron          % optical depth specified  \n')
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


class DustyWaterVapor(BaseDusty):

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
        output.write('   optical properties index = 3; cross-sections entered in file \n')
        output.write('   water_vapor_crosssec_dusty.dat\n')
        output.write(f'- temperature = {self.parameters.tdust.value} K \n')
        output.write('- density type = 1                   \n')
        output.write('- number of powers = 1              \n')
        output.write(f'- shells relative thickness = {self.parameters.shell_thickness.value}\n')
        output.write('- power = 2 \n')
        output.write('- grid type = 1                  % linear grid \n')
        output.write(f'- lambda0 = {self.parameters.tau_wavelength_microns.value} micron          % optical depth specified  \n')
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


class DustyCustomInputSpectrum(BaseDusty):
    """
    Custom input spectrum for dusty, requires custom_input_spectrum_file parameter to be specified
    with the file containing the spectrum in lambda vs F_lambda format.
    """

    def generate_input(self):

        if self.parameters.custom_input_spectrum_file is None:
            raise ValueError('custom_input_spectrum_file parameter must be specified'
                             'for DustyCustomInputSpectrum class')

        if self.parameters.blackbody.value:
            raise ValueError('Blackbody spectrum not supported for custom input spectrum')

        output = open(f'{self.file_basename}.inp', 'w')

        output.write('Spectrum = 5   \n')
        output.write(f'{self.parameters.custom_input_spectrum_file.value}\n')
        output.write('   optical properties index = 1 \n')
        output.write('   #   Sil-Ow  Sil-Oc  Sil-DL  grf-DL  amC-Hn  SiC-Pg \n')
        if self.parameters.dust_type.value == 'graphite':
            output.write('    x = 0.00    0.00   0.00    1.00    0.00    0.00 \n')
        elif self.parameters.dust_type.value == 'silicate':
            output.write('    x = 0.00    0.00   1.00    0.00    0.00    0.00 \n')
        elif self.parameters.dust_type.value == 'amorphous_carbon':
            output.write('    x = 0.00    0.00   0.00    0.00    1.00    0.00 \n')
        elif self.parameters.dust_type.value == 'silicate_carbide':
            output.write('    x = 0.00    0.00   0.00    0.00    0.00    1.00 \n')
        else:
            raise ValueError(f'Unknown dust type {self.parameters.dust_type.value}')
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
        output.write(f'- lambda0 = {self.parameters.tau_wavelength_microns.value} micron          % optical depth specified  \n')
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