import os

filelist = ['lums_AOHer_ebv0.017_nodisterr_thick10_pc10.dat','lums_ASAS-RCB-21_ebv1.0_nodisterr_thick10_pc10.dat','lums_NSV11154_ebv0.07_nodisterr_thick10_pc10.dat','lums_WISE_J175749.76-075314.9_ebv1.18_nodisterr_thick10_pc10.dat','lums_WISE_ToI_1213_ebv0.4_nodisterr_thick10_pc10.dat','lums_WISE_ToI_290_ebv0.09_nodisterr_thick10_pc10.dat']
#['lums_NSV11154_ebv0.07_nodisterr.dat','lums_AOHer_ebv0.017_nodisterr.dat','lums_ASAS-RCB-21_ebv1.0_nodisterr.dat','lums_WISE_ToI_290_ebv0.09_nodisterr.dat','lums_WISE_J175749.76-075314.9_ebv1.18_nodisterr.dat','lums_WISE_ToI_1213_ebv0.4_nodisterr.dat']

for f in filelist:
    print('Currently doing %s'%(f))
    #os.system('/scr2/viraj/anaconda3/bin/python dusty_emcee_parallel_molabs.py %s'%(f))
    os.system('/usr/local/anaconda2/bin/python plotting_parallel.py %s'%(f))
    #os.system('cp %s %s_thick10.dat'%(f,f.split('.dat')[0]))
