#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 12:25:19 2022

@author: elisa
"""

import numpy as np
import scipy as sp
import os
import sys
from pytictoc import TicToc
import subprocess
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import fileinput


from astropy.modeling.physical_models import NFW
from astropy import constants as const
from astropy.cosmology import Planck15
from astropy.wcs import WCS
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename

from scipy.special import erf, erfc, erfcinv, erfinv
from scipy import interpolate
from scipy.stats import linregress
from scipy.stats import pearsonr
from scipy.optimize import curve_fit
from scipy import stats
from numpy import genfromtxt
from scipy.signal import find_peaks, peak_widths, peak_prominences, hilbert



plt.rc('font', **{'family': "sans-serif"})
params = {'text.latex.preamble': [r'\usepackage{siunitx}', r'\usepackage{amsmath}']}
plt.rcParams['text.usetex'] = True
plt.rcParams.update(params)

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

# np.set_printoptions(suppress=True)
arcmin = 1. / 60 * np.pi / 180 # arcmin to radians

timer = TicToc()
pi = np.pi

n_fwhm = 157
n_spectrum = 3721
skip_wl_every = 1

N_g_sqrtrho = 51
N_Sflats = 21
N_r0s = 10
Msun = 1.9884e+30 * 10 ** 3  # grams
cm = 3.24078e-22  # kpc

log_g_sqrtrho = np.linspace(-15.0, -10.0, num=N_g_sqrtrho, endpoint=True) # This has to match ConstParams
log_Sflats = np.linspace(-1, 1, num=N_Sflats, endpoint=True)

log_r0Dict = { # This has to match ConstParams
    'leo': np.arange(-0.4, 1.16, 0.173),
    'eri1': np.arange(-0.35, 0.51, 0.095),
    'eri2': np.arange(-0.35, 0.51, 0.095),
    'eri': np.arange(-0.35, 0.51, 0.095),
    'gru': np.arange(-0.6, 1.1, 0.178),
    'scu': np.arange(-0.5, 0.5, 0.111),
    'hyd': np.arange(0.25, 2.1, 0.195)}

log_r0Dict_core = { # This has to match ConstParams_cored
    'leo': np.arange(-0.35, 1.17, 0.167),
    'eri1': np.arange(-0.35, 0.61, 0.106),
    'eri2': np.arange(-0.35, 0.61, 0.106),
    'eri': np.arange(-0.35, 0.61, 0.106),
    'gru': np.arange(-0.7, 1.21, 0.211),
    'scu': np.arange(-0.75, 0.9, 0.183),
    'hyd': np.arange(0.2, 2.01, 0.2)}

v_radial_dict = {
    'leo': 1.271e-4,
    'eri1': 2.540e-4,
    'eri2':  2.540e-4,
    'eri':  2.540e-4,
    'gru': -4.687e-4,
    'scu': 3.716e-4,
    'hyd': 1.011e-3}

log_g = np.arange(-13.6 , -10.5, 0.1)
N_gs = len(log_g)
Pcl = 0.95

dchan = 1.25
channel_centers = np.arange(4700, 4700 + dchan * n_spectrum, dchan)  # wavelengths of the center of the channels for all but LeoT. Angstrom
channel_centers_leo = np.arange(4699.57, 4700 + dchan * n_spectrum, dchan)# wavelengths of the center of the channels for LeoT. Angstrom
channelDict = {'leo': channel_centers_leo, 'eri2': channel_centers, 'eri3': channel_centers, 'gru': channel_centers,
               'eri1': channel_centers, 'scu': channel_centers, 'ant': channel_centers, 'hyd':channel_centers, 'eri': channel_centers,
               'leo_scu_gru_hyd_eri':channel_centers, 'eri_gru_hyd_leo_scu':channel_centers}

channel_centers_extended = np.arange(4700 - dchan * 10, 4700 + dchan * (n_spectrum + 10), dchan)  # wavelengths of the center of the channels for all but LeoT

fwhmDict = {'leo': 0.61, 'eri': 0.53, 'eri1': 0.53, 'eri2': 0.53, 'eri3': 0.53, 'gru': 0.67,  'scu': 0.50, 'hyd':0.40}
################################################################################################
# Files
################################################################################################

data_dir = 'data/'
image_dir = data_dir + 'images/'
chains_dir = 'chains/'

chi2_lists_dir = 'chi2_lists/'
chi2_null_dir = 'chi2_nulls/'
bounds_dir = 'bounds/'

# Input files
file_fwhm = data_dir + 'FWHM.txt'
file_spectrum = data_dir + '%s_spectrum_ave_ref.dat'

# Input files created by C program
file_chi2_null = 'chi2best.dat'
file_chi2_null_permanent = chi2_null_dir + 'chi2null_%s_%d.dat'
file_chi2_null_permanent_core = chi2_null_dir + 'chi2null_core_%s_%d.dat'
file_chi2_bound = 'chi2bound.dat'
file_chi2_list = 'chi2list_%s.dat'
file_chi2_list_permanent = chi2_lists_dir + 'chi2list_%s_%d.dat'
file_chi2_list_permanent_core = chi2_lists_dir + 'chi2list_core_%s_%d.dat'

# Output files
file_input_cube = 'inputfile_cube.dat'

file_bounds = bounds_dir + 'bounds_%s.dat'
file_bounds_doppler = bounds_dir + 'bounds_doppler_%s.dat'

file_bounds_cored = bounds_dir + 'bounds_cored_%s.dat'
file_bounds_doppler_cored = bounds_dir + 'bounds_doppler_cored_%s.dat'


file_evidence = bounds_dir + 'evidence_%s.dat'
file_evidence_no_peaks = bounds_dir + 'evidence_no_peaks_%s.dat'
file_evidence_no_peaks_LEE = bounds_dir + 'evidence_no_peaks_LEE_%s.dat'
file_evidence_no_peaks_doppler = bounds_dir + 'evidence_no_peaks_doppler_%s.dat'
file_evidence_no_peaks_LEE_doppler = bounds_dir + 'evidence_no_peaks_LEE_doppler_%s.dat'


file_evidence_cored = bounds_dir + 'evidence_cored_%s.dat'
file_evidence_no_peaks_cored = bounds_dir + 'evidence_no_peaks_cored_%s.dat'
file_evidence_no_peaks_LEE_cored = bounds_dir + 'evidence_no_peaks_LEE_cored_%s.dat'
file_evidence_no_peaks_doppler_cored = bounds_dir + 'evidence_no_peaks_doppler_cored_%s.dat'
file_evidence_no_peaks_LEE_doppler_cored = bounds_dir + 'evidence_no_peaks_LEE_doppler_cored_%s.dat'


file_chains = chains_dir + 'chain_%s.dat'
file_chains_npz = chains_dir + 'chain_%s.npz'
file_chains_npz_binned = chains_dir + 'chain_binned_%s.npz'
file_chi2_of_g = 'chi2_of_g/chi2_mins_%s_%d.npz'
file_chi2_of_g_doppler = 'chi2_of_g/chi2_mins_doppler_%s_%d.npz'
# file_chi2_of_g = 'chi2_of_g_old_angres/chi2_mins_%s_%d.npz'
# file_chi2_of_g_doppler = 'chi2_of_g_old_angres/chi2_mins_doppler_%s_%d.npz'


file_chains_core = chains_dir + 'chain_core_%s.dat'
file_chains_core_npz = chains_dir + 'chain_core_%s.npz'
file_chains_core_npz_binned = chains_dir + 'chain_core_binned_%s.npz'
file_chi2_of_g_core = 'chi2_of_g/chi2_mins_core_%s_%d.npz'
file_chi2_of_g_core_doppler = 'chi2_of_g/chi2_mins_doppler_core_%s_%d.npz'

################################################################################################
################################################################################################
################################################################################################


############################# HELPER FUNCTIONS ###########################################
def cleanup(name_dw):
    if os.path.exists(file_chi2_null):
        os.remove(file_chi2_null)
    
    if os.path.exists(file_chi2_bound):
        os.remove(file_chi2_bound)
    
    if os.path.exists(file_input_cube):
        os.remove(file_input_cube)
    
    if os.path.exists(file_bounds % name_dw):
        os.remove(file_bounds % name_dw)
    
    if os.path.exists('chi2disc.dat'):
        os.remove('chi2disc.dat')
    
    if os.path.exists('chi2list2D.dat'):
        os.remove('chi2list2D.dat')
    
    if os.path.exists('chi2list.dat'):
        os.remove('chi2list.dat')

def strip_line(line):
    line = line.strip()
    line = " ".join(line.split())
    return line.split(' ')

def write_mainh(name_dw, isCored):
    
    if isCored:
        string_mainh = """#ifndef MAIN_H
#define MAIN_H
#include \"ConstParams_cored_""" +name_dw+""".h\"
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <cstring>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include \"/Users/elisa/include/fitsio.h\"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_gamma.h>
#include <chrono>
#endif"""
    else:
        string_mainh = """#ifndef MAIN_H
#define MAIN_H
#include \"ConstParams_""" +name_dw+""".h\"
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <cstring>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include \"/Users/elisa/include/fitsio.h\"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_gamma.h>
#include <chrono>
#endif"""
    
    with open("main.h", "w") as fout:
        fout.write(string_mainh)

def write_makefile(isMakeEvidence):
    if isMakeEvidence:
        
        string_makefile = """# MAKEFILE
CFITSIO = /Users/elisa
FITSLIB = -L$(CFITSIO)/lib -lcfitsio -Wl,-rpath,$(CFITSIO)/lib
LIBS = `pkg-config --cflags --libs gsl`
LDFLAGS = $(FITSLIB) $(LIBS)\nrundSph: main_evidence.o dSph.o Components.o Tools.o
\tg++  main_evidence.o dSph.o Components.o Tools.o -o rundSph ${LDFLAGS}
main_evidence.o: main_evidence.cc main.h
\tg++ -c main_evidence.cc
dSph.o: dSph.cc dSph.h
\tg++ -c dSph.cc -o dSph.o
Components.o: Components.cc Components.h
\tg++ -c Components.cc
Tools.o: Tools.cc Tools.h
\tg++ -c Tools.cc
clean:
\trm *.o rundSph

# END OF MAKEFILE"""
    
    else:

        string_makefile = """# MAKEFILE
CFITSIO = /Users/elisa
FITSLIB = -L$(CFITSIO)/lib -lcfitsio -Wl,-rpath,$(CFITSIO)/lib
LIBS = `pkg-config --cflags --libs gsl`
LDFLAGS = $(FITSLIB) $(LIBS)
rundSph: main_bound.o dSph.o Components.o Tools.o
\tg++  main_bound.o dSph.o Components.o Tools.o -o rundSph ${LDFLAGS}
main_bound.o: main_bound.cc main.h
\tg++ -c main_bound.cc
dSph.o: dSph.cc dSph.h
\tg++ -c dSph.cc -o dSph.o
Components.o: Components.cc Components.h
\tg++ -c Components.cc
Tools.o: Tools.cc Tools.h
\tg++ -c Tools.cc
clean:
\trm *.o rundSph

# END OF MAKEFILE"""
    # print(string_makefile)
    with open("Makefile", "w") as fout:
        fout.write(string_makefile)

def write_makefile_spectrum():
    string_makefile = """# MAKEFILE
CFITSIO = /Users/elisa
FITSLIB = -L$(CFITSIO)/lib -lcfitsio -Wl,-rpath,$(CFITSIO)/lib
LIBS = `pkg-config --cflags --libs gsl`
LDFLAGS = $(FITSLIB) $(LIBS)\nrundSph: main_spectrum.o dSph.o Components.o Tools.o
\tg++  main_spectrum.o dSph.o Components.o Tools.o -o rundSph ${LDFLAGS}
main_spectrum.o: main_spectrum.cc main.h
\tg++ -c main_spectrum.cc
dSph.o: dSph.cc dSph.h
\tg++ -c dSph.cc -o dSph.o
Components.o: Components.cc Components.h
\tg++ -c Components.cc
Tools.o: Tools.cc Tools.h
\tg++ -c Tools.cc
clean:
\trm *.o rundSph

# END OF MAKEFILE"""
    with open("Makefile", "w") as fout:
        fout.write(string_makefile)
        
def initialize_fwhm():
    fwhm = {}
    with open(file_fwhm, "r") as fin:
        for line in fin:
            line = strip_line(line)
            fwhm[float(line[0])] = float(line[1])
    return fwhm

def read_chi2_best_permanent(name_dw, i, isCored):
    filename_chi2_null = file_chi2_null_permanent_core if isCored else file_chi2_null_permanent
    
    with open(filename_chi2_null % (name_dw, i), "r") as fin:
        # print('reading file', file_chi2_null_permanent % (name_dw, i))
        
        line = fin.readline()
        # Entries are in file are: logGamma_best, Sflat_best, r0_best, chi2best
        line = strip_line(line)
        logGamma_best = float(line[0])
        Sflat_best = float(line[1])
        r0_best = float(line[2])
        chi2null = float(line[3])
    
    return logGamma_best, Sflat_best, r0_best, chi2null

def read_chi2list_permanent(name_dw, j, isCored):
    # returns 3D array of chi2s. First index is g*sqrt(rho), second index is Sflat, third index is r0
    # j is the channel number, i_r0 is the index in log_r0Dict[name_dw]
    # NOTA BENE: this method won't work if log_r0Dict[name_dw][i_r0]>= 3
    chi2s_tmp = np.empty(N_g_sqrtrho * N_Sflats * N_r0s)
    
    filename_chi2_list = file_chi2_list_permanent_core if isCored else file_chi2_list_permanent
    with open(filename_chi2_list % (name_dw, j), "r") as fin:
        # print('Opening file ',filename_chi2_list % (name_dw, j))
        k = 0
        for line in fin:
            line = strip_line(line)
            # print('line =', k, line)
            chi2s_tmp[k] = float(line[3])
            k += 1
        
        chi2s = chi2s_tmp.reshape((N_g_sqrtrho, N_Sflats, N_r0s))
    
    return chi2s

def linint(fwhm, wl):
    wls = np.fromiter(fwhm.keys(), dtype=float)
    fwhms = np.fromiter(fwhm.values(), dtype=float)
    
    if wl < wls[0] or wl > wls[-1]:
        
        # print('linint called out of range. wavelength = %f, minimum wl = %f,  maximum wl = %f' % (wl, wls[0], wls[-1]))
        
        if wl < wls[0]:
            
            return fwhms[0]
        
        else:
            
            return fwhms[-1]
    
    i = np.where(wls < wl)[0][-1]
    a = (fwhms[i + 1] - fwhms[i]) / (wls[i + 1] - wls[i])
    b = fwhms[i] - a * wls[i]
    
    return a * wl + b

def write_inputfile_cube(ichlow, ichup, flat, nu, delta_wl, beam):
    with open(file_input_cube, "w") as fout:
        l = ['%d' % ichlow, '%d' % ichup, '%f' % flat, '%f' % nu, '%f' % delta_wl, '%f' % beam]
        s = " ".join(l) + "\n"
        fout.write(s)

def my_ceil(a, decimals=0):
    return np.true_divide(np.ceil(a * 10**decimals), 10**decimals)

def my_floor(a, decimals=0):
    return np.true_divide(np.floor(a * 10**decimals), 10**decimals)

def NFW_parameters():
    # Data from 2112.09374
    Msun = 1.9884e+30 * 10 ** 3  # grams
    cm = 3.24078e-22  # kpc
    print(Planck15.critical_density(0) / cm ** 3 / Msun, Planck15.h)
    print('\nInputs are M200 and c200')
    print('\nAnt')
    nfw = NFW(mass=10 ** 9.17, concentration=10 ** 1.26, redshift=0.0, massfactor=('critical', 200), cosmo=Planck15)
    print('rho0 =', nfw.rho_scale.value / 10 ** 8., 'x 10^8 Msun kpc^-3')
    print('rs =', nfw.r_s.value * Msun ** (1 / 3) * cm, 'kpc')
    print('rs = r200/c200 =', nfw.r_virial.value * Msun ** (1 / 3) * cm / 10 ** 1.26)
    print('log10(r200/kpc) =', np.log10(nfw.r_virial.value * Msun ** (1 / 3) * cm))
    print(nfw.mass, 10 ** 9.17)
    # print('\nLeo')
    # nfw = NFW(mass=10 ** 9.23, concentration=10 ** 1.24, redshift=0.0, massfactor=('critical', 200), cosmo=Planck15)
    # print('rho0 =', nfw.rho_scale.value / 10 ** 8., 'x 10^8 Msun kpc^-3')
    # print('rs =', nfw.r_s.value * Msun ** (1 / 3) * cm, 'kpc')
    # print('rs = r200/c200 =', nfw.r_virial.value * Msun ** (1 / 3) * cm / 10 ** 1.24)
    # print('log10(r200/kpc) =', np.log10(nfw.r_virial.value * Msun ** (1 / 3) * cm))
    #
    # print('\nEri')
    # nfw = NFW(mass=10 ** 8.92, concentration=10 ** 1.26, redshift=0.0, massfactor=('critical', 200), cosmo=Planck15)
    # print('rho0 =', nfw.rho_scale.value / 10 ** 8., 'x 10^8 Msun kpc^-3')
    # print('rs =', nfw.r_s.value * Msun ** (1 / 3) * cm, 'kpc')
    # print('rs = r200/c200 =', nfw.r_virial.value * Msun ** (1 / 3) * cm / 10 ** 1.26)
    # print('log10(r200/kpc) =', np.log10(nfw.r_virial.value * Msun ** (1 / 3) * cm))
    #
    # print('\nHya')
    # nfw = NFW(mass=10 ** 10.78, concentration=10 ** 1.15, redshift=0.0, massfactor=('critical', 200), cosmo=Planck15)
    # print('rho0 =', nfw.rho_scale.value / 10 ** 8., 'x 10^8 Msun kpc^-3')
    # print('rs =', nfw.r_s.value * Msun ** (1 / 3) * cm, 'kpc')
    # print('rs = r200/c200 =', nfw.r_virial.value * Msun ** (1 / 3) * cm / 10 ** 1.15)
    # print('log10(r200/kpc) =', np.log10(nfw.r_virial.value * Msun ** (1 / 3) * cm))
    #
    # print('\nGru')
    # nfw = NFW(mass=10 ** 8.51, concentration=10 ** 1.29, redshift=0.0, massfactor=('critical', 200), cosmo=Planck15)
    # print('rho0 =', nfw.rho_scale.value / 10 ** 8., 'x 10^8 Msun kpc^-3')
    # print('rs =', nfw.r_s.value * Msun ** (1 / 3) * cm, 'kpc')
    # print('rs = r200/c200 =', nfw.r_virial.value * Msun ** (1 / 3) * cm / 10 ** 1.29)
    # print('log10(r200/kpc) =', np.log10(nfw.r_virial.value * Msun ** (1 / 3) * cm))

def get_channels_within_resolution(wl, name_dw, nsigma):
    channels_galaxy = channelDict[name_dw]  # channel centers of the galaxy we are considering
    # Risoluzione in wavelenght. Official one from https://www.aanda.org/articles/aa/pdf/2017/12/aa30833-17.pdf
    delta_wl = 5.983 - 9.080e-4 * wl + 5.835e-8 * wl ** 2  # in Angstrom
    
    channels = np.intersect1d(np.where(channels_galaxy <= my_ceil(wl + nsigma * delta_wl / 2, decimals=1)),
                              np.where(channels_galaxy >= my_floor(wl - nsigma * delta_wl / 2, decimals=1)))
    
    ichlow = np.maximum(0, channels[0])
    ichup = np.minimum(n_spectrum, channels[-1])
    
    return ichlow, ichup, delta_wl

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]

############################# CHAIN PREPARATION ###########################################

def format_chain(name_dw):
    Msun = 1.9884e+30 * 10 ** 3  # grams
    cm = 3.24078e-22  # kpc
    
    if name_dw in ['eri1', 'eri2', 'eri3']:
        name_dw = 'eri'
    
    timer.tic()
    chains = np.loadtxt(chains_dir + name_dw + '/gs/' + name_dw + '_Boutput_chain_nfw.txt')
    timer.toc()
    print('initial shape =', chains.shape)
    
    chi2s = chains[:, 0]
    
    index = np.where(chi2s < np.min(chi2s) * 500.0)[0]  # Remove bad fits. Taken from GravSphere
    chains = chains[index]
    chi2s = chi2s[index]
    print('shape after removing big chi2 =', chains.shape, chi2s.shape)

    M200s = 10 ** chains[:, 11]
    c200s = chains[:, 12]
    u, indices = np.unique(M200s, return_index=True)
    M200s = M200s[indices]
    c200s = c200s[indices]
    chi2s = chi2s[indices]
    
    rho0s = np.empty(len(M200s))
    r0s = np.empty(len(M200s))
    with open(file_chains % name_dw, "w") as fout:
        timer.tic()
        for i in np.arange(len(M200s)):
            M200 = M200s[i]
            c200 = c200s[i]
            nfw = NFW(mass=M200, concentration=c200, redshift=0.0, massfactor=('critical', 200), cosmo=Planck15)
            rho0s[i] = nfw.rho_scale.value  # Msun kpc^-3
            r0s[i] = nfw.r_s.value * Msun ** (1 / 3) * cm  # kpc
            l = ['%.10f' % r0s[i], '%.10f' % rho0s[i], '%.10f' % chi2s[i]]
            s = " ".join(l) + "\n"
            fout.write(s)
        timer.toc()
    np.savez(file_chains_npz % name_dw, r0s=r0s, rho0s=rho0s, chi2s=chi2s)

def format_chain_scu():
    name_dw = 'scu'
    Msun = 1.9884e+30 * 10 ** 3  # grams
    cm = 3.24078e-22  # kpc
    
    timer.tic()
    chains = np.loadtxt(chains_dir + 'Sculptor_ChainsDwarfs_NFW')
    timer.toc()
    print('initial shape =', chains.shape)
    
    chi2s = -2 * chains[:, -1]
    # print(chi2s)
    
    index = np.where(chi2s < np.min(chi2s) * 500.0)[0]  # Remove bad fits. Taken from GravSphere
    chains = chains[index]
    chi2s = chi2s[index]
    print('shape after removing big chi2 =', chains.shape, chi2s.shape)
    
    rho0s = 10 ** chains[:, 0]  # Msun kpc^-3
    r0s = 10 ** chains[:, 1]  # kpc
    u, indices = np.unique(rho0s, return_index=True)
    rho0s = rho0s[indices]
    r0s = r0s[indices]
    chi2s = chi2s[indices]
    print('shape after removig duplicates =', chi2s.shape)
    
    with open(file_chains % name_dw, "w") as fout:
        timer.tic()
        for i in np.arange(len(rho0s)):
            l = ['%.10f' % r0s[i], '%.10f' % rho0s[i], '%.10f' % chi2s[i]]
            s = " ".join(l) + "\n"
            fout.write(s)
        timer.toc()
    np.savez(file_chains_npz % name_dw, r0s=r0s, rho0s=rho0s, chi2s=chi2s)

def format_chain_core(name_dw):
    
    Msun = 1.9884e+30 * 10 ** 3  # grams
    cm = 3.24078e-22  # kpc
    
    if name_dw in ['eri1', 'eri2', 'eri3']:
        name_dw = 'eri'
    
    timer.tic()
    chains = np.loadtxt(chains_dir + name_dw + '/gs/' + name_dw + '_Boutput_chain_coreNFW.txt')
    timer.toc()
    print('initial shape =', chains.shape)
    
    chi2s = chains[:, 0]
    
    index = np.where(chi2s < np.min(chi2s) * 500.0)[0]  # Remove bad fits. Taken from GravSphere
    chains = chains[index]
    chi2s = chi2s[index]
    print('shape after removing big chi2 =', chains.shape, chi2s.shape)
    
    u, indices = np.unique(chi2s, return_index=True)
    chains = chains[indices]
    chi2s = chi2s[indices]
    print('shape after unique =', chains.shape, chi2s.shape)

    M200s = 10 ** chains[:, 11] # M_sun
    c200s = chains[:, 12]
    rcs = 10**chains[:, 13] # kpc
    ns = chains[:, 14]
    rts = 10**chains[:, 15] # kpc
    deltas = chains[:, 16]
    
    rho0s = np.empty_like(M200s)
    r0s = np.empty_like(M200s)
    with open(file_chains_core % name_dw, "w") as fout:
        timer.tic()
        for i in np.arange(len(M200s)):
            M200 = M200s[i]
            c200 = c200s[i]
            nfw = NFW(mass=M200, concentration=c200, redshift=0.0, massfactor=('critical', 200), cosmo=Planck15)
            rho0s[i] = nfw.rho_scale.value  # Msun kpc^-3
            r0s[i] = nfw.r_s.value * Msun ** (1 / 3) * cm  # kpc
            l = ['%.10f' % r0s[i], '%.10f' % rho0s[i], '%.10f' % rcs[i], '%.10f' % ns[i], '%.10f' % rts[i], '%.10f' % deltas[i], '%.10f' % chi2s[i]]
            s = " ".join(l) + "\n"
            fout.write(s)
        timer.toc()
    np.savez(file_chains_core_npz % name_dw, r0s=r0s, rho0s=rho0s, chi2s=chi2s, rcs=rcs, ns=ns, rts=rts, deltas=deltas)

def format_chain_core_scu():
    name_dw = 'scu'
    Msun = 1.9884e+30 * 10 ** 3  # grams
    cm = 3.24078e-22  # kpc
    
    timer.tic()
    chains = np.loadtxt(chains_dir + 'Sculptor_ChainsCore_Dwarfs_6param.txt')
    timer.toc()
    print('initial shape =', chains.shape)
    
    chi2s = -2 * chains[:, -1]
    # print(chi2s)
    
    index = np.where(chi2s < np.min(chi2s) * 500.0)[0]  # Remove bad fits. Taken from GravSphere
    chains = chains[index]
    chi2s = chi2s[index]
    print('shape after removing big chi2 =', chains.shape, chi2s.shape)
    
    u, indices = np.unique(chi2s, return_index=True)
    chains = chains[indices]
    chi2s = chi2s[indices]
    print('shape after removig duplicates =', chi2s.shape)
    
    M200s = 10 ** chains[:, 0] # M_sun
    c200s = chains[:, 1]
    rcs = 10**chains[:, 2] # kpc
    ns = chains[:, 3]
    rts = 10**chains[:, 4] # kpc
    deltas = chains[:, 5]

    rho0s = np.empty_like(M200s)
    r0s = np.empty_like(M200s)
    with open(file_chains_core % name_dw, "w") as fout:
        timer.tic()
        for i in np.arange(len(M200s)):
            M200 = M200s[i]
            c200 = c200s[i]
            nfw = NFW(mass=M200, concentration=c200, redshift=0.0, massfactor=('critical', 200), cosmo=Planck15)
            rho0s[i] = nfw.rho_scale.value  # Msun kpc^-3
            r0s[i] = nfw.r_s.value * Msun ** (1 / 3) * cm  # kpc
            l = ['%.10f' % r0s[i], '%.10f' % rho0s[i], '%.10f' % rcs[i], '%.10f' % ns[i], '%.10f' % rts[i], '%.10f' % deltas[i], '%.10f' % chi2s[i]]
            s = " ".join(l) + "\n"
            fout.write(s)
        timer.toc()
    np.savez(file_chains_core_npz % name_dw, r0s=r0s, rho0s=rho0s, chi2s=chi2s, rcs=rcs, ns=ns, rts=rts, deltas=deltas)

def read_chain_npz(name_dw):
    
    data = np.load(file_chains_npz % name_dw)
    r0s_bas = data['r0s']
    rho0s_bas = data['rho0s']
    chi2s_bas = data['chi2s']
    return r0s_bas, rho0s_bas, chi2s_bas

def read_chain_core_npz(name_dw):
    
    data = np.load(file_chains_core_npz % name_dw)
    r0s_bas = data['r0s']
    rho0s_bas = data['rho0s']
    rcs_bas = data['rcs']
    ns_bas = data['ns']
    rts_bas = data['rts']
    deltas_bas = data['deltas']
    chi2s_bas = data['chi2s']
    return r0s_bas, rho0s_bas, rcs_bas, ns_bas, rts_bas, deltas_bas, chi2s_bas

def bin_chain(name_dw):
    
    print('Binning chain for dSph ', name_dw)
    r0s_bas, rho0s_bas, chi2s_bas = read_chain_npz(name_dw)
    print(chi2s_bas)
    log_r0s = log_r0Dict[name_dw]
    d_log_r0 = log_r0s[1] - log_r0s[0]
    bin_edges_log_r0 = np.concatenate((log_r0s - d_log_r0 / 2, [log_r0s[-1] + d_log_r0 / 2]))  # these are the 'internal' bin edges. there's one more bin at the beginning and one at the end. We will discard these data at the very edges in the loop below
    nbins_r0 = len(log_r0s) + 1
    
    nbins_rho0 = 100
    bin_edges_log_rho0 = np.linspace(np.log10(np.amin(rho0s_bas)) - 8, np.log10(np.amax(rho0s_bas)) - 8, num=nbins_rho0+1) # -8 to change units to 10^8 M_sun
    log_rho0s = (bin_edges_log_rho0[:-1] + bin_edges_log_rho0[1:]) / 2
    
    idx_chi2_superbest = np.argmin(chi2s_bas)
    print('\nbest values')
    print('log10(r0) =', np.log10(r0s_bas[idx_chi2_superbest]))
    print('r0 =', r0s_bas[idx_chi2_superbest])
    print('rho0 =', rho0s_bas[idx_chi2_superbest]*1e-8)
    print()
    chi2s_bin = np.zeros([len(log_r0s), len(log_rho0s)])
    chi2s_bin[:] = np.nan
    
    bin_indices_r0 = np.digitize(np.log10(r0s_bas), bin_edges_log_r0, right=True)
    bin_indices_rho0 = np.digitize(np.log10(rho0s_bas) - 8, bin_edges_log_rho0, right=True)
    
    for i in np.arange(0, len(log_r0s)):
        for j in np.arange(0, len(log_rho0s)):
            
            indices = np.where((bin_indices_r0 == i + 1) & (bin_indices_rho0 == j + 1))
            if indices[0].shape[0] != 0:
                chi2s_bin[i, j] = np.amin(chi2s_bas[indices])
    
    print(log_r0s.shape, log_rho0s.shape, chi2s_bin.shape)
    np.savez(file_chains_npz_binned % name_dw, log_r0s=log_r0s, log_rho0s=log_rho0s, chi2s=chi2s_bin)

def bin_chain_core(name_dw):
    np.set_printoptions(suppress=True)
    nbins_rho0 = 30
    Nbins = 5 * np.array([1, 1, 1, 1])
    
    print('\nBinning chain core for dSph ', name_dw)
    r0s_bas, rho0s_bas, rcs_bas, ns_bas, rts_bas, deltas_bas, chi2s_bas = read_chain_core_npz(name_dw)
    
    print('Range of r0s Bas full', np.amin(np.log10(r0s_bas)), np.amax(np.log10(r0s_bas)))
    
    # Prepare bins in r0. These are the ones we scanned
    log_r0s = log_r0Dict_core[name_dw]
    d_log_r0 = log_r0s[1] - log_r0s[0]
    bin_edges_log_r0 = np.concatenate((log_r0s - d_log_r0 / 2, [log_r0s[-1] + d_log_r0 / 2]))
    # print(log_r0s)
    # print(bin_edges_log_r0)
    
    # Prepare bins in rho0
    bin_edges_log_rho0 = np.linspace(np.log10(np.amin(rho0s_bas)) - 8, np.log10(np.amax(rho0s_bas)) - 8,
                                     num=nbins_rho0 + 1)  # -8 to change units to 10^8 M_sun
    log_rho0s = (bin_edges_log_rho0[:-1] + bin_edges_log_rho0[1:]) / 2
    
    
    # Find best fit values for extra parameters
    idx_chi2_superbest = np.argmin(chi2s_bas)
    print('\nbest values')
    # print('log10(rc) =', np.log10(rcs_bas[idx_chi2_superbest]))
    print('rc =', rcs_bas[idx_chi2_superbest])
    print('n = ', ns_bas[idx_chi2_superbest])
    # print('log10(rt) =', np.log10(rts_bas[idx_chi2_superbest]))
    print('rt =', rts_bas[idx_chi2_superbest])
    print('delta = ', deltas_bas[idx_chi2_superbest])
    print()
    # print('log10(r0) =', np.log10(r0s_bas[idx_chi2_superbest]))
    print('r0 =', r0s_bas[idx_chi2_superbest])
    print('rho0 =', rho0s_bas[idx_chi2_superbest]*1e-8)
    print()
    # sys.exit()
    extra_params = [np.log10(rcs_bas), ns_bas, np.log10(rts_bas), deltas_bas]
    ranges = np.empty([len(extra_params), 2])

    i = 0
    for x in extra_params:
        width = (np.amax(x) - np.amin(x))/Nbins[i]
        ranges[i, 0] = x[idx_chi2_superbest] - width / 2
        ranges[i, 1] = x[idx_chi2_superbest] + width / 2
        i += 1
        
    
    # Find the indices of the chi2 values that have the extra parameters in the best bins. These are indices in the arrays chi2s_bas, r0s_bas, rho0s_bas
    indices = np.where((extra_params[0] >= ranges[0, 0]) & (extra_params[0] <= ranges[0, 1])
                       & (extra_params[1] >= ranges[1, 0]) & (extra_params[1] <= ranges[1, 1])
                       & (extra_params[2] >= ranges[2, 0]) & (extra_params[2] <= ranges[2, 1])
                       & (extra_params[3] >= ranges[3, 0]) & (extra_params[3] <= ranges[3, 1])                       )
    # print(chi2s_bas.shape)
    # print(indices[0].shape)
    # print('Range of r0s Bas selected', np.amin(np.log10(r0s_bas[indices])), np.amax(np.log10(r0s_bas[indices])))
    
    # Make a 2D histogram of the surviving chi2s
    Z, xedges, yedges, binnumber = stats.binned_statistic_2d(np.log10(r0s_bas[indices]),
                                                             np.log10(rho0s_bas[indices]) - 8,
                                                             values=chi2s_bas[indices], statistic='min',
                                                             bins=[bin_edges_log_r0, bin_edges_log_rho0])
    # print('Z.shape =', Z.shape)
    CL90 = np.array([np.amin(chi2s_bas) + 4.61])  # chi2_best + 1.64237,
    # Z[Z>CL90] = np.nan
    CL68 = np.array([np.amin(chi2s_bas) + 2.27887])

    # fig, ax1 = plt.subplots(1, 1, figsize=(5, 5))
    # x, y = np.meshgrid(xedges, yedges)
    # im1 = ax1.pcolormesh(x, y, Z.transpose())  # , norm=matplotlib.colors.LogNorm())
    # cbar = fig.colorbar(im1, orientation='vertical', ax=ax1, label=r'$q$')
    #
    # # cbar = plt.colorbar()
    # cbar.ax1.set_ylabel('y2')
    # print('chi2_best = %f CL90 = %f' % (np.amin(chi2s_bas), np.amin(chi2s_bas) + 4.61))
    # CP = ax1.contour(log_r0s, log_rho0s, Z.transpose(), levels=CL68, linewidths=1.0, colors='r')
    # ax1.plot(np.log10(r0s_bas[idx_chi2_superbest]), np.log10(rho0s_bas[idx_chi2_superbest])-8, marker='o')
    #
    # # CP = ax1.contour(log_r0s, log_rho0s, Z.transpose(), levels=[CL90, CL68], linewidths=1.0, colors=['r', 'orange'])
    # p = CP.collections[0].get_paths()[0].vertices
    # for path in CP.collections[0].get_paths()[1:]:
    #     p = np.vstack((p, path.vertices))
    #
    # # print("\n\ndataCP\n", p)
    # # print("dataCP shape ", p.shape)
    # # print('Range log10(r0):', np.amin(p[:, 0]), np.amax(p[:, 0]))
    # # sys.exit()
    #
    #
    # ax1.set_xlabel(r'$\log_{10}(r_0 / \mathrm{kpc})$')
    # ax1.set_ylabel(r'$\log_{10}(\rho_0 / 10^8M_\odot)$')
    # # ax1.set_xlim([-1, 1.3])
    # plt.title(name_dw)
    #
    # # plt.show()
    # plt.savefig('chi2_bas_'+name_dw+'_cored_restricted.png')
    # plt.close()
    # sys.exit()
    
    # np.savez(file_chains_core_npz_binned % name_dw, log_r0s=log_r0s, log_rho0s=log_rho0s, chi2s=Z, log_rcbest=np.log10(rcs_bas[idx_chi2_superbest]), nbest=ns_bas[idx_chi2_superbest], log_rtbest=np.log10(rts_bas[idx_chi2_superbest]), deltabest=deltas_bas[idx_chi2_superbest])

def read_chain_binned(name_dw, isCored):

    if isCored:
        data = np.load(file_chains_core_npz_binned % name_dw)
    else:
        data = np.load(file_chains_npz_binned % name_dw)
        
    log_r0s = data['log_r0s']
    log_rho0s = data['log_rho0s']
    chi2s_bas = data['chi2s']
    return log_r0s, log_rho0s, chi2s_bas


############################# CALCULATE CHI2 FUNCTIONS ###########################################
def calculate_chi2list(fwhms, i, name_dw, isMakeEvidence, isCored):
    
    wl = channel_centers[i] # wavelength at which we want to perform the analysis. Needs to be always the same, even for leo
    
    ichlow, ichup, delta_wl = get_channels_within_resolution(wl, name_dw, 1)
    
    # Risoluzione angolare at the wavelength under consideration. Calculates standard deviation from FWHM.
    factor = fwhmDict[name_dw] / fwhmDict['leo'] # different angular resolution per galaxy
    beam = linint(fwhms, wl) / (2 * np.sqrt(np.log(4))) / 60  * factor # In arcmin

    # Frequency
    nu = 2.99792458e8 / (wl * 1e-10) / 1e9  # in GHz

    # prepare input file for underlying C program
    write_inputfile_cube(ichlow, ichup, 0, nu, delta_wl, beam)

    # call underlying C program
    print('Running')
    timer.tic()
    subprocess.run(['./rundSph'])
    timer.toc()

    if isCored:
        if isMakeEvidence:
            os.system('cp -f ' + file_chi2_null + ' ' + file_chi2_null_permanent_core % (name_dw, i))
        else:
            os.system('cp -f ' + file_chi2_list % name_dw + ' ' + file_chi2_list_permanent_core % (name_dw, i))
    else:
        if isMakeEvidence:
            os.system('cp -f ' + file_chi2_null + ' ' + file_chi2_null_permanent % (name_dw, i))
        else:
            os.system('cp -f ' + file_chi2_list % name_dw + ' ' + file_chi2_list_permanent % (name_dw, i))

def make_chi2list(names_dSph, i_channels, isMakeEvidence, isCored):
    
    if len(i_channels) == 0:
        i_channels = range(1, n_spectrum - 1) # exclude first and last channel
    
    # Read FWHM as a function of wavelength. FWHM is given in arcseconds (I guess)
    fwhms = initialize_fwhm()

    i = 0
    for name_dw in names_dSph:
        print("Script for bounds on dSph ", name_dw)

        # Create main.h header file including correct parameter file for dSph
        write_mainh(name_dw, isCored)
        
        # Create makefile with appropriate version of dSph.cc
        write_makefile(isMakeEvidence)
        
        # Compile C files
        os.system('make clean')
        subprocess.run(['make'])
        
        for i in i_channels:
            
            if i % skip_wl_every != 0:
                continue
            
            print('\ni =', i, end='\t')
            
            # Calls underlying C program that calculates the chi2s
            calculate_chi2list(fwhms, i, name_dw, isMakeEvidence, isCored)

    
    print('Done.')

def get_spectrum(names_dSph, i_channels):
    if len(i_channels) == 0:
        i_channels = range(1, n_spectrum - 1)  # exclude first and last channel
    
    # Read FWHM as a function of wavelength. FWHM is given in arcseconds (I guess)
    fwhms = initialize_fwhm()
    
    for name_dw in names_dSph:
        print("Script for spectrum of dSph ", name_dw)
        
        # Create main.h header file including correct parameter file for dSph
        write_mainh(name_dw, False)
        
        # Create makefile with appropriate version of dSph.cc
        write_makefile_spectrum()
        
        # Compile C files
        os.system('make clean')
        subprocess.run(['make'])
        
        spectra = np.empty(len(i_channels))
        
        idx = 0
        for i in i_channels:
            
            if i % skip_wl_every != 0:
                continue
            
            print('\ni =', i, end='\t')
            
            # Calls underlying C program that calculates the chi2s
            wl = channel_centers[i]  # wavelength at which we want to perform the analysis. Needs to be always the same, even for leo
            
            ichlow, ichup, delta_wl = get_channels_within_resolution(wl, name_dw, 1)
            
            # Risoluzione angolare at the wavelength under consideration. Calculates standard deviation from FWHM.
            beam = linint(fwhms, wl) / (2 * np.sqrt(np.log(4))) / 60  # In arcmin
            
            # Frequency
            nu = 2.99792458e8 / (wl * 1e-10) / 1e9  # in GHz
            
            # prepare input file for underlying C program
            write_inputfile_cube(ichlow, ichup, 0, nu, delta_wl, beam)
            
            # call underlying C program
            subprocess.run(['./rundSph'])
            # res = subprocess.check_output(["./rundSph"], universal_newlines=True)
            # spectra[idx] = float(res)
            idx += 1
        # np.savez('spectra_'+name_dw+'.npz', spectra=spectra, channels=i_channels)
    
    print('Done.')


############################# BOUND FUNCTIONS ###########################################
def merge_eridanus(names_dSph, i_channels, isCored):
    if len(i_channels) == 0:
        i_channels = range(1, n_spectrum - 1)
    
    print('Script for merging chi2s of Field 1-2 of Eridanus')
    for i in i_channels:
        if i % skip_wl_every != 0:
            continue
        
        # print('i =', i)
        
        chi2s = np.zeros([N_g_sqrtrho, N_Sflats, N_r0s])
        
        filename_chi2_list_permanent = file_chi2_list_permanent_core if isCored else file_chi2_list_permanent
        with open(filename_chi2_list_permanent % ('eri', i), "w") as fout:
            # print('writing to file', file_chi2_list_permanent % ('eri', i))
            for name_dw in names_dSph:
                chi2s += read_chi2list_permanent(name_dw, i, isCored)  # chi2(log_g_sqrtrho, log_Sflat, log_r0)


            r0s = 10 ** log_r0Dict_core[name_dw] if isCored else 10 ** log_r0Dict[name_dw]
            ig = 0
            for g in 10 ** log_g_sqrtrho:
                iS = 0
                for S in 10 ** log_Sflats:
                    ir = 0
                    for r in r0s:
                        l = ['%e' % g, '%f' % S, '%f' % r, '%f' % chi2s[ig, iS, ir]]
                        s = " ".join(l) + "\n"
                        fout.write(s)
                        ir += 1
                    iS += 1
                ig += 1

        chi2_null_tot = np.zeros(N_Sflats * N_r0s)
        
        filename_chi2_null_permanent = file_chi2_null_permanent_core if isCored else file_chi2_null_permanent
        with open(filename_chi2_null_permanent % ('eri', i), "w") as fout2:
            
            for name_dw in names_dSph:
                filename_chi2_null = file_chi2_null_permanent_core if isCored else file_chi2_null_permanent
                data = np.loadtxt(filename_chi2_null % (name_dw, i))
    
                logGamma_best, Sflat_best, r0_best, chi2null = data.T
                # print(logGamma_best, Sflat_best, r0_best, chi2null)
                chi2_null_tot += chi2null
            
            idx = np.argmin(chi2_null_tot)
            l = ['%e' % logGamma_best[idx], '%f' % Sflat_best[idx], '%f' % r0_best[idx], '%f' % chi2_null_tot[idx]]
            s = " ".join(l) + "\n"
            fout2.write(s)

def P_gauss(chi2, chi2_ref):
    return 1 - 0.5 * erfc(np.sqrt((chi2 - chi2_ref) / 2))

def P_chi(chi2, chi2_ref):
    return 1 - erfc(np.sqrt((chi2 - chi2_ref) / 2))

def calculate_upper_bound(chi2_mins, chi2_best, idx_chi2_best):
    i = idx_chi2_best + 1
    log_g_bound = 100
    for chi2 in chi2_mins[idx_chi2_best + 1:]:
        P0 = P_gauss(chi2_mins[i - 1], chi2_best)
        P1 = P_gauss(chi2_mins[i], chi2_best)

        if P1 >= Pcl:
            log_g_bound = log_g[i] + (log_g[i] - log_g[i - 1]) / (P1 - P0) * (Pcl - P1)
            break
        i += 1
    
    return log_g_bound

def prepare_upper_bound(chi2, chi2null, i):
    chi2_best = np.nanmin(chi2)
    idx_chi2_best = np.where(chi2 == chi2_best)[0][-1]
    chi2_worst = np.nanmax(chi2[idx_chi2_best:])
    
    # Make bound
    if chi2null < chi2_best and P_gauss(chi2_best, chi2null) > Pcl:
        # the null hypothesis (no DM) is strongly preferred
        print(
            'chi2null < chi2_best and P_gauss(chi2_best, chi2null) > Pcl. i = %d, chi2_best = %.2f, chi2null = %.2f' % (
                i, chi2_best, chi2null))
        gbound = -100

    else:
        if chi2null < chi2_best:
            print('Changing chi2_best. i = %d, chi2null = %.2f, chi2_best = %.2f, chi2null - chi2_best = %.2e' % (
            i, chi2null, chi2_best, (chi2null - chi2_best)))
            chi2_best = chi2null
            # idx_chi2_best = 0
        
        if P_gauss(chi2_worst, chi2_best) < Pcl:
            # we can't find a bound within the range of parameters we scanned
            print(
                'BOUND NOT FOUND: P_gauss(chi2_worst, chi2_best) < Pc. i = %d, chi2_best = %.2f, chi2worst = %.2f' % (
                    i, chi2_best, chi2_worst))
            gbound = 100
        else:
            gbound = calculate_upper_bound(chi2, chi2_best, idx_chi2_best)
            
    return gbound
        
def make_bounds_single(names_dSph, i_channels, isCored):
    np.set_printoptions(25)
    
    if len(i_channels) == 0:
        i_channels = range(1, n_spectrum - 1)
    
    for name_dw in names_dSph:
        print("\n\nScript for bounds on dSph %s, isCored = %s" % (name_dw, str(isCored)))
        
        log_r0s_bas, log_rho0s_bas, chi2s_bas = read_chain_binned(name_dw, isCored)

        i_bas= np.unravel_index(np.nanargmin(chi2s_bas), chi2s_bas.shape)
        # print('chi2s_bas.shape, i_bas =', chi2s_bas.shape, i_bas)
        # print('Bas best values. r0s_bas, rho0s_bas, chi2s_bas =', log_r0s_bas[i_bas[0]], log_rho0s_bas[i_bas[1]], chi2s_bas[i_bas])

        r0_return = np.empty_like(i_channels, dtype=float)
        rho0_return = np.empty_like(i_channels, dtype=float)
        filename_bounds = file_bounds_cored if isCored else file_bounds
        with open(filename_bounds % name_dw, "w") as fout:
            j = 0
            for i in i_channels:
                
                if i % skip_wl_every != 0:
                    continue
                if name_dw == 'hyd' and i == 884:
                    continue
                    
                print('\ni =', i)
                # Axion mass
                wl = channel_centers[i]
                ma = np.round(7893.1756681904 *  pi / wl, decimals=6) # in eV
                
                # Read all our chi2s
                chi2s_noi_big = read_chi2list_permanent(name_dw, i, isCored)
                print('np.nanmin(chi2s_noi_big) =', np.nanmin(chi2s_noi_big))
                chi2_tot_g_r0 = np.empty([N_gs, N_r0s])
                N_rhos = len(log_rho0s_bas)
                
                r0s = log_r0Dict_core[name_dw] if isCored else log_r0Dict[name_dw]
                rho0s_best = np.ones_like(r0s) * 1e3
                for i_r0 in np.arange(len(r0s)):
                    
                    # chi2s_noi is 2D array for fixed r0
                    chi2s_noi = chi2s_noi_big[:, :, i_r0]

                    # Interpolate chi2s_noi along the log_g_sqrtrho. Calling chi2s_noi_int returns an array of chi2s
                    # for different values of Sflat
                    chi2s_noi_int = interpolate.interp1d(log_g_sqrtrho, chi2s_noi, axis=0)
                    
                    # Initialize array to contain chi2_tot (noi+bas)
                    chi2_tot = np.empty([N_gs, N_rhos, N_Sflats])
                    
                    # For each g in our pre-set list, and for each rho0 in log_rho0s_bas,
                    # calculate g sqrt(rho0) and ask chi2s_noi_int to give the corresponding chi2_noi.
                    # Returns an array of dimension N_Sflat.
                    # Then sum this chi2_noi to the chi2_bas for the same rho0
                    for irho in np.arange(N_rhos):
                        # print('i_r0, irho =', i_r0, irho)
                        chi2_bass = chi2s_bas[i_r0, irho]
                        if np.isnan(chi2_bass):
                            chi2_tot[:, irho] = np.nan
                            continue
                        try:
                            chi2_tot[:, irho] = chi2_bass + chi2s_noi_int(log_g + 0.5 * log_rho0s_bas[irho]) # this is a 2D array. First dimension is g, second is Sflat
                        except ValueError as e:
                            print('Value Error.', irho, log_g, 0.5 * log_rho0s_bas[irho],
                                  log_g + 0.5 * log_rho0s_bas[irho])
                            sys.exit()

                    # Profile out Sflat and rho0 from chi2_tot. Now we have chi2_tot as a function of g, for our r0 under consideration
                    chi2_tot_g_r0[:, i_r0] = np.nanmin(chi2_tot, axis=(1, 2))  # [N_gs, N_rhos, N_Sflats]
                    # print(log_r0, chi2_tot_g_r0.shape, chi2_tot_g_r0[:, i_r0])
                    try:
                        idxs = np.unravel_index(np.nanargmin(chi2_tot), chi2_tot.shape)
                        rho0s_best[i_r0] = 10**log_rho0s_bas[idxs[1]]
                        # print(i_r0, 10**r0s[i_r0] ,idxs, 10**log_rho0s_bas[idxs[1]], np.nanmin(chi2_tot))
                    except ValueError as e:
                        # print('Value Error.')
                        continue
                    
                # Get chi2_null
                chi2null = read_chi2_best_permanent(name_dw, i, isCored)[3] + np.nanmin(chi2s_bas)
                
                # Profile out r0. Get chi2_tot as a function of log_g only
                chi2_mins = np.nanmin(chi2_tot_g_r0, axis=1)
                idx = np.unravel_index(np.nanargmin(chi2_tot_g_r0), chi2_tot_g_r0.shape)
                # print(chi2_mins)
                # print(idx)
                # print('r0best = %f, rho0best = %f chi2_best = %f' % (10**r0s[idx[1]], rho0s_best[idx[1]],np.nanmin(chi2_tot_g_r0)) )
                rho0_return[j] = rho0s_best[idx[1]]
                r0_return[j] = 10**r0s[idx[1]]
                # print(j, rho0_return[j], r0_return[j])
                j += 1
                gbound = prepare_upper_bound(chi2_mins, chi2null, i)
                
                print(i, gbound, 10 ** gbound)
                l = ['%e' % ma, '%e' % 10 ** gbound]
                s = " ".join(l) + "\n"
                
                fout.write(s)
                filename_chi2_of_g = file_chi2_of_g_core if isCored else file_chi2_of_g
                np.savez(filename_chi2_of_g % (name_dw, i), log_gs=log_g, chi2s=chi2_mins)
    
    print('Done.')
    return gbound, ma, r0_return, rho0_return

def make_bounds_combined(names_dSph, i_channels, isCored, isShifted):
    np.set_printoptions(25)
    if len(i_channels) == 0:
        i_channels = range(1, n_spectrum - 1)
    
    print('\n\nMaking combined bound')

    if isShifted:
        filename_bounds = file_bounds_doppler_cored if isCored else file_bounds_doppler
    else:
        filename_bounds = file_bounds_cored if isCored else file_bounds

    with open(filename_bounds % "_".join(names_dSph), "w") as fout:
        print('Writing file ', file_bounds % "_".join(names_dSph))
        
        for i in i_channels:
            chi2s = np.zeros(N_gs)
            chi2_null = 0
            
            if i % skip_wl_every != 0:
                continue

            # Frequency
            wl = channel_centers[i]
            ma = np.round(7893.1756681904 *  pi / wl, decimals=6)  # in eV
            
            i_dw = 0
            for name_dw in names_dSph:
                chi2s_bas = read_chain_binned(name_dw, isCored)[2]

                filename_chi2_of_g = file_chi2_of_g_core if isCored else file_chi2_of_g
                if isShifted:
                    # ma is the mass emitted, the actual axion mass at which I want to put the bound
                    # mass_obs is what MUSE observes
                    # For each ma, I want to find which channel I have to consider
                    masses_extended = np.round(7893.1756681904 *  pi / channel_centers_extended, decimals=6)  # in eV

                    mass_obs = ma * (1 - v_radial_dict[name_dw])
                    
                    # Observational channel corresponding to given axion mass
                    i_observed = find_nearest(masses_extended, mass_obs)[0] - 10
                    if i_observed < 1 or i_observed > n_spectrum - 2 or (883 < i_observed < 1013):
                        print('** Skipping channel %d for %s' % (i, name_dw))
                        continue
                    if name_dw == 'hyd' and i_observed == 884:
                        continue
                    else:
                        # print('Loading data for channel %d for %s' % (i, name_dw))
                        data = np.load(filename_chi2_of_g % (name_dw, i_observed))
                        chi2_null += read_chi2_best_permanent(name_dw, i_observed, isCored)[3] + np.nanmin(chi2s_bas)

                else:
                    if name_dw == 'hyd' and i == 884:
                        continue
                        
                    data = np.load(filename_chi2_of_g % (name_dw, i))
                    chi2_null += read_chi2_best_permanent(name_dw, i, isCored)[3] + np.nanmin(chi2s_bas)
                    
                chi2s += data['chi2s']

            gbound = prepare_upper_bound(chi2s, chi2_null, i)
            print(i, gbound, 10**gbound)

            l = ['%e' % ma, '%e' % 10 ** gbound]
            s = " ".join(l) + "\n"
            
            fout.write(s)

            if isShifted:
                filename_chi2_of_g = file_chi2_of_g_core_doppler if isCored else file_chi2_of_g_doppler
            else:
                filename_chi2_of_g = file_chi2_of_g_core if isCored else file_chi2_of_g
            np.savez(filename_chi2_of_g % ("_".join(names_dSph), i), log_gs=log_g, chi2s=chi2s)
    
    print('Done.')

def doppler_shift_single(names_dSph, isCored, isEvidence):
    # Works only for bounds/evidence on ALL channels
    
    for name_dw in names_dSph:
        print("\n\nScript for Doppler shifting dSph %s, isCored = %s" % (name_dw, str(isCored)))
        
        if isEvidence:
            filename = file_evidence_no_peaks_LEE_cored if isCored else file_evidence_no_peaks_LEE
        else:
            filename = file_bounds_cored if isCored else file_bounds
            
        data = np.loadtxt(filename % name_dw)
        masses_obs = data[:, 0]
        bounds = data[:, 1]
        # print(bounds)

        masses_extended =   np.round(7893.1756681904 *  pi / channel_centers_extended, decimals=6)  # in eV
        
        masses_emitted = masses_obs * (1 + v_radial_dict[name_dw]) # in eV
        masses_emitted_discrete = np.empty_like(masses_emitted)
        
        i = 0
        for mass_emitted in masses_emitted:
            
            i_emitted, masses_emitted_discrete[i] = find_nearest(masses_extended, mass_emitted)
            print('\nmass emitted, masses obs', masses_emitted_discrete[i], masses_obs[i], masses_emitted_discrete[i] - masses_obs[i])
            print('i_observed, i_emitted =', i + 1, i_emitted-10)
            i += 1
            
        if isEvidence:
            filename_shifted = file_evidence_no_peaks_LEE_doppler_cored if isCored else file_evidence_no_peaks_LEE_doppler
        else:
            filename_shifted =  file_bounds_doppler_cored if isCored else file_bounds_doppler
            
        np.savetxt(filename_shifted % name_dw, np.c_[masses_emitted_discrete, bounds], fmt=['%.6e','%.6e'], delimiter=' ', newline=os.linesep)


############################# EVIDENCE FUNCTIONS ###########################################
def make_evindence_single(names_dSph, i_channels, isCored):
    
    if len(i_channels) == 0:
        i_channels = range(1, n_spectrum - 1)
    
    for name_dw in names_dSph:  # names_dSph:#['scu']: #['leo']:
        print("\nScript for evidence on dSph ", name_dw)

        chi2s_bas = read_chain_binned(name_dw, isCored)[2]
        chi2null_bas = np.nanmin(chi2s_bas)
        
        filename_evidence = file_evidence_cored if isCored else file_evidence
        with open(filename_evidence % name_dw, "w") as fout:
            for i in i_channels:
                
                if i % skip_wl_every != 0:
                    continue
                # if name_dw == 'hyd' and i == 884:
                #     continue
                # print('i =',i)
                if 884 <= i <= 1012:
                    l = ['%e' % ma, '%f' % np.nan]
                    s = " ".join(l) + "\n"
                    fout.write(s)
                    print(i, np.nan)
                
                else:
                    # Axion mass
                    wl = channel_centers[i]
                    ma = np.round(7893.1756681904 * pi / wl, decimals=6)  # in eV
                    
                    # chi2s best file
                    chi2null_noi = read_chi2_best_permanent(name_dw, i, isCored)[3]
                    chi2null = chi2null_noi + chi2null_bas
                    # print('\nchi2null_noi, np.nanmin(chi2s_bas) =', chi2null_noi, np.nanmin(chi2s_bas))
                    
                    filename_chi2_of_g = file_chi2_of_g_core if isCored else file_chi2_of_g
                    chi2s = np.load(filename_chi2_of_g % (name_dw, i))['chi2s']
                    chi2_best = np.amin(chi2s)
                    # print('chi2_best =', chi2_best)
                    
                    if chi2null - chi2_best < 0:
                        chi2_best = chi2null
                    
                    l = ['%e' % ma, '%f' % np.sqrt(chi2null - chi2_best)]
                    s = " ".join(l) + "\n"
                
                    fout.write(s)
                    print(i, np.sqrt(chi2null - chi2_best))

    
    print('Done.')

def make_evindence_combined(names_dSph, i_channels, isCored):
    print('\n\nMaking combined evidence')
    
    if len(i_channels)  == 0 :
        i_channels = range(1, n_spectrum - 1)
    
    
    filename_evidence = file_evidence_cored if isCored else file_evidence
    with open(filename_evidence % "_".join(names_dSph), "w") as fout:
        for i in i_channels:
            
            if i % skip_wl_every != 0:
                continue
            
            # print('\ni =',i)
            
            # Axion mass
            wl = channel_centers[i]
            ma = np.round(7893.1756681904 *  pi / wl, decimals=6)  # in eV
            
            chi2null = 0
            chi2s_comb = np.zeros(N_gs)
            chi2nulls = np.empty(len(names_dSph))
            chi2_bests = np.empty([len(names_dSph), N_gs])
            
            idw = 0
            for name_dw in names_dSph:
                # print('Evidence for dSph ', name_dw)
                chi2null_noi = read_chi2_best_permanent(name_dw, i, isCored)[3]
                chi2s_bas = read_chain_binned(name_dw, isCored)[2]
                chi2null += chi2null_noi + np.nanmin(chi2s_bas)
                chi2nulls[idw] = chi2null_noi + np.nanmin(chi2s_bas)
                
                filename_chi2_of_g = file_chi2_of_g_core if isCored else file_chi2_of_g
                chi2s = np.load(filename_chi2_of_g % (name_dw, i))['chi2s']
                chi2s_comb += chi2s
                chi2_bests[idw] = chi2s
                
                idw += 1
            
            chi2superbest = np.amin(chi2s_comb)
            
            if chi2null < chi2superbest:  # no DM is favoured. No evidence.
                chi2superbest = chi2null
            
            l = ['%e' % ma, '%e' % np.sqrt(chi2null - chi2superbest)]
            s = " ".join(l) + "\n"
            fout.write(s)
    
    print('Done.')

def search_peak_boundaries(peak, array):
    # print('peak =',peak)
    for j in np.arange(peak, len(array) - 1):
        # print(j, array[j], array[j+1] )
        if array[j] > array[j + 1]:
            # print('continue')
            if j+1 == len(array)-1:
                right_ips = j+1
            continue
        if array[j] == 0:
            right_ips = j-1
        else:
            right_ips = j-1

        break
    
    for j in np.flip(np.arange(0, peak + 1)):
        # print(j, evidence[j], evidence[j-1] )
        if array[j] > array[j - 1]:
            # print('continue')
            continue
        if array[j] == 0:
            left_ips = j+1
        else:
            left_ips = j+1
        # print(left_ips)
        break

    # print('peak =',peak, left_ips, right_ips)
    return left_ips, right_ips

def fit_peaks(data, idxs, offset, degree_polyn):
    peaks, _ = find_peaks(-data)
    data = data[peaks]
    z = np.polyfit(peaks + offset, data, degree_polyn)
    p = np.poly1d(z)
    
    return p

def find_threshold(data, idxs, nloops, degree_polyn):
    
    x1 = np.arange(idxs[0])
    x2 = np.arange(idxs[-1] + 1, len(data))
    data1 = data[x1]
    data2 = data[x2]
    
    data1[data1 > 10* np.mean(data1)] = np.mean(data1)
    data2[data2 > 10* np.mean(data2)] = np.mean(data2)

    for i in np.arange(nloops[0]):
        p1 = fit_peaks(data1, idxs, 0, degree_polyn[0])
        data1 = np.minimum(data[x1], p1(x1))
    
    for i in np.arange(nloops[1]):
        p2 = fit_peaks(data2, idxs, idxs[-1] + 1, degree_polyn[1])
        data2 = np.minimum(data[x2], p2(x2))
        
    return p1, p2

def subtract_threshold_from_skyspectrum():
    
    nloops = np.array([1, 3])
    degree_polyn = np.array([5, 5])
    
    data = genfromtxt('skyspec_sum.dat', delimiter=' ')
    wls_sky_spectrum = data[:, 0]
    fluxes_sky_spectrum = np.array(data[:, 1])
    fluxes_sky_spectrum[np.where(fluxes_sky_spectrum == 0)] = np.nan
    print('fluxes_sky_spectrum.shape =', fluxes_sky_spectrum.shape)
    
    idxs = np.where(np.isnan(fluxes_sky_spectrum))[0]
    x1 = np.arange(idxs[0])
    x2 = np.arange(idxs[-1] + 1, len(fluxes_sky_spectrum))
    
    p1, p2 = find_threshold(fluxes_sky_spectrum, idxs, nloops, degree_polyn)
    threshold1 = p1(x1)
    threshold2 = p2(x2)
    threshold = np.concatenate((np.concatenate((threshold1, 1e1 * np.ones_like(idxs))), threshold2))
    print('threshold.shape =', threshold.shape)
    
    plt.plot(fluxes_sky_spectrum)
    plt.plot(threshold)
    plt.show()
    
    fluxes_sky_spectrum = fluxes_sky_spectrum - threshold
    fluxes_sky_spectrum[np.where(fluxes_sky_spectrum < 0)] = 0
    
    with open('skyspec_sum_threshold.dat', "w") as fout:
        for i in np.arange(len(wls_sky_spectrum)):
            l = ['%.10f' % wls_sky_spectrum[i], '%.10f' % fluxes_sky_spectrum[i],  '%.10f' % threshold[i]]
            s = " ".join(l) + "\n"
            # print(s)
            fout.write(s)

def find_bad_channels_using_skyspec(name_dw, isCored):
    np.set_printoptions(precision=3)
    
    print('\n*** find_bad_channels. Galaxy:', name_dw)
    
    # Load leo sky spectrum
    data = genfromtxt('skyspec_sum_threshold.dat', delimiter=' ')
    wls_sky_spectrum = data[:, 0]
    fluxes_sky_spectrum = np.array(data[:, 1])  # these are the fluxes with the baseline already subtracted
    baseline = np.array(data[:, 2])
    sigma_flux = np.std(np.concatenate((fluxes_sky_spectrum[:685], fluxes_sky_spectrum[720:880])))
    # print('Wavelength skyspectrum = ', wls_sky_spectrum[:5], wls_sky_spectrum[:-5])
    # print('Channel centers Leo T = ', channel_centers_leo[:5], channel_centers_leo[:-5])

    # print(baseline)
    # print(np.nanmean(fluxes_sky_spectrum))
    # print(sigma_flux)
    
    nsigma = 5.7  # nearby: distance <= nsigma * sigma(wl)
    peak_threshold = 5 * sigma_flux
    print('peak_threshold =', peak_threshold)
    
    peaks_fluxes, _ = find_peaks(fluxes_sky_spectrum, height=peak_threshold)
    
    # Remove peaks that don't correspond to emission lines
    peaks_fluxes = peaks_fluxes[peaks_fluxes != 2492]
    peaks_fluxes = peaks_fluxes[peaks_fluxes != 3495]
    
    bad_channels = []
    for i in range(len(peaks_fluxes)):
        peak = peaks_fluxes[i]
        wl = channel_centers_leo[peak]
        
        # Range of bad channels
        ichlow, ichup, delta_wl = get_channels_within_resolution(wl, name_dw, nsigma)
        # print(peak, np.arange(ichlow, ichup + 1))
        bad_channels.append(np.arange(ichlow, ichup + 1))
    
    bad_channels = [item for sublist in bad_channels for item in sublist]
    bad_channels = np.unique(np.array(bad_channels))
    print('Number of bad channels =', bad_channels.shape)
    # for ch in bad_channels:
    #     print(ch, end=" ")
    # print('\n')
    #
    # sys.exit()
    bad_ranges = []
    low = bad_channels[0]
    for i in np.arange(1, len(bad_channels)):
        ch = bad_channels[i]
        if ch != bad_channels[i - 1] + 1:
            bad_ranges.append([low, bad_channels[i - 1]])
            low = ch
        elif i == len(bad_channels) - 1:
            bad_ranges.append([low, ch])
            
    # print(bad_ranges)
    bad_ranges = np.array(bad_ranges)
    bad_wavelength_ranges = wls_sky_spectrum[bad_ranges]
    bad_wavelength_ranges[:, 0] = bad_wavelength_ranges[:, 0] - dchan/2
    bad_wavelength_ranges[:, 1] = bad_wavelength_ranges[:, 1] + dchan/2
    # print(bad_wavelength_ranges)
    # sys.exit()
    
    # Load evindence
    filename = file_evidence_cored if isCored else file_evidence
    data = genfromtxt(filename % name_dw, delimiter=' ')
    masses = data[:, 0]
    evidence = data[:, 1]
    print(evidence.shape)
    evidence_new = np.copy(evidence)
    print('evidence.shape =', evidence.shape)
    nanini = len(np.where(np.isnan(evidence_new))[0])
    print('initial number of nans =', nanini)
    
    evidence = np.concatenate(([evidence[1]], evidence, [evidence[-2]]))  # trick to detect peaks at boundaries of array
    peaks_evi, _ = find_peaks(evidence, height=5)
    peaks_evi = peaks_evi - 1
    evidence = evidence[1:-1]
    print('peaks_evi =', peaks_evi)
    print(evidence.shape)
    
    channels_wavelenghts_used = channel_centers[1:n_spectrum - 1]
    # print('channels_wavelenghts_used = ', channels_wavelenghts_used[:5], channels_wavelenghts_used[:-5])
    # sys.exit()

    counter = 0  # how many evidence peaks we remove
    print('\nSearching for nearby peaks')
    for i in range(len(peaks_evi)):
        peak = peaks_evi[i]
        wl = channels_wavelenghts_used[peak]

        for low, high in zip(bad_wavelength_ranges[:, 0],  bad_wavelength_ranges[:, 1]):
            if low <= wl <= high:
                # print('Throwing peak away ', low, wl, high)
                left_ips, right_ips = search_peak_boundaries(peak, evidence)
                # print('%d \t %d \t %d' % (peak, left_ips, right_ips + 1))
                evidence_new[left_ips:right_ips + 1] = np.nan
                counter += 1
                break

        

    counter_channels = len(np.where(np.isnan(evidence_new))[0])
    print('removed %d evidence peaks out of %d' % (counter, len(peaks_evi)))
    print('number of nans %d channels out of %d' % (counter_channels, n_spectrum - 2))
    print('bad_channels.shape =', bad_channels.shape)
    
    evidence_new[bad_channels - 1] = np.nan
    print('setting all bad channels to nan. Total number of removed channels =', len(np.where(np.isnan(evidence_new))[0]), nanini + len(bad_channels))
    
    
    idx = np.argpartition(-evidence_new, np.arange((~np.isnan(evidence_new)).sum()))
    evidence_tmp = evidence_new[idx]
    print('maximum evidence =', evidence_tmp[:10], idx[:10])
    
    print('\nWriting to file')
    filename = file_evidence_no_peaks_cored if isCored else file_evidence_no_peaks
    with open(filename % name_dw, "w") as fout:
        
        for i in np.arange(len(evidence_new)):
            wl = channels_wavelenghts_used[i]
            ma = np.round(7893.1756681904 * pi / wl, decimals=6)  # in eV
            l = ['%e' % ma, '%.6f' % evidence_new[i]]
            s = " ".join(l) + "\n"
            # print(s)
            fout.write(s)
    
    evidence_new = np.concatenate(
        ([evidence_new[1]], evidence_new, [evidence_new[-2]]))  # trick to detect peaks at boundaries of array
    peaks_evi_new, _ = find_peaks(evidence_new, height=5)
    peaks_evi_new = peaks_evi_new - 1
    print('peaks_evi_new =', peaks_evi_new)
    evidence_new = evidence_new[1:-1]
    
    plt.plot(1e-3 * fluxes_sky_spectrum, label=r'$10^{-3}$ * (fluxes - baseline)')
    plt.plot(peaks_fluxes, 1e-3 * fluxes_sky_spectrum[peaks_fluxes], "x", color=colors[-1])
    plt.plot(evidence, label='evidence initial')  # , color=colors[1])
    plt.plot(evidence_new, label='evidence without peaks')  # , color=colors[2])
    plt.plot(peaks_evi_new, evidence_new[peaks_evi_new], "x", color=colors[-2])
    plt.plot(5 * np.ones_like(evidence), linewidth=0.6, color='k')
    plt.title(name_dw.replace("_", " "))
    # plt.plot(3*np.ones_like(evidence),linewidth=0.5, color="gray")
    plt.legend()
    # plt.ylim([-2, 10])
    plt.xlabel('channel \#')
    plt.show()
    plt.close()
    
def apply_LEE(name_dw, isCored):
    print('\n*** apply_LEE. Galaxy:', name_dw)
    np.set_printoptions(precision=6)
    

    # Load evindence
    if name_dw == 'eri_gru_hyd_leo_scu':
        filename = file_evidence_no_peaks_doppler_cored if isCored else file_evidence_no_peaks_doppler
    else:
        filename = file_evidence_no_peaks_cored if isCored else file_evidence_no_peaks
        
    data = genfromtxt(filename % name_dw, delimiter=' ')
    evidence = data[:, 1]
    evidence_new_LEE = np.empty_like(evidence)
    print('evidence.shape =', evidence.shape)
    
    # Apply Look elsewhere effect
    counter_channels = len(np.where(np.isnan(evidence))[0])
    print('counter_channels =', counter_channels)
    Ntrials = (n_spectrum - 2 - counter_channels) / 3
    min_significance = np.sqrt(2) * erfcinv(1 / Ntrials)
    
    if name_dw == 'eri_gru_hyd_leo_scu':
        filename = file_evidence_no_peaks_LEE_doppler_cored if isCored else file_evidence_no_peaks_LEE_doppler
    else:
        filename = file_evidence_no_peaks_LEE_cored if isCored else file_evidence_no_peaks_LEE

    with open(filename % name_dw, "w") as fout:
        
        for i in np.arange(len(evidence)):
            wl = channel_centers[i + 1]
            ma = np.round(7893.1756681904 *  pi / wl, decimals=6)  # in eV
            if np.isnan(evidence[i]):
                evidence_new_LEE[i] = evidence[i]
            elif evidence[i] < min_significance:
                # print(i, evidence_new[i], min_significance)
                evidence_new_LEE[i] = 0
            else:
                if evidence[i] < 30:
                    evidence_new_LEE[i] = np.sqrt(2) * erfcinv(erfc(evidence[i] / np.sqrt(2)) * Ntrials)
                else:
                    evidence_new_LEE[i] = np.sqrt(evidence[i] ** 2 - 2 * np.log(Ntrials) + 2 * np.log(evidence[i]) -  np.log(evidence[i] ** 2 - 2 * np.log(Ntrials) + 2 * np.log(evidence[i])))
            
            l = ['%e' % ma, '%.6f' % evidence_new_LEE[i]]
            s = " ".join(l) + "\n"
            # print(s)
            fout.write(s)
    

    idx = np.argpartition(-evidence_new_LEE, np.arange((~np.isnan(evidence_new_LEE)).sum()))
    evidence_tmp = evidence_new_LEE[idx]
    print('maximum evidence LEE =', evidence_tmp[:10], idx[:10])
    
    peaks_evi_new_LEE, _ = find_peaks(evidence_new_LEE, height=5)
    print('Peaks left at 5 sigma after LEE:', peaks_evi_new_LEE)
    
    # plt.plot(evidence, label='evidence initial')  # , color=colors[2])
    plt.plot(evidence_new_LEE, label='evidence final LEE')  # , color=colors[2])
    plt.plot(peaks_evi_new_LEE, evidence_new_LEE[peaks_evi_new_LEE], "x", color=colors[-2])
    plt.title(name_dw)
    plt.plot(5*np.ones_like(evidence),linewidth=0.5, color="gray")
    plt.legend()
    # plt.ylim([-2, 10])
    # plt.show()
    plt.close()

def make_evindence_combined_nopeaks_LEE(names_dSph, i_channels, isShifted):
    
    if len(i_channels) == 0:
        i_channels = range(1, n_spectrum - 1)
        
    print('\n\nMaking combined evidence')

    # Find channels to exclude
    idxs_bad_channels_tot = np.empty(len(names_dSph), dtype=object)
    bad_masses_tot = np.empty(len(names_dSph), dtype=object)
    counter_channels = np.empty(len(names_dSph))
    
    # Load evidence for single galaxy to find channels to exclude from combined evidence
    # If isShifetd, we load the Doppler shifted evidence, so we know which actual axion
    # mass corresponds to a bad channel
    i = 0
    for name_dw in names_dSph:
        filename = file_evidence_no_peaks_LEE_doppler if isShifted else file_evidence_no_peaks_LEE
        data = genfromtxt(filename % name_dw, delimiter=' ')
        masses = data[:, 0]
        evidence = data[:, 1]
        idxs_bad_channels = np.where(np.isnan(evidence))[0] + 1 # add 1 to get the actual channel number because we throw away the first channel
        print('Number of bad channels for %s: %d ' % (name_dw, idxs_bad_channels.shape[0]))
        idxs_bad_channels_tot[i] = idxs_bad_channels
        bad_masses_tot[i] = masses[idxs_bad_channels]
        counter_channels[i] = len(idxs_bad_channels)  # number of channels removed due to atmosphere peaks
        i += 1

    idxs_bad_channels_tot = idxs_bad_channels_tot.flatten()
    idxs_bad_channels_tot = np.unique(np.array([item for sublist in idxs_bad_channels_tot for item in sublist]))

    bad_masses_tot = bad_masses_tot.flatten()
    bad_masses_tot = np.unique(np.array([item for sublist in bad_masses_tot for item in sublist]))

    print('Number of bad channels combined: %d ' % idxs_bad_channels_tot.shape[0])

    # Make combined evidence
    filename = file_evidence_no_peaks_LEE_doppler if isShifted else file_evidence_no_peaks_LEE
    with open(filename % "_".join(names_dSph), "w") as fout:
        for i in i_channels:
            if i % skip_wl_every != 0:
                continue
            # if i != 1609:
            #     continue
            # print('\ni =',i)
            
            # Axion mass
            wl = channel_centers[i]
            ma = np.round(7893.1756681904 * pi / wl, decimals=6)  # in eV. This is the actual axion mass at which we're looking for evidence
            
            # Remove channels that are excluded in at least one of the galaxies
            if ma in bad_masses_tot: # if isShifted=True, these are the actual bad axion masses, else they are the bad observed masses
                # print('Skipping channel %d. NaN.' %i)
                l = ['%e' % ma, '%e' % np.nan]
                s = " ".join(l) + "\n"
                # print(i, s)
                fout.write(s)
                continue
            
            chi2null = 0
            chi2s_comb = np.zeros(N_gs)
            chi2nulls = np.empty(len(names_dSph))
            chi2_bests = np.empty([len(names_dSph), N_gs])
            i_observeds = np.empty(len(names_dSph))
            
            # chi2s is 2D array
            idw = 0
            for name_dw in names_dSph:

                # print('Evidence for dSph no peaks with LEE', name_dw)
                chi2s_bas = read_chain_binned(name_dw, isCored=False)[2]
                chi2null_bas = np.nanmin(chi2s_bas)
                
                if isShifted:
                    # ma is the mass emitted, the actual axion mass at which I want to put the bound
                    # mass_obs is what MUSE observes
                    # For each ma, I want to find which channel I have to consider
                    masses_extended =  np.round(7893.1756681904 *  pi / channel_centers_extended, decimals=6)  # in eV
    
                    mass_obs = ma * (1 - v_radial_dict[name_dw]) # the observed frequency is the axion mass is ma
    
                    # Observational channel corresponding to given axion mass
                    i_observed = find_nearest(masses_extended, mass_obs)[0] - 10 # -10 bc channels extended has 10 extra channels at the beginning. At the end too
                    # print(name_dw, i, i_observed)
                    if i_observed < 1 or i_observed > n_spectrum - 2:
                        # print('** Skipping channel %d for %s' % (i, name_dw))
                        i_observeds[idw] = np.nan
                        continue
                    elif name_dw == 'hyd' and i_observed == 884:
                        i_observeds[idw] = np.nan
                        continue
                    else:
                        i_actual = i_observed

                else:
                    if name_dw == 'hyd' and i == 884:
                        i_observeds[idw] = np.nan
                        continue
                    i_actual = i
                    
                chi2s = np.load(file_chi2_of_g % (name_dw, i_actual))['chi2s']
                chi2null_noi = read_chi2_best_permanent(name_dw, i_actual, isCored=False)[3]
                chi2null += chi2null_noi + chi2null_bas
                chi2s_comb += chi2s

                i_observeds[idw] = i_actual
                chi2nulls[idw] = chi2null_noi + chi2null_bas
                chi2_bests[idw] = chi2s
                
                idw += 1
            
            np.set_printoptions(suppress=True)
            
            chi2superbest = np.amin(chi2s_comb)  # chi2_best = np.zeros(N_g_sqrtrho)
            idx_chi2superbest = np.argmin(chi2s_comb)

            if chi2null < chi2superbest:  # no DM is favoured. No evidence.
                chi2superbest = chi2null
                # print('Changing chi2superbest. chi2_best_noi = %.2f, chi2null = %.2f' % (chi2superbest, chi2null))
            
            evidence_beforeLEE = np.sqrt(chi2null - chi2superbest)
            # print('\nBefore LEE i, evidence =', i, evidence)
            
            # Apply LEE
            Ntrials = (n_spectrum - 2 - 2358) / 3 # I got the 2356 by looking at how many nans I have in the evidence at the end
            min_significance = np.sqrt(2) * erfcinv(1 / Ntrials)
            # print(evidence, min_significance)
            
            if evidence_beforeLEE < min_significance:
                evidence = 0
            else:
                if evidence_beforeLEE < 30:
                    evidence = np.sqrt(2) * erfcinv(erfc(evidence_beforeLEE / np.sqrt(2)) * Ntrials)
                else:
                    evidence = np.sqrt(evidence_beforeLEE ** 2 - 2 * np.log(Ntrials) + 2 * np.log(evidence_beforeLEE) - np.log(
                        evidence_beforeLEE ** 2 - 2 * np.log(Ntrials) + 2 * np.log(evidence_beforeLEE)))
            # print('After LEE i, evidence =', i, evidence)

            if evidence >= 5 or i in [1609, 3054, 3056, 3057, 3058, 3549, 3582, 3583, 3585, 3586, 3587, 3588, 3593, 3594, 3597]:
                print('\ni, evidence, evidence_beforeLEE =', i, evidence, evidence_beforeLEE)
                idw = 0
                for name_dw in names_dSph:
                    print(name_dw, i_observeds[idw],
                          # chi2nulls[idw], chi2_bests[idw][idx_chi2superbest],
                          chi2nulls[idw]-chi2_bests[idw][idx_chi2superbest],
                          chi2nulls[idw]-chi2_bests[idw][np.nanargmin(chi2_bests[idw])])
                    
                    # if name_dw == 'hyd':
                    #     print(i_observeds[idw])
                    log_gs = np.load(file_chi2_of_g % (name_dw, i_actual))['log_gs']
                    idx_chi2best = np.nanargmin(chi2_bests[idw])
                    chi2_best = chi2_bests[idw][idx_chi2best]
                    
                    plt.plot(log_gs, chi2_bests[idw] - chi2_best, label=name_dw, color=colors[idw])
                    
                    plt.plot(log_gs[idx_chi2superbest], chi2_bests[idw][idx_chi2superbest]-chi2_best, marker='.', color=colors[idw])
                    plt.plot(log_gs[idx_chi2best], 0, marker='*', color=colors[idw])

                    plt.plot(-13.6, chi2nulls[idw]-chi2_best, marker='d', color=colors[idw])


                    idw += 1
                plt.legend()
                # plt.yscale('log')
                plt.xlim([-14, -11.5])
                plt.ylim([-5,55])
                # plt.show()
            l = ['%e' % ma, '%e' % evidence]
            s = " ".join(l) + "\n"
            # print(i, s)
            fout.write(s)
    
    print('Done.')


############################# FUNCTIONS TO CALL ###########################################
def get_likelihood(channels, isCored):

    make_chi2list(['eri1'], channels, isMakeEvidence=False, isCored=isCored)
    make_chi2list(['eri1'], channels, isMakeEvidence=True, isCored=isCored)
    make_chi2list(['eri'], channels, isMakeEvidence=False, isCored=isCored)
    make_chi2list(['eri2'], channels, isMakeEvidence=True, isCored=isCored)
    merge_eridanus(['eri1', 'eri2'], channels)
    make_chi2list(['leo'], channels, isMakeEvidence=False, isCored=isCored)
    make_chi2list(['leo'], channels, isMakeEvidence=True, isCored=isCored)
    make_chi2list(['gru'], channels, isMakeEvidence=False, isCored=isCored)
    make_chi2list(['gru'], channels, isMakeEvidence=True, isCored=isCored)

def prepare_chains():
    format_chain('gru')
    format_chain('eri')
    format_chain('leo')
    format_chain('hyd')
    format_chain_scu()
    
    for name_dw in ['hyd', 'leo', 'gru', 'scu', 'eri']:
        bin_chain(name_dw)

def prepare_chains_core():
    format_chain_core('gru')
    format_chain_core('eri')
    format_chain_core('leo')
    format_chain_core('hyd')
    format_chain_core_scu()
    for name_dw in ['hyd', 'leo', 'gru', 'scu', 'eri']:
        bin_chain_core(name_dw)
        
def make_bounds(channels, isCored):

    merge_eridanus(['eri1', 'eri2'], channels, isCored=isCored)
    make_bounds_single(['eri'], channels, isCored=isCored)
    make_bounds_single(['gru'], channels, isCored=isCored)
    make_bounds_single(['hyd'], channels, isCored=isCored)
    make_bounds_single(['leo'], channels, isCored=isCored)
    make_bounds_single(['scu'], channels, isCored=isCored)
    doppler_shift_single(['leo'], isCored=isCored, isEvidence=False)
    
    doppler_shift_single(['eri', 'gru', 'hyd', 'leo', 'scu'], isCored=isCored, isEvidence=False)

    make_bounds_combined(['eri', 'gru', 'hyd', 'leo', 'scu'], channels, isCored=isCored, isShifted=True)
    
def make_evidence(channels, isCored):

    subtract_threshold_from_skyspectrum()
    for name_dw in ['eri', 'gru', 'hyd', 'leo', 'scu']:
        make_evindence_single([name_dw], channels, isCored=isCored)
        find_bad_channels_using_skyspec(name_dw, isCored=isCored)
        apply_LEE(name_dw, isCored=isCored)
        doppler_shift_single([name_dw], isCored=isCored, isEvidence=True)

    make_evindence_combined_nopeaks_LEE(['eri', 'gru', 'hyd', 'leo', 'scu'], [], isShifted=True)
    
################################################################################################
################################################################################################
# ################################################################################################
channels = [] # an empty list means use all channels
get_likelihood(channels=channels, isCored=False) # isCored=False --> NWF, isCored=True --> cored profile
# prepare_chains()
# prepare_chains_core()
# make_bounds(channels=channels, isCored=False)
# make_evidence(channels=channels, isCored=False)













