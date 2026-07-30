
import os
import re
import copy
import numpy as np
import pandas as pd
from astropy.io import fits
import matplotlib.pyplot as plt
from ifscube.io import line_fit
import matplotlib as mpl
#import figures_plot as f
import glob
from function_spectra_io import read_starlight_output, safe_float
from astropy.wcs import WCS


galaxy = '1_279073'

galaxy_directory = f'../data/Manga_{galaxy}_pca.fits'
fits_file = fits.open(galaxy_directory)

header = fits_file[0].header
data = fits_file[0].data
nz = data.shape[0]
nx, ny = data.shape[2], data.shape[1]  # salva as dimensões espaciais
sl = fits_file[0].header['crpix1']

starlight_files = glob.glob(f'../data/starlight_old/*XSL')
files_path = f'../data/starlight_old/'

print(starlight_files[0])
example = read_starlight_output(starlight_files[0].split("/")[-1], files_path)
wave = example['wavelength']
nwave = len(wave)
wave_min = wave[0]
wave_max = wave[-1]
print(f'wave min: {wave_min}, wave max: {wave_max}')

cubo     = np.full((nwave, ny, nx), np.nan)
observed = np.full((nwave, ny, nx), np.nan)
stellar  = np.full((nwave, ny, nx), np.nan)
gas      = np.full((nwave, ny, nx), np.nan)
weight   = np.full((nwave, ny, nx), np.nan)
stellar_velocity    = np.full((ny, nx), np.nan)
stellar_dispersion  = np.full((ny, nx), np.nan)
av_min              = np.full((ny, nx), np.nan)
chi2                = np.full((ny, nx), np.nan)


cnt         = np.nanmedian(data[1000:1500,:,:], axis=0)
arr2D   = cnt.copy()
arr2D   = arr2D / np.nanmean(arr2D) # para trabalhar com numeros menores, o ideal seria ser magnitudes
resultado  = np.where(arr2D == np.amax(arr2D))
print(resultado)
ycen = resultado[0][0]
xcen = resultado[1][0]
print(xcen,ycen)
ycen = round(ycen)
xcen = round(xcen)#

# calcular a distancia a partir do nucleo
rad = data[0,:,:] * 0.0
nx  = len(rad[0,:])
ny  = len(rad[:,0])
for i in range(0,nx):
    for j in range(0,ny):
        rad[j,i] = ((j - ycen)**2 + (i - xcen)**2)**0.5 * sl


rows = []
for arq in starlight_files:
    print(arq)
    nome = arq.split("/")[-1]
    partes = nome.split("_")
    y = int(partes[2])
    x = int(partes[3])

    result = read_starlight_output(nome, files_path)

    gas[:, y, x]        = result['gas']
    observed[:, y, x]   = result['observed']
    weight[:, y, x]     = result['weight']
    stellar[:, y, x]    = result['stellar']
    stellar_velocity[y, x]   = result['v0_min']
    stellar_dispersion[y, x] = result['vd_min']
    av_min[y, x]             = result['AV_min']
    chi2[y, x]               = result['chi2']

    gas[:,(rad > 2.5)] = np.nan
    observed[:,(rad > 2.5)] = np.nan
    weight[:,(rad > 2.5)] = np.nan
    stellar[:,(rad > 2.5)] = np.nan
    stellar_velocity[(rad > 2.5)] = np.nan
    stellar_dispersion[(rad > 2.5)] = np.nan
    av_min[(rad > 2.5)] = np.nan
    chi2[(rad > 2.5)] = np.nan





    rows.append((
        x, y,
        safe_float(result['chi2']),
        safe_float(result['adev']),
        safe_float(result['AV_min']),
        safe_float(result['flux_tot']),
        safe_float(result['Mini_tot']),
        safe_float(result['Mcor_tot']),
        safe_float(result['v0_min']),
        safe_float(result['vd_min']),
        safe_float(result['q_norm']),
        safe_float(result['fobs_norm']),
        safe_float(result['l_ini']),
        safe_float(result['dl'])
    ))


    #break


# converter tabela para structured array
dtype = [
    ('x', 'i4'), ('y', 'i4'),
    ('chi2', 'f8'), ('adev', 'f8'), ('AV_min', 'f8'),
    ('flux_tot', 'f8'), ('Mini_tot', 'f8'), ('Mcor_tot', 'f8'),
    ('v0_min', 'f8'), ('vd_min', 'f8'),
    ('q_norm', 'f8'), ('fobs_norm', 'f8'),
    ('l_ini', 'f8'), ('dl', 'f8')
]
table_data = np.array(rows, dtype=dtype)


cdelt = 1
crpix = 1
crval = wave[0]

header['NAXIS3']  = len(wave)
header['CRPIX3']  = crpix
header['CD3_3']  = cdelt
header['CRVAL3']  = crval
header['CTYPE3']  = 'WAVE'


# salvar em FITS com múltiplas extensões
hdu0 = fits.PrimaryHDU()  # cabeçalho vazio
hdu1 = fits.ImageHDU(observed, header=header, name="OBSERVED")
hdu2 = fits.ImageHDU(stellar, header=header, name="STELLAR")
hdu3 = fits.ImageHDU(gas, header=header, name="GAS")
hdu4 = fits.ImageHDU(weight, header=header, name="WEIGHT")
hdu5 = fits.ImageHDU(stellar_velocity, header=header, name="STELLAR_VELOCITY")
hdu6 = fits.ImageHDU(stellar_dispersion, header=header, name="STELLAR_DISPERSION")
hdu7 = fits.ImageHDU(av_min, header=header, name="AV_MIN")
hdu8 = fits.ImageHDU(chi2, header=header, name="CHI2")
hdu_table = fits.BinTableHDU(table_data, name="PARAMS")

hdul = fits.HDUList([hdu0, hdu1, hdu2, hdu3, hdu4, hdu5, hdu6, hdu7, hdu8, hdu_table])
hdul.writeto(f"../data/Manga_{galaxy}_starlight.fits", overwrite=True)








