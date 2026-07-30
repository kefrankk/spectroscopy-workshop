
##########################################################
#
# Author: Kelly F. Heckler
# Year: 2026.1
#
##########################################################

import os
import re
import numpy as np 
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


def wavelength_axis(header, axis=3):
    """
    Calculate the wavelength axis for a given FITS header.

    Parameters:
    header (astropy.io.fits.Header): FITS header containing wavelength calibration keywords.
    axis (int): The axis number for which to calculate the wavelength array. Default is 3.

    Returns:
    numpy.ndarray: Calculated wavelength values for the specified axis.
    """
    # Get reference pixel, value, and increment from the header
    crpix = header[f'CRPIX{axis}']
    crval = header[f'CRVAL{axis}']
    cdelt = header[f'CD3_{axis}']
    #cdelt = header[f'CDELT{axis}']
    naxis = header[f'NAXIS{axis}']
    
    # Create an array of pixel indices
    index = np.arange(naxis)
    
    # Calculate and return the wavelength values
    return crval + (index + 1 - crpix) * cdelt


def save_spectrum_txt(galaxy_name, wave, redshift, hdul, out_dir):
    """
    Extracts the spectrum of each spaxel from a data cube and saves it as an
    individual text file.

    The observed wavelength axis is converted to the rest frame using the
    galaxy redshift. Each spectrum is interpolated onto a linear wavelength
    grid with a sampling of 1 Å. Any NaN values resulting from the
    interpolation are replaced with zeros before saving.

    Parameters
    ----------
    galaxy_name : str
        Name of the galaxy, used as the prefix for the output filenames.
    wave : ndarray
        Observed wavelength array.
    redshift : float
        Galaxy redshift used to compute the rest-frame wavelength.
    hdul : astropy.io.fits.HDUList
        FITS file containing the spectral data cube.
    out_dir : str
        Directory where the output text files will be saved.

    Notes
    -----
    - Each output file contains two columns: rest-frame wavelength (Å) and
      interpolated flux.
    - Existing files are skipped to avoid redundant processing.
    - A reference spaxel is plotted for visual inspection of the rest-frame
      correction.
    - Spectra containing only NaN values are identified separately.
    """
    sci = hdul[0].data
    
    nz, ny, nx = sci.shape

    rest_wave = wave / (1 + redshift)

    count = 0
    for y in range(ny):
        for x in range(nx):
            # Nome do arquivo
            filename = f"{galaxy_name}_{y}_{x}.txt"
            filepath = os.path.join(out_dir, filename)

            if os.path.exists(filepath):
                #print(f"File {filename} already exists. Skipping...")
                continue

            flux = sci[:, y, x]


            # Cria novo eixo de comprimento de onda com passo de 1 Å
            wavelength_linear = np.arange(np.ceil(rest_wave.min()), np.floor(rest_wave.max()) + 1, 1)

            # Interpola o fluxo
            interp_flux = interp1d(rest_wave, flux, kind='linear', bounds_error=False, fill_value="extrapolate")

            flux_linear = interp_flux(wavelength_linear)

            onde_tem_nan = np.where(np.isnan(flux_linear))
            flux_linear[onde_tem_nan] = 0

            # Organiza em colunas
            data = np.column_stack((wavelength_linear, flux_linear))#, error_linear))

            # faz um plot de um spaxel qualquer, para conferir seo espectro está no repouso
            if y == 25 and x == 15:
                plt.plot(wavelength_linear, flux_linear)
                plt.xlabel("Wavelength (Å)")
                plt.ylabel("Flux")
                plt.title(f"Interpolated Spectrum at ({x}, {y}) for {galaxy_name}")
                plt.show()


            # Salva como texto com cabeçalho
            np.savetxt(filepath, data, fmt="%.6e", comments='')

            if np.isnan(flux).all():
                print(f"NaN values found in {filename}")
                new_file_path = os.path.join(out_dir, f"{galaxy_name}_{y}_{x}_CB19_16x5")
                np.savetxt(new_file_path, data, fmt="%.6e", comments='')
            else:
                count += 1

            #print(f"Arquivo salvo: {filepath}")

    print(f"Total de arquivos salvos: {count}")


def sort_by_pattern(file):
    name = file.stem  # sem extensão
    
    # Padrão SDSS: spec-PLATE-MJD-FIBERID
    match_sdss = re.match(r"spec-(\d+)-(\d+)-(\d+)", name)
    if match_sdss:
        plate, mjd, fiberid = map(int, match_sdss.groups())
        return (0, plate, mjd, fiberid)  # "0" força SDSS a vir antes
    
    # Padrão UUID-like
    match_uuid = re.match(r"spec-[0-9a-f]{8}-", name)
    if match_uuid:
        return (1, name)  # "1" garante que UUIDs fiquem depois
    
    # fallback: string normal
    return (2, name)



def identificar_arquivo(path):
    with open(path, "r") as f:
        primeira_linha = f.readline().strip()
    
    # caso seja o arquivo com cabeçalho
    if primeira_linha.startswith("###"):
        return True
    else:
        return False

