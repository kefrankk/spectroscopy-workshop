
##########################################################
#
# Author: Kelly F. Heckler
# Year: 2026.1
#
##########################################################

import os
import re
from pathlib import Path
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


def generate_input_list(
    name: str,
    extension: str = ".txt",
    gmos: bool = True,
    library: str = 'CB19_16x5',
    mask: str = 'mask_gmos.gm',  
    obs_dir: str = './data/', 
):
    """ 
    Generate a STARLIGHT input list.

    The function scans the directory containing the observed spectra and
    creates a STARLIGHT input list with one entry per spectrum. The stellar
    population template is selected according to the chosen library.
    
    Parameters 
    ---------- 
    name : str 
        Galaxy name (e.g., 'NGC'). 

    extension : str 
        File extension to include (default: ".txt"). 

    gmos : bool 
        If True, use GMOS configuration. 
        
    library : str 
        Template library name (e.g., 'CB19_16x5'). 
        
    mask : str 
        Mask file name (e.g., 'mask_gmos.gm'). 

    obs_dir : str
        Directory with the data.
    """

    directory = Path(obs_dir)
    output_list = Path(f"./data/{name}_input_list.txt")

    # Select template based on library
    if gmos:
        if library == 'EMILES_Padova_240':
            template = 'base_EMILES_Padova_KU_all_ages.txt'
        elif library == 'EMILES_Padova_38':
            template = 'base_EMILES_Padova_KU_SafeRanges_for_IR.txt'
        elif library == 'EMILES_Basti_120':
            template = 'base_EMILES_Basti_KU_SafeRanges_for_IR.txt'
        elif library == 'EMILES_Basti_250':
            template = 'base_EMILES_Basti_KU_all_ages_gt_0pt06.txt'
        elif library == 'CB19_16x5':
            template = "CBASE.PARSEC.chab.16x5.all"
        else:
            raise ValueError("Unknown library provided")
    else:
        raise ValueError("Only GMOS configuration is implemented")

    # Fixed components
    fixed_components = ["StCv04.C11.config", template, mask, "CCM", "0.0", "150.0"]

    processed = {f.name for f in directory.iterdir()}

    # Write output file
    with output_list.open("w") as f:

        for file in sorted(directory.glob("*" + extension)):
            file_name = file.stem
            out_name = f"{file_name}_{library}"


            # Skip spectra that already have a corresponding output placeholder.
            # This is used, for example, for spaxels with zero flux, for which a
            # dummy STARLIGHT output is created beforehand.
            if out_name in processed:
                continue

            line = f"{file.name}  {'  '.join(fixed_components)}  {out_name}"
            f.write(line + "\n")

    print(f"Output file created: {output_list}")


def read_starlight_output(file: str, filepath: str) -> dict:
    """
    Parse and extract information from STARLIGHT output files.

    This function reads one or more STARLIGHT output files, extracting both
    metadata and synthetic spectrum information. For each galaxy (identified
    by the prefix of the filename), the function stores:
    
    - Metadata key-value pairs parsed from lines in the format `value [description]`.
    - Synthetic spectrum data, including wavelength, observed flux, model flux,
      and weights.

    Parameters
    ----------
    files : list of str
        List of STARLIGHT output filenames to be read.
        Each filename is expected to start with the galaxy name
        (e.g., ``NGC1234_output.txt`` → galaxy name = ``NGC1234``).
    filepath : str
        Path to the directory containing the STARLIGHT output files.

    Returns
    -------
    dict
        A nested dictionary with the following structure:
        
        {
            galaxy_name: {
                'metadata': dict
                    Key-value pairs extracted from the header section.
                'wavelenght': list of float
                    Wavelength values of the synthetic spectrum.
                'f_obs': list of float
                    Observed flux values.
                'f_model': list of float
                    Model flux values.
                'weight': list of float
                    Weights associated with each wavelength point.
            },
            ...
        }
    """

    results_all = {}
    reading_spectrum = False

    galaxy_name = file.split('_')[0]

    #results_all = {
    #    'metadata': {},
    #    'wavelenght': [],
    #    'f_obs': [],
    #    'f_model': [],
    #    'weight': [],
    #    'residual': []
    #}

    metadata = {}
    lamb, f_obs, f_model, weight = [], [], [], []

    with open(filepath+file, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()

        if lines[0].startswith("#"):

            for idx, line in enumerate(lines):
                line = line.strip()
    
                match = re.match(r'^(.*?)\[(.*)\]', line)
                if match:
                    raw_values, description = match.groups()
                    raw_values = raw_values.strip()
                    description = description.strip()
            
                    values = raw_values.split()
                    if len(values) == 1:
                        value = values[0]  # valor único como string
                    else:
                        value = values  # lista de strings
            
                    metadata[description] = value
    
                elif line.startswith('## Synthetic spectrum '):  # Início dos espectros # elif line:  # Início dos espectros
                    reading_spectrum = True
                    continue
                
                elif reading_spectrum:
                    parts = line.split()
                    if len(parts) == 4:
                        spectrum = list(map(float, parts))
                        lamb.append(spectrum[0])
                        f_obs.append(spectrum[1]) #* 10**-17  # Convertendo para erg/s/cm²/Å
                        f_model.append(spectrum[2]) #* 10**-17  # Convertendo para erg/s/cm²/Å
                        weight.append(spectrum[3])  # is 1/error 
                    else:
                        reading_spectrum = False  # Parar de ler se a linha mudar o formato
                    continue  # Já tratou essa linha, então pula pro próximo loop

        if not lines[0].startswith("#"):
            wave_fixo = np.arange(4310, 7300 + 1, step=1)
            nwave = len(wave_fixo)
            model = np.zeros(nwave)
            #print(nwave, model)

            #for i in lines:
            #lamb.append(wave_fixo)
            #f_obs.append(model)
            #f_model.append(model)
            #weight.append(model)

    if not metadata:
        #print(f"Warning: No metadata found in file {file}. Using default values.")
        #print(min(np.array(lamb)), max(np.array(lamb)))
        results_all['chi2']     = 0
        results_all['adev']     = 0
        results_all['AV_min']   = 0
        results_all['flux_tot'] = 0
        results_all['Mini_tot'] = 0
        results_all['Mcor_tot'] = 0
        results_all['v0_min']   = 0
        results_all['vd_min']   = 0
        results_all['q_norm']   = 0
        results_all['fobs_norm']= 0
        results_all['l_ini']    = 0
        results_all['dl']       = 0

        results_all['wavelength']   = wave_fixo
        results_all['observed']     = model
        results_all['stellar']      = model
        results_all['gas']          = model
        results_all['weight']       = model

    elif metadata:
        #print(min(np.array(lamb)), max(np.array(lamb)))

        
        results_all['chi2']     = float(metadata.get('chi2/Nl_eff'))
        results_all['adev']     = float(metadata.get('adev (%)'))
        results_all['AV_min']   = float(metadata.get('AV_min  (mag)'))
        results_all['flux_tot'] = float(metadata.get('Flux_tot (units of input spectrum!)'))
        results_all['Mini_tot'] = float(metadata.get('Mini_tot (???)') )
        results_all['Mcor_tot'] = float(metadata.get('Mcor_tot (???)') )
        if metadata.get('v0_min  (km/s)').startswith('*'):
            results_all['v0_min'] = np.nan
        else:
            results_all['v0_min']   = float(metadata.get('v0_min  (km/s)'))
        results_all['vd_min']   = float(metadata.get('vd_min  (km/s)'))
        results_all['q_norm']   = float(metadata.get('q_norm = A(l_norm)/A(V)'))
        results_all['fobs_norm']= float(metadata.get('fobs_norm (in input units)'))
        results_all['l_ini']    = float(metadata.get('l_ini (A)'))
        results_all['dl']       = float(metadata.get('dl    (A)'))

        fobs_norm = metadata.get('fobs_norm (in input units)')

        results_all['wavelength']   = np.array(lamb) 
        results_all['observed']     = np.array(f_obs)  * float(fobs_norm)
        results_all['stellar']      = np.array(f_model) * float(fobs_norm)
        results_all['gas']          = np.array(np.array(results_all['observed']) - np.array(results_all['stellar']))
        results_all['weight']       = np.array(weight) * float(fobs_norm)


    return results_all

