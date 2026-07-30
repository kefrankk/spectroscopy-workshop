
##########################################################
#
# Author: Kelly F. Heckler
# Year: 2026.1
#
##########################################################


import os
from astropy.io import fits
import utils


if __name__ == "__main__":

    redshifts = {
    '1_279073': 0.0323
    }

    cwd = os.getcwd()
    galaxy_name = '1_279073'
    redshift = redshifts[galaxy_name]
    fits_file = cwd + f"/Manga_{galaxy_name}_pca.fits"

    with fits.open(fits_file) as hdul:
        data = hdul[0].data
        wave = utils.wavelength_axis(hdul[0].header)
        print(len(wave))

        
        out_dir = os.path.join(cwd, f"{galaxy_name}/")
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        utils.save_spectrum_txt(galaxy_name, wave, redshift, hdul, out_dir=out_dir)




