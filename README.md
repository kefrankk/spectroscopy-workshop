# 3D Spectroscopy of Galaxies Workshop

This project is part of the **Workshop: 3D Spectroscopy of Galaxies**, which aims to train, primarily undergraduate and graduate students, in the use of tools and methodologies to study galaxies with Integral Field Spectroscopy data.

In this project, we will carry out the complete process of data reduction and analysis.
The reduction will be performed entirely using the [IRAF](https://iraf.readthedocs.io/en/latest/index.html) (Image Reduction and Analysis Facility) software, a well-established tool for astronomical data calibration and preprocessing.
The subsequent analysis will be conducted with two complementary codes: [Starlight](http://www.starlight.ufsc.br/), employed to derive the properties of stellar populations, and [Ifscube](https://github.com/danielrd6/ifscube), used to determine the gas properties of the science objects.

## Data Reduction with IRAF

IRAF documentation is avaliable in [here](https://iraf.readthedocs.io/en/latest/index.html)



 
## Starlight

Starlight is employed to analyze the properties of stellar populations and to subtract their contribution from the observed spectrum. This spectral synthesis approach enables a more accurate characterization of galaxies by disentangling the stellar component from other physical processes present in the data.

The Starlight Manual can be found [here](https://minerva.ufsc.br/starlight/files/papers/Manual_StCv04.pdf).

To install Starlight, run the following command:

``` 
./install_starlight.sh 
```
The code has been automated to perform the installation in the current folder.


---

