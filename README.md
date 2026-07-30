# 🌌 3D Spectroscopy of Galaxies Workshop

This project is part of the **Workshop: 3D Spectroscopy of Galaxies**, which aims to train, primarily undergraduate and graduate students, in the use of tools and methodologies to study galaxies with Integral Field Spectroscopy data.

📊 The project covers the full workflow of data reduction and analysis:
    
- Reduction with [IRAF](https://iraf.readthedocs.io/en/latest/index.html)

- Stellar population analysis with [Starlight](http://www.starlight.ufsc.br/)

- Gas property determination with [Ifscube](https://github.com/danielrd6/ifscube)

👩‍🏫 Organized by: Gabriele Ilha, Kelly Heckler, Michele Bertoldo Coelho, and Angela C. Krabbe



## 📂 Repository structure

```
spectroscopy-workshop/
├── data/             # Input data and STARLIGHT files
├── notebooks/        # Jupyter notebooks
├── scripts/          # Python modules
├── install_starlight.sh
├── environment.yml
└── README.md
```

## 💻 How to Use

### 1. Clone the repository

```bash
git clone https://github.com/kefrankk/spectroscopy-workshop.git
cd spectroscopy-workshop
```

## 2. Create the Conda environment

```bash
conda env create -f environment/science.yml
```

Activate it:

```bash
conda activate science
```


## 3. Install STARLIGHT

Run

```bash
bash install_starlight.sh
```

This downloads and installs the STARLIGHT package into

```
data/STARLIGHTv04/
```
and move all not used files to 

```
data/STARLIGHTv04/examples/
```


## 🔧 Data Reduction with IRAF

IRAF documentation is avaliable in [here](https://iraf.readthedocs.io/en/latest/index.html)



 
## 🌟  Starlight

Starlight is employed to analyze the properties of stellar populations and to subtract their contribution from the observed spectrum. This spectral synthesis approach enables a more accurate characterization of galaxies by disentangling the stellar component from other physical processes present in the data.

The Starlight Manual can be found [here](https://minerva.ufsc.br/starlight/files/papers/Manual_StCv04.pdf).

To install Starlight, run the following command:

``` 
./install_starlight.sh 
```
The code has been automated to perform the installation in the current folder.


📓 “Check the notebook for a guided tutorial on Starlight fitting.”

🐍 “Use the scripts folder for modular tools to automate specific tasks.”

---

