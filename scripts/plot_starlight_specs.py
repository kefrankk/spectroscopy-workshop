
import os
import re
import gc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.ndimage import label
from pathlib import Path
import utils


font = 12
plt.rcParams.update({
    'font.family': 'DejaVu Sans',  # ou Helvetica se estiver instalada
    'text.usetex': True,
    'font.size': font,
    'axes.titlesize': font,
    'axes.labelsize': font,
    'xtick.labelsize': font,
    'ytick.labelsize': font,
    'legend.fontsize': font,
    'figure.titlesize': font,
    'font.sans-serif': ['Helvetica']
})


cluster = '1_279073'

filepath = f'../data/starlight/'

files = [f for f in os.listdir(filepath) if os.path.isfile(os.path.join(filepath, f)) and '.' not in f]


files.sort()

for file in files:
        results_all = {}

        galaxy_name = file.split('_')[0]
        library = '_'.join(file.split('_')[-2:])

        if file+'_starlight.png' in os.listdir(filepath):
            # print(f'File {file} already processed, skipping...')
            continue
        print(f'Processing {file}...')

        if galaxy_name not in results_all:
            results_all[galaxy_name] = {}

        reading_spectrum = False
        lamb, f_obs, f_model, weight = [], [], [], []
        population_data = []
        spectrum = []
        metadata = {}
        with open(filepath+file, 'r', encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()

            if not lines[0].startswith("#"):
                continue


        for idx, line in enumerate(lines):
            line = line.strip()

            if '[' in line and ']' in line:
                value, description = line.split('[')
                value = value.strip()
                description = description.strip('[]').strip()
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

        results_all[galaxy_name]['chi2'] = metadata.get('chi2/Nl_eff')
        results_all[galaxy_name]['adev'] = metadata.get('adev (%)')
        results_all[galaxy_name]['AV_min'] = metadata.get('AV_min  (mag)')
        results_all[galaxy_name]['flux_tot'] = metadata.get('Flux_tot (units of input spectrum!)')
        results_all[galaxy_name]['Mini_tot'] = metadata.get('Mini_tot (???)') 
        results_all[galaxy_name]['Mcor_tot'] = metadata.get('Mcor_tot (???)') 
        results_all[galaxy_name]['v0_min'] = metadata.get('v0_min  (km/s)')
        results_all[galaxy_name]['vd_min'] = metadata.get('vd_min  (km/s)')
        results_all[galaxy_name]['q_norm'] = metadata.get('q_norm = A(l_norm)/A(V)')
        results_all[galaxy_name]['fobs_norm'] = metadata.get('fobs_norm (in input units)')



        

        fobs_norm = metadata.get('fobs_norm (in input units)')

        lam         = np.array(lamb) 
        f_obs       = np.array(f_obs)  * float(fobs_norm)
        f_model     = np.array(f_model) * float(fobs_norm)
        gas         = np.array(np.array(f_obs) - np.array(f_model)) 
        weight     = np.array(weight) * float(fobs_norm)




        #limits_list = [(lam > 3750) & (lam < 4050), (lam > 4200) & (lam < 4500), (lam > 6450) & (lam < 6750)]
        limits_list = [(lam > 4800) & (lam < 5100), (lam > 6100) & (lam < 6400), (lam > 6450) & (lam < 6750)]

        # Identify emission line intervals
        w0      = np.array(weight) <= 0
        labels, num_features = label(w0)

        intervals = []
        for i in range(1, num_features + 1):
            indices = np.where(labels == i)[0]
            inicio = lam[indices[0]]
            fim = lam[indices[-1]]
            intervals.append((inicio, fim))


        fig = plt.figure(figsize=(14, 8))

        gd = mpl.gridspec.GridSpec(2, 1, figure=fig, height_ratios=[3, 1])
        gs_top = mpl.gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gd[0], height_ratios=[3, 1], hspace=0)

        ax_top = fig.add_subplot(gs_top[0, :])
        ax_bottom = fig.add_subplot(gs_top[1, :], sharex=ax_top)

        # Principal plot
        ax_top.plot(lam, f_obs, label='Observed', lw=1, color='k')
        ax_top.plot(lam, f_model, label='Model', lw=2, color='r')
        ax_top.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
        for  i, (inicio, fim) in enumerate(intervals):
            labels = 'Masked region' if i == 0 else None
            ax_top.axvspan(inicio, fim, color='gray', alpha=0.3, label=labels)
        ax_top.legend()
        #ax_top.set_ylim(0, 0.2)

        ax_top.set_ylabel(r'flux (10$^{-17}$ erg s$^{-1}$ cm$^{-2}$ $\AA$)')

        # Observed - stellar
        ax_bottom.plot(lam, gas, lw=1, color='gray')
        ax_bottom.set_ylabel('residue')
        ax_bottom.set_ylim(-0.05, 0.05)
        ax_bottom.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
        ax_bottom.set_xlabel(r'rest wavelength ($\AA$)')


        gs_bot = mpl.gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gd[1])

        for i, limit in enumerate(limits_list):
            
            ax = fig.add_subplot(gs_bot[i])
            ax.plot(lam[limit], f_obs[limit], label='Observed', lw=1, color='k' )
            ax.plot(lam[limit], f_model[limit], label='Model', lw=1, color='r')
            ax.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
            ax.set_xlabel(r'rest wavelength ($\AA$)')
            xmin = lam[limit].min() 
            xmax = lam[limit].max() 
            ax.set_xlim(xmin, xmax)           
            for  i, (inicio, fim) in enumerate(intervals):
                labels = 'Masked region' if i == 0 else None
                ax.axvspan(inicio, fim, color='gray', alpha=0.3, label=labels)

            if i == 0:
                # ax.legend()
                ax.set_ylabel(r'flux (10$^{-17}$ erg s$^{-1}$ cm$^{-2}$ $\AA$)')


        fig.tight_layout()
        plt.savefig(f'{filepath}{file}_starlight.png')
        #plt.show()
        plt.close(fig)
        #break


        del fig, results_all
