##########################################################
#
# Author: Kelly F. Heckler
# Year: 2026.1
#
##########################################################

import os
import utils
from pathlib import Path

name = '1_279073'

directory = Path.home() / f"1_279073"

extension = ".txt" 

output_list = Path(f"{name}_input_list.txt")

gmos = True

if gmos:
    library = 'CB19_16x5'
    mask = "mask_gmos_workshop.gm"
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
    elif library == 'XSL':
        template = 'base_XSL_Kroupa_PC_AGN_LLAMA.txt'

# Componentes fixos
fixos = ["StCv04.C11.config", template, 
         mask, "CCM", "0.0", "150.0"]


with output_list.open("w") as f:
    for file in sorted(directory.glob("*" + extension)):
        file_name = file.stem  # nome do arquivo sem extensão
        out_name = f"{file_name}_{library}"
        #print(out_name)

        if out_name in os.listdir(directory):
            continue
        linha = f"{file.name}  {'  '.join(fixos)}  {out_name}"
        f.write(linha + "\n")
