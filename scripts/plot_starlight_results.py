##########################################################
#
# Author: Michele B. Coelho
# Adapted by: Kelly Heckler
# Year: 2026.1
#
##########################################################


import os
import copy
import glob
import numpy as np
import pandas as pd
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Altere para True quando tiver poucos spaxels
spec_spaxel = False

galaxia = 'Manga_1_279073'
datacube = fits.open(f"../data/{galaxia}_pca.fits")

header = datacube[0].header
data = datacube[0].data
nz = data.shape[0]
nx, ny = data.shape[2], data.shape[1]  # salva as dimensões espaciais
sl = datacube[0].header['crpix1']


# Pegando os outputs do starlight
arqs = glob.glob(f'../data/starlight_old/*XSL')


# Header dos espectros
headerSSP = '# j     x_j(%)      Mini_j(%)     Mcor_j(%)     age_j(yr)     Z_j      (L/M)_j    YAV?  Mstars   component_j        a/Fe...       SSP_chi2r SSP_adev(%)   SSP_AV   SSP_x(%)'
header = '## Synthetic spectrum (Best Model) ##l_obs f_obs f_syn wei'

# Iniciar o data frame com as informações
df = pd.DataFrame(columns=['Galaxy', 'x', 'y', 'chi2', 'adev', 'Mini', 'Mcor',
                           'v0', 'vd', 'AV', 'xy', 'xiy', 'xio', 'xo',
                           '<logt>x', '<t>x', '<Z>x', '<[M/H]>x', 'FCx', 'HDx',
                           '<logt>m', '<t>m', '<Z>m', '<[M/H]>m'])

# Iniciar dicionários para apoiar a construção do cubo
# Para salvar como cubo de dados todos os l precisam ser iguais
#dados_l = {}
l_wave = None
dados_fobs = {}
dados_fsyn = {}
dados_wei = {}

# Definir a metalicidade utilizada nos modelos
Zsun = 0.0207

# Loop para cada arquivo
for arq in arqs:
   
   #print("Recovering the spectral fitting informations for: ", arq, '\n')
   
   # Zera o dicionário para cada spaxel
   dados_arquivo = {}

   # Recuperar a posicao do espectro no cubo
   # Depende de como foi nomeado, para exemplificar estou 
   # considerando que o nome possui formato XXXXX_1_1.spec.sc4.16x5
   arq_nm = arq.split('/')[-1]
      
   glx = arq_nm.split("_")[0:1]
   x = int(arq_nm.split("_")[3])
   y = int(arq_nm.split("_")[2])
   dados_arquivo['x'] = x
   dados_arquivo['y'] = y
   
   dados_arquivo['Galaxy'] = glx
   
   # Array para armazenar as tabelas
   dadosSSP = []
   dados = []
   
   # Ler as informações do output
   with open(arq, 'r', encoding='utf-8') as f:
      salvaSSP = False
      salva = False
      for i, linha in enumerate(f):
      
         linha_limpa = linha.strip()
         
         # Pegando chi2
         if i == 49:
            dados_arquivo['chi2'] = float(linha_limpa.split('   ')[0])
         # Pegando adev
         elif i == 50:
            dados_arquivo['adev'] = float(linha_limpa.split('   ')[0])
         # Pegando Mini
         elif i == 54:
            dados_arquivo['Mini'] = float(linha_limpa.split('   ')[0])
         # Pegando Mcor
         elif i == 55:
            dados_arquivo['Mcor'] = float(linha_limpa.split('   ')[0])
         # Pegando v0
         elif i == 57:
            if linha_limpa.startswith('*******'):
               dados_arquivo['v0'] = np.nan
            else:
               dados_arquivo['v0'] = float(linha_limpa.split('   ')[0])
         # Pegando vd
         elif i == 58:
            dados_arquivo['vd'] = float(linha_limpa.split('   ')[0])
         # Pegando AV
         elif i == 59:
            dados_arquivo['AV'] = float(linha_limpa.split('   ')[0])
         
         else:
            # Pegando a tabela das SSPs
            if linha_limpa.startswith(headerSSP):
               salvaSSP = True
               i_tab = 1
            if salvaSSP:
               if not linha_limpa:
                  salvaSSP = False
               else:
                  if i_tab > 1:
                     dadosSSP.append(linha_limpa)
                     continue
                  i_tab += 1
            
            # Pegando a tabela dos espectros
            if linha_limpa.startswith(header):
               salva = True
               i_spec = 1
            if salva:
               if not linha_limpa:
                  salva = False
                  break
               else:
                  if i_spec > 2:
                     dados.append(linha_limpa)
                  i_spec += 1
   
   # Atribuindo arrays para as colunas de interesse
   linhas_divididas = [linha.split() for linha in dadosSSP]
   (j, x_j, Mini_j, Mcor_j, age_j, Z_j, LoMj, YAV, Mstars, component_j, aoFe,
   SSP_chi2r, SSP_adev, SSP_AV, SSP_x) = zip(*linhas_divididas)
      
   j = np.array(j, dtype=int)
   x_j = np.array(x_j, dtype=float)
   Mini_j = np.array(Mini_j, dtype=float)
   Mcor_j = np.array(Mcor_j, dtype=float)
   age_j = np.array(age_j, dtype=float)
   Z_j = np.array(Z_j, dtype=float)
   LoMj = np.array(LoMj, dtype=float)
   Mstars = np.array(Mstars, dtype=float)
   component_j = np.array(component_j)
      
   # Separar SSP, FC e HD (caso haja)
   onlyFC = np.char.find(component_j, "FC") != -1
   onlyFC |= np.char.find(component_j, "Power") != -1
      
   onlyHD = np.char.find(component_j, "BB") != -1
   onlyHD |= np.char.find(component_j, "HD") != -1
      
   onlySSP = (~onlyFC) & (~onlyHD)
      
   ageSSP = age_j[onlySSP]
   Z_SSP = Z_j[onlySSP]
      
   # Calculando as idades e metalicidades médias
   # renormalizando a soma x_j para fechar 100%
   x_j = (x_j * 100)/np.sum(x_j)
   x_jSSP = x_j[onlySSP]
   logt = np.sum( x_jSSP * np.log10(ageSSP) )/100
   dados_arquivo['<logt>x'] = logt
   dados_arquivo['<t>x'] = np.power(10, logt)
   Z = np.sum( x_jSSP * Z_SSP )/100
   dados_arquivo['<Z>x'] = Z
   dados_arquivo['<[M/H]>x'] = np.log10( Z/Zsun )
      
   dados_arquivo['FCx'] = np.sum( x_j[onlyFC] )
   dados_arquivo['HDx'] = np.sum( x_j[onlyHD] )
      
   m_j = (Mcor_j * 100)/np.sum(Mcor_j)
   m_jSSP = m_j[onlySSP]
   logt = np.sum( m_jSSP * np.log10(ageSSP) )/100
   dados_arquivo['<logt>m'] = logt
   dados_arquivo['<t>m'] = np.power(10, logt)
   Z = np.sum( m_jSSP * Z_SSP )/100
   dados_arquivo['<Z>m'] = Z
   dados_arquivo['<[M/H]>m'] = np.log10( Z/Zsun )
      
   # Calculando a fração de população para os bins de idades
   # considerando xi (age<1e8), xiy (age<7e8), xio (age<2e9), xo (age>2e9)
   only_xy = ageSSP <= 1e8
   only_xiy = (ageSSP > 1e8) & (ageSSP <= 7e8)
   only_xio = (ageSSP > 7e8) & (ageSSP <= 2e9)
   only_xo = ageSSP > 2e9
      
   dados_arquivo['xy'] = np.sum( x_jSSP[only_xy] )
   dados_arquivo['xiy'] = np.sum( x_jSSP[only_xiy] )
   dados_arquivo['xio'] = np.sum( x_jSSP[only_xio] )
   dados_arquivo['xo'] = np.sum( x_jSSP[only_xo] )
      
   dados_arquivo['my'] = np.sum( m_jSSP[only_xy] )
   dados_arquivo['miy'] = np.sum( m_jSSP[only_xiy] )
   dados_arquivo['mio'] = np.sum( m_jSSP[only_xio] )
   dados_arquivo['mo'] = np.sum( m_jSSP[only_xo] )
   
   # Salvando os dados dos espectros em um dicionário
   linhas_divididas = [linha.split() for linha in dados]
   (l_obs, f_obs, f_syn, wei) = zip(*linhas_divididas)
   
   l_obs = np.array(l_obs, dtype=float)
   f_obs = np.array(f_obs, dtype=float)
   f_syn = np.array(f_syn, dtype=float)
   wei = np.array(wei, dtype=float)
   
   if spec_spaxel:
      mask = (wei <= 0)
      f_masked = np.where(mask, f_obs, np.nan)
      plt.figure(figsize=(6,2), dpi=200)
      plt.plot(l_obs, f_obs, label="observado", lw=0.6, c='k', zorder=1)
      plt.plot(l_obs, f_syn, label="ajustado", lw=0.6, c='r', zorder=5)
      plt.plot(l_obs, f_masked, label=f'mascarado', lw=0.6, color='pink', zorder=2)
      plt.title(arq_nm)
      plt.legend(frameon=False, ncol=3)
      plt.xlabel("Comprimendo de onda (Å)")
      plt.ylabel("Fluxo")
      plt.tight_layout()
      plt.savefig(f'ajuste_{arq_nm}.png')
      plt.close()
   else:
      np.savetxt(arq+".spec", np.transpose([l_obs, f_obs, f_syn, wei]))
   
#   dados_l[(x, y)] = l_obs
   dados_fobs[(x, y)] = f_obs
   dados_fsyn[(x, y)] = f_syn
   dados_wei[(x, y)] = wei
   
   if l_wave is None:
      l_wave = l_obs
   
   # Transformando as informações da população estelar em um data frame
   df = pd.concat([df, pd.DataFrame([dados_arquivo])], ignore_index=True)
#print(df)


# Salvando em um CSV os dados da população estelar
df.to_csv(f"../data/{galaxia}_STARLIGHT_output_params.csv", float_format="%.5E", index=False)



# pegando o spaxel que corresponde ao pico do continuo
cnt         = np.nanmedian(data[1000:1500,:,:], axis=0)
arr2D   = cnt.copy()
arr2D   = arr2D / np.nanmean(arr2D) # para trabalhar com numeros menores, o ideal seria ser magnitudes
resultado  = np.where(arr2D == np.amax(arr2D))
print(resultado)
ycen = resultado[0][0]
xcen = resultado[1][0]

# calcular a distancia a partir do nucleo
rad = data[0,:,:] * 0.0
nx  = len(rad[0,:])
ny  = len(rad[:,0])
for i in range(0,nx):
    for j in range(0,ny):
        rad[j,i] = ((j - ycen)**2 + (i - xcen)**2)**0.5 * sl

# definicoes das figuras:
font = {'size': 13, 'family':'serif'} # 11
mpl.rc('font', **font)      # converte todas as fontes do grafico no estilo 'font'
font = font['size']
cmap1 = copy.copy(plt.cm.viridis) #plt.cm.gist_heat
cmap1.set_bad(color='gainsboro')
ext     = [-xcen*sl*(-1),-(len(cnt[0,:])-xcen)*sl,ycen*sl*(-1),(len(cnt[:,0])-ycen)*sl]

mask_radius = 2.5  # exemplo em pixels

# Plota figura com as frações das populações estelares
colunas = ['xy', 'xiy', 'xio', 'xo',
           'my', 'miy', 'mio', 'mo']


fig, axes = plt.subplots(2, 4, figsize=(12,8), dpi=200)

for i, ax in enumerate(axes.flat):
   col = colunas[i]

   grid = np.full((ny, nx), np.nan)
   grid[df['y'].values.astype(int), df['x'].values.astype(int)] = df[col].values
   masked_grid = np.where(rad <= mask_radius, grid, np.nan)
   im = ax.imshow(masked_grid,cmap = cmap1,  
                origin = 'lower', extent = ext)

   
   # Título do subplot (em matplotlib usa-se set_title)
   ax.set_title(col, fontsize=10)
   ax.set_aspect('equal') # Mantém a proporção espacial correta da galáxia
   ax.set_xlabel(r'$\Delta$X (arcsec)')
   ax.set_ylabel(r'$\Delta$Y (arcsec)')


#plt.subplots_adjust(left=0.05, bottom=0.05, right=0.9, top=0.95, wspace=0.4, hspace=0.4)

# Adiciona UMA barra de cores global à direita compartilhada por todos os subplots
fig.canvas.draw()
pos_topo = axes[0, 3].get_position()
pos_base = axes[1, 3].get_position()
x_cbar = pos_topo.x1 + 0.02                    # Um pouco à direita da última coluna
y_cbar = pos_base.y0                           # Base = base do plot inferior
largura_cbar = 0.015                           # Espessura da barra
altura_cbar = pos_topo.y1 - pos_base.y0        # Altura exata da linha 1 até a linha 2
cax = fig.add_axes([x_cbar, y_cbar, largura_cbar, altura_cbar])
cbar = fig.colorbar(im, cax=cax)
cbar.set_label('Fração (%)', fontsize=11)

plt.subplots_adjust(top=0.932,
bottom=0.152,
left=0.051,
right=0.912,
hspace=0.525,
wspace=0.179)
plt.savefig(f'../data/{galaxia}_stellar_population.png')
plt.close()

# Plota figura com as informações da galáxia
colunas = ['<t>x', '<[M/H]>x', 'FCx', 'HDx',
           '<t>m', '<[M/H]>m', 'AV', 'adev']
fig, axes = plt.subplots(2, 4, figsize=(12,6), dpi=200)

for i, ax in enumerate(axes.flat):
   col = colunas[i]

   grid = np.full((ny, nx), np.nan)
   grid[df['y'].values.astype(int), df['x'].values.astype(int)] = df[col].values
   masked_grid = np.where(rad <= mask_radius, grid, np.nan)
   im = ax.imshow(masked_grid,cmap = cmap1,  
                   origin = 'lower', extent = ext)
   
   ax.set_aspect('equal') # Mantém a proporção espacial correta da galáxia
   ax.set_xlabel(r'$\Delta$X (arcsec)')
   ax.set_ylabel(r'$\Delta$Y (arcsec)')
   
   # Cria um divisor de eixos para acoplar a colorbar ao lado do subplot
   divider = make_axes_locatable(ax)
   cax = divider.append_axes("right", size="5%", pad=0.08) # 'size' é a largura
                                                           # 'pad', distância
   cbar = fig.colorbar(im, cax=cax)
   cbar.ax.tick_params(labelsize=8)  # Ajusta o tamanho dos números na barra
   cbar.set_label(col, fontsize=11)

plt.tight_layout()
plt.savefig(f'../data/{galaxia}_properties.png')
plt.close()