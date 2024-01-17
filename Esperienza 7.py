#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt


# Sostituisci con il nome del tuo file Excel
nome_colonna_1 = 'V_d'
nome_colonna_2 = 'ΔV_d'
nome_colonna_3 = 'I_d (mA)'
nome_colonna_4 = 'ΔI_d'
nome_colonna_5 = 'V_dn'
nome_colonna_6 = 'ΔV_dn'
nome_colonna_7 = 'I_d (microA)n'
nome_colonna_8 = 'ΔI_dn'
# Carica il file Excel in un DataFrame pandas
df = pd.read_excel(io=os.path.abspath('Esperienza 7.xlsx'), sheet_name='Circuito1 errori - new')
# Converti le colonne specificate in float, gestendo i valori non convertibili come NaN
df[nome_colonna_1] = pd.to_numeric(df[nome_colonna_1], errors='coerce')
df[nome_colonna_2] = pd.to_numeric(df[nome_colonna_2], errors='coerce')
df[nome_colonna_3] = pd.to_numeric(df[nome_colonna_3], errors='coerce')
df[nome_colonna_4] = pd.to_numeric(df[nome_colonna_4], errors='coerce')
df[nome_colonna_5] = pd.to_numeric(df[nome_colonna_5], errors='coerce')
df[nome_colonna_6] = pd.to_numeric(df[nome_colonna_6], errors='coerce')
df[nome_colonna_7] = pd.to_numeric(df[nome_colonna_7], errors='coerce')
df[nome_colonna_8] = pd.to_numeric(df[nome_colonna_8], errors='coerce')

# Estrai i valori dalle due colonne specificate
V_d = df[nome_colonna_1].values
ΔV_d = df[nome_colonna_2].values
I_d = df[nome_colonna_3].values*10**3
ΔI_d = df[nome_colonna_4].values
V_dn = df[nome_colonna_5].values
ΔV_dn = df[nome_colonna_6].values
I_dn = df[nome_colonna_7].values*10**3
ΔI_dn = df[nome_colonna_8].values



plt.errorbar(V_d, I_d, xerr=ΔV_d, yerr=ΔI_d,fmt='s', markersize=2., capsize=2.5, label='Dati polarizzazione diretta')
plt.errorbar(V_dn, I_dn, xerr=ΔV_dn, yerr=ΔI_dn,fmt='s', markersize=2., capsize=2.5, label='Dati polarizzazione inversa')
plt.xlabel('$V_d$ [V]')
plt.ylabel('$I_d/mA$ [mA]')
plt.title('Curva caratteristica diodo a giunzione')
plt.grid(True)  # Abilita la griglia
plt.axhline(0, color='black')
plt.legend()
plt.show()


# In[3]:


# Secondo subplot

# Imposta la frazione come etichetta
label_text = r'$\frac{I_d}{mA}$'

filtrati = df[df['V_d'] >= 0.6]
def fit_lineare_pesato(x,y,w) :
    """
    questa funzione calcola i parametri della retta di best fit y = a + b*x 
    usando le formule dei minimi quadrati pesati.
    N.B. w è il vettore dei pesi, dove il peso w_i è l'inverso del quadrato dell'incertezza 
    sulla coordinata y_i dell'i-esimo punto.    
    """
    S_0 = np.sum(w)
    S_x = np.sum(x*w)     
    S_xx = np.sum((x*x)*w)     
    S_y = np.sum(y*w)
    S_xy = np.sum(x*y*w)
    D = S_0 * S_xx - S_x**2
    a = (S_xx * S_y - S_x * S_xy) / D
    b = (S_0 * S_xy - S_x * S_y) / D
    var_a = S_xx / D
    var_b = S_0 / D
    cov_ab = -S_x / D
    sigma_a = np.sqrt(var_a)
    sigma_b = np.sqrt(var_b)
    #Compute chi^2 = \sum w_i (y_i - (a + b * x_i))^2
    chi2 = np.sum (w * (y-(a+b*x))**2)
    print(f"a = {a}+/-{sigma_a}")
    print(f"b = {b}+/-{sigma_b}")
    print(f"cov(a,b) = {cov_ab}")
    print(f"chi/ndof= {chi2}/{len(x)-2} = {chi2/(len(x)-2)}")

    return a,b,sigma_a,sigma_b,cov_ab,chi2
a, b, sigma_a, sigma_b, cov_ab, chi2 = fit_lineare_pesato(filtrati['V_d'], np.log(filtrati['I_d (mA)']), 1/(filtrati['ΔI_d']) )
text_box = {
    'boxstyle': 'round',
    'facecolor': 'lightblue',
    'alpha': 1
}

# Testo con andate a capo
testo = " log($I_d$) = b$\cdot V_d$ + a \n a = -20,7 $\pm$ 0.7 \n b = 19.9 $\pm$ 0.6 $V^{-1}$ "
# Posiziona il testo in alto a destra
plt.text(0.71, 0.17, testo, ha='left', va='top', transform=plt.gca().transAxes, bbox=text_box)
# Filtra i dati dove I_d è maggiore di 0.6 e prendi i valori corrispondenti di V_d
plt.errorbar(filtrati['V_d'], np.log(filtrati['I_d (mA)']), xerr=filtrati['ΔV_d'], yerr=filtrati['ΔI_d'],
             fmt='s', markersize=2., capsize=2.5, label='Dati polarizzazione diretta > 0.6')
plt.plot(filtrati['V_d'],filtrati['V_d']*b+a, color='red', ls='-.', lw=1, label='bestfit')
plt.xlabel('$V_d$ [V]')
plt.ylabel('log($I_d/mA$)')
plt.title('Curva caratteristica polarizzazione diretta(logaritmo)')
plt.grid(True)  # Abilita la griglia
plt.legend()
plt.show()


# In[4]:


plt.errorbar(V_d, I_d, xerr=ΔV_d, yerr=ΔI_d,fmt='s', markersize=2., capsize=2.5, label='Dati polarizzazione diretta')
plt.xlabel('$V_d$ [V]')
plt.ylabel('$I_d/mA$ [mA]')
plt.title('Curva caratteristica polarizzazione diretta')
plt.grid(True)  # Abilita la griglia
plt.legend()
plt.show()


# In[5]:


V_dN = V_dn[~np.isnan(V_dn)]
ΔV_dN = ΔV_dn[~np.isnan(ΔV_dn)]
I_dN = V_dn[~np.isnan(I_dn)]
ΔI_dN=ΔI_dn[~np.isnan(ΔI_dn)]  
def fit_lineare_pesato(x,y,w) :
    """
    questa funzione calcola i parametri della retta di best fit y = a + b*x 
    usando le formule dei minimi quadrati pesati.
    N.B. w è il vettore dei pesi, dove il peso w_i è l'inverso del quadrato dell'incertezza 
    sulla coordinata y_i dell'i-esimo punto.    
    """
    S_0 = np.sum(w)
    S_x = np.sum(x*w)     
    S_xx = np.sum((x*x)*w)     
    S_y = np.sum(y*w)
    S_xy = np.sum(x*y*w)
    D = S_0 * S_xx - S_x**2
    a = (S_xx * S_y - S_x * S_xy) / D
    b = (S_0 * S_xy - S_x * S_y) / D
    var_a = S_xx / D
    var_b = S_0 / D
    cov_ab = -S_x / D
    sigma_a = np.sqrt(var_a)
    sigma_b = np.sqrt(var_b)
    #Compute chi^2 = \sum w_i (y_i - (a + b * x_i))^2
    chi2 = np.sum (w * (y-(a+b*x))**2)
    print(f"a = {a}+/-{sigma_a}")
    print(f"b = {b}+/-{sigma_b}")
    print(f"cov(a,b) = {cov_ab}")
    print(f"chi/ndof= {chi2}/{len(x)-2} = {chi2/(len(x)-2)}")

    return a,b,sigma_a,sigma_b,cov_ab,chi2
a, b, sigma_a, sigma_b, cov_ab, chi2 = fit_lineare_pesato(V_dN, I_dN, ΔI_dN )
text_box = {
    'boxstyle': 'round',
    'facecolor': 'lightblue',
    'alpha': 1
}

# Testo con andate a capo
testo = " $I_d$ = b$\cdot V_d$ + a \n a = 0.0 $\pm$ 0.3 mA\n b = 1.0 $\pm$ 0.3 $mA/V$ "
# Posiziona il testo in alto a destra
plt.text(0.7, 0.16, testo, ha='left', va='top', transform=plt.gca().transAxes, bbox=text_box)



plt.errorbar(V_dN, I_dN, xerr=ΔV_dN, yerr=ΔI_dN,fmt='none', markersize=2.,ecolor='orange', capsize=2.5, label='Dati polarizzazione inversa')
plt.plot(V_dN, (V_dN)*b+a, color='purple', ls='-.', lw=1, label='bestfit')
plt.xlabel('$V_d$ [V]')
plt.ylabel('$I_d/mA$ [mA]')
plt.title('Curva caratteristica polarizzazione inversa')
plt.grid(True)  # Abilita la griglia
plt.legend()
plt.show()





# In[6]:


nome_colonna_1 = 'V_0max'
nome_colonna_2 = 'σ_{V_0^{max}}'
nome_colonna_3 = 'V_Rmax'
nome_colonna_4 = 'σ_{V_R^{max}}'
nome_colonna_5 = 'V0'
nome_colonna_6 = 'sigmaV0'
nome_colonna_7 = 'VR'
nome_colonna_8 = 'sigmaVR'
df = pd.read_excel(io=os.path.abspath('Esperienza 7.xlsx'), sheet_name='Circuito2 Clipper')
# Converti le colonne specificate in float, gestendo i valori non convertibili come NaN
df[nome_colonna_1] = pd.to_numeric(df[nome_colonna_1], errors='coerce')
df[nome_colonna_2] = pd.to_numeric(df[nome_colonna_2], errors='coerce')
df[nome_colonna_3] = pd.to_numeric(df[nome_colonna_3], errors='coerce')
df[nome_colonna_4] = pd.to_numeric(df[nome_colonna_4], errors='coerce')
df[nome_colonna_5] = pd.to_numeric(df[nome_colonna_5], errors='coerce')
df[nome_colonna_6] = pd.to_numeric(df[nome_colonna_6], errors='coerce')
df[nome_colonna_7] = pd.to_numeric(df[nome_colonna_7], errors='coerce')
df[nome_colonna_8] = pd.to_numeric(df[nome_colonna_8], errors='coerce')
V0 = df[nome_colonna_1].values
ΔV0 = df[nome_colonna_2].values
VR = df[nome_colonna_3].values
ΔVR = df[nome_colonna_4].values
V0F = df[nome_colonna_5].values
ΔV0F = df[nome_colonna_6].values
VRF = df[nome_colonna_7].values
ΔVRF = df[nome_colonna_8].values


plt.errorbar(V0, VR, xerr=ΔV0, yerr=ΔVR,fmt='s', markersize=2., capsize=2.5, label='Con diodo')
plt.errorbar(V0F, VRF, xerr=ΔV0F, yerr=ΔVRF,fmt='s', markersize=2., capsize=2.5, label='Senza diodo')
plt.xlabel('$V_0$ [V]')
plt.ylabel('$V_R$ [V]')
plt.title('Limitatore di tensione con montaggio clipper')
plt.grid(True)  # Abilita la griglia
plt.legend()
plt.show()
print(V0)
print(ΔV0)
print(VR)
print(ΔVR)
print(V0F)
print(ΔV0F)
print(VRF)
print(ΔVRF)

