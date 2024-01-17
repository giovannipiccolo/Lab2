#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy as sp
import math
from matplotlib import rc

from uncertainties import ufloat
from scipy.optimize import curve_fit



# In[2]:


import numpy as np
import matplotlib.pyplot as plt

V0 = 9.81
Vc = np.array([2.39, 4.14, 5.47, 6.47, 7.25, 7.83, 8.29, 8.64, 8.9, 9.11, 9.27, 9.39, 9.49, 9.56, 9.62, 9.66, 9.7, 9.73, 9.75, 9.77])
T = np.array([30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360, 390, 420, 450, 480, 510, 540, 570, 600])
R = 1 - (Vc / V0)
err_R = 0.002 * np.ones(len(Vc))

plt.figure(figsize=(12, 5))

# Primo subplot
plt.subplot(1, 2, 1)
plt.errorbar(T, 1 - (Vc / V0), yerr=err_R, fmt='none', ecolor='b', capsize=5, label='Dati RC')
plt.xlabel('t (s)')
plt.ylabel('$R_C (t)$')
plt.legend()
plt.title('Carica lenta del condensatore')
plt.grid()

# Secondo subplot
plt.subplot(1, 2, 2)
logR = np.log(R)
err_log = err_R/R
text_box = {
    'boxstyle': 'round',
    'facecolor': 'lightblue',
    'alpha': 1
}

# Testo con andate a capo
testo = " a =0.009 $\pm$ 0.004 \n b = (8.89 $\pm$ 0.02)e-03 Hz \n $log(R_C (t))= a + bt$ "
# Posiziona il testo in alto a destra
plt.text(0.05, 0.15, testo, ha='left', va='top', transform=plt.gca().transAxes, bbox=text_box)

def fit_lineare_pesato(x, y, w):
    S_0 = np.sum(w)
    S_x = np.sum(x * w)
    S_xx = np.sum(x * x * w)
    S_y = np.sum(y * w)
    S_xy = np.sum(x * y * w)
    D = S_0 * S_xx - S_x**2
    a = (S_xx * S_y - S_x * S_xy) / D
    b = (S_0 * S_xy - S_x * S_y) / D
    var_a = S_xx / D
    var_b = S_0 / D
    cov_ab = -S_x / D
    sigma_a = np.sqrt(var_a)
    sigma_b = np.sqrt(var_b)
    chi2 = np.sum(w * (y - (a + b * x))**2)
    
    print(f"a = {a} ± {sigma_a}")
    print(f"b = {b} ± {sigma_b}")
    print(f"cov(a,b) = {cov_ab}")
    print(f"chi/ndof = {chi2}/{len(x)-2} = {chi2/(len(x)-2)}")

    return a, b, sigma_a, sigma_b, cov_ab, chi2

a, b, sigma_a, sigma_b, cov_ab, chi2 = fit_lineare_pesato(T, logR, 1 / err_R**2)

plt.errorbar(T, logR, yerr=err_R, fmt='none', ecolor='orange', capsize=5, label='Dati RC')
plt.plot(T, a + b * T, label='Carica bestfit')
plt.xlabel('t (s)')
plt.ylabel('log($R_C (t))$')
plt.legend()
plt.title('Carica lenta del condensatore (logaritmo)')
plt.grid()

plt.tight_layout()
plt.show()


# In[3]:


import numpy as np
import matplotlib.pyplot as plt

V0 = 9.81
Vc = np.array([7.52, 5.78, 4.38, 3.402, 2.626, 2.019, 1.563, 1.208, 0.931, 0.728, 0.563, 0.44, 0.344, 0.271, 0.213, 0.169, 0.134, 0.107, 0.086, 0.069])
T = np.array([30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360, 390, 420, 450, 480, 510, 540, 570, 600])
R1 = (Vc / V0)
err_R = 0.002 * np.ones(len(Vc))

plt.figure(figsize=(12, 5))

# Primo subplot
plt.subplot(1, 2, 1)
plt.errorbar(T, R1, yerr=err_R, fmt='.', ecolor='b', capsize=5, label='Dati CR')
plt.xlabel('t (s)')
plt.ylabel('$R_s (t)$')
plt.legend()
plt.title('Scarica lenta del condensatore')
plt.grid()

# Secondo subplot
plt.subplot(1, 2, 2)
logR1 = np.log(R1)
err_log = err_R/R1
text_box = {
    'boxstyle': 'round',
    'facecolor': 'lightblue',
    'alpha': 1
}

# Testo con andate a capo
testo = " a =0.083 $\pm$ 0.005 \n b = (8.27 $\pm$ 0.02)e-03 Hz \n $log(R_C (t))= a + bt$ "
# Posiziona il testo in alto a destra
plt.text(0.05, 0.15, testo, ha='left', va='top', transform=plt.gca().transAxes, bbox=text_box)


def fit_lineare_pesato(x, y, w):
    S_0 = np.sum(w)
    S_x = np.sum(x * w)
    S_xx = np.sum(x * x * w)
    S_y = np.sum(y * w)
    S_xy = np.sum(x * y * w)
    D = S_0 * S_xx - S_x**2
    a = (S_xx * S_y - S_x * S_xy) / D
    b = (S_0 * S_xy - S_x * S_y) / D
    var_a = S_xx / D
    var_b = S_0 / D
    cov_ab = -S_x / D
    sigma_a = np.sqrt(var_a)
    sigma_b = np.sqrt(var_b)
    chi2 = np.sum(w * (y - (a + b * x))**2)
    
    print(f"a = {a} ± {sigma_a}")
    print(f"b = {b} ± {sigma_b}")
    print(f"cov(a,b) = {cov_ab}")
    print(f"chi/ndof = {chi2}/{len(x)-2} = {chi2/(len(x)-2)}")

    return a, b, sigma_a, sigma_b, cov_ab, chi2

a, b, sigma_a, sigma_b, cov_ab, chi2 = fit_lineare_pesato(T, logR1, 1 / err_R**2)

plt.errorbar(T, logR1, yerr=err_R, fmt='.', ecolor='blue', capsize=5, label='Dati CR')
plt.plot(T, a + b * T, label='Scarica bestfit')
plt.xlabel('t (s)')
plt.ylabel('log($R_s (t)$)')
plt.legend()
plt.title('Scarica lenta del condensatore(logaritmo)')
plt.grid()

plt.tight_layout()
plt.show()


# In[4]:


import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

# Sostituisci con il nome del tuo file Excel
nome_colonna_1 = 'Vc/V0'
nome_colonna_2 = 't/T'

# Carica il file Excel in un DataFrame pandas
df = pd.read_excel(io=os.path.abspath('Misure oscilloscopio VEREABEM.xlsx'), sheet_name='500')
# Converti le colonne specificate in float, gestendo i valori non convertibili come NaN
df[nome_colonna_1] = pd.to_numeric(df[nome_colonna_1], errors='coerce')
df[nome_colonna_2] = pd.to_numeric(df[nome_colonna_2], errors='coerce')

# Estrai i valori dalle due colonne specificate
valori_colonna_1 = df[nome_colonna_1].values
valori_colonna_2 = df[nome_colonna_2].values

# Filtra i valori della colonna_x compresi tra -1000 e 1000
indici_filtrati = np.where((valori_colonna_2 > 0) & (valori_colonna_2 <= 0.0009))
x_filtrati = valori_colonna_2[indici_filtrati]
y_corrispondenti = valori_colonna_1[indici_filtrati]

# Moltiplica i valori di x per 1000
x_filtrati_ms = x_filtrati * 1000  # Converti in millisecondi

# Creazione del subplot con asse x logaritmico
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Primo subplot
ax1.scatter(x_filtrati_ms, y_corrispondenti)
ax1.set_xlabel('t [ms]')
ax1.set_ylabel('R$_C$(t)')
ax1.set_title('Circuito RC in regime impulsivo con  $\\tau <<$ T')
ax1.grid(True)  # Abilita la griglia

# Secondo subplot
ax2.scatter(x_filtrati_ms, y_corrispondenti)
ax2.set_xlabel('log(t)')
ax2.set_ylabel('R$_C$(t)')
ax2.set_title('Circuito RC in regime impulsivo con  $\\tau <<$ T (logaritmo)')
ax2.set_xscale('log')
ax2.grid(True)  # Abilita la griglia


plt.tight_layout()
plt.show()


# In[6]:


import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

# Sostituisci con il nome del tuo file Excel
nome_colonna_1 = 'Vc/V0'
nome_colonna_2 = 't/T'
nome_colonna_3= 'VF1'
nome_colonna_4= 'VF2'
# Carica il file Excel in un DataFrame pandas
df = pd.read_excel(io=os.path.abspath('Misure oscilloscopio VEREABEM.xlsx'), sheet_name='9750')
# Converti le colonne specificate in float, gestendo i valori non convertibili come NaN
df[nome_colonna_1] = pd.to_numeric(df[nome_colonna_1], errors='coerce')
df[nome_colonna_2] = pd.to_numeric(df[nome_colonna_2], errors='coerce')
df[nome_colonna_3] = pd.to_numeric(df[nome_colonna_3], errors='coerce')
df[nome_colonna_4] = pd.to_numeric(df[nome_colonna_4], errors='coerce')
# Estrai i valori dalle due colonne specificate
valori_colonna_1 = df[nome_colonna_1].values
valori_colonna_2 = df[nome_colonna_2].values
V0 = df[nome_colonna_3].values
Vt = df[nome_colonna_4].values
err_y=np.sqrt(((1/(-0.728-V0))**2)+(((Vt+0.728)/(V0+0.728)**2)**2) + ((Vt-V0)/((V0+0.728)**2)**2))*(np.sqrt((((0.01*3*Vt)/(2.58))**2)+((0.1*0.1)/np.sqrt(3))**2) + ((0.001)/(np.sqrt(3)))**2)
#io = err_y/valori_colonna_2
# Filtra i valori della colonna_x compresi tra -1000 e 1000
indici_filtrati = np.where((valori_colonna_2 > 0.005) & (valori_colonna_2 <= 0.049))
x_filtrati = valori_colonna_2[indici_filtrati]
y_corrispondenti = valori_colonna_1[indici_filtrati]

# Moltiplica i valori di x per 1000
x_filtrati_ms = x_filtrati * 1000  # Converti in millisecondi

# Creazione del subplot con asse x logaritmico
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Primo subplot
ax1.scatter(x_filtrati_ms, y_corrispondenti)
ax1.plot(x_filtrati_ms, x_filtrati_ms*b +a )
ax1.set_xlabel('t [ms]')
ax1.set_ylabel('R$_C$(t)')
ax1.set_title('Circuito RC in regime impulsivo con  $\\tau \\approx$ T')
ax1.grid(True)  # Abilita la griglia

# Secondo subplot
ax2.scatter(x_filtrati_ms, y_corrispondenti)
ax2.set_xlabel('log(t)')
ax2.set_ylabel('R$_C$(t)')
ax2.set_title('Circuito RC in regime impulsivo con  $\\tau \\approx$ T (logaritmo)')
ax2.set_yscale('log')
ax2.grid(True)  # Abilita la griglia


plt.tight_layout()
plt.show()

err_y=np.sqrt(((1/(-0.728-V0))**2)+(((Vt+0.728)/(V0+0.728)**2)**2) + ((Vt-V0)/((V0+0.728)**2)**2))*(np.sqrt((((0.01*3*Vt)/(2.58))**2)+((0.1*0.1)/np.sqrt(3))**2) + ((0.001)/(np.sqrt(3)))**2)

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
a, b, sigma_a, sigma_b, cov_ab, chi2 = fit_lineare_pesato(x_filtrati_ms,y_corrispondenti , np.ones(len(y_corrispondenti))*1000)

ax1.plot(x_filtrati_ms, x_filtrati_ms*b +a )

print(err_y)


# In[8]:


print(V0)
print(Vc)


# In[ ]:





# In[1]:


import matplotlib.pyplot as plt

# Creazione di un set di subplot
fig, (ax1, ax2) = plt.subplots(1, 2)

# Esempio di tracciamento di qualche cosa in ciascun subplot
ax1.plot([1, 2, 3], [4, 5, 6])
ax2.plot([3, 2, 1], [4, 5, 6])

# Aggiunta di testo al subplot ax2
ax2.text(0.5, 0.5, 'Testo nel Subplot 2', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes)

plt.show()


# In[ ]:




