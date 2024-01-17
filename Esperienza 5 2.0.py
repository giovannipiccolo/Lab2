#!/usr/bin/env python
# coding: utf-8

# In[48]:


import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt


# In[28]:


import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import random
import random

# Sostituisci con il nome del tuo file Excel
nome_colonna_1 = 'Vc/V0'
nome_colonna_2 = 't/T'
nome_colonna_3= 'VF1'
nome_colonna_4= 'VF2'
# Carica il file Excel in un DataFrame pandas
df = pd.read_excel(io=os.path.abspath('Misure oscilloscopio VEREABEM.xlsx'), sheet_name='500')
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

# Filtra i valori della colonna_x compresi tra -1000 e 1000
indici_filtrati = np.where((valori_colonna_2 > 0) & (valori_colonna_2 <= 0.0009))
x_filtrati = valori_colonna_2[indici_filtrati]
y_corrispondenti = valori_colonna_1[indici_filtrati]

# Moltiplica i valori di x per 1000
x_filtrati_ms = x_filtrati * 1000  # Converti in millisecondi
err_y=np.sqrt(((1/(-0.728-V0))**2)+(((Vt+0.728)/(V0+0.728)**2)**2) + ((Vt-V0)/((V0+0.728)**2)**2))*(np.sqrt((((0.01*3*Vt)/(2.58))**2)+((0.1*0.1)/np.sqrt(3))**2) + ((0.001)/(np.sqrt(3)))**2)
err_y_filtered = err_y[indici_filtrati]
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
a, b, sigma_a, sigma_b, cov_ab, chi2 = fit_lineare_pesato(x_filtrati_ms,np.log(y_corrispondenti) , err_y_filtered*10)
# Creazione del subplot con asse x logaritmico
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
text_box = {
    'boxstyle': 'round',
    'facecolor': 'lightblue',
    'alpha': 1
}

# Testo con andate a capo
testo = " $log(R$_C$(t))$= b$\cdot t$ + a \n a = -0.817 $\pm$ 0.009 \n b = -11.43 $\pm$ 0.06 $ms^{-1}$ "
# Posiziona il testo in alto a destra
ax2.text(0.01, 0.14, testo, ha='left', va='top', transform=ax2.transAxes, bbox=text_box)



# Primo subplot
ax1.errorbar(x_filtrati_ms, y_corrispondenti, yerr=err_y_filtered, fmt='s', markersize=2., capsize=2.5,label='Dati RC')
ax1.set_xlabel('t [ms]')
ax1.set_ylabel('R$_C$(t)')
ax1.set_title('Circuito RC in regime impulsivo con  $\\tau <<$ T')
ax1.legend()
ax1.grid(True)  # Abilita la griglia

# Secondo subplot
ax2.errorbar(x_filtrati_ms, np.log(y_corrispondenti),xerr=np.ones(len(x_filtrati_ms))*0.000006, yerr=err_y_filtered*10, fmt='s', markersize=2., capsize=2.5, label='Dati RC')
ax2.plot(x_filtrati_ms, x_filtrati_ms*b+ a,color='red', ls='-.', lw=1, label='Bestfit')
ax2.set_xlabel('t [ms]')
ax2.set_ylabel('log(R$_C$(t))')
ax2.set_title('Circuito RC in regime impulsivo con  $\\tau <<$ T (logaritmo)')
ax2.legend()
ax2.grid(True)  # Abilita la griglia


plt.tight_layout()
plt.show()


# In[29]:


import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

# Sostituisci con il nome del tuo file Excel
nome_colonna_1 = 'Vc/V0'
nome_colonna_2 = 't/T'

# Carica il file Excel in un DataFrame pandas
df = pd.read_excel(io=os.path.abspath('Misure oscilloscopio VEREABEM.xlsx'), sheet_name='9750')
# Converti le colonne specificate in float, gestendo i valori non convertibili come NaN
df[nome_colonna_1] = pd.to_numeric(df[nome_colonna_1], errors='coerce')
df[nome_colonna_2] = pd.to_numeric(df[nome_colonna_2], errors='coerce')

# Estrai i valori dalle due colonne specificate
valori_colonna_1 = df[nome_colonna_1].values
valori_colonna_2 = df[nome_colonna_2].values

# Filtra i valori della colonna_x compresi tra -1000 e 1000
indici_filtrati = np.where((valori_colonna_2 > 0.001) & (valori_colonna_2 <= 0.0175))
x_filtrati = valori_colonna_2[indici_filtrati]
y_corrispondenti = valori_colonna_1[indici_filtrati]
err_y_filtered = err_y[indici_filtrati]
# Moltiplica i valori di x per 1000
x_filtrati_ms = x_filtrati * 1000  # Converti in millisecondi
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
a, b, sigma_a, sigma_b, cov_ab, chi2 = fit_lineare_pesato(x_filtrati_ms,np.log(y_corrispondenti) , err_y_filtered*10)
# Creazione del subplot con asse x logaritmico
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
text_box = {
    'boxstyle': 'round',
    'facecolor': 'lightblue',
    'alpha': 1
}

# Testo con andate a capo
testo = " $log(R$_C$(t))$= b$\cdot t$ + a \n a = -4.42 $\pm$ 0.03 \n b = -0.066 $\pm$ 0.041 $ms^{-1}$ "
# Posiziona il testo in alto a destra
ax2.text(0.01, 0.14, testo, ha='left', va='top', transform=ax2.transAxes, bbox=text_box)

# Primo subplot
ax1.errorbar(x_filtrati_ms, y_corrispondenti,xerr=np.ones(len(x_filtrati_ms))*0.000006, yerr=err_y_filtered/100, fmt='s', markersize=2., capsize=2.5,label='Dati RC')
ax1.set_xlabel('t [ms]')
ax1.set_ylabel('R$_C$(t)')
ax1.set_title('Circuito RC in regime impulsivo con  $\\tau \\approx$ T')
ax1.legend()
ax1.grid(True)  # Abilita la griglia

# Secondo subplot
ax2.errorbar(x_filtrati_ms, np.log(y_corrispondenti), yerr=err_y_filtered, fmt='s', markersize=2., capsize=2.5,label='Dati RC')
ax2.plot(x_filtrati_ms, x_filtrati_ms*b+ a,color='red', ls='-.', lw=1, label='Bestfit')
ax2.set_xlabel('t [ms]')
ax2.set_ylabel('log(R$_C$(t))')
ax2.set_title('Circuito RC in regime impulsivo con  $\\tau \\approx$ T (logaritmo)')
ax2.legend()
ax2.grid(True)  # Abilita la griglia


plt.tight_layout()
plt.show()


# In[30]:


import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

# Sostituisci con il nome del tuo file Excel
nome_colonna_1 = 'Vc/V0'
nome_colonna_2 = 't/T'

# Carica il file Excel in un DataFrame pandas
df = pd.read_excel(io=os.path.abspath('Misure oscilloscopio VEREABEM.xlsx'), sheet_name='40000')
# Converti le colonne specificate in float, gestendo i valori non convertibili come NaN
df[nome_colonna_1] = pd.to_numeric(df[nome_colonna_1], errors='coerce')
df[nome_colonna_2] = pd.to_numeric(df[nome_colonna_2], errors='coerce')

# Estrai i valori dalle due colonne specificate
valori_colonna_1 = df[nome_colonna_1].values
valori_colonna_2 = df[nome_colonna_2].values

# Filtra i valori della colonna_x compresi tra -1000 e 1000
indici_filtrati = np.where((valori_colonna_2 > 0.000037) & (valori_colonna_2 <= 0.00005))
x_filtrati = valori_colonna_2[indici_filtrati]
y_corrispondenti = valori_colonna_1[indici_filtrati]
err_y_filtered = err_y[indici_filtrati]
# Moltiplica i valori di x per 1000
x_filtrati_ms = x_filtrati * 1000  # Converti in millisecondi
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
a, b, sigma_a, sigma_b, cov_ab, chi2 = fit_lineare_pesato(x_filtrati_ms,np.log(y_corrispondenti) , err_y_filtered*10)
# Creazione del subplot con asse x logaritmico
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
text_box = {
    'boxstyle': 'round',
    'facecolor': 'lightblue',
    'alpha': 1
}

# Testo con andate a capo
testo = " $log(R$_C$(t))$= b$\cdot t$ + a \n a = 0.414 $\pm$ 0.003 \n b = -9.47 $\pm$ 0.02 $ms^{-1}$ "
# Posiziona il testo in alto a destra
ax2.text(0.01, 0.14, testo, ha='left', va='top', transform=ax2.transAxes, bbox=text_box)

# Primo subplot
ax1.errorbar(x_filtrati_ms, y_corrispondenti,xerr=np.ones(len(x_filtrati_ms))*0.000006, yerr=err_y_filtered/5, fmt='s', markersize=2., capsize=2.5,label='Dati RC')
ax1.set_xlabel('t [ms]')
ax1.set_ylabel('R$_C$(t)')
ax1.set_title('Circuito RC in regime impulsivo con  $\\tau >>$ T')
ax1.legend()
ax1.grid(True)  # Abilita la griglia
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
a, b, sigma_a, sigma_b, cov_ab, chi2 = fit_lineare_pesato(x_filtrati_ms,np.log(y_corrispondenti) , err_y_filtered*10)


# Secondo subplot
ax2.errorbar(x_filtrati_ms, np.log(y_corrispondenti), yerr=err_y_filtered/10, fmt='s', markersize=2., capsize=2.5,label='Dati RC')
ax2.plot(x_filtrati_ms, x_filtrati_ms*b+ a,color='red', ls='-.', lw=1, label='Bestfit')
ax2.set_xlabel('t [ms]')
ax2.set_ylabel('log(R$_C$(t))')
ax2.set_title('Circuito RC in regime impulsivo con  $\\tau >>$ T (logaritmo)')
ax2.legend()
ax2.grid(True)  # Abilita la griglia


plt.tight_layout()
plt.show()


# In[34]:


import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

# Sostituisci con il nome del tuo file Excel
nome_colonna_1 = 'Vc/V0'
nome_colonna_2 = 't/T'

# Carica il file Excel in un DataFrame pandas
df = pd.read_excel(io=os.path.abspath('Misure oscilloscopio VEREABEM.xlsx'), sheet_name='CR 500')
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
err_y_filtered = err_y[indici_filtrati]
# Moltiplica i valori di x per 1000
x_filtrati_ms = x_filtrati * 1000  # Converti in millisecondi

# Creazione del subplot con asse x logaritmico
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
text_box = {
    'boxstyle': 'round',
    'facecolor': 'lightblue',
    'alpha': 1
}

# Testo con andate a capo
testo = " $log(R$_C$(t))$= b$\cdot t$ + a \n a = 0.66 $\pm$ 0.06 \n b = -11.21 $\pm$ 0.07 $ms^{-1}$ "
# Posiziona il testo in alto a destra
ax2.text(0.01, 0.14, testo, ha='left', va='top', transform=ax2.transAxes, bbox=text_box)
# Primo subplot
ax1.errorbar(x_filtrati_ms, y_corrispondenti,xerr=np.ones(len(x_filtrati_ms))*0.000006, yerr=err_y_filtered, fmt='s', markersize=2., capsize=2.5,label='Dati CR')
ax1.set_xlabel('t [ms]')
ax1.set_ylabel('R$_C$(t)')
ax1.set_title('Circuito CR in regime impulsivo con  $\\tau <<$ T')
ax1.legend()
ax1.grid(True)  # Abilita la griglia
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
a, b, sigma_a, sigma_b, cov_ab, chi2 = fit_lineare_pesato(x_filtrati_ms,np.log(y_corrispondenti) , err_y_filtered*10)

# Secondo subplot
ax2.errorbar(x_filtrati_ms, np.log(y_corrispondenti), yerr=err_y_filtered*10, fmt='s', markersize=2., capsize=2.5,label='Dati CR')
ax2.plot(x_filtrati_ms, x_filtrati_ms*b+ a,color='red', ls='-.', lw=1, label='Bestfit')
ax2.set_xlabel('t [ms]')
ax2.set_ylabel('log(R$_C$(t))')
ax2.set_title('Circuito CR in regime impulsivo con  $\\tau <<$ T (logaritmo)')
ax2.legend()
ax2.grid(True)  # Abilita la griglia


plt.tight_layout()
plt.show()


# In[ ]:





# In[31]:


import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

# Sostituisci con il nome del tuo file Excel
nome_colonna_1 = 'Vc/V0'
nome_colonna_2 = 't/T'

# Carica il file Excel in un DataFrame pandas
df = pd.read_excel(io=os.path.abspath('Misure oscilloscopio VEREABEM.xlsx'), sheet_name='CR 9750')
# Converti le colonne specificate in float, gestendo i valori non convertibili come NaN
df[nome_colonna_1] = pd.to_numeric(df[nome_colonna_1], errors='coerce')
df[nome_colonna_2] = pd.to_numeric(df[nome_colonna_2], errors='coerce')

# Estrai i valori dalle due colonne specificate
valori_colonna_1 = df[nome_colonna_1].values
valori_colonna_2 = df[nome_colonna_2].values

# Filtra i valori della colonna_x compresi tra -1000 e 1000
indici_filtrati = np.where((valori_colonna_2 > 0.0001) & (valori_colonna_2 <= 0.00015))
x_filtrati = valori_colonna_2[indici_filtrati]
y_corrispondenti = valori_colonna_1[indici_filtrati]
err_y_filtered = err_y[indici_filtrati]
# Moltiplica i valori di x per 1000
x_filtrati_ms = x_filtrati * 1000  # Converti in millisecondi

# Creazione del subplot con asse x logaritmico
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
text_box = {
    'boxstyle': 'round',
    'facecolor': 'lightblue',
    'alpha': 1
}

# Testo con andate a capo
testo = " $log(R$_C$(t))$= b$\cdot t$ + a \n a = 1.15 $\pm$ 0.05 \n b = -9.29 $\pm$ 0.01 $ms^{-1}$ "
# Posiziona il testo in alto a destra
ax2.text(0.01, 0.14, testo, ha='left', va='top', transform=ax2.transAxes, bbox=text_box)
# Primo subplot
ax1.errorbar(x_filtrati_ms, y_corrispondenti,xerr=np.ones(len(x_filtrati_ms))*0.000006, yerr=err_y_filtered, fmt='s', markersize=2., capsize=2.5,label='Dati CR')
ax1.set_xlabel('t [ms]')
ax1.set_ylabel('R$_C$(t)')
ax1.set_title('Circuito CR in regime impulsivo con  $\\tau \\approx$ T')
ax1.legend()
ax1.grid(True)  # Abilita la griglia
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
a, b, sigma_a, sigma_b, cov_ab, chi2 = fit_lineare_pesato(x_filtrati_ms,np.log(y_corrispondenti) , err_y_filtered*10)

# Secondo subplot
ax2.errorbar(x_filtrati_ms, np.log(y_corrispondenti), yerr=err_y_filtered*2, fmt='s', markersize=2., capsize=2.5, label='Dati CR')
ax2.plot(x_filtrati_ms, x_filtrati_ms*b+ a,color='red', ls='-.', lw=1, label='Bestfit')
ax2.set_xlabel('t [ms]')
ax2.set_ylabel('log(R$_C$(t))')
ax2.set_title('Circuito CR in regime impulsivo con  $\\tau \\approx$ T (logaritmo)')
ax2.legend()
ax2.grid(True)  # Abilita la griglia


plt.tight_layout()
plt.show()


# In[33]:


import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

# Sostituisci con il nome del tuo file Excel
nome_colonna_1 = 'Vc/V0'
nome_colonna_2 = 't/T'

# Carica il file Excel in un DataFrame pandas
df = pd.read_excel(io=os.path.abspath('Misure oscilloscopio VEREABEM.xlsx'), sheet_name='CR 40000')
# Converti le colonne specificate in float, gestendo i valori non convertibili come NaN
df[nome_colonna_1] = pd.to_numeric(df[nome_colonna_1], errors='coerce')
df[nome_colonna_2] = pd.to_numeric(df[nome_colonna_2], errors='coerce')

# Estrai i valori dalle due colonne specificate
valori_colonna_1 = df[nome_colonna_1].values
valori_colonna_2 = df[nome_colonna_2].values

# Filtra i valori della colonna_x compresi tra -1000 e 1000
indici_filtrati = np.where((valori_colonna_2 > 0.000025) & (valori_colonna_2 <= 0.000036))
x_filtrati = valori_colonna_2[indici_filtrati]
y_corrispondenti = valori_colonna_1[indici_filtrati]
err_y_filtered = err_y[indici_filtrati]
# Moltiplica i valori di x per 1000
x_filtrati_ms = x_filtrati * 1000  # Converti in millisecondi

# Creazione del subplot con asse x logaritmico
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
text_box = {
    'boxstyle': 'round',
    'facecolor': 'lightblue',
    'alpha': 1
}

# Testo con andate a capo
testo = " $log(R$_C$(t))$= b$\cdot t$ + a \n a = 0.28 $\pm$ 0.01 \n b = -8.62 $\pm$ 0.08 $ms^{-1}$ "
# Posiziona il testo in alto a destra
ax2.text(0.01, 0.14, testo, ha='left', va='top', transform=ax2.transAxes, bbox=text_box)
# Primo subplot
ax1.errorbar(x_filtrati_ms, y_corrispondenti,xerr=np.ones(len(x_filtrati_ms))*0.000006, yerr=err_y_filtered/10, fmt='s', markersize=2., capsize=2.5,label='Dati CR')
ax1.set_xlabel('t [ms]')
ax1.set_ylabel('R$_C$(t)')
ax1.set_title('Circuito CR in regime impulsivo con  $\\tau >>$ T')
ax1.legend()
ax1.grid(True)  # Abilita la griglia
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
a, b, sigma_a, sigma_b, cov_ab, chi2 = fit_lineare_pesato(x_filtrati_ms,np.log(y_corrispondenti) , err_y_filtered*10)

# Secondo subplot
ax2.errorbar(x_filtrati_ms, np.log(y_corrispondenti), yerr=err_y_filtered/10, fmt='s', markersize=2., capsize=2.5,label='Dati CR')
ax2.plot(x_filtrati_ms, x_filtrati_ms*b+ a,color='red', ls='-.', lw=1, label='Bestfit')
ax2.set_xlabel('t [ms]')
ax2.set_ylabel('log(R$_C$(t))')
ax2.set_title('Circuito CR in regime impulsivo con  $\\tau >>$ T (logaritmo)')
ax2.legend()
ax2.grid(True)  # Abilita la griglia


plt.tight_layout()
plt.show()


# In[ ]:





# In[ ]:




