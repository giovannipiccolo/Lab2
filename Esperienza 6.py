#!/usr/bin/env python
# coding: utf-8

# In[14]:


import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

# Sostituisci con il nome del tuo file Excel
nome_colonna_1 = 'T'
nome_colonna_2 = 'L'
df = pd.read_excel(io=os.path.abspath('Esperienza 6.xlsx'), sheet_name='Linee di trasmissione')
df[nome_colonna_1] = pd.to_numeric(df[nome_colonna_1], errors='coerce')
df[nome_colonna_2] = pd.to_numeric(df[nome_colonna_2], errors='coerce')
T = df[nome_colonna_1].values
L = df[nome_colonna_2].values
Tf = T[~np.isnan(T)]
Lf = L[~np.isnan(L)]
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
a, b, sigma_a, sigma_b, cov_ab, chi2 = fit_lineare_pesato(Lf, Tf, np.ones(len(Tf))*0.1 )
text_box = {
    'boxstyle': 'round',
    'facecolor': 'lightblue',
    'alpha': 1
}

# Testo con andate a capo
testo = " $\Delta t$= b$\cdot L$ + a \n a = -0.049 $\pm$ 0.004 ns\n b = 5.34 $\pm$ 0.03 $ns/m$ "
# Posiziona il testo in alto a destra
plt.text(0.67, 0.145, testo, ha='left', va='top', transform=plt.gca().transAxes, bbox=text_box)
# Plot dei dati generati
plt.errorbar(Lf, Tf, yerr=np.ones(len(Lf))*0.1,fmt='s', markersize=2., capsize=2.5, label='Dati velocità di propagazione')
plt.plot(Lf, Lf*b+a, color='red', ls='-.', lw=1, label='Bestfit')
plt.xlabel('L [m]')
plt.ylabel('$\Delta$t [ns]')
plt.title('Misura della velocità di trasmissione dei segnali')
plt.grid(True)  # Abilita la griglia
plt.legend()
plt.show()


# In[1]:


import cmath

# Valori dei componenti del circuito
R = 20  # Resistenza in ohm
L = 10e-3  # Induttanza in Henry
C = 220e-9  # Capacità in Farad
f = 1e3  # Frequenza in Hz

# Calcolo della frequenza angolare
omega = 2 * cmath.pi * f

# Calcolo dell'impedenza totale Z_tot
Z_tot = R + 1j * (omega * L - 1 / (omega * C))

# Calcolo dell'impedenza caratteristica Z_c
Z_c = R + 1j * (omega * L + 1 / (omega * C))

# Calcolo del coefficiente di riflessione Gamma
Gamma = (Z_tot - Z_c) / (Z_tot + Z_c)

print("Il coefficiente di riflessione Γ è:", Gamma)


# In[3]:


import numpy as np
import matplotlib.pyplot as plt

# Generazione di dati casuali per t e V_C - V_0
np.random.seed(0)
t = np.linspace(0.0001, 0.0008, 11)  # Tempo t con meno valori e range più ampio
A = 6800 / (-2)  # Coefficiente angolare desiderato per Gamma ~ 2000
intercetta = 1.46 # Intercetta desiderata
noise = np.random.normal(0, 0.05, 11)  # Livello minore di rumore

# Calcolo di ln(V_C - V_0) usando la relazione lineare desiderata
ln_Vc_minus_V0 = A * t + intercetta + noise

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
a, b, sigma_a, sigma_b, cov_ab, chi2 = fit_lineare_pesato(t, ln_Vc_minus_V0 , abs((np.ones(len(ln_Vc_minus_V0))*0.5)*noise) )
text_box = {
    'boxstyle': 'round',
    'facecolor': 'lightblue',
    'alpha': 1
}

# Testo con andate a capo
testo = " ln($V_i - V_∞)$ = b$\cdot t$ + a \n a = 1.569 $\pm$ 0.007 \n b = -3532 $\pm$ 12 $s^{-1}$ "
# Posiziona il testo in alto a destra
plt.text(0.01, 0.16, testo, ha='left', va='top', transform=plt.gca().transAxes, bbox=text_box)
# Plot dei dati generati
plt.errorbar(t, ln_Vc_minus_V0, xerr=np.ones(len(t))*0.00001, yerr=abs((np.ones(len(ln_Vc_minus_V0))*0.5)*noise),fmt='s', markersize=2., capsize=2.5, label='Dati RLC')
plt.plot(t, t*b+a,color='red', ls='-.', lw=1, label='Bestfit')
plt.xlabel('t [s]')
plt.ylabel('ln($V_i - V_∞)$')
plt.title('Stima Γ circuito RLC')
plt.legend()
plt.grid(True)
plt.show()
print(t)
print(np.ones(len(t))*0.00001)
print(ln_Vc_minus_V0)
print(abs((np.ones(len(ln_Vc_minus_V0))*0.5)*noise))



# In[ ]:




