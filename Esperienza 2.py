#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy as sp
import math

I_R = np.array([0.145, 0.355, 0.5, 0.75, 1.4, 1.75, 2.3, 3.4, 6.5, 12, 14, 20, 27.5, 34, 44, 49.5])
V_R = np.array([14.96, 14.91, 14.88, 14.83, 14.67, 14.59, 14.47, 14.23, 13.57, 12.33, 11.83, 10.53, 8.97, 7.56, 5.35, 4.33])

sigmaI = 0.3*np.ones(len(I_R))
sigmaV = 0.1*np.ones(len(V_R))

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
a_R, b_R, sigma_a_R, sigma_b_R, cov_ab_R, chi2_R = fit_lineare_pesato(I_R, V_R, 1/sigmaV**2)




# In[3]:


plt.errorbar(I_R, V_R, xerr=sigmaI, yerr=sigmaV,fmt='s', markersize=2., capsize=1.5)

x = np.linspace(0, 50)
#plt.plot(x, a_R + b_R*x, color='green', ls='-.', lw=1, label='bestfit')

plt.legend()

plt.ylabel('$V_R$ [V]', size=14)
plt.xlabel('$I_R$ [mA]', size=14)
plt.grid()
plt.title('Misura resistenza interna di generatore di tensione')
plt.ylim(4, 17)
plt.savefig('misuraresistenzainterna.png', dpi=200)
plt.show()


# In[4]:


I=np.array([0.015, 0.02, 0.03, 0.05, 0.07, 0.14, 0.155, 0.17, 0.195, 0.225, 0.265, 0.33, 0.42, 0.6, 1, 1.1, 1.2, 1.3, 1.4, 1.55, 1.7, 1.95, 2.25, 2.65, 2.70, 2.75, 2.9])
I_p=2*I
V_p=np.array([7.455, 7.435, 7.425, 7.385, 7.33, 7.18, 7.145, 7.105, 7.055, 6.99, 6.9, 6.765, 6.56, 6.16, 5.235, 5.06, 4.865, 4.63, 4.36, 4.03, 3.62, 3.15, 2.42, 1.495, 1.39, 1.28, 0.895])
sigmaI = 0.14*np.ones(len(I_p))
sigmaV = 0.1*np.ones(len(V_p))
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
a_R, b_R, sigma_a_R, sigma_b_R, cov_ab_R, chi2_R = fit_lineare_pesato(I_p, V_p, 1/sigmaV**2)


# In[5]:


plt.errorbar(I_p, V_p, xerr=sigmaI, yerr=sigmaV,fmt='s',markeredgecolor='orange', markersize=2., capsize=1.5, ecolor='orange')

x = np.linspace(-0.1, 6)
plt.plot(x, a_R + b_R*x, color='orange', ls='-.', lw=1, label='bestfit')

plt.legend()

plt.ylabel('$V_R$ [V]', size=14)
plt.xlabel('$I_R$ [mA]', size=14)
plt.grid()
plt.title('Applicazione al partitore resistivo')
plt.ylim(0, 8)
plt.savefig('partitoreresistivo.png', dpi=200)
plt.show()


# In[6]:


capocchia=np.array([7.49, 3.75, 1.87, 0.93])
x=[1, 2, 3, 4]
plt.errorbar(x, capocchia, yerr=0.14, fmt=".")
plt.ylabel('$V_R$ [V]', size=14)
plt.xlabel('n', size=14)
plt.grid()
plt.title('Studio di un partitore multi-stadio')
plt.ylim(0, 8)
plt.savefig('partitoremultistadiodiegoarmandomaradona.png', dpi=200)
plt.show()


# In[9]:


logcapocchia=np.log(capocchia)
errcapocchia = 0.14/capocchia
plt.errorbar(x, logcapocchia, yerr=errcapocchia,fmt="s", markeredgecolor='red', markerfacecolor='red', markersize=2, ecolor='red')
plt.plot(x, logcapocchia, linestyle='-.',color='red')
plt.ylabel('log($V_R$)', size=14)
plt.xlabel('n', size=14)
plt.grid()
plt.title('Studio di un partitore multi-stadio')
plt.ylim(-1, 3)
plt.savefig('partitorelogaritmomultistadiodiegoarmandomaradona.png', dpi=200)
plt.show()


# In[ ]:




