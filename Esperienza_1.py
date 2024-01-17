#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy as sp
import math


# In[4]:


# configurazione a monte

I_m = np.array([0.02, 0.06, 0.11, 0.15, 0.2, 0.25, 0.29, 0.34, 0.39, 0.43]) #in \A
V_m = np.array([0.06,0.24,0.44,0.64,0.8,1,1.2,1.4,1.56, 1.76]) # in V

# configurazione a valle

I_v = np.array([0.02, 0.06, 0.1, 0.14, 0.18, 0.22, 0.26, 0.31, 0.35, 0.39]) # in \muA
V_v = np.array([0.08, 0.28, 0.5, 0.7, 0.92, 1.14, 1.36, 1.56, 1.76, 1.96]) # in V

fs_V = 2
dV = (fs_V/50)
err_V = (dV/2) * np.ones(len(V_m))
print(err_V)
fs_I = 500
dI = (fs_I/50)
err_I = (dI/2000)
err_I_m = err_I * np.ones(len(I_m)) # in \muA
print(err_I)
sigmaV = err_V*2/np.sqrt(12)
sigmaI = err_I*2/np.sqrt(12)

def fit_lineare_pesato(x,y,w) :
    """
    questa funzione calcola i parametri della retta di best fit y = a + b*x 
    usando le formule dei minimi quadrati pesati.
    N.B. w è il vettore dei pesi, dove il peso w_i è l'inverso del quadrato dell'incertezza 
    sulla coordinata y_i dell'i-esimo punto.    
    """
    S_0 = np.sum(w)
    S_x = np.sum(x*w)     
    S_xx = np.sum(x*x*w)     
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
a_m, b_m, sigma_a_m, sigma_b_m, cov_ab_m, chi2_m = fit_lineare_pesato(I_m, V_m, 1/sigmaV**2)
a_v, b_v, sigma_a_v, sigma_b_v, cov_ab_v, chi2_v = fit_lineare_pesato(I_v, V_v, 1/sigmaV**2)


# In[4]:


plt.errorbar(I_m, V_m, xerr=err_I_m, yerr=err_V,fmt='s', markersize=2., capsize=1.5, label='Dati conf. a monte')
plt.errorbar(I_v, V_v, xerr=err_I_m, yerr=err_V,fmt='s', markersize=2., capsize=1.5, label='Dati conf. a valle')

x = np.linspace(0, 0.45)
plt.plot(x, a_m + b_m*x, color='teal', ls='-.', lw=1, label='bestfit conf. a monte')
plt.plot(x, a_v + b_v*x, color='orange', ls='-.', lw=1, label='bestfit conf. a valle')
plt.legend()

plt.ylabel('V [V]', size=14)
plt.xlabel('I [A]', size=14)
plt.grid()
plt.title('Stima di R con metodo voltamperometrico')
plt.ylim(0, 2)
plt.savefig('metodo_voltamperometrico.png', dpi=200)
plt.show()


# In[5]:


I_p = np.array([280, 260, 240, 230, 210, 190, 180, 170, 160, 150, 135, 125, 115, 100, 95, 85, 75, 70, 60, 50, 45, 30, 25, 20, 15, 5, 0, 0, -5, -10, -15, -20, -25, -30, -35, -40, -50, -55, -55, -60, -65, -70, -75, -80, -80, -85, -90, -95, -100, -110, -115])
R_p = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390,	400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500])
dI = (fs_I/50)
err_I = (dI/2000000)
err_I_p = err_I * np.ones(len(I_p))
err_R_p = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5])


# In[6]:


plt.plot(R_p, I_p/100000, linestyle='none', color='purple')
plt.errorbar(R_p, I_p/1000000,xerr= err_R_p, yerr = err_I_p, linestyle= 'none', ecolor='purple')
plt.title('I vs $R_3$ Ponte di Wheatstone')
plt.xlabel('R [$\Omega$]', size=14)
plt.ylabel('I [A]', size=14)
plt.ylim(-0.0002, 0.00029)
plt.savefig('Ponte di Wheatstone', dpi=200)
plt.grid( color="white")
plt.axhline(0, color='black')
plt.axvline(0, color='black')
plt.show()


# In[7]:


I_p1 = np.array( [15, 5, 0, 0, -5, -10])
R_p1 = np.array([240, 250, 260, 270, 280, 290])
a_p= -0.0000004429
sigma_a_p= 0.06281
b_p= 0.0001182
sigma_b_p=16.68


# In[8]:


plt.plot(R_p, I_p/100000, linestyle='none', color='pink')
plt.errorbar(R_p, I_p/1000000,xerr= err_R_p, yerr = err_I_p, linestyle= 'none', ecolor='pink')
x1 = np.linspace(0, 239)
x = np.linspace(240, 290)
x2 = np.linspace(291, 500)
plt.plot(x1, x1*a_p+b_p, ls='dotted', color='black')
plt.plot(x, x*a_p+b_p, ls='solid', color='black', label='bestfit')
plt.plot(x2, x2*a_p+b_p, ls='dotted', color='black')
plt.legend()
plt.savefig('Ponte di wheatstone bestfit.png')
plt.errorbar(I_v, V_v, xerr=sigmaI, yerr=sigmaV, fmt='s', markersize=2., capsize=1, color='black', label='punti sperimentali (err. stat)')
plt.plot(I_v, a_v + b_v*I_v, color='red', ls='--', label='retta di bestfit')
plt.title('I vs $R_3$ Ponte di Wheatstone')
plt.xlabel('R [$\Omega]$', size=14)
plt.ylabel('I [A]', size=14)
plt.ylim(-0.0002, 0.00029)
plt.grid( color="white")
plt.axhline(0, color='black')
plt.axvline(0, color='black')
plt.show()


# In[19]:


import plotly.graph_objects as go

fig = go.Figure(data=go.Scatter(
        x=R_p,
        y=I_p,
        error_y=dict(
            type='data',
            array=err_I_p,
            visible=True),
        error_x=dict(
            type='data', 
            array=err_R_p,
            visible=True)
    ))
fig.show()


# In[17]:


fig = go.Figure(data=go.Scatter(
        x=[1, 2, 3, 4],
        y=[2, 1, 3, 4],
        error_x=dict(
            type='percent',
            value=10)
    ))
fig.show()


# In[26]:


import plotly.graph_objects as go
import numpy as np

x = R_p
y = I_p


fig = go.Figure()

fig.add_trace(go.Scatter(
    x=x, y=y,
    mode='markers',
    name='measured',
    error_y=dict(
        type='data',
        array=err_I_p,
        color='purple',
        thickness=1.5,
        width=3,
    ),
    error_x=dict(
        type='data',
        array=err_R_p,
        color='purple',
        thickness=1.5,
        width=3,
    ),
    marker=dict(color='purple', size=1)
))
fig.show()


# In[ ]:




