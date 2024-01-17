#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy as sp
import math


# In[25]:


A_RC=np.array([1.000, 0.992, 0.992, 0.985, 0.985, 0.969, 0.939, 0.878, 0.847, 0.817, 0.779, 0.737, 0.714, 0.676, 0.650, 0.620, 0.599, 0.534, 0.359, 0.187, 0.114, 0.0390, 0.0252])
err_RC=np.array([0.016, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.014, 0.013, 0.013, 0.012, 0.011, 0.011, 0.011, 0.010, 0.010, 0.009, 0.008, 0.006, 0.003, 0.002, 0.0010, 0.0004])
f_RC=np.array([50, 75.19, 100, 299.4, 598.8, 990.1, 2976, 5952, 6944, 7874, 8929, 10000, 11013, 12077, 12887, 13966, 14793, 17986, 29940, 59880, 100000, 299400, 502510])
err_f_RC=np.array([0.02, 0.03, 0.06, 0.1, 0.2, 0.6, 1, 2, 3, 4, 5, 6, 2, 3, 3, 4, 4, 6, 11, 20, 60, 100, 150])
LOGerr_f_RC= err_f_RC/ f_RC
A_CR=np.array([0.0435, 0.00641, 0.00908, 0.0259, 0.0519, 0.0870, 0.246,0.454, 0.504, 0.549, 0.610, 0.646, 0.690, 0.714, 0.733, 0.767, 0.809, 0.839, 0.901, 0.931, 0.947, 0.969, 0.992])
err_CR=np.array([0.00013, 0.00015, 0.00018, 0.0004, 0.0008, 0.0014, 0.004, 0.007, 0.008, 0.009, 0.009, 0.010, 0.010, 0.011, 0.011, 0.012, 0.012, 0.013, 0.014, 0.014, 0.015, 0.015, 0.015])
f_CR=np.array([49.500, 74.07, 100.00, 297.62, 606.1, 1010.1, 3012.0, 5988, 6896, 7874, 9009, 9804, 10965, 11972, 12626, 13812, 14793, 17730, 29940, 59980, 101010, 297600, 500000])
err_f_CR=np.array([0.010, 0.03, 0.06, 0.09, 0.2, 0.6, 1.0, 2, 3, 4, 5, 6, 2, 3, 3, 4, 4, 6, 11, 20, 60, 100, 150])

LOGerr_f_CR= err_f_CR/f_CR
plt.errorbar(np.log(f_RC), A_RC, xerr=LOGerr_f_RC, yerr=err_RC,fmt='s', markersize=2., capsize=2.5, label='Dati RC')
plt.errorbar(np.log(f_CR), A_CR, xerr=LOGerr_f_CR, yerr=err_CR,fmt='s', markersize=2., capsize=2.5, label='Dati CR')
plt.plot(np.log(f_RC),A_RC, color='red', ls='-.', lw=1, label=' RC')
plt.plot(np.log(f_CR),A_CR, color='green', ls='-.', lw=1, label=' CR')
plt.legend()
plt.ylabel('A', size=14)
plt.xlabel('log($\\nu$)', size=14)
plt.grid()
plt.title('Studio dei circuiti RC e CR in regime sinusoidale')
plt.savefig('capocchia dolce', dpi=200)
plt.show()


# In[23]:


A_RC_B=np.array([0.847, 0.817, 0.779, 0.737, 0.714, 0.676, 0.650, 0.620, 0.599])
err_RC_B=np.array([0.013, 0.013, 0.012, 0.011, 0.011, 0.011, 0.010, 0.010, 0.009,])
f_RC_B=np.array([5952, 6944, 7874, 8929, 10000, 11013, 12077, 12887, 13966])
err_f_RC_B=np.array([3, 4, 5, 6, 2, 3, 3, 4, 4])
LOGerr_f_RC_B= err_f_RC_B/ f_RC_B
A_CR_B=np.array([0.454, 0.504, 0.549, 0.610, 0.646, 0.690, 0.714, 0.733, 0.767])
err_CR_B=np.array([ 0.007, 0.008, 0.009, 0.009, 0.010, 0.010, 0.011, 0.011, 0.012])
f_CR_B=np.array([ 5988, 6896, 7874, 9009, 9804, 10965, 11972, 12626, 13812])
err_f_CR_B=np.array([2, 3, 4, 5, 6, 2, 3, 3, 4])
LOGerr_f_CR_B= err_f_CR_B/f_CR_B

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
a_RC, b_RC, sigma_a_RC, sigma_b_RC, cov_ab_RC, chi2_RC = fit_lineare_pesato(np.log(f_RC_B), A_RC_B,  1/err_RC_B**2)
a_CR, b_CR, sigma_a_CR, sigma_b_CR, cov_ab_CR, chi2_CR = fit_lineare_pesato(np.log(f_CR_B), A_CR_B,  1/err_CR_B**2)

x=np.linspace(8.6, 9.6)
plt.errorbar(np.log(f_RC_B), A_RC_B, xerr=LOGerr_f_RC_B, yerr=err_RC_B,fmt='s', markersize=2., capsize=2.5, label='Dati RC')
plt.errorbar(np.log(f_CR_B), A_CR_B, xerr=LOGerr_f_CR_B, yerr=err_CR_B,fmt='s', markersize=2., capsize=2.5, label='Dati CR')

plt.plot(x,a_RC+b_RC*x, color='red', ls='-.', lw=1, label='bestfit RC')
plt.plot(x,a_CR+b_CR*x, color='green', ls='-.', lw=1, label='bestfit CR')
plt.legend()
plt.ylabel('A', size=14)
plt.xlabel('log($\\nu$)', size=14)
plt.grid()
plt.title('Studio dei circuiti RC e CR bestfit')
plt.savefig('mAgIkArP', dpi=200)
plt.show()


# In[29]:


ni_RC=np.array([50, 75.19,100, 299.40, 598.8, 990.1, 2976, 5952, 6944, 7874, 8929, 10000, 11013, 12077, 12887, 13966, 14793, 17986, 29940, 59880, 100000, 299400, 502510])
err_ni_RC=np.array([0.02, 0.03, 0.06, 0.11, 0.2, 0.6, 1.0, 2, 3, 4, 5, 6, 2, 3, 3, 4, 4, 5, 11, 20, 60, 110, 150])
LOGerr_ni_RC= err_ni_RC/ ni_RC
phirc=np.array([0.0003769, 0.0007180, -0.00130, -0.00319, -0.00545, -0.00852, -0.2520, -0.470, -0.540, -0.604, -0.667, -0.735, -0.789, -0.835, -0.8664, -0.8951, -0.9109, -0.9944, -1.1810, -1.3846, -1.4891, -1.4972, -1.4990])
err_phirc=np.array([0.0000006, 0.0000009, 0.00010, 0.00002, 0.00003, 0.00006, 0.0010, 0.002, 0.003, 0.004, 0.005, 0.006, 0.005, 0.005, 0.0005, 0.0006, 0.0006, 0.0008, 0.0006, 0.0007, 0.0009, 0.008, 0.0010])
phicr=np.array([1.4560, 1.5359, 1.5460, 1.5259, 1.5156, 1.4980, 1.3096, 1.0760, 1.0053, 0.9549, 0.9057,0.8501, 0.8129, 0.7557, 0.7044, 0.6839, 0.6581, 0.5659, 0.3857, 0.2047, 0.12760, 0.03550, 0.0207])
err_phicr=np.array([0.0008, 0.0009, 0.0010, 0.0005, 0.0007, 0.0010, 0.0006, 0.0006, 0.0005, 0.0005, 0.0005, 0.0006, 0.0006, 0.0005, 0.0005, 0.0002, 0.0003, 0.0003, 0.0003, 0.0002, 0.00010, 0.00010, 0.0009])
ni_CR=np.array([49.500, 74.07,100.00, 297.62, 606.1, 1010.1, 3012.0, 5988, 6896, 7874, 9009, 9804, 10965, 11792, 12626, 13812, 14793, 17730, 29940, 59880, 101000, 297620, 500000])
err_ni_CR=np.array([0.010, 0.03, 0.06, 0.09, 0.2, 0.6, 1.0, 2, 3, 4, 5, 6, 2, 3, 3, 4, 4, 6, 11, 20, 60, 110, 150])
LOGerr_ni_CR= err_ni_CR/ ni_CR
plt.errorbar(np.log10(ni_RC),phirc,xerr=LOGerr_ni_RC*100,yerr=err_phirc, fmt='.',label='Dati RC')
plt.errorbar(np.log10(ni_CR),phicr,xerr=LOGerr_ni_CR*100,yerr=err_phicr, fmt='.',label='Dati CR')
#plt.plot(np.log10(ni_RC),phirc, ls='-.',c='red',label=' RC')
#plt.plot(np.log10(ni_CR),phicr,ls='-.', c='green', label=' CR')
plt.legend()
plt.ylabel('$\Phi$', size=14)
plt.xlabel('log($\\nu$)')
plt.title('Sfasamento dei circuiti RC e CR')
plt.savefig('mAgIkArP', dpi=200)
plt.grid()
plt.show()


# In[28]:


ni_RC=np.array([50, 75.19,100, 299.40, 598.8, 990.1, 2976, 5952, 6944, 7874, 8929, 10000, 11013, 12077, 12887, 13966, 14793, 17986, 29940, 59880, 100000, 299400, 502510])
err_ni_RC=np.array([0.02, 0.03, 0.06, 0.11, 0.2, 0.6, 1.0, 2, 3, 4, 5, 6, 2, 3, 3, 4, 4, 5, 11, 20, 60, 110, 150])
LOGerr_ni_RC= err_ni_RC/ ni_RC
phirc=np.array([0.0003769, 0.0007180, -0.00130, -0.00319, -0.00545, -0.00852, -0.2520, -0.470, -0.540, -0.604, -0.667, -0.735, -0.789, -0.835, -0.8664, -0.8951, -0.9109, -0.9944, -1.1810, -1.3846, -1.4891, -1.4972, -1.4990])
err_phirc=np.array([0.0000006, 0.0000009, 0.00010, 0.00002, 0.00003, 0.00006, 0.0010, 0.002, 0.003, 0.004, 0.005, 0.006, 0.005, 0.005, 0.0005, 0.0006, 0.0006, 0.0008, 0.0006, 0.0007, 0.0009, 0.008, 0.0010])
phicr=np.array([1.4560, 1.5359, 1.5460, 1.5259, 1.5156, 1.4980, 1.3096, 1.0760, 1.0053, 0.9549, 0.9057,0.8501, 0.8129, 0.7557, 0.7044, 0.6839, 0.6581, 0.5659, 0.3857, 0.2047, 0.12760, 0.03550, 0.0207])
err_phicr=np.array([0.0008, 0.0009, 0.0010, 0.0005, 0.0007, 0.0010, 0.0006, 0.0006, 0.0005, 0.0005, 0.0005, 0.0006, 0.0006, 0.0005, 0.0005, 0.0002, 0.0003, 0.0003, 0.0003, 0.0002, 0.00010, 0.00010, 0.0009])
ni_CR=np.array([49.500, 74.07,100.00, 297.62, 606.1, 1010.1, 3012.0, 5988, 6896, 7874, 9009, 9804, 10965, 11792, 12626, 13812, 14793, 17730, 29940, 59880, 101000, 297620, 500000])
err_ni_CR=np.array([0.010, 0.03, 0.06, 0.09, 0.2, 0.6, 1.0, 2, 3, 4, 5, 6, 2, 3, 3, 4, 4, 6, 11, 20, 60, 110, 150])
LOGerr_ni_CR= err_ni_CR/ ni_CR
plt.errorbar(np.log10(ni_RC),np.arctan(phirc),xerr=LOGerr_ni_RC*100,yerr=err_phirc, fmt='.',label='Dati RC')
plt.errorbar(np.log10(ni_CR),np.arctan(phicr),xerr=LOGerr_ni_CR*100,yerr=err_phicr, fmt='.',label='Dati CR')
#plt.plot(np.log10(ni_RC),np.arctan(phirc), ls='-.', label=' RC', c='red')
#plt.plot(np.log10(ni_CR),np.arctan(phicr), ls='-.', label=' CR', c='green')
plt.title('Sfasamento dei circuiti RC e CR(arcotangente)')
plt.ylabel('arctan($\Phi$)', size=14)
plt.xlabel('log($\\nu$)')
plt.legend()
plt.grid()
plt.show()


# In[190]:


x=np.linspace(3.5,13, num=23)
A_RC=np.array([1.000, 0.992, 0.992, 0.985, 0.985, 0.969, 0.939, 0.878, 0.847, 0.817, 0.779, 0.737, 0.714, 0.676, 0.650, 0.620, 0.599, 0.534, 0.359, 0.187, 0.114, 0.0390, 0.0252])
err_RC=np.array([0.016, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.014, 0.013, 0.013, 0.012, 0.011, 0.011, 0.011, 0.010, 0.010, 0.009, 0.008, 0.006, 0.003, 0.002, 0.0010, 0.0004])
f_RC=np.array([50, 75.19, 100, 299.4, 598.8, 990.1, 2976, 5952, 6944, 7874, 8929, 10000, 11013, 12077, 12887, 13966, 14793, 17986, 29940, 59880, 100000, 299400, 502510])
err_f_RC=np.array([0.02, 0.03, 0.06, 0.1, 0.2, 0.6, 1, 2, 3, 4, 5, 6, 2, 3, 3, 4, 4, 6, 11, 20, 60, 100, 150])
LOGerr_f_RC= err_f_RC/ f_RC
A_CR=np.array([0.0435, 0.00641, 0.00908, 0.0259, 0.0519, 0.0870, 0.246,0.454, 0.504, 0.549, 0.610, 0.646, 0.690, 0.714, 0.733, 0.767, 0.809, 0.839, 0.901, 0.931, 0.947, 0.969, 0.992])
err_CR=np.array([0.00013, 0.00015, 0.00018, 0.0004, 0.0008, 0.0014, 0.004, 0.007, 0.008, 0.009, 0.009, 0.010, 0.010, 0.011, 0.011, 0.012, 0.012, 0.013, 0.014, 0.014, 0.015, 0.015, 0.015])
f_CR=np.array([49.500, 74.07, 100.00, 297.62, 606.1, 1010.1, 3012.0, 5988, 6896, 7874, 9009, 9804, 10965, 11972, 12626, 13812, 14793, 17730, 29940, 59980, 101010, 297600, 500000])
err_f_CR=np.array([0.010, 0.03, 0.06, 0.09, 0.2, 0.6, 1.0, 2, 3, 4, 5, 6, 2, 3, 3, 4, 4, 6, 11, 20, 60, 100, 150])

LOGerr_f_CR= err_f_CR/f_CR
plt.errorbar(x, A1, xerr=LOGerr_f_RC, yerr=err_RC,fmt='s', markersize=2., capsize=2.5, label='Dati RC')
plt.errorbar(x, A2, xerr=LOGerr_f_CR, yerr=err_CR,fmt='s', markersize=2., capsize=2.5, label='Dati CR')

a=np.exp(x)*2*math.pi
A1=1/(1+np.sqrt((a**2)*((2.2*10**-9)**2)*(8200**2)))
A2=a*(2.2*10**-9)*8200/(1+np.sqrt((a**2)*((2.2*10**-9)**2)*(8200**2)))
plt.plot(x,A1, color='teal', ls='-.', lw=1, label=' CR')
plt.plot(x,A2, color='orange', ls='-.', lw=1, label=' CR')
plt.legend()
plt.ylabel('A', size=14)
plt.xlabel('log($\\nu$)', size=14)
plt.grid()
plt.title('Studio dei circuiti RC e CR in regime sinusoidale')
plt.savefig('capocchia dolce', dpi=200)
plt.show()
print(A)
#A_RC_B=np.array([0.847, 0.817, 0.779, 0.737, 0.714, 0.676, 0.650, 0.620, 0.599])
err_RC_B=np.array([0.012, 0.011, 0.011, 0.011, 0.010])
f_RC_B=np.array([5952, 6944, 7874, 8929, 10000, 11013, 12077, 12887, 13966])
err_f_RC_B=np.array([ 5, 6, 2, 3, 3])
LOGerr_f_RC_B= err_f_RC_B/ XX
#A_CR_B=np.array([0.454, 0.504, 0.549, 0.610, 0.646, 0.690, 0.714, 0.733, 0.767])
err_CR_B=np.array([0.009, 0.009, 0.010, 0.010, 0.011])
f_CR_B=np.array([ 5988, 6896, 7874, 9009, 9804, 10965, 11972, 12626, 13812])
err_f_CR_B=np.array([4, 5, 6, 2, 3])
LOGerr_f_CR_B= err_f_CR_B/XX

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
a_RC, b_RC, sigma_a_RC, sigma_b_RC, cov_ab_RC, chi2_RC = fit_lineare_pesato(XX, A_1,  1/err_RC_B**2)
a_CR, b_CR, sigma_a_CR, sigma_b_CR, cov_ab_CR, chi2_CR = fit_lineare_pesato(XX, A_2,  1/err_CR_B**2)

x=np.linspace(8.6, 9.6)
#plt.errorbar(np.log(f_RC_B), A_RC_B, xerr=LOGerr_f_RC_B, yerr=err_RC_B,fmt='s', markersize=2., capsize=2.5, label='Dati RC')
#plt.errorbar(np.log(f_CR_B), A_CR_B, xerr=LOGerr_f_CR_B, yerr=err_CR_B,fmt='s', markersize=2., capsize=2.5, label='Dati CR')
plt.errorbar(XX, A_1, xerr=LOGerr_f_RC_B/10, yerr=err_RC_B,fmt='s', markersize=2., capsize=2.5, label='Dati RC')
plt.errorbar(XX, A_2, xerr=LOGerr_f_CR_B/10, yerr=err_CR_B,fmt='s', markersize=2., capsize=2.5, label='Dati CR')
XX=np.linspace(8,10, num=5)
a=np.exp(XX)*2*math.pi
A_2=a*(2.2*10**-9)*8200/(1+np.sqrt((a**2)*((2.2*10**-9)**2)*(8200**2)))
A_1=1/(1+np.sqrt((a**2)*((2.2*10**-9)**2)*(8200**2)))
plt.plot(XX,a_RC+b_RC*XX, color='teal', ls='-.', lw=1, label='bestfit RC')
plt.plot(XX,a_CR+b_CR*XX, color='orange', ls='-.', lw=1, label='bestfit CR')
plt.legend()
plt.ylabel('A', size=14)
plt.xlabel('log($\\nu$)', size=14)
plt.grid()
plt.title('Studio dei circuiti RC e CR bestfit')
plt.savefig('mAgIkArP', dpi=200)
plt.show()


# In[227]:


import numpy as np
x=np.linspace(1.5,5.5, num=23)
a=np.exp(x)*2*math.pi
phirc=(-a*(2.2*10**-9)*8200)
phirc1=np.arctan(phirc)
phicr=np.arctan(1/(a*(2.2*10**-9)*8200))
err_ni_RC=np.array([0.02, 0.03, 0.06, 0.11, 0.2, 0.6, 1.0, 2, 3, 4, 5, 6, 2, 3, 3, 4, 4, 5, 11, 20, 60, 110, 150])
LOGerr_ni_RC= err_ni_RC/ ni_RC
err_phirc=np.array([0.0000006, 0.0000009, 0.00010, 0.00002, 0.00003, 0.00006, 0.0010, 0.002, 0.003, 0.004, 0.005, 0.006, 0.005, 0.005, 0.0005, 0.0006, 0.0006, 0.0008, 0.0006, 0.0007, 0.0009, 0.008, 0.0010])
#phicr=np.array([1.4560, 1.5359, 1.5460, 1.5259, 1.5156, 1.4980, 1.3096, 1.0760, 1.0053, 0.9549, 0.9057,0.8501, 0.8129, 0.7557, 0.7044, 0.6839, 0.6581, 0.5659, 0.3857, 0.2047, 0.12760, 0.03550, 0.0207])
err_phicr=np.array([0.0008, 0.0009, 0.0010, 0.0005, 0.0007, 0.0010, 0.0006, 0.0006, 0.0005, 0.0005, 0.0005, 0.0006, 0.0006, 0.0005, 0.0005, 0.0002, 0.0003, 0.0003, 0.0003, 0.0002, 0.00010, 0.00010, 0.0009])
ni_CR=np.array([49.500, 74.07,100.00, 297.62, 606.1, 1010.1, 3012.0, 5988, 6896, 7874, 9009, 9804, 10965, 11792, 12626, 13812, 14793, 17730, 29940, 59880, 101000, 297620, 500000])
err_ni_CR=np.array([0.010, 0.03, 0.06, 0.09, 0.2, 0.6, 1.0, 2, 3, 4, 5, 6, 2, 3, 3, 4, 4, 6, 11, 20, 60, 110, 150])
LOGerr_ni_CR= err_ni_CR/ ni_CR
plt.errorbar(x,phirc,xerr=LOGerr_ni_RC*100,yerr=err_phirc, fmt='.')
plt.errorbar(x,phicr,xerr=LOGerr_ni_CR*100,yerr=err_phicr, fmt='.')
plt.plot(x,phirc, ls='-.', c='teal')
plt.plot(x,phicr,ls='-.', c='orange')
plt.grid()
plt.show()
print(phirc)



# In[ ]:




