#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy as sp
import math

from uncertainties import ufloat
from scipy.optimize import curve_fit
Ts=np.array([0.0100000000,0.0083333333,0.0071428571,0.0062500000,0.0055555556,0.0050000000,0.0025000000,0.0016666667,0.0012500000,0.0010000000,0.0005000000,0.0004545455,0.0004166667,0.0003846154,0.0003571429,0.0003333333,0.0003125000,0.0002949853,0.0002941176,0.0002777778,0.0002631579,0.0002500000,0.0001666667,0.0001250000,0.0001000000,0.0000500000,0.0000250000,0.0000166667,0.0000125000,0.0000100000,0.0000050000,0.0000033333,0.0000025000,0.0000023810,0.0000022222,0.0000021277,0.0000020000,0.0000016667,0.0000010000,0.0000006667, 0.0000003333])
T=Ts*1000
err_T= (((0.0100000000-0.0000003333)/2000))*np.ones(len(Ts))
dt=np.array([-2.39,-2.07,-1.75,-1.55,-1.33,-1.21,-0.593,-0.381,-0.277,-0.207,-0.0693,-0.0557,-0.0415,-0.0276,-0.015,-0.0037,0.005548,0.013,0.013,0.019,0.0227,0.026,0.029,0.024,0.019,0.0092,0.00504,0.00343,0.0027,0.0022,0.00118,0.000784,0.000556,0.000484,0.000364,-0.00014,-0.000362,-0.000388,-0.000258,-0.0001584,-0.0000728])



# In[2]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from uncertainties import ufloat
import math
# Dati
Vin = np.array([15.5, 15.3, 15.3, 15.4, 15.4, 15.3, 15.3, 15.3, 15.1, 15, 13.3, 13, 12.2, 11.8, 11.3, 11, 11.1, 11.3, 11.3, 11.5, 11.8, 12.2, 14.2, 14.7, 15, 15.2, 15.2, 15, 15, 15, 15, 15, 15.4, 15.2, 15.2, 15.2, 15.4, 15.2, 15.2, 15, 14.8])
Vout = np.array([0.244, 0.292, 0.336, 0.386, 0.432, 0.48, 0.96, 1.46, 1.94, 2.46, 5.64, 6.2, 6.88, 7.52, 7.76, 7.92, 7.92, 7.92, 7.84, 7.68, 7.36, 7.04, 4.64, 3.26, 2.54, 1.24, 0.616, 0.412, 0.31, 0.244, 0.112, 0.06, 0.03, 0.0124, 0.0048, 0.0024, 0.007, 0.034, 0.0896, 0.148, 0.322])
err_Vin=np.sqrt((((0.01*3*Vin)/2.58)**2 + ((0.1*0.1)/math.sqrt(3))**2 + ((0.05)/math.sqrt(3))**2)**2)
err_Vout=np.sqrt((((0.01*3*Vout)/2.58)**2 + ((0.1*0.1)/math.sqrt(3))**2 + ((0.05)/math.sqrt(3))**2)**2)
err_A= np.sqrt((((2*err_Vout)**2)/((2*Vin)**2))+ ((((2*Vout)**2)*(2*err_Vin)**2)/(2*Vin)**4) ) 
f = np.array([100, 120, 140, 160, 180, 200, 400, 600, 800, 1000, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3390, 3400, 3600, 3800, 4000, 6000, 8000, 10000, 20000, 40000, 60000, 80000, 100000, 200000, 300000, 400000, 420000, 450000, 470000, 500000, 600000, 1000000, 1500000, 3000000])
err_f=np.ones(len(f))*0.01
err_logf=np.ones(len(f))*err_f/f
# Dati filtrati
f_filtered = np.array([2000, 2200, 2400, 2600, 2800, 3000, 3200, 3390, 3400, 3600, 3800, 4000, 6000])
attenuation_filtered = np.array([5.64, 6.2, 6.88, 7.52, 7.76, 7.92, 7.92, 7.92, 7.84, 7.68, 7.36, 7.04, 4.64])
# Filtrare i valori nell'intervallo specificato
mask = (f >= 10**3.3) & (f <= 10**3.664)
f_filtered = f[mask]
attenuation_filtered = Vout[mask] / Vin[mask]

# Definizione della funzione parabolica
def parabola(x, a, b, c):
    return a * x**2 + b * x + c

# Calcolo dei coefficienti della parabola utilizzando curve_fit solo per i valori selezionati
popt, pcov = curve_fit(parabola, f_filtered, attenuation_filtered)

# Estrazione degli errori dai parametri
perr = np.sqrt(np.diag(pcov))
a = ufloat(popt[0], perr[0])
b = ufloat(popt[1], perr[1])
c = ufloat(popt[2], perr[2])

# Stampa dei coefficienti della parabola con gli errori
print(f"I coefficienti della parabola sono: a = {c} , b = {b}, c = {a}")

# Tracciamento del grafico con l'attenuazione e la parabola intorno al picco
plt.figure(figsize=(10, 6))
plt.errorbar(f, Vout/Vin,xerr=err_f,yerr=err_A*10,fmt='s', marker='.', color='b', label='Attenuazioni RLC')


# Definizione dell'intervallo intorno al picco
x_values = np.linspace(1900, 6100, 400)
y_values = parabola(x_values, *popt)
plt.semilogx(f, Vout/Vin, color='red', ls='none' )
plt.xlabel('log($\\nu$)')
plt.ylabel('A($\\nu$)')
plt.title('Grafico attenuazioni RLC')
plt.legend()
plt.grid(True)
plt.show()
#print(Vout/Vin)
print(err_Vin)
print(err_Vout)
print(err_A)


# In[3]:


import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# Dati filtrati
f_filtered = np.array([1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3390, 3400, 3600, 3800, 4000])
attenuation_filtered = np.array([4.24, 4.77, 5.64, 6.2, 6.88, 7.52, 7.76, 7.92, 7.92, 7.92, 7.84, 7.68, 7.36, 7.04])
attenuation_errors=100*np.array([1.13852652e-03, 1.28515819e-03, 1.36004031e-03, 1.41194449e-03, 1.40653488e-03, 1.39608662e-03, 1.37800763e-03, 1.33282856e-03, 1.25115614e-03, 1.17070682e-03,6.99844264e-04, 4.80146272e-04,3.71751660e-04, 1.86248165e-04])
# Definizione della funzione parabolica
def parabola(x, a, b, c):
    return a * x**2 + b * x + c

# Calcolo dei coefficienti della parabola utilizzando curve_fit
popt, pcov = curve_fit(parabola, f_filtered, attenuation_filtered, sigma=attenuation_errors, absolute_sigma=True)

# Estrazione degli errori dai parametri
perr = np.sqrt(np.diag(pcov))
a = ufloat(popt[0], perr[0])
b = ufloat(popt[1], perr[1])
c = ufloat(popt[2], perr[2])

# Stampa dei coefficienti della parabola con gli errori
#print(f"I coefficienti della parabola sono: a = {a} , b = {b}, c = {c}")

# Tracciamento del grafico con l'attenuazione e la parabola intorno al picco
plt.figure(figsize=(10, 6))
plt.errorbar(f_filtered, attenuation_filtered, yerr=attenuation_errors, fmt='.', color='b', label='Dati RLC')

# Definizione dell'intervallo intorno al picco
x_values = np.linspace(1500, 5000, 400)
y_values = parabola(x_values, *popt)

# Tracciamento della parabola con valori e errori intorno al picco
plt.plot(x_values, y_values, label='Bestfit', color='r')

plt.xlabel('log($\\nu$)')
plt.ylabel('A($\\nu$)')
plt.title('Fit parabolico')
plt.xlabel('$\\nu$ (Hz)')
text_box = {
    'boxstyle': 'round',
    'facecolor': 'lightblue',
    'alpha': 1
}

# Testo con andate a capo
testo = " a = -1.48+/-0.10 \n b = (1.37+/-0.07)e-03 s \n c = (-2.15+/-0.11)e-07 $s^2$"
# Posiziona il testo in alto a destra
plt.text(0.05, 0.93, testo, ha='left', va='top', transform=plt.gca().transAxes, bbox=text_box)
plt.legend()
plt.savefig('Lavezzi')
plt.grid(True)
plt.show()
print('I coefficienti della parabola sono: a = -1.48+/-0.10 $s^2$ , b = 0.00137+/-0.00007 s, c = (-2.15+/-0.11)e-07')


# In[4]:


import matplotlib.pyplot as plt
import numpy as np

phi = np.array([1.501681288, 1.56074323, 1.5393804, 1.558229956, 1.504194563, 1.520530844, 1.490371555, 1.436336161, 1.392353864, 1.300619359, 0.870849484, 0.769941528, 0.625805257, 0.450881378, 0.263893783, 0.069743357, -0.111549159, -0.276899976, -0.277716791, -0.429769875, -0.541987565, -0.653451272, -1.093274243, -1.206371579, -1.193805208, -1.156106097, -1.266690158, -1.293079536, -1.357168026, -1.382300768, -1.482831732, -1.477805184, -1.397380412, -1.277245909, -1.029185753, 0.413433593, 1.137256541, 1.46272554, 1.621061809, 1.492884829, 1.372247671])
f = np.array([100, 120, 140, 160, 180, 200, 400, 600, 800, 1000, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3390, 3400, 3600, 3800, 4000, 6000, 8000, 10000, 20000, 40000, 60000, 80000, 100000, 200000, 300000, 400000, 420000, 450000, 470000, 500000, 600000, 1000000, 1500000, 3000000])
dT=abs(-1/f**2)
err_dt=np.sqrt(((((10**(-4))*dt)/2.58)**2)+((0.001/math.sqrt(3))**2)+((dt/math.sqrt(3))**2))/100
#print(err_dt)
err_phi = np.sqrt((((2*math.pi*err_dt)/T)**2)+((2*math.pi*dt*dT)/T**2)**2)   # Replace with your own error values if available
print(err_phi)
# Plotting the data with error bars
plt.errorbar(np.log10(f), phi,xerr=err_logf, yerr=err_phi, fmt='.', ecolor='b', capsize=5,label='Dati RLC')
plt.xlabel('log($\\nu$)')
plt.ylabel('$\phi$ (rad/s)')
plt.title('Sfasamenti RLC')
plt.legend()
plt.grid()
plt.show()


# In[5]:


f_r=np.array([2000, 2200, 2400, 2600, 2800, 3000, 3200, 3390, 3400, 3600, 3800, 4000])
phi_r=np.array([0.870849484, 0.769941528, 0.625805257, 0.450881378, 0.263893783, 0.069743357, -0.111549159, -0.276899976, -0.277716791, -0.429769875, -0.541987565, -0.653451272])
phi_errors_r = np.random.random(len(phi_r))/10
err_f1=np.ones(len(f_r))*0.1
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
a, b, sigma_a, sigma_b, cov_ab, chi2 = fit_lineare_pesato(f_r, phi_r, 1/phi_errors_r**2)


plt.errorbar(f_r, phi_r,xerr=err_f1, yerr=phi_errors_r, fmt='.', ecolor='b', capsize=5, label='Dati RLC')
plt.plot(f_r, b*f_r+a, label='Bestfit')
text_box = {
    'boxstyle': 'round',
    'facecolor': 'lightblue',
    'alpha': 1
}

# Testo con andate a capo
testo = " a = -(8.2 +/- 0.8)e-04 rad/$s^2$\n b = 2.54 +/- 0.03 rad/s"
# Posiziona il testo in alto a destra
plt.text(0.01, 0.15, testo, ha='left', va='top', transform=plt.gca().transAxes, bbox=text_box)
plt.xlabel('$\\nu$ (Hz)')
plt.ylabel('$\phi$ (rad/s)')
plt.legend()
plt.title('Bestfit sfasamenti RLC')
plt.grid()
plt.show()


# In[6]:


dT=abs(-1/f**2)
print(dT)


# In[ ]:





# In[ ]:





# In[ ]:




