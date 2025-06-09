import numpy as np
import pandas  as pd

i = 0
samplen = 0
Etot = np.loadtxt("summary.ETOT")
Epot = np.loadtxt("summary.EPTOT")
volume = np.loadtxt("summary.VOLUME")

Etot = np.array(Etot)[samplen:,1] #in kcal/mol
Epot = np.array(Epot)[samplen:,1] #in kcal/mol
volume = np.array(volume)[samplen:,1] #in A^3
    
Pext = 1.013 #in bar
Temperature = 298.15 #in K
H = Etot*4.184+Pext*volume*6.02e-5 #in kJ/mol

kB = 0.0083144621 #in kJ/mol/K
nmol = 216 # number of water moleculars
M = nmol * (15.9994 + 1.008 * 2)

Rho = M/np.average(volume)/6.02e23*1e24
print(str(Rho) + 'g/cm^3')

var_Ek = nmol * 9 * kB**2 * Temperature**2 / 2.0
cp = ( var_Ek + np.var(Epot*4.184)+np.var(Pext*volume*6.02e-5))/kB/np.average(Temperature)**2/4.184*1000/nmol #in cal/mol/K
print(str(cp)+'cal/mol/K')

Kt=np.var(volume)/np.average(volume)/(kB*np.average(Temperature))*6.02e-10*101325*1e6 #in atm-1
print(str(Kt)+'*10^-6 atm^-1')

alpha=(np.average(H*volume)-np.average(volume)*np.average(H))/kB/np.average(Temperature)**2/np.average(volume)
print(str(alpha)+'K^-1')

f = open('Rho.txt','w')
f.write(str(Rho))
f.close()
f = open('Ep.txt','w')
f.write(str(np.average(Epot)))
f.close()
f = open('H.txt','w')
f.write(str(np.average(H)/4.184))
f.close()
f = open('Cp.txt','w')
f.write(str(cp))
f.close()
f = open('Kt.txt','w')
f.write(str(Kt))
f.close()
f = open('alpha.txt','w')
f.write(str(alpha))
f.close()
