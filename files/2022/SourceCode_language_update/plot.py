import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

# plot figure 1

mpl.rcParams['font.family'] = 'Arial'
F = 3

def plot_one(fname, i, N1, N2, n, fmt):
    a = pd.read_csv(fname, header=None, sep='\s+').values
    t = a[:,0]
    theta = a[:,1:F]
    varphi = a[:,F:2*F-1]
    gp = a[:,2*F-1]
    c_SW = a[:,2*F:3*F] + 1j*a[:,3*F:4*F]
    c_CR = a[:,5*F:6*F] + 1j*a[:,6*F:7*F]
    c_EX = a[:,8*F:9*F] + 1j*a[:,9*F:10*F]

    # plt.subplot(N1, N2, n)
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.plot(t, abs(c_SW[:,i])**2, fmt, label=r'SW')
    plt.plot(t[::10], abs(c_EX[::10,i])**2, 'k.', label=r'Exact')
    plt.xlim(0,t[-1])

plt.figure(figsize=(4,3), dpi=300)
plt.tick_params(axis='both', which='major', labelsize=10)
plt.ticklabel_format(style='sci', scilimits=(-1,2), axis='y')
plot_one("weak-SW-L.out", 0, 1, 1, 1, 'r-')
plot_one("weak-SW-Huo1-H1.out", 0, 1, 1, 1, 'b-')
plot_one("weak-SW-Huo2-H2.out", 0, 1, 1, 1, 'g--')
plt.ylim(0.9996,1)
plt.xticks([0,5,10], fontsize=10)
plt.yticks([1.0], rotation=0, fontsize=10)
plt.xlabel('Time (au)', fontsize=12)
plt.ylabel(r'$|c_1|^2$', fontsize=12, position=(1000,0.5))
plt.savefig('fig1a.png',bbox_inches = 'tight',dpi=300)

plt.figure(figsize=(4,3), dpi=300)
plt.tick_params(axis='both', which='major', labelsize=10)
plot_one("strong-SW-L.out", 0, 1, 1, 1, 'r-')
plot_one("strong-SW-Huo1-H1.out", 0, 1, 1, 1, 'b-')
plot_one("strong-SW-Huo2-H2.out", 0, 1, 1, 1, 'g--')
plt.ylim(0,0.8)
plt.yticks([0,0.4,0.8], fontsize=10)
plt.xticks([0,0.25,0.5], fontsize=10)
plt.xlabel('Time (au)', fontsize=12)
plt.ylabel(r'$|c_1|^2$', fontsize=12, position=(1000,0.5))
plt.savefig('fig1b.png',bbox_inches = 'tight',dpi=300)

# plot figure 2

def plot_two(fname, i, N1, N2, n, fmt):
    a = pd.read_csv(fname, header=None, sep='\s+').values
    t = a[:,0]
    theta = a[:,1:F]
    varphi = a[:,F:2*F-1]
    gp = a[:,2*F-1]
    c_SW = a[:,2*F:3*F] + 1j*a[:,3*F:4*F]
    c_CR = a[:,5*F:6*F] + 1j*a[:,6*F:7*F]
    c_EX = a[:,8*F:9*F] + 1j*a[:,9*F:10*F]

    # plt.subplot(N1, N2, n)
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.plot(t, np.real(c_SW[:,i]), fmt, label=r'SW')
    if(fname=='weak-SW-L.out' or fname=='strong-SW-L.out' ):
        plt.plot(t, np.real(c_CR[:,i]), 'violet', label=r'SW+gp')
    if(fname[0:4]=='weak'):
        eve = 2
    if(fname[0:4]=='stro'):
        eve = 8
    plt.plot(t[::eve], np.real(c_EX[::eve,i]), 'k.', ms=4, label=r'Exact')
    plt.xlim(0,t[-1])

plt.figure(figsize=(4,3), dpi=300)
plt.tick_params(axis='both', which='major', labelsize=10)
plt.ticklabel_format(style='sci', scilimits=(-1,2), axis='y')
plot_two("weak-SW-L.out", 0, 1, 1, 1, 'r-')
plt.xticks([0,5,10], fontsize=10)
plt.yticks(fontsize=10)
plt.xlabel('Time (au)', fontsize=12)
plt.ylabel(r'$\Re c_1$', fontsize=12, position=(1000,0.5))
plt.xlim(0,5)
plt.savefig('fig2a.png',bbox_inches = 'tight',dpi=300)

plt.figure(figsize=(4,3), dpi=300)
plt.tick_params(axis='both', which='major', labelsize=10)
plot_two("strong-SW-L.out", 0, 1, 1, 1, 'r-')
plt.yticks(fontsize=10)
plt.xticks([0,0.25,0.5], fontsize=10)
plt.xlabel('Time (au)', fontsize=12)
plt.ylabel(r'$\Re c_1$', fontsize=12, position=(1000,0.5))
plt.ylim(-0.8,1)
plt.xlim(-0.00,0.5)
plt.savefig('fig2b.png',bbox_inches = 'tight',dpi=300)

plt.figure(figsize=(4,3), dpi=300)
plt.tick_params(axis='both', which='major', labelsize=10)
plt.ticklabel_format(style='sci', scilimits=(-1,2), axis='y')
plot_two("weak-SW-Huo1-H1.out", 0, 1, 1, 1, 'b-')
plot_two("weak-SW-Huo2-H2.out", 0, 1, 1, 1, 'g-')
plt.xticks([0,5,10], fontsize=10)
plt.yticks(fontsize=10)
plt.xlabel('Time (au)', fontsize=12)
plt.ylabel(r'$\Re c_1$', fontsize=12, position=(1000,0.5))
plt.xlim(0,5)
plt.savefig('fig2c.png',bbox_inches = 'tight',dpi=300)

plt.figure(figsize=(4,3), dpi=300)
plt.tick_params(axis='both', which='major', labelsize=10)
plot_two("strong-SW-Huo1-H1.out", 0, 1, 1, 1, 'b-')
plot_two("strong-SW-Huo2-H2.out", 0, 1, 1, 1, 'g-')
plt.yticks(fontsize=10)
plt.xticks([0,0.25,0.5], fontsize=10)
plt.xlabel('Time (au)', fontsize=12)
plt.ylabel(r'$\Re c_1$', fontsize=12, position=(1000,0.5))
plt.ylim(-0.8,1)
plt.xlim(-0.00,0.5)
plt.savefig('fig2d.png',bbox_inches = 'tight',dpi=300)

# plot figure 3

F=2

plt.figure(figsize=(5,3), dpi=300)
plt.tick_params(axis='both', which='major', labelsize=10)
a = pd.read_csv("Berry-SW-L.out", header=None, sep='\s+').values
t = a[:,0]
theta = a[:,1:F]
varphi = a[:,F:2*F-1]
gp = a[:,2*F-1]
c_SW = a[:,2*F:3*F] + 1j*a[:,3*F:4*F]
c_CR = a[:,5*F:6*F] + 1j*a[:,6*F:7*F]
c_EX = a[:,8*F:9*F] + 1j*a[:,9*F:10*F]
alpha = t[:]/t[-1]*(np.pi*4)

plt.plot(alpha/np.pi, np.real(c_EX[:,1]), 'k-', label=r'Meyer-Miller $\Re b$')
plt.plot(alpha/np.pi, np.real(c_SW[:,1]), 'r-', label=r'SW $\Re b$')
plt.plot(alpha/np.pi, np.real(c_CR[:,1]), 'c--', label=r'SW+gp $\Re b$')
plt.plot(alpha/np.pi, np.imag(c_EX[:,1]), 'b-', label=r'Meyer-Miller $\Im b$')
plt.plot(alpha/np.pi, np.imag(c_SW[:,1]), 'r--', label=r'SW $\Im b$')
plt.plot(alpha/np.pi, np.imag(c_CR[:,1]), 'violet', ls='--', label=r'SW+gp $\Im b$')

plt.xlim(0,2)
plt.xlabel(r'$\Phi(t) / \pi$')
plt.ylabel('Value')
plt.legend(fontsize=7, ncol=2, loc=1, bbox_to_anchor=(1,0.94))
plt.savefig('fig3.png',bbox_inches = 'tight',dpi=300)
