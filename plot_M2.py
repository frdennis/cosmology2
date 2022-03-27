import numpy as np
import matplotlib.pyplot as plt 

data = np.loadtxt("recombination.txt")
path = 'figures/milestone2/'

def saha(x): 
    OmegaB0      = 0.05
    pi = 3.1415926535

    km          = 1e3
    N           = 1
    J           = N 
    W           = 1
    Mpc         = 3.08567758e22
    eV          = 1.60217653e-19
    H0          = 0.67*100 * km/Mpc
  
    k_b         = 1.38064852e-23
    m_e         = 9.10938356e-31
    m_H         = 1.6735575e-27 
    c           = 2.99792458e8
    G           = 6.67430e-11 
    hbar        = 1.054571817e-34

    rho_c0 = 3*H0**2/(8*pi*G)
    dimension_Xe_saha = hbar**(-3)* k_b**(3./2)
    epsilon_0 = 13.605693122994 * eV
    TCMB0       = 2.7255

    a           = np.exp(x)
    n_H         = (OmegaB0*rho_c0/(m_H*a**3))
    A  = dimension_Xe_saha*np.exp(-(epsilon_0)/(k_b*TCMB0/a))/n_H * pow(m_e*(TCMB0/a)/(2*pi),3./2.)
    Xe = 2/(1+ np.sqrt(1 + 4/A))
    return Xe

x           = data[:,0]
Xe          = data[:,1]
ne          = data[:,2]
tau         = data[:,3]
dtaudx      = data[:,4]
ddtauddx    = data[:,5]
g_tilde     = data[:,6]
dgdx        = data[:,7]
ddgddx      = data[:,8]
Tb          = data[:,9]

a = np.exp(x)
def myindex(inarray, val):
    # returns value in inarray that is closest to val
    minarray = np.abs(inarray - val)
    return np.where(np.min(minarray) == minarray)[0][0]
i = myindex(tau, 1) # where tau is 1
print(i)
j = np.argmin(np.abs(Xe - 0.5)) # where Xe is 0.5
print(j)
print(x[j])

print("\n")
print(r"Time of decoupling (x=", x[i], ", z=", 1/np.exp(x[i]) - 1, ")")
print(r"Time of Xe = 0.5 (x =", x[j], ", z =", 1/np.exp(x[j]) - 1)
print("Saha at Xe = 0.5 (x =", x[myindex(saha(x), 0.5)], ", z =", 1/np.exp(x[myindex(saha(x),0.5)]) - 1)
print(r"Freeze out abundance at Xe(x=0) = %.3e" % Xe[myindex(a,1)])


#plt.plot(x, ne)
#plt.yscale("log")
#plt.show()

figsize = (7,5)

fig, ax = plt.subplots(figsize=figsize)
ax.plot(a, Xe)#, label=r'$X_e$')
#ax.vlines(a[i], a.max()*0.00005, a.max()*1.01, 'k', ls='--')
ax.plot(a[j], Xe[j], 'ko', label=r'Recombination (a = %.3e, $X_e$=%.3e)'%(a[j],Xe[j]))
ax.set_ylabel(r'$X_e$')
ax.set_xlabel(r'Scalefactor $a$')
ax.set_yscale("log")
ax.set_xscale("log")
ax.set_xlim(np.exp(-12),np.exp(0))
ax.grid()
plt.legend()
plt.savefig(path+'electron_fraction.pdf')
plt.show()

"""
fig, ax = plt.subplots(figsize=figsize)
ax.set_yscale("log")
ax.set_xscale("log")
ax.plot(a, ne)
plt.show()
"""

fig, ax = plt.subplots(figsize=figsize)
ax.set_ylabel(r"Optical depth ($\tau$, $\tau'$, $\tau''$)")
ax.set_xlabel(r'Scalefactor $a$')
ax.set_yscale("log")
ax.set_xscale("log")
ax.set_ylim(1e-8, 1e8)
ax.set_xlim(np.exp(-12),np.exp(0))
ax.plot(a, tau, label=r"$\tau$")
ax.plot(a, -dtaudx, label=r"-$\tau'$")
ax.plot(a, ddtauddx, label=r"$\tau''$")
ax.vlines(a[i], tau.min(), tau.max(), 'k', ls='--')
ax.plot(a[i], tau[i], 'ko', label=r'Decoupling (a = %.3e, $\tau$=%.3e)'%(a[i],tau[i]))
ax.legend()
ax.grid()
plt.savefig(path+'optical_depth.pdf')
plt.show()

fig, ax = plt.subplots(figsize=figsize, nrows=3, sharex=True)
ax[0].set_ylabel(r"Visibility function $g$")
ax[1].set_ylabel(r"$g'$")
ax[2].set_ylabel(r"$g''$")
ax[2].set_xlabel(r'Scalefactor $a$')
#ax.set_yscale("log")
g =  g_tilde#/np.sum(np.abs(g_tilde))
ax[0].plot(a, g)
ax[0].plot(a[i], g[i], 'ko')
ax[0].vlines(a[i], g.min(), g.max(), 'k', ls='--', label='Decoupling')
ax[0].legend()
dgdx = dgdx#/np.sum(np.abs(dgdx))
ax[1].plot(a, dgdx)
ax[1].plot(a[i], dgdx[i], 'ko')
ax[1].vlines(a[i], dgdx.min(), dgdx.max(), 'k', ls='--')
ddgddx = ddgddx#/np.sum(np.abs(ddgddx))
ax[2].plot(a, ddgddx)
ax[2].plot(a[i], ddgddx[i], 'ko')
ax[2].vlines(a[i], ddgddx.min(), ddgddx.max(), 'k', ls='--')
[[axi.set_xlim(4e-4, 2e-3), axi.grid(), axi.set_xscale("log")] for axi in ax]
plt.tight_layout()
plt.savefig(path+'visibility_function.pdf')
plt.show()

def T(a):
    return 2.7255/a

fig, ax = plt.subplots(figsize=figsize)
ax.plot(a, Tb, label=r'$T_b$')
ax.plot(a, T(a), label=r'$T$', ls='--')
ax.vlines(a[i], Tb.min(), Tb.max(), 'k', ls='--', label='Decoupling')
ax.set_ylabel(r"Baryon Temperature [K]")
ax.set_xlabel(r'Scalefactor $a$')
ax.set_yscale("log")
ax.set_xscale("log")
#ax.set_xlim(np.exp(-12),np.exp(0))
ax.grid()
ax.legend()
plt.savefig(path+'baryon_temperature.pdf')
plt.show()

fig, ax = plt.subplots(figsize=figsize)
ax.plot(a, Xe, label=r'$X_e$')
ax.plot(a, saha(x), label=r'$X_e$ from Saha only', ls='--')
ax.hlines(Xe[j], a.min(), a.max(), 'k', ls='--', label=r'Recombination ($X_e = 0.5$)')
ax.grid()
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlim(1e-4,1e-2); ax.set_ylim(1e-5, 0.5e1)
ax.set_xlabel(r'Scalefactor $a$'); ax.set_ylabel(r'Electron fraction $X_e$')
plt.legend()
plt.savefig(path+"saha.pdf")
plt.show()