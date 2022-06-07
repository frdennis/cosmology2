from logging import NullHandler
from xml.etree.ElementTree import PI

from pyrsistent import v
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("perturbations_k0.01.txt")
data2 = np.loadtxt("perturbations_k0.1.txt")
data3 = np.loadtxt("perturbations_k0.001.txt")
data4 = np.loadtxt("perturbations_k0.3.txt")

fig, ax = plt.subplots()
ax.plot(data[:,0], data[:,-1])
ax.plot(data2[:,0], data2[:,-1])
ax.plot(data3[:,0], data3[:,-1])
ax.plot(data4[:,0], data4[:,-1])
ax.set_xlim(-8,0); ax.set_ylim(-1,3)
plt.show()

#quit()

path = "figures/milestone3/"

x           = data[:,0]
a           = np.exp(x)

OmegaB0       = 0.05
OmegaCDM0     = 0.267
OmegaLambda0  = 0.682945
OmegaK0       = 0
OmegaNu0      = 0
OmegaR0       = 5.50896e-05

a_mat_rad = (OmegaR0 + OmegaNu0)/(OmegaB0 + OmegaCDM0)
a_mat_lambda = ((OmegaCDM0+OmegaB0)/(OmegaLambda0))**(1/3)
a_acc = ((OmegaCDM0+OmegaB0)/(2*OmegaLambda0))**(1/3)
a_rec = (1.0 / (1.0 + 2000.0))
a_dec = (1.0 / (1.0 + 1100.0))
a_today = 1
print("Value of x at recombination: ", np.log(a_rec))

""" Density contrast """
deltacdm    = data[:,1]
deltacdm2   = data2[:,1]
deltacdm3   = data3[:,1]

deltab      = data[:,2]
deltab2     = data2[:,2]
deltab3      = data3[:,2]

Theta0      = data[:,5]
Theta0_2    = data2[:,5]
Theta0_3    = data3[:,5]

fig, ax = plt.subplots(figsize=(6,5)) # (6,5) 
ax.plot(a, deltacdm2, label='k = 0.1 h/Mpc', color='black')
ax.plot(a, abs(deltab2), ls='--', color='black')
ax.plot(a, 3*abs(Theta0_2), ls=':', color='black')
#ax.plot(a, 3*abs(data2[:,11]), ls='-.', color='darkred')

ax.plot(a, deltacdm, label='k = 0.01 h/Mpc', color='blue')
ax.plot(a, abs(deltab), ls='--', color='blue')
ax.plot(a, 3*abs(Theta0), ls=':', color='blue')
#ax.plot(a, 3*abs(data[:,11]), ls='-.', color='lime')

ax.plot(a, deltacdm3, label='k = 0.001 h/Mpc', color='red')
ax.plot(a, abs(deltab3), ls='--', color='red')
ax.plot(a, 3*abs(Theta0_3), ls=':', color='red')
#ax.plot(a, 3*abs(data3[:,11]), ls='-.', color='dimgray')
# Vlines 
#ax.vlines(a_rec, 1e-3, 1e4, color="k", ls="--", label='Recombination')
ax.vlines(a_dec, 1e-4, 1e5, color="gray", ls="--", label='Photon Decoupling')
#ax.vlines(a_mat_rad, 1e-3, 1e4, color="k", ls="--", label='Matter-Radiation eq.')
ax.vlines(a_mat_lambda, 1e-4, 1e5, color="pink", ls="--", label='Matter-DE eq.')

ax.set_xlabel(r"Scalefactor $a$")
ax.set_ylabel(r"$\delta_{CDM}$ (solid), $|\delta_b|$ (dashed), $|\delta_\gamma|$ (dotted)")
ax.set_yscale("log"); ax.set_xscale("log")
ax.set_xlim(1e-6, 1e1); ax.set_ylim(3e-4, 3e4)
ax.grid(); ax.legend()
plt.savefig(path+"densities.pdf")
plt.show()

""" Neutrino, photon and baryon density contrasts"""

Theta0      = 3*abs(data[:,5])
Theta0_2    = 3*abs(data2[:,5])
Theta0_3    = 3*abs(data3[:,5])

Nu0     = 3*abs(data[:,11])
Nu0_2   = 3*abs(data2[:,11])
Nu0_3   = 3*abs(data3[:,11])

fig, ax = plt.subplots(figsize=(6,5)) # (6,5)
ax.plot(a, abs(deltab2), label='k = 0.1 h/Mpc', color='black')
ax.plot(a, Nu0_2, ls='--', color='black')

ax.plot(a, abs(deltab), label='k = 0.01 h/Mpc', color='blue')
#ax.plot(a, Theta0, ls='--', color='blue')
ax.plot(a, Nu0, ls='--', color='blue')

ax.plot(a, abs(deltab3), label='k = 0.001 h/Mpc', color='red')
#ax.plot(a, Theta0_3, ls='--', color='red')
ax.plot(a, Nu0_3, ls='--', color='red')

ax.vlines(a_dec, 1e-4, 1e5, color="gray", ls="--", label='Photon Decoupling')
#ax.vlines(a_mat_rad, 1e-3, 1e4, color="k", ls="--", label='Matter-Radiation eq.')
ax.vlines(a_mat_lambda, 1e-4, 1e5, color="pink", ls="--", label='Matter-DE eq.')

ax.set_xlabel(r"Scalefactor $a$")
ax.set_ylabel(r"$\delta_{b}$ (solid), $|\delta_\nu|$ (dashed)")
ax.set_xlim(1e-6, 1e1); ax.set_ylim(5e-5, 4e4)
ax.set_yscale("log"); ax.set_xscale("log")
ax.grid(); ax.legend()
plt.savefig(path+"densities_tc.pdf")
plt.show()

""" Velocities """
v_cdm   = data[:,3]
v_cdm2  = data2[:,3]
v_cdm3  = data3[:,3]

v_b     = data[:,4]
v_b2    = data2[:,4]
v_b3    = data3[:,4]

v_photons  = - 3*data[:,6] #3*abs(data[:,6])
v_photons2 = - 3*data2[:,6] #3*abs(data2[:,6])
v_photons3 = - 3*data3[:,6] #3*abs(data3[:,6])

v_nu     = - 3*data[:,12]
v_nu2    = - 3*data2[:,12]
v_nu3    = - 3*data3[:,12]

fig, ax = plt.subplots(figsize=(6,5))
ax.plot(a, v_cdm2, label='k = 0.1 h/Mpc', color='black')
ax.plot(a, abs(v_b2), ls='--', color='black')
ax.plot(a, abs(v_photons2), ls=':', color='black')
#ax.plot(a, abs(v_nu2), ls='-.', color='darkred')

ax.plot(a, v_cdm, label='k = 0.01 h/Mpc', color='blue')
ax.plot(a, abs(v_b), ls='--', color='blue')
ax.plot(a, abs(v_photons), ls=':', color='blue')
#ax.plot(a, abs(v_nu), ls='-.', color='lime')

ax.plot(a, v_cdm3, label='k = 0.001 h/Mpc', color='red')
ax.plot(a, abs(v_b3), ls='--', color='red')
ax.plot(a, abs(v_photons3), ls=':', color='red')
#ax.plot(a, abs(v_nu3), ls='-.', color='dimgray')

ax.vlines(a_dec, 1e-7, 1e5, color="gray", ls="--", label='Photon Decoupling')
ax.vlines(a_mat_rad, 1e-7, 1e5, color="orange", ls="--", label='Matter-Radiation eq.')
ax.vlines(a_mat_lambda, 1e-7, 1e5, color="pink", ls="--", label='Matter-DE eq.')

ax.set_ylabel(r"$v_{CDM}$ (solid), $v_b$ (dashed), $v_\gamma$ (dotted)")
ax.set_xlabel(r"Scalefactor $a$")
ax.set_xlim(1e-7, 1e1); ax.set_ylim(1e-7, 3e1)
ax.set_yscale("log"); ax.set_xscale("log")
ax.grid(); ax.legend()
plt.savefig(path+"velocities.pdf")
plt.show()

""" v photons and neutrinos """
fig, ax = plt.subplots(figsize=(6,5))
ax.plot(a, abs(v_b2), label='k = 0.1 h/Mpc', color='black')
#ax.plot(a, abs(v_photons2), ls='--', color='black')
ax.plot(a, abs(v_nu2), ls='--', color='black')

ax.plot(a, abs(v_b), label='k = 0.01 h/Mpc', color='blue')
#ax.plot(a, abs(v_photons), ls='--', color='blue')
ax.plot(a, abs(v_nu), ls='--', color='blue')
#ax.plot(a, abs(v_nu), ls='-.', color='lime')

ax.plot(a, abs(v_b3), label='k = 0.001 h/Mpc', color='red')
#ax.plot(a, abs(v_photons3), ls='--', color='red')
ax.plot(a, abs(v_nu3), ls='--', color='red')

ax.vlines(a_dec, 1e-7, 1e5, color="gray", ls="--", label='Photon Decoupling')
ax.vlines(a_mat_rad, 1e-7, 1e5, color="orange", ls="--", label='Matter-Radiation eq.')
ax.vlines(a_mat_lambda, 1e-7, 1e5, color="pink", ls="--", label='Matter-DE eq.')

ax.set_ylabel(r"$v_{b}$ (solid), $v_\nu$ (dashed)")
ax.set_xlabel(r"Scalefactor $a$")
ax.set_xlim(1e-7, 1e1); ax.set_ylim(1e-7, 3e1)
ax.set_yscale("log"); ax.set_xscale("log")
ax.grid(); ax.legend()
plt.savefig(path+"velocities_relativistic.pdf")
plt.show()

""" Anisotropisc stress """
Psi         = data[:,15]
Psi2        = data2[:,15]
Psi3        = data3[:,15]

Phi         = data[:,14]
Phi2        = data2[:,14]
Phi3        = data3[:,14]

fig, ax = plt.subplots(figsize=(6,5))
ax.plot(a, Psi2 + Phi2, label='k = 0.1 h/Mpc', color='red')
ax.plot(a, Psi + Phi, label='k = 0.01 h/Mpc', color='g')
ax.plot(a, Psi3 + Phi3, label='k = 0.001 h/Mpc', color='k')

ax.vlines(a_dec, -1e-1, 1e1, color="gray", ls="--", label='Photon Decoupling')
ax.vlines(a_mat_rad, -1e-1, 1e1, color="orange", ls="--", label='Matter-Radiation eq.')
ax.vlines(a_mat_lambda, -1e-1, 1e1, color="pink", ls="--", label='Matter-DE eq.')

ax.set_xscale("log"); 
ax.set_ylim(-2e-2, 1e-1); ax.set_xlim(1e-6, 1e1)
ax.grid(); ax.legend()
ax.set_xlabel(r"Scalefactor $a$"); 
ax.set_ylabel(r"Anisotropic stress $(\Psi + \Phi)$")
plt.savefig(path+'stress.pdf')
plt.show()

Thetap0     = data[:,8]
Thetap0_2   = data2[:,8]
Thetap0_3   = data3[:,8]

Thetap1     = data[:,9]
Thetap1_2   = data2[:,9]
Thetap1_3   = data3[:,9]

Thetap2     = data[:,10]
Thetap2_2   = data2[:,10]
Thetap2_3   = data3[:,10]

fig, ax = plt.subplots(figsize=(12,4), ncols=3)
ax[0].plot(a, Thetap0, label='k = 0.01 h/Mpc', color='g')
ax[0].plot(a, Thetap0_2, label='k = 0.1 h/Mpc', color='red')
ax[0].plot(a, Thetap0_3, label='k = 0.001 h/Mpc', color='k')
ax[0].set_ylabel(r'$\Theta^P_0$')
ax[0].set_ylim(-0.02, 0.035)

ax[1].plot(a, Thetap1, label='k = 0.01 h/Mpc', color='g')
ax[1].plot(a, Thetap1_2, label='k = 0.1 h/Mpc', color='red')
ax[1].plot(a, Thetap1_3, label='k = 0.001 h/Mpc', color='k')
ax[1].set_ylabel(r'$\Theta^P_1$')
ax[1].set_xlabel(r'Scalefactor $a$')
ax[1].set_ylim(-0.007, 0.015)

ax[2].plot(a, Thetap2, label='k = 0.01 h/Mpc', color='g')
ax[2].plot(a, Thetap2_2, label='k = 0.1 h/Mpc', color='red')
ax[2].plot(a, Thetap2_3, label='k = 0.001 h/Mpc', color='k')
ax[2].set_ylabel(r'$\Theta^P_2$')
ax[2].set_ylim(-0.007, 0.015)

plt.tight_layout()
[[axi.grid(), axi.set_xscale("log"), axi.set_xlim(5e-5, 2e-1), axi.vlines(a_dec, -1e1, 1e1, color="gray", ls="--", label='Photon Decoupling'), 
axi.vlines(a_mat_rad, -1e1, 1e1, color="orange", ls="--", label='Matter-Radiation eq.'), 
axi.vlines(a_mat_lambda, -1e1, 1e1, color="pink", ls="--", label='Matter-DE eq.')] for axi in ax]
ax[0].legend()
plt.savefig(path+'polarization.pdf')
plt.show()

""" Curvature """

Phi         = data[:,14]
Phi2        = data2[:,14]
Phi3        = data3[:,14]

fig, ax = plt.subplots(figsize=(6,5))
ax.plot(a, Phi, label='k = 0.01/Mpc', color='g')
ax.plot(a, Phi2, label='k = 0.1/Mpc', color='red')
ax.plot(a, Phi3, label='k = 0.001/Mpc', color='k')

ax.vlines(a_dec, -1e-1, 1e1, color="gray", ls="--", label='Photon Decoupling')
ax.vlines(a_mat_rad, -1e-1, 1e1, color="orange", ls="--", label='Matter-Radiation eq.')
ax.vlines(a_mat_lambda, -1e-1, 1e1, color="pink", ls="--", label='Matter-DE eq.')

ax.set_xscale("log");
ax.set_xlim(1e-6, 1e1); ax.set_ylim(-0.01, 0.8)
ax.set_ylabel(r"$\Phi$"); ax.set_xlabel(r"Scalefactor $a$")
ax.grid(); ax.legend()
ax.set_title(r'Spatial curvature perturbation $\Phi$')
plt.savefig(path+"curvature.pdf")
plt.show()

Psi         = data[:,15]
Psi2        = data2[:,15]
Psi3        = data3[:,15]

fig, ax = plt.subplots(figsize=(6,5))
ax.plot(a, Psi, label='k = 0.01/Mpc', color='g')
ax.plot(a, Psi2, label='k = 0.1/Mpc', color='red')
ax.plot(a, Psi3, label='k = 0.001/Mpc', color='k')

ax.vlines(a_dec, -1e1, 1e1, color="gray", ls="--", label='Photon Decoupling')
ax.vlines(a_mat_rad, -1e1, 1e1, color="orange", ls="--", label='Matter-Radiation eq.')
ax.vlines(a_mat_lambda, -1e1, 1e1, color="pink", ls="--", label='Matter-DE eq.')

ax.set_xscale("log");
ax.set_xlim(1e-6, 1e1); ax.set_ylim(-0.7, 0.1)
ax.set_ylabel(r"$\Psi$"); ax.set_xlabel(r"Scalefactor $a$")
ax.grid(); ax.legend()
ax.set_title(r'Newtonian potential perturbation $\Psi$')
plt.savefig(path+"potential.pdf")
plt.show()
quit()

Theta0      = data[:,5]
Theta0_2    = data2[:,5]
Theta0_3    = data3[:,5]

Psi         = data[:,15]
Psi2        = data2[:,15]
Psi3        = data3[:,15]

fig, ax = plt.subplots(figsize=(6,5))
ax.plot(x, Theta0 + Psi, label='k = 0.01/Mpc', color='g')
ax.plot(x, Theta0_2 + Psi2, label='k = 0.1/Mpc', color='red')
ax.plot(x, Theta0_3 + Psi3, label='k = 0.001/Mpc', color='k')
ax.grid()
ax.set_title(r'$\Theta_0$')
plt.show()


Theta1      = data[:,6]
Theta1_2    = data2[:,6]
Theta1_3    = data3[:,6]

fig, ax = plt.subplots(figsize=(6,5))
ax.plot(x, Theta1, label='k = 0.01/Mpc', color='g')
ax.plot(x, Theta1_2, label='k = 0.1/Mpc', color='red')
ax.plot(x, Theta1_3, label='k = 0.001/Mpc', color='k')
ax.grid()
ax.set_title(r'$\Theta_1$')
plt.show()



Theta2      = data[:,7]
Theta2_2    = data2[:,7]
Theta2_3    = data3[:,7]

fig, ax = plt.subplots(figsize=(6,5))
ax.plot(x, Theta2, label='k = 0.01/Mpc', color='g')
ax.plot(x, Theta2_2, label='k = 0.1/Mpc', color='red')
ax.plot(x, Theta2_3, label='k = 0.001/Mpc', color='k')
ax.grid()
ax.set_title(r'$\Theta_2$')
plt.show()
quit()
"""
Thetap0     = data[:,8]
Thetap0_2   = data2[:,8]
Thetap0_3   = data3[:,8]

Thetap1     = data[:,9]
Thetap1_2   = data2[:,9]
Thetap1_3   = data3[:,9]

Thetap2     = data[:,10]
Thetap2_2   = data2[:,10]
Thetap2_3   = data3[:,10]

fig, ax = plt.subplots(figsize=(10,5), ncols=3)
ax[0].plot(x, Thetap0, label='k = 0.01/Mpc', color='g')
ax[0].plot(x, Thetap0_2, label='k = 0.1/Mpc', color='red')
ax[0].plot(x, Thetap0_3, label='k = 0.001/Mpc', color='k')
ax[0].set_title(r'$\Theta^P_0$')

ax[1].plot(x, Thetap1, label='k = 0.01/Mpc', color='g')
ax[1].plot(x, Thetap1_2, label='k = 0.1/Mpc', color='red')
ax[1].plot(x, Thetap1_3, label='k = 0.001/Mpc', color='k')
ax[1].set_title(r'$\Theta^P_1$')

ax[2].plot(x, Thetap2, label='k = 0.01/Mpc', color='g')
ax[2].plot(x, Thetap2_2, label='k = 0.1/Mpc', color='red')
ax[2].plot(x, Thetap2_3, label='k = 0.001/Mpc', color='k')
ax[2].set_title(r'$\Theta^P_2$')
[[axi.grid(), axi.set_yscale("log")] for axi in ax]
plt.show()

Nu0     = data[:,11]
Nu0_2   = data2[:,11]
Nu0_3   = data3[:,11]

Nu1     = data[:,12]
Nu1_2   = data2[:,12]
Nu1_3   = data3[:,12]

Nu2     = data[:,13]
Nu2_2   = data2[:,13]
Nu2_3   = data3[:,13]

fig, ax = plt.subplots(figsize=(10,5), ncols=3)
ax[0].plot(x, Nu0, label='k = 0.01/Mpc', color='g')
ax[0].plot(x, Nu0_2, label='k = 0.1/Mpc', color='red')
ax[0].plot(x, Nu0_3, label='k = 0.001/Mpc', color='k')
ax[0].set_title(r'$\mathcal{N}_0$')

ax[1].plot(x, Nu1, label='k = 0.01/Mpc', color='g')
ax[1].plot(x, Nu1_2, label='k = 0.1/Mpc', color='red')
ax[1].plot(x, Nu1_3, label='k = 0.001/Mpc', color='k')
ax[1].set_title(r'$\mathcal{N}_1$')

ax[2].plot(x, Nu2, label='k = 0.01/Mpc', color='g')
ax[2].plot(x, Nu2_2, label='k = 0.1/Mpc', color='red')
ax[2].plot(x, Nu2_3, label='k = 0.001/Mpc', color='k')
ax[2].set_title(r'$\mathcal{N}_2$')
[[axi.grid(), axi.set_yscale("log")] for axi in ax]
plt.show()
"""

Phi         = data[:,14]
Phi2        = data2[:,14]
Phi3        = data3[:,14]

fig, ax = plt.subplots(figsize=(6,5))
ax.plot(x, Phi, label='k = 0.01/Mpc', color='g')
ax.plot(x, Phi2, label='k = 0.1/Mpc', color='red')
ax.plot(x, Phi3, label='k = 0.001/Mpc', color='k')
ax.grid()
ax.set_title(r'$\Phi$')
plt.show()

Psi         = data[:,15]
Psi2        = data2[:,15]
Psi3        = data3[:,15]

fig, ax = plt.subplots(figsize=(6,5))
ax.plot(x, Psi, label='k = 0.01/Mpc', color='g')
ax.plot(x, Psi2, label='k = 0.1/Mpc', color='red')
ax.plot(x, Psi3, label='k = 0.001/Mpc', color='k')
ax.grid()
ax.set_title(r'$\Psi$')
plt.show()