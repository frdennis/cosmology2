import numpy as np
import matplotlib.pyplot as plt 

data = np.loadtxt("cosmology.txt")

m_to_Mpc = 3.24077929e-23
H0 = 3e8/67 # km/s/Mpc
c = 2.99792458e8 # m/s
Hp_convert = 3.086e19

x           = data[:,0]
eta         = data[:,1]*m_to_Mpc
Hp          = data[:,2]*Hp_convert
dHp_dx      = data[:,3]
ddHpddx     = data[:,4]
OmegaB      = data[:,5]
OmegaCDM    = data[:,6]
OmegaLambda = data[:,7]
OmegaR      = data[:,8]
OmegaNu     = data[:,9]
OmegaK      = data[:,10]
H           = data[:,11]
t           = data[:,12]
r           = data[:,13]*m_to_Mpc


a = np.e**x

OmegaB0       = 0.05
OmegaCDM0     = 0.267
OmegaLambda0  = 0.682945
OmegaK0       = 0
OmegaR0       = 5.50896e-05
OmegaNu0      = 3.046*7.0/8.0*(4/11)**(4/3)*OmegaR0

a_mat_rad = (OmegaR0 + OmegaNu0)/(OmegaB0 + OmegaCDM0)
a_mat_lambda = ((OmegaCDM0+OmegaB0)/(OmegaLambda0))**(1/3)
a_acc = ((OmegaCDM0+OmegaB0)/(2*OmegaLambda0))**(1/3)
a_today = 1

z = 1/a - 1

dA = a*r
dL = dA/a**2

lumdata = np.loadtxt("supernovadata.txt")
z_L = lumdata[:,0]
d_L = lumdata[:,1]
error = lumdata[:,2]

""" Redshift data """
fig, ax = plt.subplots(figsize=(6,5))
ax.plot(z, dL/1000, label="Theoretical prediction")
ax.errorbar(z_L, d_L, yerr=error, label="Supernova data", color='k', fmt='o',)
ax.set_yscale("log"); ax.set_xscale("log")
ax.set_xlim(0.5e-2, 1e1); ax.set_ylim(1e-2, 1e2)
ax.legend()
ax.set_ylabel("Distance [Gpc]")
ax.set_xlabel("Redshift z")
ax.grid()
plt.savefig('figures/milestone1/'+'redshift_data.pdf')
plt.show()

""" Luminosity data """
fig, ax = plt.subplots(figsize=(6,5))
ax.plot(z, dL, label="Luminosity distance")
ax.plot(z, r, label="Comoving distance")
ax.plot(z, dA, label="Angular distance")
ax.plot(z, 4400*z, label=r"Naive Hubble Distance ($d = (c/H_0)z$)", color="k", ls="--")

#ax.vlines(1/a_today-1, H.min(), H.max()*1.01, color='k', ls='--', label='Today')
#ax.vlines(1/a_mat_rad-1, H.min(), H.max()*1.01, color='b', ls='--', label='Radiation-Matter eq')
#ax.vlines(1/a_mat_lambda-1, H.min(), H.max()*1.01, color='r', ls='--', label='Matter-Dark energy eq')
#ax.vlines(1/a_acc-1, H.min(), H.max()*1.01, color='g', ls='--', label='Accelerated universe')

ax.legend(); ax.grid()
ax.set_xlabel("Redshift z"); ax.set_ylabel("Distance [Mpc]")
ax.set_xscale("log"); ax.set_yscale("log")
plt.savefig('figures/milestone1/'+'luminosity_data.pdf')
plt.show()

""" a(t) """

def myindex(inarray, val):
    # returns value in inarray that is closest to val
    minarray = np.abs(inarray - val)
    return np.where(np.min(minarray) == minarray)[0][0]
#print("Age of the universe: %.3f Gyrs" % (t[0]/10**9/365/24/60/60))
s_to_Gyrs = 1/10**9/365/24/60/60
fig, ax = plt.subplots(figsize=(6,5))
ax.plot(t*s_to_Gyrs, a, label='a(t)')
ax.grid()
ax.hlines(1,t.min()*s_to_Gyrs, t.max()*s_to_Gyrs, color='k', ls='--')
minindx = np.where(abs(a-1) == np.min(abs(a-1)))[0][0]
print(minindx)
# vline at a(t) = 1
ax.vlines(t[myindex(a, a_today)]*s_to_Gyrs, a.min(), a.max()*1.01, color="k", ls="--", label='Today')
ax.vlines(t[myindex(a, a_mat_rad)]*s_to_Gyrs, a.min(), a.max()*1.01, color='b', ls='--', label='Radiation-Matter eq')
ax.vlines(t[myindex(a, a_mat_lambda)]*s_to_Gyrs, a.min(), a.max()*1.01, color='r', ls='--', label='Matter-Dark energy eq')
ax.vlines(t[myindex(a, a_acc)]*s_to_Gyrs, a.min(), a.max()*1.01, color='g', ls='--', label='Accelerated universe')
ax.set_yscale("log"); ax.set_xscale("log")
ax.set_xlabel("$t [Gyrs]$")
ax.set_ylabel("Scalefactor $a$")
ax.legend()
plt.savefig('figures/milestone1/'+'scalefactor.pdf')
plt.show()

""" Reduced Hubble parameter """
fig, ax = plt.subplots(figsize=(6,5))
ax.plot(a, H)
ax.grid()
ax.set_yscale("log"); ax.set_xscale("log")
ax.set_ylabel("$H/H_0$")
ax.set_xlabel("Scalefactor $a$")
ax.vlines(a_today, H.min(), H.max()*1.01, color='k', ls='--', label='Today')
ax.vlines(a_mat_rad, H.min(), H.max()*1.01, color='b', ls='--', label='Radiation-Matter eq')
ax.vlines(a_mat_lambda, H.min(), H.max()*1.01, color='r', ls='--', label='Matter-Dark energy eq')
ax.vlines(a_acc, H.min(), H.max()*1.01, color='g', ls='--', label='Accelerated universe')
ax.legend()
plt.savefig('figures/milestone1/'+'reduced_hubble.pdf')
plt.show()

""" Hp """
fig, ax = plt.subplots(figsize=(6,5))
ax.plot(a, Hp)
ax.grid()
ax.set_yscale("log"); ax.set_xscale("log")
ax.set_ylabel("$\mathcal{H}$ [km/s/Mpc]")
ax.set_xlabel("Scalefactor $a$")
ax.vlines(a_today, Hp.min(), Hp.max()*1.01, color='k', ls='--', label='Today')
ax.vlines(a_mat_rad, Hp.min(), Hp.max()*1.01, color='b', ls='--', label='Radiation-Matter eq')
ax.vlines(a_mat_lambda, Hp.min(), Hp.max()*1.01, color='r', ls='--', label='Matter-Dark energy eq')
ax.vlines(a_acc, Hp.min(), Hp.max()*1.01, color='g', ls='--', label='Accelerated universe')
ax.legend()
plt.savefig('figures/milestone1/'+'Hp.pdf')
plt.show()

""" Conformal time eta(x) """
fig, ax = plt.subplots(figsize=(6,5))
ax.plot(a, eta); ax.grid()
ax.set_ylabel(r'$\eta$ [Mpc]')
ax.set_xlabel(r'Scalefactor $a$')
ax.set_xscale("log"); ax.set_yscale("log")
ax.vlines(a_today, eta.min(), 1.01*eta.max(), color='k', ls='--', label='Today')
ax.vlines(a_mat_rad, eta.min(), 1.01*eta.max(), color='b', ls='--', label='Radiation-Matter eq')
ax.vlines(a_mat_lambda, eta.min(), 1.01*eta.max(), color='r', ls='--', label='Matter-Dark energy eq')
ax.vlines(a_acc, eta.min(), 1.01*eta.max(), color='g', ls='--', label='Accelerated universe')
ax.legend()
plt.savefig('figures/milestone1/'+'eta.pdf')
plt.show()

""" Density parameters """
fig, ax = plt.subplots(figsize=(10,5))
ax.plot(a, (OmegaR+OmegaNu), label="Relativistic")
ax.plot(a, (OmegaB+OmegaCDM), label="Matter")
ax.plot(a, OmegaLambda, label="Lambda")
# Vlines: 
mat_rad_tmp = abs((OmegaB+OmegaCDM) - (OmegaR + OmegaNu))
indx_mat_rad = np.where(mat_rad_tmp == np.min(mat_rad_tmp[:(len(mat_rad_tmp)//2)]))[0][0]
print(indx_mat_rad)
mat_darke_tmp = abs((OmegaB + OmegaCDM) - OmegaLambda)
indx_mat_darke = np.where(mat_darke_tmp == np.min(mat_darke_tmp))[0][0]

ax.vlines(a_today, OmegaLambda.min(), OmegaLambda.max()*1.01, color='k', ls='--', label='Today')
ax.vlines(a_mat_rad, OmegaLambda.min(), OmegaLambda.max()*1.01, color='b', ls='--', label='Radiation-Matter eq')
ax.vlines(a_mat_lambda, OmegaLambda.min(), OmegaLambda.max()*1.01, color='r', ls='--', label='Matter-Dark energy eq')
ax.vlines(a_acc, OmegaLambda.min(), OmegaLambda.max()*1.01, color='g', ls='--', label='Accelerated universe')
ax.legend(); ax.grid()
ax.set_xscale("log")
ax.set_ylabel(r'$\Omega_i$'); ax.set_xlabel('Scalefactor a')
plt.savefig('figures/milestone1/'+'density_parameters.pdf')
plt.show()

fig, ax = plt.subplots(figsize=(6,5))
ax.plot(a, (OmegaR+OmegaNu)+(OmegaB+OmegaCDM)+OmegaLambda - 1)
ax.set_ylabel(r'Devianse in $\Sigma \Omega$ from 1'); ax.set_xlabel('Scalefactor a')
ax.set_xscale('log'); ax.grid()
plt.savefig('figures/milestone1/'+'error.pdf')
plt.show()

""" eta H """
fig, ax = plt.subplots(figsize=(10,5))
ax.set_yscale("log")
y = data[:,1]*data[:,2]/c #eta*Hp/c/m_to_Mpc
ax.plot(a, y); ax.grid()
ax.set_xscale("log")
ax.set_xlabel('Scalefactor a'); ax.set_ylabel(r"$\frac{\eta(x)\mathcal{H}(x)}{c}$")
ax.vlines(a_today, y.min(), y.max()*1.01, color='k', ls='--', label='Today')
ax.vlines(a_mat_rad, y.min(), y.max()*1.01, color='b', ls='--', label='Radiation-Matter eq')
ax.vlines(a_mat_lambda, y.min(), y.max()*1.01, color='r', ls='--', label='Matter-Dark energy eq')
ax.vlines(a_acc, y.min(), y.max()*1.01, color='g', ls='--', label='Accelerated universe')
ax.legend()
plt.savefig('figures/milestone1/'+'eta_H.pdf')
plt.show()

""" Derivatives of H_p """
fig, ax = plt.subplots(figsize=(10,5))
ax.plot(a, dHp_dx/Hp, label="first derivative")
ax.plot(a, ddHpddx/Hp, label="second derivative")
ax.set_xscale('log')
ax.set_ylabel(r'$\frac{1}{\mathcal{H}}\frac{d\mathcal{H}}{dx}, \frac{1}{\mathcal{H}}\frac{d^2 \mathcal{H}}{dx^2}$')
ax.set_xlabel('Scalefactor $a$')
ax.vlines(a_today, (dHp_dx/Hp).min(), (ddHpddx/Hp).max()*1.01, color='k', ls='--', label='Today')
ax.vlines(a_mat_rad, (dHp_dx/Hp).min(), (ddHpddx/Hp).max()*1.01, color='b', ls='--', label='Radiation-Matter eq')
ax.vlines(a_mat_lambda, (dHp_dx/Hp).min(), (ddHpddx/Hp).max()*1.01, color='r', ls='--', label='Matter-Dark energy eq')
ax.vlines(a_acc, (dHp_dx/Hp).min(), (ddHpddx/Hp).max()*1.01, color='g', ls='--', label='Accelerated universe')
ax.legend()
ax.grid()
plt.savefig('figures/milestone1/'+"derivatives.pdf")
plt.show()