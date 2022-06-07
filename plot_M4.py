from re import A
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp

data = np.loadtxt("cells.txt")
path = "figures/milestone4/"

def plot_power_spectrum(filename):
    """ Power Spectrum """
    # load x,y data
    x = data[:,0]
    y = data[:,1]

    # load observation data 
    low_obs_data = np.loadtxt("planck_cell_low.txt")
    high_obs_data = np.loadtxt("planck_cell_high.txt")

    # plot
    # double normfactor  = (ell * (ell+1)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(), 2);
    fig, ax = plt.subplots()
    ax.plot(x,y,label='Our model')
    ax.errorbar(low_obs_data[:,0], low_obs_data[:,1], [low_obs_data[:,3], low_obs_data[:,2]], fmt='o', markersize=3, color='k')
    ax.errorbar(high_obs_data[:,0], high_obs_data[:,1], [high_obs_data[:,3], high_obs_data[:,2]], fmt='o', markersize=3, color='k', label='Planck 2018')
    ax.set_ylabel(r"$\frac{C_\ell^{TT} \; (\ell (\ell + 1))}{2 \pi} (\mu K)^2$")
    ax.set_xlabel(r"Multipole $\ell$")
    ax.set_xscale("log")
    ax.grid(); ax.legend()
    plt.savefig(path+filename)
    plt.show()

def plot_polarization(filename):
    """ Polarization """

    high_obs_data = np.loadtxt("planck_cell_high_EE.txt")
    
    x       = data[:,0]
    TCMB0   = 2.7255
    y       = data[:,2] * 1e5 * (1e6 * TCMB0)**2

    fig, ax = plt.subplots()
    ax.plot(x,y, label="Our model")

    l = high_obs_data[:,0]
    norm = 1/(l*(l+1))*1e5*2*3.14159
    Dl = high_obs_data[:,1]*norm
    Dl_pluss = high_obs_data[:,3]*norm
    Dl_minus = high_obs_data[:,2]*norm
    ax.errorbar(l, Dl, [Dl_pluss, Dl_minus], fmt='o', markersize=3, color='k', label='Planck 2018')
    ax.set_xlabel(r"Multipole $\ell$")
    ax.set_ylabel(r"$C_\ell^{EE} \; (10^{-5} \mu K^2)$")
    ax.grid(); ax.legend()
    ax.set_xscale("log")
    plt.savefig(path+filename)
    plt.show()

def plot_pol_temp(filename):
    x = data[:,0]
    y = data[:,3]

    obs_data = np.loadtxt("planck_cell_high_TE.txt")
    l = obs_data[:,0]
    Dl = obs_data[:,1]
    Dl_pluss = obs_data[:,3]
    Dl_minus = obs_data[:,2]
    fig, ax = plt.subplots()
    ax.plot(x,y,label="Our model")
    ax.errorbar(l, Dl, [Dl_pluss, Dl_minus], fmt='o', markersize=3, color='k', label='Planck 2018')
    ax.set_xlabel(r"Multipole $\ell$")
    ax.set_ylabel(r"$\frac{C_\ell^{TE} \; (\ell (\ell + 1))}{2 \pi} (\mu K)^2$")
    ax.grid(); ax.legend()
    ax.set_xscale("log")
    plt.savefig(path+filename)
    plt.show()


def plot_theta(filename):
    """ Theta """

    data = np.loadtxt("theta.txt")

    x = data[:,0]

    y = data[:,1]
    y2 = data[:,2]
    y3 = data[:,3]
    y4 = data[:,4]
    y5 = data[:,5]

    fig, ax = plt.subplots()
    ax.plot(x,y, label = r'$\ell = 6$')
    ax.plot(x,y2, label=r'$\ell = 100$')
    ax.plot(x,y3, label=r'$\ell = 200$')
    ax.plot(x,y4, label=r'$\ell = 500$')
    ax.plot(x,y5, label=r'$\ell = 1000$')
    ax.set_xlim(0, 1200); ax.set_ylim(-0.0065, 0.02)
    ax.set_ylabel(r'$\Theta_\ell$'); ax.set_xlabel(r'$\eta_0 k$')
    ax.legend()
    ax.grid()
    plt.tight_layout()
    plt.savefig(path+filename)
    plt.show()

def plot_neutrino(filename):
    # load x,y data
    x = data[:,0]
    y = data[:,4]

    # plot

    fig, ax = plt.subplots()
    ax.plot(x,y,label='Our model')
    ax.set_ylabel(r"$\frac{C_\ell^{\nu} \; (\ell (\ell + 1))}{2 \pi} (\mu K)^2$")
    ax.set_xlabel(r"Multipole $\ell$")
    ax.set_xscale("log")
    ax.grid()
    plt.tight_layout()
    plt.savefig(path+filename)
    plt.show()

def plot_map():
    x = data[:,0]
    TCMB0   = 2.7255
    inv_norm = 2*3.1415/(x*(x+1)) / (1e6 * TCMB0)**2
    y = data[:,1]*inv_norm
    resolution = 256

    a = hp.sphtfunc.synfast(y, resolution)
    wmap_map_I_smoothed = hp.smoothing(a, fwhm=np.radians(.1))
    hp.mollview(a, title="CMB map")
    plt.savefig(path+"cmb_map.pdf")
    plt.show()


    y = data[:,2]

    a = hp.sphtfunc.synfast(y, resolution)
    wmap_map_I_smoothed = hp.smoothing(a, fwhm=np.radians(.1))
    hp.mollview(a, title="E-mode polarization map")
    plt.savefig(path+"polarization_map.pdf")
    plt.show()


    y = data[:,4]*inv_norm

    a = hp.sphtfunc.synfast(y, resolution)
    wmap_map_I_smoothed = hp.smoothing(a, fwhm=np.radians(.1))
    hp.mollview(a, title="Neutrino map")
    plt.savefig(path+"neutrino_map.pdf")
    plt.show()

def plot_theta_2(filename):
    """ Theta """

    data = np.loadtxt("theta.txt")

    x = data[:,0]

    y = data[:,6]
    y2 = data[:,7]
    y3 = data[:,8]
    y4 = data[:,9]
    y5 = data[:,10]

    fig, ax = plt.subplots()
    ax.plot(x,y, label = r'$\ell = 6$')
    ax.plot(x,y2, label=r'$\ell = 100$')
    ax.plot(x,y3, label=r'$\ell = 200$')
    ax.plot(x,y4, label=r'$\ell = 500$')
    ax.plot(x,y5, label=r'$\ell = 1000$')
    #ax.set_yscale("log")
    ax.set_xlim(0, 1100); ax.set_ylim(0, 1.25e-6)
    ax.set_ylabel(r'$\Theta_\ell^2H_0/ (kc)$'); ax.set_xlabel(r'$\eta_0 k$')
    ax.legend(); ax.grid()
    plt.tight_layout()
    plt.savefig(path+filename)
    plt.show()

def plot_matter_power_spectrum(filename):

    data = np.loadtxt("theta.txt")
    dr7  = np.loadtxt("reid_DR7.txt")
    wmap = np.loadtxt("wmap_act.txt")

    x = data[:,11]
    y = data[:,12]
    k_eq =  0.0155824

    fig, ax = plt.subplots()
    ax.plot(x,y,label='Our model')
    ax.vlines(k_eq, 0.1*y.min(), 10*y.max(), ls='--', color='k', label=r'$k_{eq}$')
    ax.errorbar(dr7[:,0], dr7[:,1], dr7[:,2], label='SDSS Galaxies (DR7 LRG)', fmt='o', markersize=3, color='k')
    ax.errorbar(wmap[:,0], wmap[:,1], wmap[:,2]-wmap[:,1], label='CMB (WMAP + ACT)', fmt='o', markersize=3, color='r')

    ax.set_yscale("log"); ax.set_xscale("log")
    ax.set_xlim(1e-3, 4e-1); ax.set_ylim(1e2, 6e4)
    ax.set_xlabel(r"Wavenumber $k$ (h/Mpc)"); ax.set_ylabel(r"P(k) (Mpc/h)$^3$")
    ax.grid(); ax.legend()
    plt.savefig(path+filename)
    plt.show()

### Call functions
#plot_power_spectrum("power_spectrum.pdf")
#plot_polarization("polarization_power_spectrum.pdf")
#plot_pol_temp("pol_temp.pdf")
#plot_neutrino("neutrino.pdf") # ??
#plot_theta("theta.pdf")
#plot_map()
#plot_theta_2("theta2.pdf")
#plot_matter_power_spectrum("matter_pk.pdf")