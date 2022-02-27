#ifndef _BACKGROUNDCOSMOLOGY_HEADER
#define _BACKGROUNDCOSMOLOGY_HEADER
#include <iostream>
#include <fstream>
#include "Utils.h"

using Vector = std::vector<double>;

class BackgroundCosmology{
  private:
   
    // Cosmological parameters
    double h;                       // Little h = H0/(100km/s/Mpc)
    double OmegaB;                  // Baryon density today
    double OmegaCDM;                // CDM density today
    double OmegaLambda;             // Dark energy density today
    double Neff;                    // Effective number of relativistic species (3.046 or 0 if ignoring neutrinos)
    double TCMB;                    // Temperature of the CMB today in Kelvin
   
    // Derived parameters
    double OmegaR;                  // Photon density today (follows from TCMB)
    double OmegaNu;                 // Neutrino density today (follows from TCMB and Neff)
    double OmegaK;                  // Curvature density = 1 - OmegaM - OmegaR - OmegaNu - OmegaLambda
    double OmegaM;                  // Combined mass (CDM + Baryons)
    double OmegaRad;                // Combined radiation (R + Nu)
    double H0;                      // The Hubble parameter today H0 = 100h km/s/Mpc

    // Useful parameter
    int mytruth;

    // Start and end of x-integration (can be changed)
    double x_start = log(1e-20); // Constants.x_start;
    double x_end   = 2.0; // Constants.x_end;

    // Vectors to be made
    Vector x_array;
    Vector Hp_array;

    // Splines to be made
    Spline eta_of_x_spline{"eta"};
    Spline Hp_of_x_spline{"Hp"};
    Spline t_of_x_spline{"t"};
    Spline x_of_t_spline{"x"};
 
  public:

    // Constructors 
    BackgroundCosmology() = delete;
    BackgroundCosmology(
        double h, 
        double OmegaB, 
        double OmegaCDM, 
        double OmegaK,
        double Neff, 
        double TCMB
        );

    // Print some useful info about the class
    void info() const;

    // Do all the solving
    void solve();

    // Output some results to file
    void output(const std::string filename) const;

    // Get functions that we must implement
    double eta_of_x(double x) const;
    double H_of_x(double x) const;
    double Hp_of_x(double x) const;
    double dHpdx_of_x(double x) const;
    double ddHpddx_of_x(double x) const;
    double t_of_x(double x) const;
    double x_of_t(double t) const;
    double comoving_distance(double eta) const;

    double get_OmegaB(double x = 0.0) const; 
    double get_OmegaM(double x = 0.0) const; 
    double get_OmegaR(double x = 0.0) const;
    double get_OmegaRtot(double x = 0.0) const; 
    double get_OmegaNu(double x = 0.0) const;
    double get_OmegaCDM(double x = 0.0) const; 
    double get_OmegaLambda(double x = 0.0) const; 
    double get_OmegaK(double x = 0.0) const; 
    double get_OmegaMnu(double x = 0.0) const; 
    double get_H0() const;
    double get_h() const;
    double get_Neff() const;
    double get_TCMB(double x = 0.0) const;
};

#endif
