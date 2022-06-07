#include"PowerSpectrum.h"

//====================================================
// Constructors
//====================================================

PowerSpectrum::PowerSpectrum(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec, 
    Perturbations *pert,
    double A_s,
    double n_s,
    double kpivot_mpc) : 
  cosmo(cosmo), 
  rec(rec), 
  pert(pert),
  A_s(A_s),
  n_s(n_s),
  kpivot_mpc(kpivot_mpc)
{}

//====================================================
// Do all the solving
//====================================================
void PowerSpectrum::solve(){

  //=========================================================================
  // TODO: Choose the range of k's and the resolution to compute Theta_ell(k)
  //=========================================================================

  //=========================================================================
  // TODO: Make splines for j_ell. 
  // Implement generate_bessel_function_splines
  //=========================================================================
  generate_bessel_function_splines();
  std::cout << "Bessel functions generated!" << '\n'; 

  //=========================================================================
  // TODO: Line of sight integration to get Theta_ell(k)
  // Implement line_of_sight_integration
  //=========================================================================
  line_of_sight_integration(k_array);
  std::cout << "Line of sight integration complete!" << '\n'; 

  //=========================================================================
  // TODO: Integration to get Cell by solving dCell^f/dlogk = Delta(k) * f_ell(k)^2
  // Implement solve_for_cell
  //=========================================================================
  auto cell_TT = solve_for_cell(log_k_array, thetaT_ell_of_k_spline, thetaT_ell_of_k_spline);
  cell_TT_spline.create(ells, cell_TT, "Cell_TT_of_ell");
  
  auto cell_EE = solve_for_cell(log_k_array, thetaE_ell_of_k_spline, thetaE_ell_of_k_spline);
  cell_EE_spline.create(ells, cell_EE, "Cell_EE_of_ell");

  auto cell_TE = solve_for_cell(log_k_array, thetaT_ell_of_k_spline, thetaE_ell_of_k_spline);
  cell_TE_spline.create(ells, cell_TE, "Cell_TE_of_ell");

  auto cell_Nu = solve_for_cell(log_k_array, NuT_ell_of_k_spline, NuT_ell_of_k_spline);
  cell_Nu_spline.create(ells, cell_TE, "Cell_TE_of_ell");
  std::cout << "C_ell spline created!" << '\n'; 
  
  //=========================================================================
  // TODO: Do the same for polarization...
  //=========================================================================
  // ...
  // ...
  // ...
  // ...
}

//====================================================
// Generate splines of j_ell(z) needed for LOS integration
//====================================================

void PowerSpectrum::generate_bessel_function_splines(){
  Utils::StartTiming("besselspline");
  
  // Make storage for the splines
  j_ell_splines = std::vector<Spline>(ells.size()); 
  double z_min = 0.0;
  double z_max = k_max * eta0;
  int n_z = z_max / (2. * M_PI / 16); //10000; //  10 000
  Vector z = Utils::linspace(z_min, z_max, n_z); 
    
  //=============================================================================
  // TODO: Compute splines for bessel functions j_ell(z)
  // Choose a suitable range for each ell
  // NB: you don't want to go larger than z ~ 40000, then the bessel routines
  // might break down. Use j_ell(z) = Utils::j_ell(ell, z)
  //=============================================================================

  for(size_t i = 0; i < ells.size(); i++){
    const int ell = ells[i];
    Vector j_tmp(n_z);

    for(int j = 0; j < n_z; j++){
      j_tmp[j] = Utils::j_ell(ell, z[j]);
    }

    j_ell_splines[i].create(z, j_tmp);

    // Make the j_ell_splines[i] spline
  }

  Utils::EndTiming("besselspline");
}

//====================================================
// Do the line of sight integration for a single
// source function
//====================================================

Vector2D PowerSpectrum::line_of_sight_integration_single(
    Vector & k_array, 
    std::function<double(double,double)> &source_function){
  Utils::StartTiming("lineofsight");

  // Make storage for the results
  Vector2D result = Vector2D(ells.size(), Vector(k_array.size()));

  // #pragma omp parallel
  for (int i = 0; i < ells.size(); i++){
    const int ell = ells[i];
    for(int ik = 0; ik < k_array.size(); ik++){
      //=============================================================================
      // TODO: Implement to solve for the general line of sight integral 
      // F_ell(k) = Int dx jell(k(eta-eta0)) * S(x,k) for all the ell values for the 
      // given value of k
      //=============================================================================

      double sum = 0;
      for(int ix = 0; ix < x_array.size()-1; ix++){
        double arg1    = k_array[ik] * (eta0 - cosmo->eta_of_x(x_array[ix]));
        double arg2    = k_array[ik] * (eta0 - cosmo->eta_of_x(x_array[ix+1]));
        double Theta1 = source_function(x_array[ix],k_array[ik]) * j_ell_splines[i](arg1);
        double Theta2 = source_function(x_array[ix+1],k_array[ik]) * j_ell_splines[i](arg2);
        double dx    = x_array[ix+1] - x_array[ix];

        sum += (Theta1 + Theta2)/2. * dx;
      }
      result[i][ik] = sum;
      // Store the result for Theta_ell(k) in results[ell][ik]
    }
  }

  Utils::EndTiming("lineofsight");
  return result;
}

//====================================================
// Do the line of sight integration
//====================================================
void PowerSpectrum::line_of_sight_integration(Vector & k_array){
  const int nells      = ells.size();
  
  // Make storage for the splines we are to create
  thetaT_ell_of_k_spline = std::vector<Spline>(nells);
  
  //============================================================================
  // TODO: Solve for Theta_ell(k) and spline the result
  //============================================================================

  // Make a function returning the source function
  std::function<double(double,double)> source_function_T = [&](double x, double k){
    return pert->get_Source_T(x,k);
  };

  // Do the line of sight integration
  Vector2D thetaT_ell_of_k = line_of_sight_integration_single(k_array, source_function_T);

  for(size_t i = 0; i < nells; i++){
    thetaT_ell_of_k_spline[i].create(k_array, thetaT_ell_of_k[i]);
  }

  //============================================================================
  // TODO: Solve for ThetaE_ell(k) and spline
  //============================================================================
  if(Constants.polarization){
    thetaE_ell_of_k_spline = std::vector<Spline>(nells);

    std::function<double(double,double)> source_function_E = [&](double x, double k){
      return pert->get_Source_E(x,k);
    };

    Vector2D thetaE_ell_of_k = line_of_sight_integration_single(k_array, source_function_E);

    for(size_t i = 0; i < nells; i++){
      thetaE_ell_of_k_spline[i].create(k_array, thetaE_ell_of_k[i]);
    }
  }
  if(Constants.neutrinos){
    NuT_ell_of_k_spline = std::vector<Spline>(nells);

    std::function<double(double,double)> source_function_Nu = [&](double x, double k){
      return pert->get_Source_Nu(x,k);
    };

    Vector2D Nu_ell_of_k = line_of_sight_integration_single(k_array, source_function_Nu);

    for(size_t i = 0; i < nells; i++){
      NuT_ell_of_k_spline[i].create(k_array, Nu_ell_of_k[i]);
    }
  }
}

//====================================================
// Compute Cell (could be TT or TE or EE) 
// Cell = Int_0^inf 4 * pi * P(k) f_ell g_ell dk/k
//====================================================
Vector PowerSpectrum::solve_for_cell(
    Vector & k_array,
    std::vector<Spline> & f_ell_spline,
    std::vector<Spline> & g_ell_spline){
  const int nells      = ells.size();
  Vector result(nells);

  //============================================================================
  // TODO: Integrate Cell = Int 4 * pi * P(k) f_ell g_ell dk/k
  // or equivalently solve the ODE system dCell/dlogk = 4 * pi * P(k) * f_ell * g_ell
  //============================================================================
  //#pragma omp parallel for
  for(int i = 0; i < nells; i++){
    result[i] = 0;
    for(int ik = 0; ik < k_array.size()-1; ik++){
      double k1 = k_array[ik];
      double k2 = k_array[ik+1];
      double dk = k2 - k1;
      
      double Cell1 = 4.*M_PI*primordial_power_spectrum(k1)*f_ell_spline[i](k1)*g_ell_spline[i](k1)/k1;
      double Cell2 = 4.*M_PI*primordial_power_spectrum(k2)*f_ell_spline[i](k2)*g_ell_spline[i](k2)/k2;
      result[i] += (Cell1 + Cell2)/2.*dk;
    }
  }

  return result;
}

//====================================================
// The primordial power-spectrum
//====================================================

double PowerSpectrum::primordial_power_spectrum(const double k) const{
  return A_s * pow( Constants.Mpc * k / kpivot_mpc , n_s - 1.0);
}

//====================================================
// P(k) in units of (Mpc)^3
//====================================================

double PowerSpectrum::get_matter_power_spectrum(const double x, const double k_mpc) const{
  // Compute the matter power spectrum
  double P_primordial = primordial_power_spectrum(k_mpc)*2*M_PI*M_PI/pow(k_mpc,3.);
  double OmegaM0      = (cosmo->get_OmegaB() + cosmo->get_OmegaCDM());

  double A = pow(Constants.c*k_mpc/cosmo->get_H0(),2.0);
  double B = exp(x)*pert->get_Phi(x,k_mpc)/OmegaM0*2./3.;

  double Delta_M = A*B;

  double pofk = abs(Delta_M*Delta_M)*P_primordial;
  return pofk * pow(cosmo->get_h() / Constants.Mpc, 3.);
}

//====================================================
// Get methods
//====================================================
double PowerSpectrum::get_cell_TT(const double ell) const{
  return cell_TT_spline(ell);
}
double PowerSpectrum::get_cell_TE(const double ell) const{
  return cell_TE_spline(ell);
}
double PowerSpectrum::get_cell_EE(const double ell) const{
  return cell_EE_spline(ell);
}

//====================================================
// Output the cells to file
//====================================================

void PowerSpectrum::output(std::string filename) const{
  // Output in standard units of muK^2
  std::ofstream fp(filename.c_str());
  const int ellmax = int(ells[ells.size()-1]);
  auto ellvalues = Utils::linspace(2, ellmax, ellmax-1);

  auto print_data = [&] (const double ell) {
    double normfactor  = (ell * (ell+1)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(), 2);
    double normfactorN = (ell * (ell+1)) / (2.0 * M_PI) 
      * pow(1e6 * cosmo->get_TCMB() *  pow(4.0/11.0, 1.0/3.0), 2);
    double normfactorL = (ell * (ell+1)) * (ell * (ell+1)); // / (2.0 * M_PI);

    fp << ell                                 << " ";
    fp << cell_TT_spline( ell ) * normfactor  << " ";

    if(Constants.polarization){
      fp << cell_EE_spline( ell ) * normfactorL  << " ";
      fp << cell_TE_spline( ell ) * normfactor*sqrt(normfactorL)  << " ";
    }
    if(Constants.neutrinos){
      fp << cell_Nu_spline( ell ) * normfactor << " ";
    }
    fp << "\n";
  };
  std::for_each(ellvalues.begin(), ellvalues.end(), print_data);
}

void PowerSpectrum::output_x(std::string filename) const{
  // Output in standard units of muK^2
  double OmegaR = cosmo->get_OmegaR(0.0);
  double OmegaM0      = (cosmo->get_OmegaB() + cosmo->get_OmegaCDM());
  double omeganu = 3.046*7.0/8.0*pow((4.0/11.0), (4.0/3.0))*OmegaR;
  double a_eq = (omeganu + OmegaR)/OmegaM0;

  std::cout << "k_eq : " << a_eq*cosmo->H_of_x(log(a_eq))/Constants.c*Constants.Mpc / cosmo->get_h() << '\n';

  std::ofstream fp(filename.c_str());
  // Vector k_array = Utils::linspace(k_min, k_max, 500);

  auto print_data = [&] (const double x) {
    fp << x*cosmo->eta_of_x(0.0) << " "; //*Constants.c/cosmo->H_of_x(0.0)  << " ";
    // Theta_ells 
    fp << thetaT_ell_of_k_spline[4](x)      << " "; // ell = 6
    fp << thetaT_ell_of_k_spline[19](x)     << " "; // ell = 100
    fp << thetaT_ell_of_k_spline[24](x)     << " "; // ell = 200
    fp << thetaT_ell_of_k_spline[32](x)     << " "; // ell = 500
    fp << thetaT_ell_of_k_spline[42](x)     << " "; // ell = 1000

    double H0 = cosmo->get_H0();
    fp << H0/Constants.c*thetaT_ell_of_k_spline[4](x)*thetaT_ell_of_k_spline[4](x)/x      << " "; // ell = 6
    fp << H0/Constants.c*thetaT_ell_of_k_spline[19](x)*thetaT_ell_of_k_spline[19](x)/x    << " "; // ell = 100
    fp << H0/Constants.c*thetaT_ell_of_k_spline[24](x)*thetaT_ell_of_k_spline[24](x)/x    << " "; // ell = 200
    fp << H0/Constants.c*thetaT_ell_of_k_spline[32](x)*thetaT_ell_of_k_spline[32](x)/x    << " "; // ell = 500
    fp << H0/Constants.c*thetaT_ell_of_k_spline[42](x)*thetaT_ell_of_k_spline[42](x)/x    << " "; // ell = 1000
    
    fp << x*Constants.Mpc / cosmo->get_h() << " ";
    fp << get_matter_power_spectrum(0.0, x) << " "; // Mpc ??
    fp << "\n";
  };
  std::for_each(k_array.begin(), k_array.end(), print_data);
}