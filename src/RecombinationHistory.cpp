#include"RecombinationHistory.h"

//====================================================
// Constructors
//====================================================
   
RecombinationHistory::RecombinationHistory(
    BackgroundCosmology *cosmo, 
    double Yp) :
  cosmo(cosmo),
  Yp(Yp)
{
  /*
  Constructor: Constructs
  */
}

//====================================================
// Do all the solving we need to do
//====================================================

void RecombinationHistory::solve(){
  /*
  Initiates the functions that solves tau and g + their derivatives + Tb
  */
    
  // Compute and spline Xe, ne, Tb
  solve_number_density_electrons();
   
  // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx, ...
  solve_for_optical_depth_tau();
}

void RecombinationHistory::solve_number_density_electrons(){
  /*
  Solves and splines the number density and fraction of electrons and Tb.
  */
  Utils::StartTiming("Xe");
  
  // Set up arrays to store our values to be splined
  Vector x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);
  Vector Xe_arr(npts_rec_arrays);
  Vector ne_arr(npts_rec_arrays);
  Vector Tb_arr(npts_rec_arrays);

  // Set up current values to be changed each step
  double Xe_current;
  double ne_current;
  bool saha_regime = true; // turn false when outside regime

  int i = 0; // keeps track of step number
  while(saha_regime){
    // Get X_e from solving the Saha equation 
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);
    Xe_current = Xe_ne_data.first;
    ne_current = Xe_ne_data.second;

    // Electron fraction and number density
    saha_regime = Xe_current > Xe_saha_limit; // check saha limit
    if(not saha_regime){
      // if we are not in saha limit we break the loop and 
      // move on the peebles algorithm
      break;
    }
    // if the above does not pass (no break) we record
    Xe_arr[i] = Xe_current;
    ne_arr[i] = log(ne_current);
    Tb_arr[i] = TCMB0*exp(-x_array[i]); // baryon temp is approx photon temp
    i += 1; // udate to next time step
  }
  // The Peebles algorithm: 
  // Once outside saha regime we make a temporary spline of 
  // the solved peebles equation that we can extract values
  // of...

  // Define ODE
  ODESolver peebles_Xe_ode;
  ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
    return rhs_peebles_ode(x, Xe, dXedx);
  };

  // Define temporary x-array to solve over
  Vector x_array_tmp = Utils::linspace(x_array[i], x_end, 1000);

  // Set the IC using last Xe step... y(a -> 0) = 0
  double Xeini = Xe_current;
  double y_ini = 0.0;
  Vector Xe_ic{Xeini, y_ini};

  // Solve the ODE
  peebles_Xe_ode.solve(dXedx, x_array_tmp, Xe_ic);

  // Get the solution S(0) = Xe, S(1) = y...
  auto Xe_tmp = peebles_Xe_ode.get_data_by_component(0);
  auto y_tmp = peebles_Xe_ode.get_data_by_component(1);

  // Go back to Tb using y...
  Vector Tb_tmp(y_tmp.size());
  for(int k=0; k<y_tmp.size();k++){
    Tb_tmp[k] = TCMB0*exp(-x_array_tmp[k])*(y_tmp[k] + 1);
  }

  // Make temporary splines of solutions
  Spline Xe_tmp_spline;
  Spline Tb_tmp_spline;
  Xe_tmp_spline.create(x_array_tmp, Xe_tmp);
  Tb_tmp_spline.create(x_array_tmp, Tb_tmp);
  Xe_tmp_spline.set_out_of_bounds_warning(true);
  Tb_tmp_spline.set_out_of_bounds_warning(true);

  // Set up values to be changed
  double OmegaB;
  double n_H;
  double ne;

  // Now we are able to get the values in the Peebles regime
  while(not saha_regime and i < npts_rec_arrays){
    n_H         = (OmegaB0*rho_c0/(m_H*pow(exp(x_array[i]), 3)));

    Xe_current = Xe_tmp_spline(x_array[i]); 
    ne_current = n_H*Xe_current;
    // Record values
    Xe_arr[i] = Xe_current;
    ne_arr[i] = log(ne_current);
    Tb_arr[i] = Tb_tmp_spline(x_array[i]);
    // Update timestep
    i += 1;
  }
  // Make final splines...
  Xe_of_x_spline.create(x_array, Xe_arr);
  Xe_of_x_spline.set_out_of_bounds_warning(true);

  ne_of_x_spline.create(x_array, ne_arr);
  ne_of_x_spline.set_out_of_bounds_warning(true);

  Tb_of_x_spline.create(x_array, Tb_arr);
  Tb_of_x_spline.set_out_of_bounds_warning(true);

  Utils::EndTiming("Xe");
}

std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
  // Solve the Saha equation to get ne and Xe
  const double a           = exp(x);
  const double n_H         = (OmegaB0*rho_c0/(m_H*pow(a,3)));
  
  // The RHS in Saha equation
  double A  = dimension_Xe_saha*exp(-(epsilon_0)/(k_b*TCMB0/a))/n_H * pow(m_e*(TCMB0/a)/(2*pi),3./2.);
  
  double Xe = sqrt(A*A + 4*A)/(A + 2);
  double ne = n_H*Xe;

  return std::pair<double,double>(Xe, ne);
}

int RecombinationHistory::rhs_peebles_ode(double x, const double *Xe, double *dXedx){
  // The right hand side of the dXedx Peebles ODE
  // Current value of a, X_e and y
  const double X_e         = Xe[0]; // select from pointer
  const double y           = Xe[1];
  const double a           = exp(x);

  // Calculate the temperatures
  double T        = TCMB0/a; 
  double Tb       = T*(y+1);
  
  // Cosmological parameter
  const double H           = cosmo->H_of_x(x); // unit H0 = 100h km/s/Mpc

  // Peebles equations
  double n_H      = (1 - Yp)*3*pow(cosmo->H_of_x(0),2)*OmegaB0/(8*pi*G*m_H*pow(a,3.0)); // 1/m3
  double n1s      = (1. - X_e)*n_H; // 1/m3
  double lambda_a = (pow(hbar,-3)*pow(c,-3)) * H*pow(3*epsilon_0,3)/(pow(8*pi, 2)*n1s); // 1/s
  double phi2     = 0.448*log(epsilon_0/(Tb*k_b)); // dim. less
  double alpha2   = (c * pow(k_b,-1./2)) * 64*pi/sqrt(27*pi)*pow(alpha/m_e,2)*sqrt(epsilon_0/(Tb))*phi2; // m3/s
  double beta     = (pow(hbar,-3)*pow(k_b, 3./2)) * alpha2*pow(m_e*Tb/(2*pi), 3./2) * exp(-(epsilon_0)/(Tb*k_b)); // 1/s
  double beta2    = (pow(hbar,-3)*pow(k_b, 3./2)) * alpha2*pow(m_e*Tb/(2*pi), 3./2) * exp(-(1./4.)*(epsilon_0)/(Tb*k_b)); // 1/s

  // RHS of peebles
  double Cr = (lambda_2s1s + lambda_a)/(lambda_2s1s + lambda_a + beta2); // dim. less
  dXedx[0] = Cr/H * (beta*(1. - X_e) - n_H*alpha2 * pow(X_e,2));

  // RHS of baryon temperature using y
  double A  = 8./3.*m_H/m_e*OmegaR0/OmegaB0/(a*cosmo->H_of_x(x))*n_H*X_e*sigma_T;
  dXedx[1] = - y - 1 - A*y;

  return GSL_SUCCESS;
}

void RecombinationHistory::solve_for_optical_depth_tau(){
  /*
    Solve for the optical depth tau, compute the 
    visibility function and spline the result
  */
  Utils::StartTiming("opticaldepth");

  // Set up x-arrays to integrate over
  const int npts = 1000;
  Vector x_array = Utils::linspace(x_end, x_start, npts);

  // The ODE system dtau/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){
    // Set the derivative for photon optical depth
    dtaudx[0] = - c*ne_of_x(x)*Constants.sigma_T/(cosmo->H_of_x(x));
    return GSL_SUCCESS;
  };
  // initial conditions when we begin at x_end
  double tauini = 0;
  Vector tau_ic{tauini};

  // Solve the ODE
  ODESolver tau_ode;
  tau_ode.solve(dtaudx, x_array, tau_ic);
  auto tau = tau_ode.get_data_by_component(0);
  
  // Flip tau 
  Vector tau_norm(npts);
  Vector x_array_norm(npts);
  for (int i = 0; i<npts; i++){
    tau_norm[npts-i-1] = tau[i];
    x_array_norm[npts-i-1] = x_array[i];
  }

  // Make spline of solution
  tau_of_x_spline.create(x_array_norm, tau_norm);
  tau_of_x_spline.set_out_of_bounds_warning(true);

  // Make seperate spline for tau derivative
  Vector deriv_tau(npts);
  for(int i = 0; i<deriv_tau.size(); i++){
    deriv_tau[i] = - c*ne_of_x(x_array_norm[i])*Constants.sigma_T/(cosmo->H_of_x(x_array_norm[i]));
  }
  // Create spline for tau derivative
  deriv_tau_of_x_spline.create(x_array_norm, deriv_tau);
  deriv_tau_of_x_spline.set_out_of_bounds_warning(true);

  // Calculate the visibility function
  Vector g_tilde_arr(npts);
  double sum = 0.0; // Use sum to normalize g
  for (int i = 0; i<npts; i++){
    g_tilde_arr[i] = -dtaudx_of_x(x_array_norm[i])*exp(-tau_of_x(x_array_norm[i]));
    sum += g_tilde_arr[i]; // add every element to sum 
  }
  // Normalize g by dividing every element by sum
  for (int i = 0; i<npts; i++){
    g_tilde_arr[i] = g_tilde_arr[i]/sum; 
  }

  // Spline g
  g_tilde_of_x_spline.create(x_array_norm, g_tilde_arr);
  g_tilde_of_x_spline.set_out_of_bounds_warning(true);

  Utils::EndTiming("opticaldepth");
}

//====================================================
// Get methods
//====================================================

double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const{
  return deriv_tau_of_x_spline(x);
}

double RecombinationHistory::ddtauddx_of_x(double x) const{
  return deriv_tau_of_x_spline.deriv_x(x);
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{
  return g_tilde_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{

  return g_tilde_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::Xe_of_x(double x) const{
  return Xe_of_x_spline(x);
}

double RecombinationHistory::ne_of_x(double x) const{

  return exp(ne_of_x_spline(x));
}

double RecombinationHistory::Tb_of_x(double x) const{

  return Tb_of_x_spline(x);
}

double RecombinationHistory::get_Yp() const{
  return Yp;
}

//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const{
  std::cout << "\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "Yp:          " << Yp          << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts       = 6500;
  const double x_min   = x_start;
  const double x_max   = x_end;

  Vector x_array = Utils::linspace(x_min, x_max, npts);
  auto print_data = [&] (const double x) {
    fp << x                    << " ";
    fp << Xe_of_x(x)           << " ";
    fp << ne_of_x(x)           << " ";
    fp << tau_of_x(x)          << " ";
    fp << dtaudx_of_x(x)       << " ";
    fp << ddtauddx_of_x(x)     << " ";
    fp << g_tilde_of_x(x)      << " ";
    fp << dgdx_tilde_of_x(x)   << " ";
    fp << ddgddx_tilde_of_x(x) << " ";
    fp << Tb_of_x(x)           << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

