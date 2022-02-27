#include "BackgroundCosmology.h"

//====================================================
// Constructors
//====================================================
    
BackgroundCosmology::BackgroundCosmology(
    double h, 
    double OmegaB, 
    double OmegaCDM, 
    double OmegaK,
    double Neff, 
    double TCMB) :
  h(h),
  OmegaB(OmegaB),
  OmegaCDM(OmegaCDM),
  OmegaK(OmegaK),
  Neff(Neff), 
  TCMB(TCMB)
{
  // We compute the missing constants:
  H0      = Constants.H0_over_h*h;            // Hubble constant
  OmegaNu = 0.0;                              // Neutrino density today (master students ignore this)
  OmegaK  = 0.0;                              // Flat curvature (Curvature density = 1 - OmegaM - OmegaR - OmegaNu - OmegaLambda)
  double pi = 3.14159265358979;
  OmegaR  = 2*std::pow(pi,2)/30 * std::pow(Constants.k_b*TCMB,4)/pow(Constants.hbar, 3)/pow(Constants.c,5)
  * 8*pi*Constants.G/(3*std::pow(H0,2));  // Photon density today (follows from TCMB):
  // and lastly, dark energy parameter which follows from a(today) = 1 ...
  OmegaLambda = 1 - (OmegaB + OmegaCDM + OmegaR + OmegaNu + OmegaK);

  // Combining the parameters simplifies some calculations:
  OmegaM = OmegaB + OmegaCDM;
  OmegaRad = OmegaR + OmegaNu;

  // It is useful to collect some variables here
  mytruth = (OmegaK < 0)*1 + (OmegaK == 0)*2 + (OmegaK > 0)*3;
  int npts = 300;

  // We create our vectors
  x_array = Utils::linspace(x_start, x_end, npts);   // Range of x and number of points for splines
  
  // We can calculate Hp and its derivatives using the following method:
  /*
  Vector Hp_array_tmp(npts); // Make temporary array to later overwrite Hp_array

  // Fill Hp_array_tmp using x_array
  for(size_t i = 0; i < x_array.size(); i++) 
    Hp_array_tmp[i] = Hp_of_x(x_array[i]); 

  Hp_array = Hp_array_tmp; // Overwrite Hp_array to Hp_array_tmp

  // We create some useful splines
  Hp_of_x_spline.create(x_array, Hp_array);
  Hp_of_x_spline.set_out_of_bounds_warning(true);
  */
}

void BackgroundCosmology::solve(){
  /*
  Solve does exactly that, we solve differential equations and set up splines to be used
  in functions.
  */
  // First we want to solve eta
  Utils::StartTiming("Eta");
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){
    // Set the rhs of the detax ODE
    detadx[0] = Constants.c/Hp_of_x(x);
    return GSL_SUCCESS;
  };

  // Set the IC. For eta should be large negative number
  double etaini = Constants.c/Hp_of_x(x_start);
  Vector eta_ic{etaini};

  // Solve the ODE
  ODESolver ode;
  ode.solve(detadx, x_array, eta_ic);

  // Get the solution (we only have one component so index 0 holds the solution)
  auto eta_array = ode.get_data_by_component(0);

  // Make eta_of_x_spline of solution
  eta_of_x_spline.create(x_array, eta_array);
  eta_of_x_spline.set_out_of_bounds_warning(true);

  Utils::EndTiming("Eta");

  // Now we solve t(x)
  Utils::StartTiming("t");
  // Define ODE
  ODEFunction dtdx = [&](double x, const double *t, double *dtdx){
    dtdx[0] = 1/H_of_x(x);
    return GSL_SUCCESS;
  };
  double tinit = 1/(2*H_of_x(x_start));
  Vector t_ic{tinit};

  // Now we solve
  ode.solve(dtdx, x_array, t_ic);

  auto t_array = ode.get_data_by_component(0);

  t_of_x_spline.create(x_array, t_array);
  t_of_x_spline.set_out_of_bounds_warning(true);

  x_of_t_spline.create(t_array, x_array);
  x_of_t_spline.set_out_of_bounds_warning(true);

  Utils::EndTiming("t");
}

//====================================================
// Get methods
//====================================================

double BackgroundCosmology::H_of_x(double x) const{
  double H_of_x = H0*std::sqrt(OmegaM*std::exp(-3*x) 
  + OmegaRad*std::exp(-4*x) + OmegaK*std::exp(-2*x) + OmegaLambda);

  return H_of_x;
}

double BackgroundCosmology::Hp_of_x(double x) const{
  double Hp = H0*std::sqrt(OmegaM*std::exp(-x) + OmegaRad*std::exp(-2*x) 
  + OmegaK + OmegaLambda*std::exp(2*x));
  return Hp;
}

double BackgroundCosmology::dHpdx_of_x(double x) const{
  // Analytical derivative
  double v = OmegaM*std::exp(-x) + OmegaRad*std::exp(-2*x) + OmegaK + OmegaLambda*std::exp(2*x);
  double u = -OmegaM*std::exp(-x) - 2*OmegaRad*std::exp(-2*x) + 2*OmegaLambda*std::exp(2*x);

  return 0.5*H0*u/std::sqrt(v);
  // return Hp_of_x_spline.deriv_x(x); // Numerical derivative
}

double BackgroundCosmology::ddHpddx_of_x(double x) const{
  // Product in square root in Hp(x)
  double v = OmegaM*std::exp(-x) + OmegaRad*std::exp(-2*x) + OmegaK + OmegaLambda*std::exp(2*x); 
  // derivative of v
  double u = -OmegaM*std::exp(-x) - 2*OmegaRad*std::exp(-2*x) + 2*OmegaLambda*std::exp(2*x);
  // s = e^x*u
  double s = 0.5*u/std::sqrt(v); //-3*OmegaM*std::exp(-2*x) + -4*OmegaRad*std::exp(-3*x) + -2*OmegaK*std::exp(-x); 
  // useful simplification
  double w = (OmegaM*std::exp(-x) + 4*OmegaRad*std::exp(-2*x) + 4*OmegaLambda*std::exp(2*x)); //6*OmegaM*std::exp(-2*x) + 12*OmegaRad*std::exp(-3*x) + 2*OmegaK*std::exp(-x); 

  return H0/(2*std::sqrt(v)) * (w - u*u/(2*v)); //H0/(2*v) * (std::sqrt(v)*w - u*s); 
  //return Hp_of_x_spline.deriv_xx(x);; // Numerical double derivative
}

double BackgroundCosmology::t_of_x(double x) const{
  /*
  Returns the time of the universe as a function of x=log(a)
  ARGS:
    - x (double)  : Input conformal time
  RETURNS:
    - t_of_x_spline(x)  (double)  : The time given x using spline
  
  */
  return t_of_x_spline(x);
}

double BackgroundCosmology::x_of_t(double t) const{
  /*
  Returns the conformal time as a function of t
  ARGS: 
    - t (double)  : Time
  RETURNS:
    - x_of_t_spline(x)  (double)  : The conformal time given t 
  */
 return x_of_t_spline(t);
}

double BackgroundCosmology::comoving_distance(double eta) const{
  // define chi 
  double chi = eta_of_x(0) - eta;
  // use cases instead of if statemenets

  switch(mytruth){
    case(1):
      return chi*sin(sqrt(abs(OmegaK))*H0*chi/Constants.c)/(sqrt(abs(OmegaK))*H0*chi/Constants.c);
    case(2):
      return chi;
    case(3):
      return chi*sinh(sqrt(abs(OmegaK))*H0*chi/Constants.c)/sqrt(abs(OmegaK))*H0*chi/Constants.c;
  }
  return 0.0;
}

double BackgroundCosmology::get_OmegaB(double x) const{ 
  if(x == 0.0) return OmegaB;
  else return OmegaB/(std::exp(3*x)*pow(H_of_x(x)/H0,2));
}

double BackgroundCosmology::get_OmegaR(double x) const{ 
  if(x == 0.0) return OmegaR;
  else return OmegaR/(std::exp(4*x)*std::pow(H_of_x(x),2)/pow(H0,2));
}

double BackgroundCosmology::get_OmegaNu(double x) const{ 
  if(x == 0.0) return OmegaNu;
  else return OmegaNu/(std::exp(4*x)*std::pow(H_of_x(x),2)/pow(H0,2));
}

double BackgroundCosmology::get_OmegaCDM(double x) const{
  if(x == 0.0) return OmegaCDM;
  else return OmegaCDM/(std::exp(3*x)*std::pow(H_of_x(x),2)/pow(H0,2));
}

double BackgroundCosmology::get_OmegaLambda(double x) const{ 
  if(x == 0.0) return OmegaLambda;
  else return OmegaLambda/(pow(H_of_x(x),2)/pow(H0,2));
}

double BackgroundCosmology::get_OmegaK(double x) const{ 
  if(x == 0.0) return OmegaK;
  else return OmegaK/(std::exp(2*x)*pow(H_of_x(x),2)/pow(H0,2));
}

double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
}

double BackgroundCosmology::get_H0() const{ 
  return H0; 
}

double BackgroundCosmology::get_h() const{ 
  return h; 
}

double BackgroundCosmology::get_Neff() const{ 
  return Neff; 
}

double BackgroundCosmology::get_TCMB(double x) const{ 
  if(x == 0.0) return TCMB;
  return TCMB * exp(-x); 
}

//====================================================
// Print out info about the class
//====================================================
void BackgroundCosmology::info() const{ 
  std::cout << "\n";
  std::cout << "Info about cosmology class:\n";
  std::cout << "OmegaB:      " << OmegaB      << "\n";
  std::cout << "OmegaCDM:    " << OmegaCDM    << "\n";
  std::cout << "OmegaLambda: " << OmegaLambda << "\n";
  std::cout << "OmegaK:      " << OmegaK      << "\n";
  std::cout << "OmegaNu:     " << OmegaNu     << "\n";
  std::cout << "OmegaR:      " << OmegaR      << "\n";
  std::cout << "Neff:        " << Neff        << "\n";
  std::cout << "h:           " << h           << "\n";
  std::cout << "TCMB:        " << TCMB        << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string filename) const{
  const double x_min = log(1e-5);
  const double x_max =  1.0;
  const int    n_pts =  200;
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  std::cout << "Age of the universe is " << t_of_x(0)/pow(10,9)/365/24/60/60 << " Gyrs\n";
  std::cout << "Conformal time today is " << eta_of_x(0)/Constants.c/pow(10,9)/365/24/60/60 << " Gyrs\n";

  double omeganu = 3.046*7.0/8.0*pow((4.0/11.0), (4.0/3.0))*OmegaR;
  double a_m_rad_eq = (omeganu + OmegaR)/OmegaM; // Including neutrinos (number from calculator)
  //double a_m_rad_eq = OmegaRad/OmegaM;
  double  m_rad_eq = log(a_m_rad_eq);

  double a_l_m_eq = pow(OmegaM/OmegaLambda, 1.0/3.0);
  double l_m_eq = log(a_l_m_eq);

  double a_acc = pow((OmegaCDM+OmegaB)/(2*OmegaLambda), 1.0/3.0);
  double l_acc = log(a_acc);

  std::cout << "Radiation-Matter Equality:" << "\n";
  std::cout << "t = " << t_of_x(m_rad_eq)/365/24/60/60 << " yrs\n";
  std::cout << "x = " << log(a_m_rad_eq) << "\n";
  std::cout << "a = " << a_m_rad_eq << "\n";
  std::cout << "z = " << 1/a_m_rad_eq - 1 << "\n\n";
  std::cout << "Matter-Dark Energy Equality:" << "\n";
  std::cout << "t = " << t_of_x(l_m_eq)/pow(10,9)/365/24/60/60 << " Gyrs\n";
  std::cout << "x = " << log(a_l_m_eq) << "\n";
  std::cout << "a = " << a_l_m_eq << "\n";
  std::cout << "z = " << 1/a_l_m_eq - 1 << "\n\n";
  std::cout << "Accelerated universe:" << "\n";
  std::cout << "t = " << t_of_x(l_acc)/pow(10,9)/365/24/60/60 << " Gyrs\n";
  std::cout << "x = " << l_acc << "\n";
  std::cout << "a = " << a_acc << "\n";
  std::cout << "z = " << 1/a_acc - 1 << "\n\n";

  //std::cout << pow((omeganu + OmegaR)/OmegaM,2)/(2*H0*sqrt(omeganu + OmegaR))/365/24/60/60 << "\n";
  //std::cout << sqrt(2*sqrt(omeganu + OmegaR))/pow(2.0/3.0*sqrt(OmegaM),3.0/2.0)/H0/365/24/60/60 << "\n";
  std::cout << pow(((omeganu + OmegaR)/OmegaM*pow(2*H0*sqrt(OmegaM), 3.0/2.0)/pow(2.0/3.0*H0*sqrt(omeganu + OmegaR), 6.0)),2.0/9.0)/365/24/60/60  << "\n";

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                               << " ";
    fp << eta_of_x(x)                     << " ";
    fp << Hp_of_x(x)                      << " ";
    fp << dHpdx_of_x(x)                   << " ";
    fp << ddHpddx_of_x(x)                 << " "; 
    fp << get_OmegaB(x)                   << " ";
    fp << get_OmegaCDM(x)                 << " ";
    fp << get_OmegaLambda(x)              << " ";
    fp << get_OmegaR(x)                   << " ";
    fp << get_OmegaNu(x)                  << " ";
    fp << get_OmegaK(x)                   << " ";
    fp << H_of_x(x)/H0                    << " "; 
    fp << t_of_x(x)                       << " ";
    fp << comoving_distance(eta_of_x(x))  << " ";
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

