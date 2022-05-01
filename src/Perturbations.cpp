#include"Perturbations.h"

//====================================================
// Constructors
//====================================================

Perturbations::Perturbations(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec) : 
  cosmo(cosmo), 
  rec(rec)
{}

//====================================================
// Do all the solving
//====================================================

void Perturbations::solve(){

  // Integrate all the perturbation equation and spline the result
  integrate_perturbations();
  std::cout << "Integration complete!" << '\n';

  // Compute source functions and spline the result
  //compute_source_functions();
  //std::cout << "Source functions complete!" << '\n';
}

//====================================================
// The main work: integrate all the perturbations
// and spline the results
//====================================================

void Perturbations::integrate_perturbations(){
  Utils::StartTiming("integrateperturbation");

  //===================================================================
  // TODO: Set up the k-array for the k's we are going to integrate over
  // Start at k_min end at k_max with n_k points with either a
  // quadratic or a logarithmic spacing
  //===================================================================
  Vector k_array(n_k);

  // Logarithmic spacing
  Vector k_tmp = Utils::linspace(log(k_min), log(k_max), n_k);
  // Quadratic spacing
  // Vector k_tmp = Utils::linspace(sqrt(k_min), sqrt(k_max), n_k);

  // Fill k_array
  for(int i = 0; i < n_k; i++){
    k_array[i] = exp(k_tmp[i]);
    // k_array[i] = pow(k_tmp[i], 2);
  }

  int N = n_k*(n_x-1);
  Vector delta_cdm_arr(N);
  Vector delta_b_arr(N);
  Vector v_cdm_arr(N);
  Vector v_b_arr(N);
  Vector Phi_arr(N);
  Vector Pi_arr(N);
  Vector Psi_arr(N);

  Vector2D Theta_arr(Constants.n_ell_theta, Vector(N));
  Vector2D Thetap_arr(Constants.n_ell_thetap, Vector(N));
  Vector2D Nu_arr(Constants.n_ell_neutrinos, Vector(N));

  // Loop over all wavenumbers
  //#pragma omp parallel for schedule(dynamic, 1)
  for(int ik = 0; ik < n_k; ik++){

    // Progress bar...
    if( (10*ik) / n_k != (10*ik+10) / n_k ) {
      std::cout << (100*ik+100)/n_k << "% " << std::flush;
      if(ik == n_k-1) std::cout << std::endl;
    }

    // Current value of k
    double k = k_array[ik];

    // Find value to integrate to (x_value, x_index)
    auto x_end_tight = get_tight_coupling_time(k);

    //===================================================================
    // TODO: Tight coupling integration
    // Remember to implement the routines:
    // set_ic : The IC at the start
    // rhs_tight_coupling_ode : The dydx for our coupled ODE system
    //===================================================================

    // Set up initial conditions in the tight coupling regime
    auto y_tight_coupling_ini = set_ic(x_start, k);

    // The tight coupling ODE system
    ODEFunction dydx_tight_coupling = [&](double x, const double *y, double *dydx){
      return rhs_tight_coupling_ode(x, k, y, dydx);
    };

    // Integrate from x_start -> x_end_tight
    Vector x_tc = Utils::linspace(x_start, x_end_tight.first, x_end_tight.second);
    ODESolver tight_coupling;//(1e-6, 1e-12, 1e-12);
    tight_coupling.solve(dydx_tight_coupling, x_tc, y_tight_coupling_ini);
    Vector2D y_tc = tight_coupling.get_data_transpose();
    auto y_tc_final = tight_coupling.get_final_data();
    
    //===================================================================
    // TODO: Full equation integration
    // Remember to implement the routines:
    // set_ic_after_tight_coupling : The IC after tight coupling ends
    // rhs_full_ode : The dydx for our coupled ODE system
    //===================================================================

    // Set up initial conditions (y_tight_coupling is the solution at the end of tight coupling)
    auto y_full_ini = set_ic_after_tight_coupling(y_tc_final, x_end_tight.first, k);

    // The full ODE system
    ODEFunction dydx_full = [&](double x, const double *y, double *dydx){
      return rhs_full_ode(x, k, y, dydx);
    };

    // Integrate from x_end_tight -> x_end
    Vector x_full = Utils::linspace(x_end_tight.first, x_end, n_x - x_end_tight.second);
    ODESolver full_ode;//(1e-6, 1e-12, 1e-12);
    full_ode.solve(dydx_full, x_full, y_full_ini);
    Vector2D y_full = full_ode.get_data_transpose();

   // We loop over ix and record values
   // First loop records the tight coupling regime
   for (int ix=0;ix<x_end_tight.second;ix++){
      int i = ix + (n_x-1)*ik;
      delta_cdm_arr[i] = y_tc[Constants.ind_deltacdm_tc][ix];
      delta_b_arr[i] = y_tc[Constants.ind_deltab_tc][ix];
      v_cdm_arr[i] = y_tc[Constants.ind_vcdm_tc][ix];
      v_b_arr[i] = y_tc[Constants.ind_vb_tc][ix];
      Phi_arr[i] = y_tc[Constants.ind_Phi_tc][ix];

      double x = x_tc[ix];
      double Hp         = cosmo->Hp_of_x(x);
      double ck_over_hp = Constants.c*k/Hp;
      double dtaudx     = rec->dtaudx_of_x(x);

      // Record Theta values in tight coupling regime
      double Theta0 = y_tc[Constants.ind_start_theta + 0][ix];
      double Theta1 = y_tc[Constants.ind_start_theta + 1][ix];
      double Theta2;
      if(Constants.polarization){
        Theta2 = -8.0*ck_over_hp/(15.0*dtaudx)*Theta1;
      }
      else{
        Theta2 = -20.0*ck_over_hp/(45.0*dtaudx)*Theta1;
      }

      Theta_arr[0][i] = Theta0;
      Theta_arr[1][i] = Theta1;
      Theta_arr[2][i] = Theta2;

      for(int l = 3; l < Constants.n_ell_theta; l++){
        Theta_arr[l][i] = - l/(2.0*l+1.0)*ck_over_hp/dtaudx*Theta_arr[l-1][i];
      }

      if(Constants.polarization){
        // Record polarization values in tight coupling regime
        Thetap_arr[0][i] = 5.0/4.0*Theta2;
        Thetap_arr[1][i] = -ck_over_hp/(4.0*dtaudx)*Theta2;
        Thetap_arr[2][i] = 0.25*Theta2;
        for(int l = 3; l < Constants.n_ell_thetap; l++){
          Thetap_arr[l][i] = -l/(2.0*l+1.0)*ck_over_hp/dtaudx*Thetap_arr[l-1][i];
        }
        Pi_arr[i] = Theta_arr[2][i] + Thetap_arr[0][i] + Thetap_arr[2][i];
      }
      else{
        Pi_arr[i] = Theta_arr[2][i];
      }

      if(Constants.neutrinos){
        // Record neutrino values in tight coupling regime
        for(int l = 0; l < Constants.n_ell_neutrinos; l++){
          Nu_arr[l][i] = y_tc[Constants.ind_start_nu_tc + l][ix];
        }
        Psi_arr[i] = - Phi_arr[i] - 12.0*pow(H0/(Constants.c*k*exp(x)),2)*(OmegaR0*Theta_arr[2][i] + OmegaNu0*Nu_arr[2][i]);
      }
      else{
        Psi_arr[i] = - Phi_arr[i] - 12.0*pow(H0/(Constants.c*k*exp(x)),2)*(OmegaR0*Theta_arr[2][i]);
      }
   }
   // Record values in the "full ODE" regime
   for (int ix=1;ix<n_x - x_end_tight.second;ix++){
      int i = (ix + x_end_tight.second-1) + (n_x-1)*ik;
      // Record scalar values
      delta_cdm_arr[i] = y_full[Constants.ind_deltacdm][ix];
      delta_b_arr[i]   = y_full[Constants.ind_deltab][ix];
      v_cdm_arr[i]     = y_full[Constants.ind_vcdm][ix];
      v_b_arr[i]       = y_full[Constants.ind_vb][ix];
      Phi_arr[i]       = y_full[Constants.ind_Phi][ix];
      
      // Record Theta values
      for(int l = 0; l < Constants.n_ell_theta; l++){
        Theta_arr[l][i]   = y_full[Constants.ind_start_theta+l][ix];
      }
      
      double x = x_full[ix];

      // If polarization then record these 
      if(Constants.polarization){
        for(int l = 0; l < Constants.n_ell_theta; l++){
          Thetap_arr[l][i]  = y_full[Constants.ind_start_thetap+l][ix];
        }
        // Also record Pi as follows
        Pi_arr[i]  = Theta_arr[2][i] + Thetap_arr[0][i] + Thetap_arr[2][i];
      }
      else{
        // If we ignore polarization record Pi as:
        Pi_arr[i]  = Theta_arr[2][i];
      }

      // If neutrinos then record these
      if(Constants.neutrinos){
        for(int l = 0; l < Constants.n_ell_theta; l++){
          Nu_arr[l][i]      = y_full[Constants.ind_start_nu+l][ix];
        }
        // Also record Psi as 
        Psi_arr[i] = - Phi_arr[i] - 12.0*pow(H0/(Constants.c*k*exp(x)),2.0)*(OmegaR0*Theta_arr[2][i] + OmegaNu0*Nu_arr[2][i]);
      }
      else{
        // If no neutrinos then record Psi as
        Psi_arr[i] = - Phi_arr[i] - 12.0*pow(H0/(Constants.c*k*exp(x)),2.0)*(OmegaR0*Theta_arr[2][i]);
      }
   }
   
    //===================================================================
    // TODO: remember to store the data found from integrating so we can
    // spline it below
    //
    // To compute a 2D spline of a function f(x,k) the data must be given 
    // to the spline routine as a 1D array f_array with the points f(ix, ik) 
    // stored as f_array[ix + n_x * ik]
    // Example:
    // Vector x_array(n_x);
    // Vector k_array(n_k);
    // Vector f(n_x * n_k);
    // Spline2D y_spline;
    // f_spline.create(x_array, k_array, f_array);
    // We can now use the spline as f_spline(x, k)
    //
    // NB: If you use Theta_spline then you have to allocate it first,
    // before using it e.g.
    // Theta_spline = std::vector<Spline2D>(n_ell_theta);
    //
    //===================================================================
    //...
    //...

  }
  Utils::EndTiming("integrateperturbation");

  //=============================================================================
  // TODO: Make all splines needed: Theta0,Theta1,Theta2,Phi,Psi,...
  //=============================================================================
  // ...
  // ...

  Vector x_array = Utils::linspace(x_start, x_end, n_x-1);

  // Spline scalar values
  delta_cdm_spline.create(x_array, k_array, delta_cdm_arr);
  delta_b_spline.create(x_array, k_array, delta_b_arr);
  v_cdm_spline.create(x_array, k_array, v_cdm_arr);
  v_b_spline.create(x_array, k_array, v_b_arr); 
  Psi_spline.create(x_array, k_array, Psi_arr);
  Phi_spline.create(x_array, k_array, Phi_arr);
  Pi_spline.create(x_array, k_array, Pi_arr);

  // Spline the big boys
  Theta_spline  = std::vector<Spline2D>(Constants.n_ell_theta);
  Theta_p_spline = std::vector<Spline2D>(Constants.n_ell_thetap);
  Nu_spline     = std::vector<Spline2D>(Constants.n_ell_neutrinos);

  // Spline theta
  for(int l = 0; l<Constants.n_ell_theta; l++){
    Theta_spline[l].create(x_array, k_array, Theta_arr[l]);
  }
  // If polarization is enabled spline this
  if(Constants.polarization){
    for(int l = 0; l<Constants.n_ell_theta; l++){
      Theta_p_spline[l].create(x_array, k_array, Thetap_arr[l]);
    }
  }
  // If neutrinos are enabled spline this
  if(Constants.neutrinos){
    for(int l = 0; l<Constants.n_ell_theta; l++){
      Nu_spline[l].create(x_array, k_array, Nu_arr[l]);
    }
  }
}

//====================================================
// Set IC at the start of the run (this is in the
// tight coupling regime)
//====================================================
Vector Perturbations::set_ic(const double x, const double k) const{

  // The vector we are going to fill
  Vector y_tc(Constants.n_ell_tot_tc);

  //=============================================================================
  // Compute where in the y_tc array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const int n_ell_tot_tc        = Constants.n_ell_tot_tc;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // References to the tight coupling quantities
  double &delta_cdm    =  y_tc[Constants.ind_deltacdm_tc];
  double &delta_b      =  y_tc[Constants.ind_deltab_tc];
  double &v_cdm        =  y_tc[Constants.ind_vcdm_tc];
  double &v_b          =  y_tc[Constants.ind_vb_tc];
  double &Phi          =  y_tc[Constants.ind_Phi_tc];
  double *Theta        = &y_tc[Constants.ind_start_theta_tc];
  double *Nu           = &y_tc[Constants.ind_start_nu_tc];

  //=============================================================================
  // TODO: Set the initial conditions in the tight coupling regime
  //=============================================================================
  // ...
  // ...

  double fnu    = OmegaNu0/(OmegaR0 + OmegaNu0);
  double Psi;
  if(neutrinos){
    Psi = - 1.0 / (3.0/2.0 + 2.0*fnu/5.0);
  }
  else{
    Psi = - 2.0/3.0;
  }
  double ckhp   = Constants.c*k/cosmo->Hp_of_x(x);
  double dtaudx = rec->dtaudx_of_x(x);

  // SET: Scalar quantities (Gravitational potential, baryons and CDM)
  if(neutrinos){
    Phi       = - (1.0 + 2.0/5.0*fnu*neutrinos)*Psi;
  }
  else{
    Phi       = - Psi;
  }
  delta_cdm = - 3.0/2.0 * Psi;
  v_cdm     = - ckhp/2.0*Psi;
  delta_b   = - 3.0/2.0 * Psi;
  v_b       = - ckhp/2.0*Psi;

  // SET: Photon temperature perturbations (Theta_ell)
  Theta[0] = - Psi/2.0;
  Theta[1] =  ckhp/6.0*Psi;

  // SET: Neutrino perturbations (N_ell)
  if(neutrinos){
    Nu[0] = -0.5*Psi;
    Nu[1] = ckhp/6.0*Psi;
    Nu[2] = pow(Constants.c*k*exp(x)/H0,2)*(Phi+Psi)/(12*OmegaNu0);
    for (int l = 3; l < n_ell_neutrinos_tc; l++){
      Nu[l] = ckhp / (2.0 * l + 1.0) * Nu[l-1];
    }
  }
  return y_tc;
}

//====================================================
// Set IC for the full ODE system after tight coupling 
// regime ends
//====================================================

Vector Perturbations::set_ic_after_tight_coupling(
    const Vector &y_tc, 
    const double x, 
    const double k) const{

  // Make the vector we are going to fill
  Vector y(Constants.n_ell_tot_full);
  
  //=============================================================================
  // Compute where in the y array each component belongs and where corresponding
  // components are located in the y_tc array
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Number of multipoles we have in the full regime
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // Number of multipoles we have in the tight coupling regime
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;

  // References to the tight coupling quantities
  const double &delta_cdm_tc    =  y_tc[Constants.ind_deltacdm_tc];
  const double &delta_b_tc      =  y_tc[Constants.ind_deltab_tc];
  const double &v_cdm_tc        =  y_tc[Constants.ind_vcdm_tc];
  const double &v_b_tc          =  y_tc[Constants.ind_vb_tc];
  const double &Phi_tc          =  y_tc[Constants.ind_Phi_tc];
  const double *Theta_tc        = &y_tc[Constants.ind_start_theta_tc];
  const double *Nu_tc           = &y_tc[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set
  double &delta_cdm       =  y[Constants.ind_deltacdm];
  double &delta_b         =  y[Constants.ind_deltab];
  double &v_cdm           =  y[Constants.ind_vcdm];
  double &v_b             =  y[Constants.ind_vb];
  double &Phi             =  y[Constants.ind_Phi];
  double *Theta           = &y[Constants.ind_start_theta];
  double *Theta_p         = &y[Constants.ind_start_thetap];
  double *Nu              = &y[Constants.ind_start_nu];

  //=============================================================================
  // TODO: fill in the initial conditions for the full equation system below
  // NB: remember that we have different number of multipoles in the two
  // regimes so be careful when assigning from the tc array
  //=============================================================================
  // ...
  // ...
  // ...
  double ckhp   = Constants.c*k/cosmo->Hp_of_x(x);
  double dtaudx = rec->dtaudx_of_x(x);

  // SET: Scalar quantities (Gravitational potental, baryons and CDM)
  Phi         = Phi_tc;
  delta_cdm   = delta_cdm_tc;
  v_cdm       = v_cdm_tc;
  delta_b     = delta_b_tc;
  v_b         = v_b_tc;

  // SET: Photon temperature perturbations (Theta_ell)
  Theta[0] = Theta_tc[0];
  Theta[1] = Theta_tc[1];
  if (polarization){
    Theta[2] = -8.0*ckhp/(15.0*dtaudx)*Theta[1]; 
  }
  else{
    Theta[2] = -20.0*ckhp/(45.0*dtaudx)*Theta[1];
  }
  for (int l=3;l<n_ell_theta;l++){
    Theta[l] = -l/(2.0*l + 1.0)*ckhp/dtaudx*Theta[l-1];
  }

  // SET: Photon polarization perturbations (Theta_p_ell)
  if(polarization){
    Theta_p[0] = 5.0*Theta[2]/4.0;
    Theta_p[1] = -ckhp/(4.0*dtaudx)*Theta[2];
    Theta_p[2] = 0.25*Theta[2];
    for (int l = 3; l<n_ell_thetap; l++){
      Theta_p[l] = -l/(2.0*l + 1.0)*ckhp/dtaudx*Theta_p[l-1];
    }
  }

  // SET: Neutrino perturbations (N_ell)
  if(neutrinos){
    for (int l = 0; l<n_ell_neutrinos; l++){
      Nu[l] = Nu_tc[l];
    }
  }
  return y;
}

//====================================================
// The time when tight coupling end
//====================================================

std::pair<double,int> Perturbations::get_tight_coupling_time(const double k) const{
  /* 
  Compute and return x for when tight coupling ends
  Remember all the three conditions in Callin
  */
  double x_tight_coupling_end = 0.0;
  double x_rec = log(1.0 / (1.0 + 2000.0)); //-7.1625; // -8.0; // A little before recombination
  Vector x = Utils::linspace(x_start, x_end, n_x);

  for(int i=0; i<n_x; i++){
    /*
    if (abs(rec->dtaudx_of_x(x[i])) < 10*std::min(1.0, Constants.c*k/cosmo->Hp_of_x(x[i])))
      return x[i];
    */
    if(abs(Constants.c*k/(cosmo->Hp_of_x(x[i]) * rec->dtaudx_of_x(x[i]))) > 0.1){
      return std::pair<double,int>(x[i-1], i-1);
    }
    else if(abs(rec->dtaudx_of_x(x[i])) < 10){
      return std::pair<double,int>(x[i-1], i-1);
    }
    else if (x[i] > x_rec)
      return std::pair<double,int>(x[i-1], i-1);
  }

  std::cout << "Could not find x at end of tight coupling" << std::endl;
  return std::pair<double,int>(x_tight_coupling_end,0);
}

//====================================================
// After integrsating the perturbation compute the
// source function(s)
//====================================================
void Perturbations::compute_source_functions(){
  Utils::StartTiming("source");

  //=============================================================================
  // TODO: Make the x and k arrays to evaluate over and use to make the splines
  //=============================================================================
  // ...
  // ...
  Vector k_array;
  Vector x_array;

  // Make storage for the source functions (in 1D array to be able to pass it to the spline)
  Vector ST_array(k_array.size() * x_array.size());
  Vector SE_array(k_array.size() * x_array.size());

  // Compute source functions
  for(auto ix = 0; ix < x_array.size(); ix++){
    const double x = x_array[ix];
    for(auto ik = 0; ik < k_array.size(); ik++){
      const double k = k_array[ik];

      // NB: This is the format the data needs to be stored 
      // in a 1D array for the 2D spline routine source(ix,ik) -> S_array[ix + nx * ik]
      const int index = ix + n_x * ik;

      //=============================================================================
      // TODO: Compute the source functions
      //=============================================================================
      // Fetch all the things we need...
      // const double Hp       = cosmo->Hp_of_x(x);
      // const double tau      = rec->tau_of_x(x);
      // ...
      // ...

      // Temperatur source
      ST_array[index] = 0.0;

      // Polarization source
      if(Constants.polarization){
        SE_array[index] = 0.0;
      }
    }
  }

  // Spline the source functions
  ST_spline.create (x_array, k_array, ST_array, "Source_Temp_x_k");
  if(Constants.polarization){
    SE_spline.create (x_array, k_array, SE_array, "Source_Pol_x_k");
  }

  Utils::EndTiming("source");
}

//====================================================
// The right hand side of the perturbations ODE
// in the tight coupling regime
//====================================================

// Derivatives in the tight coupling regime
int Perturbations::rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx){

  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  const double &delta_b         =  y[Constants.ind_deltab_tc];
  const double &v_cdm           =  y[Constants.ind_vcdm_tc];
  const double &v_b             =  y[Constants.ind_vb_tc];
  const double &Phi             =  y[Constants.ind_Phi_tc];
  const double *Theta           = &y[Constants.ind_start_theta_tc];
  const double *Nu              = &y[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm_tc];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab_tc];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm_tc];
  double &dv_bdx          =  dydx[Constants.ind_vb_tc];
  double &dPhidx          =  dydx[Constants.ind_Phi_tc];
  double *dThetadx        = &dydx[Constants.ind_start_theta_tc];
  double *dNudx           = &dydx[Constants.ind_start_nu_tc];

  //=============================================================================
  // TODO: fill in the expressions for all the derivatives
  //=============================================================================

  // Some values given x
  double Hp       = cosmo->Hp_of_x(x);
  double ckhp     = Constants.c*k/Hp;
  double R        = 4.0*OmegaR0/(3.0*OmegaB0)*exp(-x);
  double dtaudx   = rec->dtaudx_of_x(x);

  double Theta2;
  if(polarization){
    Theta2   = -8.0*ckhp/(15.0*dtaudx)*Theta[1];
  }
  else{
    Theta2 = -20.0*ckhp/(45.0*dtaudx)*Theta[1];
  }

  double Psi;
  if(neutrinos){
    Psi      = - Phi - 12.0*pow(H0/(Constants.c*k*exp(x)),2.0)*(OmegaR0*Theta2 + OmegaNu0*Nu[2]);
    dPhidx   = Psi - pow(ckhp,2)/3.0*Phi + pow(H0/Hp, 2.0)/2.0*(OmegaCDM0*exp(-x)*delta_cdm 
    + OmegaB0*exp(-x)*delta_b + 4.0*OmegaR0*exp(-2*x)*Theta[0] + 4.0*OmegaNu0*exp(-2*x)*Nu[0]);
  }
  else{
    Psi      = - Phi - 12.0*pow(H0/(Constants.c*k*exp(x)),2.0)*(OmegaR0*Theta2);
    dPhidx   = Psi - pow(ckhp,2)/3.0*Phi + pow(H0/Hp, 2.0)/2.0*(OmegaCDM0*exp(-x)*delta_cdm 
    + OmegaB0*exp(-x)*delta_b + 4.0*OmegaR0*exp(-2*x)*Theta[0]);
  }
  dThetadx[0] = - ckhp*Theta[1] - dPhidx;

  // SET: Scalar quantities (Phi, delta, v, ...)
  double q = (-((1.0 - R ) * dtaudx + (1.0 + R) * rec->ddtauddx_of_x(x)) * (3.0 * Theta[1] + v_b)
    - ckhp * Psi + (1.0 - cosmo->dHpdx_of_x(x) / Hp) * ckhp * (-Theta[0] + 2.0 * Theta2)
    - ckhp * dThetadx[0]) / ((1.0 + R) * dtaudx + cosmo->dHpdx_of_x(x) / Hp - 1.0);

  ddelta_cdmdx  = ckhp * v_cdm - 3.0 * dPhidx;
  dv_cdmdx      = - v_cdm - ckhp * Psi;
  ddelta_bdx    = ckhp * v_b - 3.0 * dPhidx;
  dv_bdx        = 1.0 / (1.0 + R) * (- v_b - ckhp * Psi + R*(q+ckhp*(-Theta[0]+2.0*Theta2)-ckhp*Psi));

  // SET: Photon multipoles (Theta_ell)
  dThetadx[1]   = 1.0/3.0*(q-dv_bdx);

  // SET: Neutrino mutlipoles (Nu_ell)
  if(neutrinos){
    dNudx[0] = -ckhp*Nu[1] - dPhidx;
    dNudx[1] = ckhp/3.0*Nu[0] - 2.0*ckhp/3*Nu[2] + ckhp/3.0*Psi;
    for (int l = 2; l<n_ell_neutrinos_tc-1; l++){
      dNudx[l] = l*ckhp/(2.0*l + 1.0)*Nu[l-1] - (l + 1.0) * ckhp/(2.0 * l + 1.0)*Nu[l+1];
    }
    int l = n_ell_neutrinos_tc - 1;
    dNudx[l] = ckhp*Nu[l-1] - Constants.c*(l+1.0)/(Hp*cosmo->eta_of_x(x))*Nu[l];
  }

  return GSL_SUCCESS;
}

//====================================================
// The right hand side of the full ODE
//====================================================

int Perturbations::rhs_full_ode(double x, double k, const double *y, double *dydx){
  
  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Index and number of the different quantities
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm];
  const double &delta_b         =  y[Constants.ind_deltab];
  const double &v_cdm           =  y[Constants.ind_vcdm];
  const double &v_b             =  y[Constants.ind_vb];
  const double &Phi             =  y[Constants.ind_Phi];
  const double *Theta           = &y[Constants.ind_start_theta];
  const double *Theta_p         = &y[Constants.ind_start_thetap];
  const double *Nu              = &y[Constants.ind_start_nu];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm];
  double &dv_bdx          =  dydx[Constants.ind_vb];
  double &dPhidx          =  dydx[Constants.ind_Phi];
  double *dThetadx        = &dydx[Constants.ind_start_theta];
  double *dTheta_pdx      = &dydx[Constants.ind_start_thetap];
  double *dNudx           = &dydx[Constants.ind_start_nu];

  // Cosmological parameters and variables
  double Hp       = cosmo->Hp_of_x(x);
  double R        = 4.0 * OmegaR0 / (3.0 * OmegaB0) * exp(-x);
  double dtaudx   = rec->dtaudx_of_x(x);
  double ckhp     = Constants.c * k / Hp;
  double Pi; 
  if(polarization){
    Pi       = Theta[2] + (Theta_p[0] + Theta_p[2]);
  }
  else{
    Pi       = Theta[2];
  }
  //=============================================================================
  // TODO: fill in the expressions for all the derivatives
  //=============================================================================
  double Psi;
  if(neutrinos){
    Psi = - Phi - 12.0*pow(H0/(Constants.c*k*exp(x)),2)*(OmegaR0*Theta[2] + OmegaNu0*Nu[2]);
    dPhidx        = Psi - pow(ckhp,2.0)/3.0*Phi 
      + pow(H0/Hp, 2.0)/2*(OmegaCDM0*exp(-x)*delta_cdm 
      + OmegaB0*exp(-x)*delta_b + 4*OmegaR0*exp(-2*x)*Theta[0] 
      + 4*OmegaNu0*exp(-2*x)*Nu[0]);
  }
  else{
    Psi = - Phi - 12.0*pow(H0/(Constants.c*k*exp(x)),2)*(OmegaR0*Theta[2]);
    dPhidx        = Psi - pow(ckhp,2.0)/3*Phi 
      + pow(H0/Hp, 2.0)/2*(OmegaCDM0*exp(-x)*delta_cdm 
      + OmegaB0*exp(-x)*delta_b + 4*OmegaR0*exp(-2*x)*Theta[0]);
  }

  // SET: Scalar quantities (Phi, delta, v, ...)
  ddelta_cdmdx  = ckhp * v_cdm - 3.0 * dPhidx;
  dv_cdmdx      = - v_cdm - ckhp * Psi;
  ddelta_bdx    = ckhp * v_b - 3.0 * dPhidx;
  dv_bdx        = - v_b - ckhp * Psi + dtaudx * R * (3.0 * Theta[1] + v_b);
 
  // SET: Photon multipoles (Theta_ell)
  dThetadx[0] = - ckhp*Theta[1] - dPhidx;
  dThetadx[1] = ckhp/3*Theta[0] - 2.0*ckhp/3*Theta[2] + ckhp/3.0*Psi + dtaudx*(Theta[1] + v_b/3.0);
  for (int l=2;l<n_ell_theta-1;l++){
    dThetadx[l] = l*ckhp/(2.0*l + 1.0)*Theta[l-1] - (l + 1.0)*ckhp/(2.0*l + 1.0)*Theta[l+1] 
    + dtaudx*(Theta[l] - 0.1*Pi*(l == 2)); 
  }
  int l = n_ell_theta - 1;
  dThetadx[l] = ckhp*Theta[l-1] - Constants.c*(l+1)/(Hp*cosmo->eta_of_x(x))*Theta[l] + dtaudx*Theta[l];

  // SET: Photon polarization multipoles (Theta_p_ell)
  if(polarization){
    dTheta_pdx[0] = -ckhp*Theta_p[1] + dtaudx*(Theta_p[0] - 0.5*Pi);
    for (int l=1;l<n_ell_thetap-1;l++){
      dTheta_pdx[l] = l*ckhp/(2.0*l+1.0)*Theta_p[l-1]
        - (l+1.0)*ckhp/(2.0*l+1.0)*Theta_p[l+1] 
        + dtaudx*(Theta_p[l] - 0.1*Pi*(l == 2));
    }
    int l = n_ell_thetap - 1; 
    dTheta_pdx[l] = ckhp*Theta_p[l-1] - Constants.c*(l+1)/(Hp*cosmo->eta_of_x(x))*Theta_p[l] + dtaudx*Theta_p[l];
  }

  // SET: Neutrino mutlipoles (Nu_ell)
  if(neutrinos){
    dNudx[0] = -ckhp*Nu[1] - dPhidx;
    dNudx[1] = ckhp/3.0*Nu[0] - 2.0*ckhp/3*Nu[2] + ckhp/3.0*Psi;
    for (int l=2;l<n_ell_neutrinos-1;l++){
      dNudx[l] = l*ckhp/(2.0*l+1.0)*Nu[l-1] - (l+1.0)*ckhp/(2.0*l+1.0)*Nu[l+1];
    }
    int l = n_ell_neutrinos -1;
    dNudx[l] = ckhp*Nu[l-1] - Constants.c*(l+1)/(Hp*cosmo->eta_of_x(x))*Nu[l];
  }

  return GSL_SUCCESS;
}

//====================================================
// Get methods
//====================================================

double Perturbations::get_delta_cdm(const double x, const double k) const{
  return delta_cdm_spline(x,k);
}
double Perturbations::get_delta_b(const double x, const double k) const{
  return delta_b_spline(x,k);
}
double Perturbations::get_v_cdm(const double x, const double k) const{
  return v_cdm_spline(x,k);
}
double Perturbations::get_v_b(const double x, const double k) const{
  return v_b_spline(x,k);
}
double Perturbations::get_Phi(const double x, const double k) const{
  return Phi_spline(x,k);
}
double Perturbations::get_Psi(const double x, const double k) const{
  return Psi_spline(x,k);
}
double Perturbations::get_Pi(const double x, const double k) const{
  return Pi_spline(x,k);
}
double Perturbations::get_Source_T(const double x, const double k) const{
  return ST_spline(x,k);
}
double Perturbations::get_Source_E(const double x, const double k) const{
  return SE_spline(x,k);
}
double Perturbations::get_Theta(const double x, const double k, const int ell) const{
  return Theta_spline[ell](x,k);
}
double Perturbations::get_Theta_p(const double x, const double k, const int ell) const{
  return Theta_p_spline[ell](x,k);
}
double Perturbations::get_Nu(const double x, const double k, const int ell) const{
  return Nu_spline[ell](x,k);
}

//====================================================
// Print some useful info about the class
//====================================================

void Perturbations::info() const{
  std::cout << "\n";
  std::cout << "Info about perturbations class:\n";
  std::cout << "x_start:       " << x_start                << "\n";
  std::cout << "x_end:         " << x_end                  << "\n";
  std::cout << "n_x:     " << n_x              << "\n";
  std::cout << "k_min (1/Mpc): " << k_min * Constants.Mpc  << "\n";
  std::cout << "k_max (1/Mpc): " << k_max * Constants.Mpc  << "\n";
  std::cout << "n_k:     " << n_k              << "\n";
  if(Constants.polarization)
    std::cout << "We include polarization\n";
  else
    std::cout << "We do not include polarization\n";
  if(Constants.neutrinos)
    std::cout << "We include neutrinos\n";
  else
    std::cout << "We do not include neutrinos\n";

  std::cout << "Information about the perturbation system:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm         << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab           << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm             << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb               << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi              << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta      << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta          << "\n";
  if(Constants.polarization){
    std::cout << "ind_start_thetap:   " << Constants.ind_start_thetap   << "\n";
    std::cout << "n_ell_thetap:       " << Constants.n_ell_thetap       << "\n";
  }
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu       << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos    << "\n";
  }
  std::cout << "n_ell_tot_full:     " << Constants.n_ell_tot_full       << "\n";

  std::cout << "Information about the perturbation system in tight coupling:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm_tc      << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab_tc        << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm_tc          << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb_tc            << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi_tc           << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta_tc   << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta_tc       << "\n";
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu_tc    << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos_tc << "\n";
  }
  std::cout << "n_ell_tot_tc:       " << Constants.n_ell_tot_tc         << "\n";
  std::cout << std::endl;
}

//====================================================
// Output some results to file for a given value of k
//====================================================

void Perturbations::output(const double k, const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts = 5000;
  auto x_array = Utils::linspace(x_start, x_end, npts);
  auto print_data = [&] (const double x) {
    double arg = k * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));
    fp << x                  << " ";
    fp << get_delta_cdm(x,k) << " ";
    fp << get_delta_b(x,k)   << " ";
    fp << get_v_cdm(x,k)     << " ";
    fp << get_v_b(x,k)       << " ";
    
    fp << get_Theta(x,k,0)   << " ";
    fp << get_Theta(x,k,1)   << " ";
    fp << get_Theta(x,k,2)   << " ";

    if(Constants.polarization){
      fp << get_Theta_p(x,k,0) << " ";
      fp << get_Theta_p(x,k,1) << " ";
      fp << get_Theta_p(x,k,2) << " ";
    }
    else{
      fp << 0 << " ";
      fp << 0 << " ";
      fp << 0 << " ";
    }

    if(Constants.neutrinos){
      fp << get_Nu(x,k,0)   << " ";
      fp << get_Nu(x,k,1)   << " ";
      fp << get_Nu(x,k,2)   << " ";
    }
    else{
      fp << 0 << " ";
      fp << 0 << " ";
      fp << 0 << " ";
    }

    fp << get_Phi(x,k)       << " ";
    fp << get_Psi(x,k)       << " ";
    fp << get_Pi(x,k)        << " ";

    /*
    fp << get_Source_T(x,k)  << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(5,   arg)           << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(50,  arg)           << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(500, arg)           << " ";
    */
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

