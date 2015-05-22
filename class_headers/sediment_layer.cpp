//sediment_layer.cpp
// implementation of the dParticle object
#include <fstream>
#include <iostream>
#include <math.h>
#include "sediment_layer.hpp"
#include "depo_particle.hpp"
#include "depo_particle_info.hpp"
using namespace std;


/**************************************************************************\
**
**  sediment_layer::sediment_layer:  Constructor for the layers.
**			   The intial age of the particle is zero
**			   The location is assigned
\**************************************************************************/
void sediment_layer::create()
 {

 }


// this initializes a sediment layer with empty particles
void sediment_layer::create(int npt,depo_particle_info& DPI)
 {
  dpi = DPI;
  n_p_types = npt;
  vector<depo_particle> dp_temp(npt);
  particle_vector = dp_temp;
  for (int i = 0; i<n_p_types; i++)
   {
    mass_of_particle_types.push_back(0.0);
   }
 }

// this initializes a sediment layer with known particles
void sediment_layer::create(vector<depo_particle>& dp_vec, depo_particle_info& DPI)
 {
  dpi = DPI;
  n_p_types = dp_vec.size();
  particle_vector = dp_vec;
  vector<double> mpv;
  for (int i = 0; i<n_p_types; i++)
   {
    mpv.push_back( particle_vector[i].getMass() );
   }
  mass_of_particle_types = mpv;
 }

// this constructor takes a depo_particle_info object and a vector of particle masses,
// then creates the particle_vec by assigning the masses to particles of each type
void sediment_layer::create(depo_particle_info& DPI, vector<double> dpart_masses)
 {

  dpi = DPI;
  // get the number of particles
  n_p_types = dpi.get_n_types();
  //cout <<"sediment_layer line 59 n_p_types is: " << n_p_types << endl;

  // this vector will be filled with the inital particles and then
  // assigned to particel_vector
  vector<depo_particle> temp_dpv;

  // check to see if the mass vec is the correct size
  int n_masses = dpart_masses.size();
  if (n_masses != n_p_types)
   {
    cout << "your mass vector does not match the number of particles\n" <<
             " AN EMPTY VECTOR IS USED INSTEAD \n";
    vector<double> temp(n_p_types,0.0);
   }
  else
   {
    // insert the particles
    for (int i = 0; i<n_p_types; i++)
     {
      depo_particle dp(i,dpart_masses[i]);
      temp_dpv.push_back(dp);
     }
   }

  // now assign the layer variables
  particle_vector = temp_dpv;
  mass_of_particle_types= dpart_masses;
 }


/**************************************************************************\
**
**  sediment_layer::calcualte_buoyant_weight:  this function calcualtes
** 			the total buoyant weight of the layer. This number is
**			then used for the porosity calculation
\**************************************************************************/
double sediment_layer::calculate_buoyant_weight()
 {
  double bouy_wgt =0;
  double rho_w = 1000;
  double g = 9.80;
  // loop through the particle types
  // the buoyant weight is (rho_p - rho_w)*V_p*g or
  //			   (rho_p - rho_w)*m/rho_p*g
  for (int i = 0; i<n_p_types; i++)
   {
    bouy_wgt += particle_vector[i].getMass()*g*(1-rho_w/dpi.rho_for_individual_type(i));
   }
  return bouy_wgt;
 }



/**************************************************************************\
**
**  sediment_layer::decay_and_mass_loss
** 			this function calculates decay (primarily of
**				organic carbon)
**		NOTE THE k VALUES MUST HAVE THE INVERSE UNITS OF
**		THE UNITS OF THE TIMESTEP, dt
\**************************************************************************/
void sediment_layer::decay_and_mass_loss(double dt, double mp_depth, double temperatureincrease, double kfactor)
 {
  vector<double> k_0_vec = dpi.get_k_0_vec();
  double k;
//  double kfactor=.25; //MK- This is the factor by which k increases per degree warming. Now trying to pass it from column_model()

  k_0_vec[3]=k_0_vec[3]+(kfactor*temperatureincrease*k_0_vec[3]);  // MK- I added this line to manipulate labile k0
  k_0_vec[2]=k_0_vec[2]+(kfactor*temperatureincrease*k_0_vec[2]);  // MK- I added this line to manipulate refractory k0

  //cout << "\n\nLine 74 sediment_layer.cpp, midpoint_depth: " << mp_depth << endl;
  for (int i = 0; i<n_p_types; i++)
   {
    if (k_0_vec[i] > 0)
     {

      vector<double> gamma_vec = dpi.get_gamma_vec();
      if (gamma_vec[i] == 0)
       k = k_0_vec[i];
      else
       k = k_0_vec[i]*exp(-mp_depth/gamma_vec[i]);
      //cout << "k_0 is: " << k_0_vec[i] << " and gamma is: " << gamma_vec[i] << " and k is "<< k << endl;
      double multiplier = 1-k*dt+0.5*dt*dt*k*k;  		// second order accurate decay
      //cout << "multiplier is: " << multiplier << endl;
      particle_vector[i].linear_decay_multiplier(multiplier);
      //cout << "mass before: " << mass_of_particle_types[i]
      //     << " and mass after: " << particle_vector[i].getMass() << endl;
      mass_of_particle_types[i] = particle_vector[i].getMass();
     }
   }

 }


/**************************************************************************\
**
**  sediment_layer::root_growth_and_death
** 			this function grows roots and transfers dead roots
**		to carbon pools. This function uses mass change explicity,
**		no rates are involved.
\**************************************************************************/
void sediment_layer::root_growth_and_death(int root_type, vector<int> carbon_types,
                               double new_root_mass, vector<double> delta_carbon_masses)
 {
  //depo_particle dp;
  //cout << "LINE 161 sed_layer.cpp initiating depo particle: " << endl;
  //dp.dp_properties_to_screen();
  //cout << endl;


  int n_carbon_types = carbon_types.size();
  //cout <<"LINE 167 sed_layer.cpp root type: " << root_type << endl;
  //for (int i = 0; i< n_carbon_types ; i++)
  // cout << "ctype["<<i<<"]: " << carbon_types[i] << endl;

  // change the root mass (if delta_root_mass is negative, this is mass loss
  //cout << "LINE 172 sediment_layer.cpp new root mass is: " << new_root_mass << endl;
  //cout << "LINE 173 particle_vector.size(): " << particle_vector.size() << endl;

  //cout << "LINE 175 sediment_layer.cpp particle data: " << endl;
  //print_particle_info_to_screen();

  particle_vector[root_type].updateMass(new_root_mass);
  //cout << "LINE 179 mass_of_particle_types.size(): " << mass_of_particle_types.size() << endl;
  //cout << "LINE 180 particle_vector.size(): " << particle_vector.size() << endl;



  //cout << "LINE 171 particle_vector[root_type].getMass(): " << particle_vector[root_type].getMass() << endl;

  //cout << particle_vector[0].getMass() << endl;
  mass_of_particle_types[root_type] = particle_vector[root_type].getMass();
  //cout << "LINE 171 sediment_layer.cpp mass of roots: " << mass_of_particle_types[root_type] << endl;


  int n_carbon_pools = carbon_types.size();
  for (int i = 0; i< n_carbon_pools; i++)
   {
    //cout << "LINE 175 sediment_layer.cpp pool number: " << i << endl;
    particle_vector[ carbon_types[i] ].changeMass( delta_carbon_masses[i] );
    mass_of_particle_types[carbon_types[i]] = particle_vector[carbon_types[i]].getMass();
   }
  //cout << "LINE 176 sediment_later.cpp done!"<< endl;
 }



/**************************************************************************\
**
**  sediment_layer::get_total_layer_mass
** 			gets the total mass of the layer
**
\**************************************************************************/
double sediment_layer::get_total_layer_solid_mass()
 {
  double tlm = 0;
  //cout << "n_p_t is: " << n_p_types << endl;
  for (int i = 0; i<n_p_types; i++)
   {
    //cout << " type: " <<  i << " mass: " << mass_of_particle_types[i] << endl;
    tlm += mass_of_particle_types[i];
   }
  return tlm;
 }

double sediment_layer::get_mass_of_particle_type(int type)
 {
  if (type > n_p_types)
   {
    cout << "LINE 85 sediment layer.cpp, you have entered an incorrect particle type\n";
   }
  return mass_of_particle_types[type];
 }


/**************************************************************************\
**
**  sediment_layer::instant_compaction_compression_index(double overburden,
                                                         double hydrostatic)
** 			compacts the layer using the compression index
**			(e.g., Lamb and Whitman 1930) relationship
** 			assumes excess pore pressures diffuse instantaneously
**			returns a vector with three values:
**			[0]: porosity
**			[1]: layer_thicknes
**
\**************************************************************************/
vector<double> sediment_layer::instant_compaction_compression_index
                                  (double eff_stress)
 {
  vector<double> p_solid_volumes(n_p_types, 0.0);
  vector<double> p_pore_volumes(n_p_types,0.0);
  vector<double> output_data(2,0.0);
  double n;
  double e;
  double layer_thickness = 0;
  double total_mass = 0;
  double total_pore_volume = 0.0;
  double total_solid_volume = 0.0;
  for(int i = 0; i<n_p_types; i++)
   {
    p_solid_volumes[i] = mass_of_particle_types[i]/dpi.rho_for_individual_type(i);
    if(eff_stress == 0)
     {
      e = dpi.e_0_for_individual_type(i)- dpi.compression_index_for_individual_type(i)*
                                log10(0.01/dpi.sigma_0_for_individual_type(i));
     }
    else
     {
      e = dpi.e_0_for_individual_type(i)- dpi.compression_index_for_individual_type(i)*
                                log10(eff_stress/dpi.sigma_0_for_individual_type(i));
     }

    n = e/(e+1);

    p_pore_volumes[i] = n*p_solid_volumes[i]/(1-n);
    total_pore_volume += p_pore_volumes[i];
    total_solid_volume += p_solid_volumes[i];
    layer_thickness += p_pore_volumes[i]+p_solid_volumes[i];
    total_mass += mass_of_particle_types[i];
    //cout << "n: " << n << endl;
   }
  //cout << "layer total mass is: " << total_mass << endl;
  output_data[0] = total_pore_volume/ (total_pore_volume+total_solid_volume);
  output_data[1] = layer_thickness;
  return output_data;
 }


/**************************************************************************\
**
**  sediment_layer::print_particle_info_to_screen();
** 			loops through teh particles in the layer giving
**			the properties of each
**
\**************************************************************************/
void sediment_layer::print_particle_info_to_screen()
 {
  int n_p = particle_vector.size();
  for (int i = 0; i< n_p; i++)
   {
	cout <<"particle number " << i << endl;
	particle_vector[i].dp_properties_to_screen();
	cout << endl;
   }
 }
