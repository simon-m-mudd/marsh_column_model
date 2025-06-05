//dParticle.cpp
// implementation of the dParticle object

#include "depo_particle_info.hpp"
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

/**************************************************************************\
**
**  dParticle::dParticle:  Constructor for particles.
**			   The intial age of the particle is zero
**			   The location is assigned
\**************************************************************************/
void depo_particle_info::create()
 {
  //cout << "depo_particle_info.cpp LINE 18 error, you need to have a filename\n";
 }

// this function loads the parameters from the file with the name fname
void depo_particle_info::create( const char* fname )
 {
  ifstream in;
  in.open(fname);


  vector<int> t_type_index;
  vector<string> t_type_name;
  vector<double> t_e_0_vec;		// holds the void ratio at reference stress
  vector<double> t_compression_index_vec;
  					// compression index
  vector<double> t_rho_vec;		// solid density
  vector<double> t_mass_vec;		// mass of the particle
  vector<double> t_k_0_vec;		// decay rates of the material
  vector<double> t_gamma_vec;		// length scales of decay
  vector<double> t_sigma_0_vec;		// reference stress

  string t_ptype;
  int temp_ti;
  int n_t;
  double t_e0,t_compression_index,t_sigma_0,t_rho,t_mass,t_k_0,t_gamma;

  in >> n_t;

  n_types = n_t;
  for (int i = 1; i<=n_t; i++)
   {
    in >> temp_ti >> t_ptype >> t_e0 >> t_sigma_0 >> t_compression_index >> t_rho >> t_mass >> t_k_0 >> t_gamma;
    t_type_index.push_back(temp_ti);
    t_type_name.push_back(t_ptype);
    t_e_0_vec.push_back(t_e0);
    t_sigma_0_vec.push_back(t_sigma_0);
    t_compression_index_vec.push_back(t_compression_index);
    t_rho_vec.push_back(t_rho);
    t_mass_vec.push_back(t_mass);
    t_k_0_vec.push_back(t_k_0);
    t_gamma_vec.push_back(t_gamma);

   }
  type_index = t_type_index;
  type_name = t_type_name;
  e_0_vec = t_e_0_vec;
  compression_index_vec = t_compression_index_vec;
  sigma_0_vec = t_sigma_0_vec;
  rho_vec = t_rho_vec;
  initial_mass_vec = t_mass_vec;
  k_0_vec = t_k_0_vec;
  gamma_vec = t_gamma_vec;
 }


string depo_particle_info::name_for_individual_type(int type)
 {
  if (type > n_types)
     {
      cout << "LINE 74 depo_particle_info.cpp: the type you requested doesn't exist\n";
     }
  return type_name[type];
 }

double depo_particle_info::rho_for_individual_type(int type)
 {
  if (type > n_types)
   {
    cout << "LINE 74 depo_particle_info.cpp: the type you requested doesn't exist\n";
   }
  return rho_vec[type];
 }

double depo_particle_info::initial_mass_for_individual_type(int type)
 {
  if (type > n_types)
   {
    cout <<"LINE 74 depo_particle_info.cpp: the type you requested doesn't exist\n";
   }
  return initial_mass_vec[type];
 }


double depo_particle_info::sigma_0_for_individual_type(int type)
 {
  if (type > n_types)
   {
    cout <<"LINE 74 depo_particle_info.cpp: the type you requested doesn't exist\n";
   }
  return sigma_0_vec[type];
 }

double depo_particle_info::e_0_for_individual_type(int type)
 {
  if (type > n_types)
   {
    cout <<"LINE 74 depo_particle_info.cpp: the type you requested doesn't exist\n";
   }
  return e_0_vec[type];
 }

double depo_particle_info::compression_index_for_individual_type(int type)
 {
  if (type > n_types)
   {
    cout <<"LINE 74 depo_particle_info.cpp: the type you requested doesn't exist\n";
   }
  return compression_index_vec[type];
 }

