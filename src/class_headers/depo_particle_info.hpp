//depo_particle_info.h
// header file for the discrete particle object
// this version of the dicrete particle keeps track
// of two integers, the time and the position

#include<iostream>
#include<vector>
#include <string>
using namespace std;


#ifndef depo_particle_info_H
#define depo_particle_info_H

class depo_particle_info
{
 public:
    depo_particle_info()				{ create(); }
    depo_particle_info( const char* fname)          	{ create(fname); }

    int get_n_types()					{ return n_types; }
    vector<int> get_type_index()			{ return type_index; }
    vector<string> get_type_name()			{ return type_name; }
    vector<double> get_e_0_vec()			{ return e_0_vec; }
    vector<double> get_compression_index_vec()		{ return compression_index_vec; }
    vector<double> get_sigma_0_vec()			{ return sigma_0_vec; }
    vector<double> get_rho_vec()			{ return rho_vec; }
    vector<double> get_initial_mass_vec	()		{ return initial_mass_vec; }
    vector<double> get_k_0_vec()			{ return k_0_vec; }
    vector<double> get_gamma_vec()			{ return gamma_vec; }

    string name_for_individual_type(int type);
    double rho_for_individual_type(int type);
    double initial_mass_for_individual_type(int type);
    double sigma_0_for_individual_type(int type);
    double e_0_for_individual_type(int type);
    double compression_index_for_individual_type(int type);

 private:
    int n_types;
    vector<int> type_index;
    vector<string> type_name;
    vector<double> e_0_vec;		// holds the porosity at zero overburden
    vector<double> compression_index_vec;
    					// compression index
    vector<double> sigma_0_vec;		// reference stress
    vector<double> rho_vec;		// solid density
    vector<double> initial_mass_vec;	// mass of the particle
    vector<double> k_0_vec;		// decay rates of the material
    vector<double> gamma_vec;		// length scales of decay


    void create();
    void create(const char*);

};

#endif
