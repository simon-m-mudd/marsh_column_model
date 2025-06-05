//sediment_layer2.h
// header file for the sediment layer object
// a sediment layer is deposited in a single model timestep. It contains:
// a) inorganic sediment (varios particle sizes)
// b) roots (live biomass)
// c) labile organic carbon
// d) refractory organic carbon
// e) 210Pb
// f) 137Cs
// other components can be added later

#include<iostream>
#include "depo_particle.hpp"
#include "depo_particle_info.hpp"
#include<vector>
using namespace std;

#ifndef SEDIMENT_LAYER_H
#define SEDIMENT_LAYER_H

class sediment_layer
{
 public:
    sediment_layer()					{ create(); }
    sediment_layer(int n_p_types, depo_particle_info& DPI)
    							{ create(n_p_types, DPI); }
    sediment_layer(vector<depo_particle>& dp_vec, depo_particle_info& DPI)
    							{ create(dp_vec, DPI); }
    sediment_layer(depo_particle_info& DPI, vector<double> dpart_masses)
     							{ create(DPI, dpart_masses); }


    vector<double> get_mass_of_particle_types()		{ return mass_of_particle_types; }
    int get_n_p_types()					{ return n_p_types; }

    double calculate_buoyant_weight();			// get the buoyant weight for effective stress calculation
    double get_total_layer_solid_mass();		// gets the total solid mass of the layer
    double get_mass_of_particle_type(int type);		// returns the mass of particle of type 'type'
    void decay_and_mass_loss(double dt, double depth, double temperatureincrease, double kfactor); // decay of carbon and radioisotopes. MK- I added temperatureincrease and kfactor
    void root_growth_and_death(int root_type, vector<int> carbon_types,
                               double new_root_mass, vector<double> delta_carbon_masses);
                               				// this function updates root and carbon pools. This function deals explicitly
                               				// with masses and not rates. If the rate is known, the mass change must
                               				// be computed extenally of this function.
                               				// the function is designed to be flexible
    vector<double> instant_compaction_compression_index(double eff_stress);
    							// compacts layer, returns a vector with:
    							// element[0] is the porosity, n
    							// element[1] is the layer thickness

    void print_particle_info_to_screen();	// prints the information about each individual
    										// particle within the layer


 private:
    int n_p_types;				// the number of particel types in the layer. This wille
    							// depend on the number of inoganic particle size fractions
    							// it is also flexible such that later iterations of the model can
    							// add additional carbon pools or radioisotopes
    vector <  depo_particle >   particle_vector;
    							// this vecor contains the particles that make up the layer. The vector
    							// is of dimension n_p_types
    vector<double> mass_of_particle_types;
    							// this vector stores the mass of each particle type
    depo_particle_info dpi;		// the particle info object, this gives the properties of each of the
    							// particles contained within the layer

    void create();
    void create(int,depo_particle_info&);
    void create(vector<depo_particle>&, depo_particle_info&);
    void create(depo_particle_info&, vector<double>);


};

#endif
