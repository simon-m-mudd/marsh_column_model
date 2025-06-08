#include "sediment_layer.hpp"
#include "depo_particle_info.hpp"
#include<vector>
#include<fstream>
using namespace std;

#ifndef SEDIMENT_STACK_H
#define SEDIMENT_STACK_H

class sediment_stack
{
 public:
    sediment_stack()					{ create(); }
    sediment_stack(double MSL)				{create(MSL); }

    void deposit_surface_sediment(depo_particle_info DPI, vector<double> dpart_masses);
    void deposit_surface_sediment(sediment_layer sl);
    vector<double> create_growth_index(double gamma);
    vector<double> create_growth_index_deep(double gamma, double efolding_multiplier,
                                                        double exp_multiplier);
    vector<double> create_growth_index_inf(double gamma);
    void root_growth_and_death(vector<double> growth_index, int root_type, vector<int> carbon_types,
                                           vector<double> chi_org, double BG_Bmass, double tot_mort_mass);

    void decay_and_mass_loss(double dt, double temperatureincrease, double kfactor); //MK- I added temperatureincrease and kfactor
    double get_marsh_surface_elevation();
    vector<double> get_total_stack_mass_of_each_type(int n_p_types);
    vector<double> get_stack_slice_mass_of_each_type(int n_p_types,
						double slice_top_depth, double slice_bottom_depth);
	vector<double> get_stack_slice_mass_of_each_type(int n_p_types,
 										int bottom_layer, int top_layer);
 	double get_layer_thickness(int layer);
 	double get_layer_top_elevation(int layer);
   int get_n_layers();
   void calculate_overburden_and_thickness();



   void print_layer_thickness_to_screen(double t_ime);
   void print_layer_mass_to_screen(double t_ime);
   void print_layer_porosity_to_screen(double t_ime);
   void print_layer_top_elevations_to_screen(double t_ime);
   void print_particle_masses_to_screen(double t_ime);
   void print_ind_particles_of_layer_to_screen(int layer_number);

   void print_layer_thickness(double t_ime, ofstream& fout);
   void print_layer_mass(double t_ime, ofstream& fout);
   void print_layer_porosity(double t_ime, ofstream& fout);
   void print_layer_top_elevations(double t_ime, ofstream& fout);
   void print_layer_particle_proportions(double t_ime, ofstream& fout);
   void print_layer_particle_masses(double t_ime, ofstream& fout);
   void print_single_layer_particle_masses(double t_ime, ofstream& fout, int layer_number);
   void print_band_layer_particle_masses(double t_ime, ofstream& fout,
 										int bottom_layer, int top_layer);
 	void print_stack_stats(double t_ime, ofstream& fout);

 private:
    vector<sediment_layer> layers;
    double basement_elevation;
    double mean_sea_level_elevation;

    vector<double>  total_layer_mass;
    vector<double>  layer_top_elevations;
    vector<double>  layer_thickness;
    vector<double>  layer_porosity;


    void create();
    void create(double);

};

#endif
