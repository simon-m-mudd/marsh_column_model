//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Marshdriver
// The main driver function for the MarshColumn model
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Copyright (C) 22025 Simon M. Mudd 2025
//
// Developer can be contacted by simon.m.mudd _at_ ed.ac.uk
//
//    Simon Mudd
//    University of Edinburgh
//    School of GeoSciences
//    Drummond Street
//    Edinburgh, EH8 9XP
//    Scotland
//    United Kingdom
//
// This program is free software;
// you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation;
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY;
// without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
//
// You should have received a copy of the
// GNU General Public License along with this program;
// if not, write to:
// Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor,
// Boston, MA 02110-1301
// USA
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


#define _USE_MATH_DEFINES

#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>
#include <math.h>
#include <string>
#include <time.h>


//unix directory routing
#include "./class_headers/depo_particle.hpp"
#include "./class_headers/depo_particle_info.hpp"
#include "./class_headers/sediment_layer.hpp"
#include "./class_headers/sediment_stack.hpp"
#include "./class_headers/marsh_util_fxns.hpp"
#include "./class_headers/LSDParameterParser.hpp"
#include "./class_headers/LSDStatsTools.hpp"
using namespace std;

int main (int nNumberofArgs,char *argv[])
{

    string version_number = "1.1d";
    string citation = "https://www.doi.org/10.5281/zenodo.997388";


    cout << "================================================================" << endl;
    cout << "|| Welcome to the MarshColumn   !                             ||" << endl;
    cout << "|| This program drives the MuddPILE lanscape evolution model. ||" << endl;
    cout << "|| One day this model will have documentation.                ||" << endl;
    cout << "|| This program was developed by Simon M. Mudd                ||" << endl;
    cout << "||  at the University of Edinburgh                            ||" << endl;
    cout << "================================================================" << endl;
    cout << "|| If you use these routines please cite:                     ||" << endl;

    cout << "================================================================" << endl;
    cout << "|| Documentation can be found at:                             ||" << endl;

    cout << "================================================================" << endl;
    cout << "|| This is MuddPILE version                                   ||" << endl;
    cout << "|| " << version_number  << endl;
    cout << "|| If the version number has a d at the end it is a           ||" << endl;
    cout << "||  development version.                                      ||" << endl;
    cout << "================================================================" << endl;


    // Get the arguments
    vector<string> path_and_file = DriverIngestor(nNumberofArgs,argv);
    string path_name = path_and_file[0];
    string f_name = path_and_file[1];


    if(f_name == "lsdtt_version.txt")
    {
        ofstream ofs;
        ofs.open("./muddpiledriver-version.txt");
        ofs << version_number << endl;
        ofs.close();

        exit(0);
    }

    // load parameter parser object
    LSDParameterParser LSDPP(path_name,f_name);

    // for the chi tools we need georeferencing so make sure we are using bil format
    LSDPP.force_bil_extension();

    // maps for setting default parameters
    map<string,int> int_default_map;
    map<string,float> float_default_map;
    map<string,bool> bool_default_map;
    map<string,string> string_default_map;

    // this will contain the help file
    map< string, vector<string> > help_map;


    float_default_map["root_efold"] = 0.11;
    help_map["root_efold"] = {  "float","0.11","E-folding deth in metres of the root growth profile.","Default set for North inlet."};

    float_default_map["labile_frac"] = 0.842;
    help_map["labile_frac"] = {  "float","0.842","Fraction of root matter that is labile.","Default set for North inlet."};

    float_default_map["silt_frac"] = 1.0;
    help_map["silt_frac"] = {  "float","1.0","Fraction of suspended sediment that is silt.","Default set for a silt dominated marsh."};

    float_default_map["tidal_amplitude"] = 0.5;
    help_map["tidal_amplitude"] = {  "float","0.5","The tidal amplitude.","Default set for 1m tidal range."};

    float_default_map["total_ssc_conc"] = 0.01;
    help_map["tidal_amplitude"] = {  "float","0.01","Total suspended sediment concentration in kg/m^3.","Default set for 1m tidal range."}; 





    float_default_map["silt_fraction"] = 0.8;
    help_map["silt_conc"] = {  "float","0.8","Total suspended sediment concentration in kg/m^3.","Default set for 1m tidal range."};  
    
    float_default_map["fine_silt_fraction"] = 0.2;
    help_map["fine_silt_conc"] = {  "float","0.2","Total suspended sediment concentration in kg/m^3.","Default set for 1m tidal range."};     
    
    float_default_map["labile_frac"] = 0.85;
    help_map["labile_frac"] = {  "float","0.85","Labile carbon fraction.","Default set for 1m tidal range."};  

    float_default_map["refrac_frac"] = 0.15;
    help_map["refrac_frac"] = {  "float","0.15","Refractory carbon fraction.","Default set for 1m tidal range."};  

    float_default_map["peak_biomass"] = 2500;
    help_map["peak_biomass"] = {  "float","0.15","Peak biomass.","Default set for 1m tidal range."};  

    float_default_map["RSLR"] = 0.005;
    help_map["RSLR"] = {  "float","0.005","Relative sea level ris in m/yr.","Default set for 1m tidal range."};  

    float_default_map["theta_root_efolding"] = 0.2;
    help_map["theta_root_efolding"] = {  "float","0.2","The efolding depth of root production as a fration of the tidal amplitude.","Default set for 1m tidal range."}; 






    string_default_map["depo_particle_info_file"] = "dpinfo.param";
    help_map["depo_particle_info_file"] =  {  "string","dpinfo.param","The name of the paramter file for holding depositional particle information.","You must include the extension."}; 

    bool_default_map["test_column_initiation"] = false;
    help_map["test_column_initiation"] = {  "bool","false","This test the column initiation and then exits.","For testing the column and also testing the layer combination function."}; 



    //=========================================================================
    //
    //.#####....####...#####....####...##...##..######..######..######..#####..
    //.##..##..##..##..##..##..##..##..###.###..##........##....##......##..##.
    //.#####...######..#####...######..##.#.##..####......##....####....#####..
    //.##......##..##..##..##..##..##..##...##..##........##....##......##..##.
    //.##......##..##..##..##..##..##..##...##..######....##....######..##..##.
    //
    //..####...##..##..######...####...##..##...####..                         
    //.##..##..##..##..##......##..##..##.##...##.....                         
    //.##......######..####....##......####.....####..                         
    //.##..##..##..##..##......##..##..##.##.......##.                         
    //..####...##..##..######...####...##..##...####..                         
    //============================================================================

    // Use the parameter parser to get the maps of the parameters required for the analysis
    LSDPP.parse_all_parameters(float_default_map, int_default_map, bool_default_map,string_default_map);
    map<string,float> this_float_map = LSDPP.get_float_parameters();
    map<string,int> this_int_map = LSDPP.get_int_parameters();
    map<string,bool> this_bool_map = LSDPP.get_bool_parameters();
    map<string,string> this_string_map = LSDPP.get_string_parameters();

    if(f_name == "cry_for_help.txt")
    {
        cout << "I am going to print the help and exit." << endl;
        cout << "You can find the help in the file:" << endl;
        cout << "./MarshColumn-README.csv" << endl;
        string help_prefix = "MarshColumn-README";
        LSDPP.print_help(help_map, help_prefix, version_number, citation);
        exit(0);
    }

    // location of the files
    string OUT_DIR = LSDPP.get_write_path();
    string OUT_ID = LSDPP.get_write_fname();
    
    cout << "Write prefix is: " << OUT_DIR+OUT_ID << endl;

    // read the depo_particle file
    cout << "Let me read the depo partile info file" << endl;
    depo_particle_info dpi(  this_string_map["depo_particle_info_file"].c_str());
	int np_types = dpi.get_n_types();			// number of particle types




	double RSLR = double(this_float_map["RSLR"]);                         // sea level rise rate (m/yr)
	double conc_fs;                     // concentration of fine silt
	double conc_silt;                   // concentration of silt
	double tot_conc;                    // total SSC
	double labile_frac;
	double refrac_frac;
	double root_efold;
	double fine_silt_frac;
	double silt_frac;
	double gi_multiplier;
	double effective_svel;
	double tA;				// the tidal amplitude
	double peak_Bmass  = double(this_float_map["peak_biomass"]);		// this is returned from the model loop
    double theta_root_efold = double(this_float_map["theta_root_efolding"]);





    cout << "Let me read the parameter file" << endl;
    double thickness_iteration_threshold = 1e-6;
                            // a convergence threshold for
                            // used in compression calculations
                            // see sediment_stack::calculate_overburden_and_thickness()

    cout << "np_types is: " << np_types << endl;

    // parameters for the sedimentation component
    double Tidal_Period;			// tidal period in hours
    double Tidal_Amplitude;		// tidal amplitude in meters
    double marsh_surface_elevation;	// surface elevation in meters
    double marsh_surface_elevation_old;
    double marsh_surface_elevation_day0;
                                    // used for differencing
    double Mean_Tide;			// the mean tide above datum to start
    double flow_velocity = 0.01;		// flow velocity at column
    double alpha_T = 1.0e6;		// parameter for trapping
    double beta_T = 0.04;			// parameter for trapping
    double epsilon = 2.08;		// parameter for trapping
    double gamma = 0.718;			// parameter for trapping
    double t_ime;
    double tidal_periods_per_day;
    int number_of_particle_types;


    // some parameters for settling and trapping
    vector<double> particle_concentrations_0;
    vector<double> particle_diameters;
    vector<double> particle_settling_velocities;


    vector<double> mass_trap;

    // parameters for the biomass_component
    double MHT;				// mean high tide in m
    double max_depth;			// max depth of macrophytes (m)
    double min_depth;			// minimum depth of macrophytes (m)
    double max_bmass;			// maximum biomass in g/m^2
    double min_bmass;			// minimum biomass during the year
    double theta_bmin;		// ratio between the peak biomass and the minimum biomass
    double max_growth;		// maximum growth rate
    double nu_Gp;				// ratio between the max growth rate and
                                // the peak biomass in 1/days
    double min_growth;		// minimum growth rate
    double nu_Gmin;			// the ratio between the minimum growth rate and
                                // teh peak biomass in 1/days
    double phase_shift;		// the delay between the peak grwoth season and the peak biomass
                                // in days
    double theta_bm;			// slope of relationship between B_ag/B_bg and
                                // the depth below MHT
    double D_mbm;				// intercept of the B_ag/B_bg curve

    double gamma_bg_biomass;	// the e-folding length of
                                // belowground biomass production
    double top_bg_biomass_rate;
                                // the rate of belowground biomass
                                // production at the marsh surface
    double gamma_bg_biomass_ratio;
                                // ratio between the belowground biomass
                                // e-folding length and the tidal amplitude

    double day;
    double day_peak_season;	// the day when biomass peaks


    double Bmass;				// the aboveground biomass
    double Delta_Bmass;		// change in the aboveground biomass
    double BG_Bmass;			// belowground bmass
    double Delta_BG_Bmass;	// change in belowground biomass
    double BG_mort;			// mortality of belowground biomass
    int root_type;			// the type of the roots (see dpart.dlist)
    double Bmass_old;

    int n_carbon_types;		// the number of carbon types
    vector<int> carbon_types; // indices into carbon tyes, these should
                                // match the indeces in dpart.dlist
    vector<double> chi_carbon_types;
                                // the mass fractions of the carbon types
                                // the mass fractions should add to 1

    vector<double> bm_and_mort(2,0.0);
                                // vector containing the aboveground biomass
                                // and the death of belowground roots
    double BG_to_AG_ratio;

    double accretion_rate;
    double accretion_ratio;

    // parameters for radioisotopes
    int Pb_type;				// index in the d_particle_info of Pb210
    double CRS_Pb_supply;		// this is in g/yr

    // the following algorthm loads the parameter file....................//
    ifstream paramfile;													//
    paramfile.open("N_inlet_final.param");									//
    string temp_string;													//
    double temp_double;													//
    int temp_int;															//
    paramfile >> temp_string >> Tidal_Period								//
                >> temp_string >> Tidal_Amplitude							//
                >> temp_string >> Mean_Tide									//
                >> temp_string >> tidal_periods_per_day						//
                >> temp_string >> max_depth									//
                >> temp_string >> min_depth									//
                >> temp_string >> max_bmass									//
                >> temp_string >> theta_bmin								//
                >> temp_string >> nu_Gp										//
                >> temp_string >> nu_Gmin									//
                >> temp_string >> phase_shift
                >> temp_string >> gamma_bg_biomass_ratio					//
                >> temp_string >> RSLR										//
                >> temp_string >> theta_bm									//
                >> temp_string >> D_mbm										//
                >> temp_string >> day_peak_season							//
                >> temp_string >> number_of_particle_types					//
                >> temp_string >> Pb_type
                >> temp_string >> CRS_Pb_supply
                >> temp_string >> root_type									//
                >> temp_string >> n_carbon_types							//
                >> temp_string;												//
    for (int i = 0; i< n_carbon_types; i++)								//
    {																	//
        paramfile >> temp_int;												//
        carbon_types.push_back(temp_int);									//
    }																	//
    paramfile >> temp_string;												//
    for (int i = 0; i< n_carbon_types; i++)
    {
        paramfile >> temp_double;
        chi_carbon_types.push_back(temp_double);
    }
    paramfile >> temp_string;
    for (int i = 0; i<number_of_particle_types; i++)
    {
        paramfile >> temp_double;
        particle_concentrations_0.push_back(temp_double);
    }
    paramfile >> temp_string;
    for (int i = 0; i<number_of_particle_types; i++)
    {
        paramfile >> temp_double;
        particle_diameters.push_back(temp_double);
    }
    paramfile >> temp_string;
    for (int i = 0; i<number_of_particle_types; i++)
    {
        paramfile >> temp_double;
        particle_settling_velocities.push_back(temp_double);
    }

	MHT = Mean_Tide+Tidal_Amplitude;		// mean high tide in m
	//gamma_bg_biomass = gamma_bg_biomass_ratio*Tidal_Amplitude;

	// reset the sea level rise and concnetrations for this run
    double SLR =RSLR;
	particle_concentrations_0[0] = conc_silt;
	particle_concentrations_0[1] = conc_fs;

	// set the particle settling velocity
	particle_settling_velocities[0] = effective_svel;

	chi_carbon_types[0] = refrac_frac;
	chi_carbon_types[1] = labile_frac;

	gamma_bg_biomass = Tidal_Amplitude*theta_root_efold;


    if(this_bool_map["test_column_initiation"])
    {
        cout << "You asked me to test initiation the column with sand" << endl;
    }

    if(this_bool_map["test_TKE_settling"])
    {
        cout << "You asked me to test the particle settling" << endl;
    }

    cout << "I've done everything you asked and am now going back to sleep. Good night." <<endl;

}
