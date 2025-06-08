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


    float_default_map["initial_layer_mass"] = 25;
    help_map["initial_layer_mass"] = {  "float","25","The intial mass of the sediment layers..","In kg (column is always assumed 1 m^2) assumed to be 0 index sediment type from d_particle_info."}; 


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

    int_default_map["end_year"] = 2;
    help_map["end_year"] = {  "int","2","The end year of the simulation.","Default is for a short simulation"};     




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
    cout << "This is: " << this_string_map["depo_particle_info_file"] << endl;
    depo_particle_info dpi(  this_string_map["depo_particle_info_file"].c_str());
	int np_types = dpi.get_n_types();			// number of particle types
    cout << "np_types is: " << np_types << endl;



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





    //==================================================================================================
    // parameters for TKE module see Mudd et al 2010
    double alpha = 0.55;
    double beta = 0.40;
    double mu = 0.00066;
    double phi = 0.29;
    double kappa = 0.224;
    double nu = 10e-6;                // kinematic viscosity of water a 20C
    double gamma_tr = 0.718;        // exponent from Nepf (1999)
    double epsilon = 2.08;        // exponent from Nepf (1999)
    double max_flow_velocity;                // flow velocity at column
    double g = 9.80;
    double reference_flow_velocity = 0.05; //default=0.05m/s.
    double limiting_flow_velocity = 0.1;
    double reference_biomass = 200;
    //MK- For comparsion Christiansen measured velocity 0.0003-.003 m/s 45 m from channel. EOSL biomass at UPC is 500-1000 g/m2.
    double alpha_T = alpha*kappa*pow(mu,gamma_tr-epsilon)/pow(nu,gamma_tr);
    double beta_T = phi*(gamma_tr-epsilon)+beta;
    double aa = 0.46;                // coeficient in equation of alhpa_1
                                                    // as fxn of solid fraction; Tanino and Nepf
                                                    // eq. 13
    double bb = 3.8;                // coeficient in equation of alhpa_1
                                                    // as fxn of solid fraction; Tanino and Nepf
                                                    // eq. 13
    double alpha_turb = 0.9;// turbulence coefficient Nepf (1999)
    double alpha_zero = 11;        // from tanino and nepf (2008)

    double reference_slope;
    reference_slope = get_slope(aa, bb, alpha_turb, alpha_zero,
                                alpha, beta, mu, phi, nu, g, reference_biomass, reference_flow_velocity);
    //
    //==================================================================================================

    // these are vectors for settling and trapping mass
	vector<double> tmass;
    vector<double> smass;

	double start_depth_frac;

	int i_start_frac=1;	//MK- I changed this from 8 to 4. Want the marsh to start low in the tide frame and build up.
	start_depth_frac = double(i_start_frac)*0.1;



    cout << "Let me read the parameter file" << endl;
    double thickness_iteration_threshold = 1e-6;
                            // a convergence threshold for
                            // used in compression calculations
                            // see sediment_stack::calculate_overburden_and_thickness()


    // parameters for the sedimentation component
    double Tidal_Period;			// tidal period in hours
    double Tidal_Amplitude;		// tidal amplitude in meters
    double marsh_surface_elevation;	// surface elevation in meters
    double marsh_surface_elevation_old;
    double marsh_surface_elevation_day0;
                                    // used for differencing
    double Mean_Tide;			// the mean tide above datum to start
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
                                // the peak biomass in 1/days
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


    cout << "I've read the parameter file. here are some parameters:" << endl;
    cout << "nu_gp: " << nu_Gp << ", mean tide: " << Mean_Tide << ", root type: " << root_type << endl;

	// calcualte the particle settling velocities explicitly
	effective_svel = calculate_w_s(particle_diameters[0]);

    // some parametrs for time series
	double SS_acc_ratio_trigger = 0.0001;
	double cycles_per_month  = 1.0;
	double months_per_year = 12.0;
	int timesteps_per_year = int(double(months_per_year*cycles_per_month+0.5));
	double dt = 365.0/(cycles_per_month*months_per_year);			// time spacing
	double yr_time;




    if(this_bool_map["test_column_initiation"])
    {
        cout << "You asked me to test initiation the column" << endl;
        sediment_stack sed_stack = initiate_column(dpi,double(this_float_map["initial_layer_mass"]), Mean_Tide);


        cout << "The number of layers is: " << sed_stack.get_n_layers() << endl;
        double t_ime = 0;
        sed_stack.print_layer_top_elevations_to_screen(t_ime);

    }

    if(this_bool_map["test_TKE_settling"])
    {
        cout << "You asked me to test the particle settling" << endl;
    }

    cout << "I've done everything you asked and am now going back to sleep. Good night." <<endl;


    if(this_bool_map["run_column_timeseries"])
    {
        string fname = OUT_DIR+OUT_ID+"_data.csv";
        ofstream data_out;
        data_out.open(fname.c_str());

        data_out << "SLRR,Tidal_amplitude,start_depth_frac,yr,conc_silt,effective_svel,";
		data_out <<	"MHT,MHT-marsh_surface_elevation,max_depth,marsh_surface_elevation-marsh_surface_elevation_day0,peak_Bmass" << endl;       

        string col_out_fname_prefix = "column_out_";
        string col_out_fname;




        // vectors for holding timeseries information
        vector<double> TimeSeries_MHT;	
        vector<double> TimeSeries_Biomass;
        vector<double> TimeSeries_Elevation;
        vector<double> TimeSeries_Accretion;
        vector<double> TimeSeries_Yr;
        vector<double> TimeSeries_Ref;
        vector<double> TimeSeries_Lab;
        vector<double> TimeSeries_Silt;



        // initiate the stack
        sediment_stack sed_stack = initiate_column(dpi,double(this_float_map["initial_layer_mass"]), Mean_Tide);


        cout << "The number of layers is: " << sed_stack.get_n_layers() << endl;
        double t_ime = 0;
        sed_stack.print_layer_top_elevations_to_screen(t_ime);


        // set up elevations
        marsh_surface_elevation_old =sed_stack.get_marsh_surface_elevation();
        marsh_surface_elevation = marsh_surface_elevation_old;

        // two variables used for printing the stats_out file
        vector<double> particle_masses = sed_stack.get_total_stack_mass_of_each_type(np_types);
        vector<double> particle_masses_old = particle_masses;
        vector<double> start_masses = particle_masses;
        vector<double> particle_masses_day0 = particle_masses;

        cout  << "marsh surf elev: " << marsh_surface_elevation << endl;

        // calucalte the max depth based on the McKee and Patrick data
        max_depth = -(0.0406-0.6663*tA*2);

        // place mean high tide at intermdiate depth in biomass curve
        double intermediate_depth = (max_depth-min_depth)*start_depth_frac+min_depth;
        double starting_depth = .06; // MK- I added this line, and the next. Usually doing .06 
        MHT = marsh_surface_elevation + starting_depth;
        //MHT = marsh_surface_elevation + intermediate_depth; //MK- I commented this one out
        Mean_Tide = MHT-Tidal_Amplitude;


        cout << "MHT: " << MHT << "S tarting elevation: " << marsh_surface_elevation << endl;



        double SS_acc_ratio_trigger = 0.0001;
        double cycles_per_month  = 1.0;
        double months_per_year = 12.0;
        int timesteps_per_year = int(double(months_per_year*cycles_per_month+0.5));
        double dt = 365.0/(cycles_per_month*months_per_year);			// time spacing
        double yr_time;


	    vector<int> annual_band_layer_top;		// gives the location of the top
											    // layer of an annual band
        t_ime = 0;						// initilaize the time


        // first, get the top layer of the initial pile of sediment
        annual_band_layer_top.push_back(sed_stack.get_n_layers()-1);
        cout << "LINE 359 starting n_layers: " << sed_stack.get_n_layers() << endl;

        // start a loop
        accretion_ratio = 10;                       // SSM need to figure out what this does
        int dead_biomass_counter = 0;
        int low_acc_ratio_counter = 0;
        int yr = 0;
        while (yr < this_int_map["end_year"])
        {
            yr++;
            t_ime = yr*365;
            yr_time = t_ime/365;

            SLR= RSLR; 
            

            // reset the biomass
            Bmass = 0;

            // recalcualte the growth index. see notes in sediment_stack.cpp
            vector<double> growth_index = sed_stack.create_growth_index_inf(gamma_bg_biomass);


            marsh_surface_elevation = sed_stack.get_marsh_surface_elevation();
            marsh_surface_elevation_day0 = marsh_surface_elevation;
            cout << "I am starting year " << yr << " and the main marsh_surface elevation is: " << marsh_surface_elevation << endl;

            // calcualte t./column.out
            //the peak biomass for the year
            double temperatureincrease = 0;
            double bfactor = 0;
            peak_Bmass =  get_peak_biomass(marsh_surface_elevation, MHT,
            max_depth, min_depth, max_bmass,temperatureincrease,yr,bfactor);	//MK- I added temperature,yr,boriginal to this line


            min_bmass = peak_Bmass*theta_bmin;
            max_growth = nu_Gp*peak_Bmass;
            min_growth = nu_Gmin*peak_Bmass;
            particle_masses = sed_stack.get_total_stack_mass_of_each_type(np_types);
            particle_masses_day0 = sed_stack.get_total_stack_mass_of_each_type(np_types);
            //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
            // now go day by day, calculating the biomass, root growth and death
            // and sedimentaion
            //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            day = 0;
            int count = 0;
            for (int yo = 1; yo <= timesteps_per_year; yo++)
            //for (int yo = 1; yo <= 2; yo++)
            {
                marsh_surface_elevation_old = marsh_surface_elevation;
                particle_masses_old = particle_masses;

                day +=dt;
                count++;
                t_ime = t_ime + dt;
                yr_time = t_ime/365;
                //cout << "LINE 363 main day is: " << day << endl;

                // get the rate of growth and decay
                // this returns a vector that is the change in the aboveground and
                // beolowground biomass during the time period dt
                //
                //cout << "LINE 418 starting growth and mortality function\n";

                //bm_and_mort =  biomass_and_mortality(day, peak_Bmass, min_bmass,
                //          day_peak_season, max_growth, min_growth, dt);
                bm_and_mort =  biomass_and_mortality2(day, peak_Bmass, min_bmass,
                        day_peak_season, max_growth, min_growth, dt, phase_shift);
                Bmass_old = Bmass;
                Bmass = bm_and_mort[0];
                Delta_Bmass = Bmass - Bmass_old;

                // calculate ratio between BG_and aboverground biomass

                // in this model, the tidal amplitude can cahnge and with it the range of
                // plant growth changes. For large amplitudes, large depths below MHT can
                // occur, so we need to adjust the empirically determined
                // BG_to_AG ratio
                // option 1: limit the ratio to a minimum value
                BG_to_AG_ratio = theta_bm*(MHT-marsh_surface_elevation)+D_mbm;
                if (BG_to_AG_ratio <0.1)
                {
                    BG_to_AG_ratio = 0.1;
                }


                //cout << "LINE 444 depth: " << MHT-marsh_surface_elevation << " and B_mass is: " << Bmass << endl;
                //cout << " ratio is: " << BG_to_AG_ratio << endl;

                // now get the belowground biomass and its change
                // NOTE we divide by 1000 becuase sedimentation is in kg whereas biomass
                // parameters are in g
                BG_Bmass = Bmass*BG_to_AG_ratio/1000;
                Delta_BG_Bmass = Delta_Bmass*BG_to_AG_ratio/1000;
                BG_mort = bm_and_mort[1]*BG_to_AG_ratio/1000;


                //*****Added these 2 lines to make Simon tke code more stable ////
                if (Bmass == 0)
                {Bmass = 1;}
                /////////////////////////////


                //************ Simon- new function call
                // calculate the maximum flow velocity using a reference water surface slope
                            // that is defined by a reference velocity set at a reverence biomass (see line 289)
                max_flow_velocity =  get_u(aa, bb, alpha_turb, alpha_zero,
                                    alpha, beta, mu, phi,
                                    nu, g, Bmass, reference_slope);

                if (max_flow_velocity > limiting_flow_velocity)
                {
                    max_flow_velocity =limiting_flow_velocity ;
                }
                //**************




                //cout << "LINE 471 main: BG_Bmass: " << BG_Bmass << " BG_mort: " << BG_mort <<endl;

                //cout << "LINE 459 masses before " << endl;
                //sed_stack.print_particle_masses_to_screen(t_ime);


                sed_stack.root_growth_and_death(growth_index, root_type,
                                carbon_types, chi_carbon_types, BG_Bmass, BG_mort);

                //cout << "LINE 466 masses after root growth death " << endl;
                //sed_stack.print_particle_masses_to_screen(t_ime);

                // compress the stack
                //cout << "LINE 453 starting compression" << endl;
                sed_stack.calculate_overburden_and_thickness();
                //cout << "LINE 455 compressed" << endl;

                // now start building the sedimentation vector
                // the sedimentation routine used here takes a depo_particle_info object
                // and a vector. The vector is the same size as the number of particles
                vector<double> depo_masses(np_types, 0.0);	// the vector that stores
                                                            // the deposition masses

                // now get the masses due to trapping and settling. This uses particle
                // concentrations within the flow column, as well as settling velocities
                //cout <<"LINE 465 trapping" << endl;
                //cout <<"Mean_Tide: " << Mean_Tide << " surf elev: " << marsh_surface_elevation << endl;

                //Simon- this new line replaced the original


                mass_trap = get_trap_TKE_eff_sett(Tidal_Period, Tidal_Amplitude,
                                        Mean_Tide, marsh_surface_elevation, particle_concentrations_0,
                                        particle_diameters, particle_settling_velocities, max_flow_velocity,
                                        Bmass, alpha_T, beta_T, epsilon, gamma_tr, aa, bb, alpha_turb,
                                        alpha_zero, alpha, beta, mu, phi, nu, g, smass,tmass);
                // add the trapped mass into the depo_masses vector
                for (int i = 0; i<np_types; i++)
                {
                    depo_masses[i] += dt*mass_trap[i]*tidal_periods_per_day;
                }
                //cout << "LINE 476 trapped" << endl;

                // now add the organogenic deposition from the surface
                for (int i = 0; i< n_carbon_types; i++)
                {
                    // the equation below is divided by 1000 because mortality is calculated in
                    // g/m^2 and deposition is in kg
                    depo_masses[ carbon_types[i] ] += chi_carbon_types[i]*bm_and_mort[1]/1000;
                }

                // now deposit the radiogenic species.
                depo_masses[ Pb_type ] += dt*CRS_Pb_supply/365.0;

                // now deposit the sediment
                sed_stack.deposit_surface_sediment(dpi,depo_masses);
                // and compress it
                sed_stack.calculate_overburden_and_thickness();


                //cout << "LINE 517 masses after surface depo " << endl;
                //sed_stack.print_particle_masses_to_screen(t_ime);

                // print to screen the masses of the individual particle types in the layers
                //cout << "LINE 490 main printing particle masses before decay " << endl;
                //sed_stack.print_particle_masses_to_screen(t_ime);
                //cout << endl;

                // now do some decay
                double kfactor = 1;            // no temp dependent decay here
                sed_stack.decay_and_mass_loss(dt,temperatureincrease,kfactor);	//MK- I added temperatureincrease and kfactor

                //cout << "LINE 528 masses after decay " << endl;
                //sed_stack.print_particle_masses_to_screen(t_ime);


                // print to screen the masses of the individual particle types in the layers
                //cout << "LINE 445 main printing particle masses after decay " << endl;
                //sed_stack.print_particle_masses_to_screen(t_ime);
                //cout << endl;

                // compress again
                sed_stack.calculate_overburden_and_thickness();

                // increase sea level
                //cout << "MHT: " << MHT << " dt/365: " << dt/365.0 << endl;
                MHT += dt*SLR/365.0;
                Mean_Tide+= dt*SLR/365.0;

                // get data for stats file
                particle_masses = sed_stack.get_total_stack_mass_of_each_type(np_types);
                marsh_surface_elevation = sed_stack.get_marsh_surface_elevation();
                //cout << endl << endl << "LINE 533, surf elev: " << marsh_surface_elevation << endl << endl;

            }	// !! end sub annual loop

            //cout << "main, LINE 538, printing stack masses" << endl;

            // now get the top layer of the annual band
            annual_band_layer_top.push_back(sed_stack.get_n_layers()-1);
            //cout << "n_layers: " << sed_stack.get_n_layers() << endl;

            accretion_rate = marsh_surface_elevation- marsh_surface_elevation_day0;
            accretion_ratio = accretion_rate/SLR;
            if (yr%1 == 0)
            {
                cout << "year is: " << yr;
                cout << " accretion rate = " << accretion_rate << " MHT= "<< MHT
                << " depth bel MHT: " << MHT-marsh_surface_elevation << " biomass= " << peak_Bmass << " SLRR= " << SLR*1000
                << endl;
            }

            int ab_pointer_sz = annual_band_layer_top.size();

            if (peak_Bmass == 0)
            {
                dead_biomass_counter++;
            }
            if ( (1-accretion_ratio)*(1-accretion_ratio) < SS_acc_ratio_trigger)
            {
                low_acc_ratio_counter++;
            }

            //MK- add some lines storing MHT, Biomass,elevation, and accretion rate to a vector
            //temperaturevector.push_back(temperature);
            TimeSeries_Yr.push_back(yr);
            TimeSeries_MHT.push_back(MHT);
            TimeSeries_Biomass.push_back(peak_Bmass);
            TimeSeries_Elevation.push_back(marsh_surface_elevation);
            TimeSeries_Accretion.push_back(accretion_rate);


            // SMM: adding the function that gets the total mass of all the kinds of particles in the
            // column
            // note: the only particles that matter are in index 0 (silt) and 2,3,4 (labile, refractory and roots)
            particle_masses = sed_stack.get_total_stack_mass_of_each_type(np_types);
            TimeSeries_Silt.push_back(particle_masses[0]);
            TimeSeries_Ref.push_back(particle_masses[2]);
            TimeSeries_Lab.push_back(particle_masses[3]);


            // now print the data
            data_out << SLR << "," << "," << Tidal_Amplitude << ","
                        << start_depth_frac << "," << yr << ","
                        << conc_silt << "," << effective_svel << ","
                        << MHT << "," << MHT-marsh_surface_elevation << ","
                        << max_depth << ","
                        << marsh_surface_elevation-marsh_surface_elevation_day0
                        << "," <<peak_Bmass;

            // print the accretion rate for each particle type
            for (int i = 0; i< np_types; i++)
            {
                data_out << "," << particle_masses[i]-particle_masses_day0[i];
            }
            data_out << endl;

            if (yr%2 == 0)
            {
                ofstream col_out;
                col_out_fname = "col_out_" + itoa(yr)+".csv";
                col_out.open(col_out_fname.c_str());
                col_out << "layer,masses"<< endl;
                
                
                int ab_pointer_sz = annual_band_layer_top.size();
                //cout << "n annual bands: " << ab_pointer_sz << endl;
                vector<double> pm_band;
                for (int i = 1; i<ab_pointer_sz; i++)
                {
                    col_out << i << "," << sed_stack.get_layer_top_elevation(annual_band_layer_top[i]);
                    pm_band = sed_stack.get_stack_slice_mass_of_each_type(np_types,
                                                        annual_band_layer_top[i-1]+1,
                                                annual_band_layer_top[i]);
                    for (int type = 0; type < np_types; type++)
                    {
                        col_out << "," << pm_band[type];
                    }
                    col_out << endl;
                }

                col_out.close();
            }



        }		// !! end annual loop

        cout << " ENDING THE LOOP, SLR: " << SLR << endl;
        data_out.close();

 


        // MK- New function to save time series data
        ofstream series_out;
        string fname2_prefix="series.";
        string fname2_suffix=".txt";
        string fname2_k=itoa(0);
        string fname2_b=itoa(0);
        //    if (kfactor!=0)
        //    {fname2_k=itoa(1);}
        //    if (bfactor!=0)
        //    {fname2_b=itoa(1);}
        string fname2=fname2_prefix+fname2_k+fname2_b+fname2_suffix;
        series_out.open(fname2.c_str());
        //series_out.open("series.txt", ios::out); // This simpler version didn't allow multiple output files for multiple parameters
        for (int q=0; q<TimeSeries_Yr.size(); q++)
        {
            series_out << TimeSeries_Yr[q] << " " << TimeSeries_MHT[q] << " " << TimeSeries_Elevation[q] << " " << TimeSeries_Biomass[q] << " "
                    << TimeSeries_Silt[q] << " " << TimeSeries_Ref[q] << " " << TimeSeries_Lab[q] << " " << endl;
        }
        series_out.close();

        //filename code: 00= no enhanced decay or productivity, 10= enhanced decay only, 01=enhanced prod only, 11=enhanced decay and productivity,
        ////////////////////

    }

}
