// crit_SLR.cpp
// this model uses the OIMAS-N model to caluclate the
// critical rate of sea level rise required to
// drown a marsh

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

// windows directory routing: use if compiling on DOS system (bcc32)
//#include ".\class_headers\depo_particle.cpp"
//#include ".\class_headers\depo_particle_info.cpp"
//#include ".\class_headers\sediment_layer.cpp"
//#include ".\class_headers\sediment_stack.cpp"

//unix directory routing
#include "./class_headers/depo_particle.cpp"
#include "./class_headers/depo_particle_info.cpp"
#include "./class_headers/sediment_layer.cpp"
#include "./class_headers/sediment_stack.cpp"
using namespace std;

// this is the function used to calculate deposition
vector<double> get_trapping_settling(double Tidal_Period, double Tidal_Amplitude,
                           double Mean_Tide, double marsh_surface_elevation,
                           vector<double> particle_concentrations_0,
                           vector<double> particle_diameters,
                           vector<double> particle_settling_velocities,
                           double flow_velocity,
                           double Biomass, double alpha_T, double beta_T,
                           double epsilon, double gamma);

vector<double> get_trapping_settling_mudd(double Tidal_Period, double Tidal_Amplitude,
                           double Mean_Tide, double marsh_surface_elevation,
                           vector<double> particle_concentrations_0,
                           vector<double> particle_diameters,
                           vector<double> particle_settling_velocities,
                           double flow_velocity,
                           double Biomass, double alpha_T, double beta_T,
                           double epsilon, double gamma_tr);


//****************************** This section is new, from Simon //
vector<double> get_trap_TKE_eff_sett(double Tidal_Period, double Tidal_Amplitude,
                           double Mean_Tide, double marsh_surface_elevation,
                           vector<double> particle_concentrations_0,
                           vector<double> particle_diameters,
                           vector<double> particle_settling_velocities,
                           double peak_flow_velocity,
                           double Biomass, double alpha_T, double beta_T,
                           double epsilon, double gamma_tr,
                           double aa, double bb, double alpha_turb, double alpha_zero,
                           double alpha, double beta, double mu, double phi,
                           double nu, double g, vector<double>& smass, vector<double>& tmass);

double calculate_w_s(double d_p);
double get_k(double aa, double bb, double alpha_turb, double alpha_zero,
                                 double alpha, double beta, double mu, double phi,
                                 double nu, double g, double B, double u);
double get_u(double aa, double bb, double alpha_turb, double alpha_zero,
                                 double alpha, double beta, double mu, double phi,
                                 double nu, double g, double B, double S);
double get_slope(double aa, double bb, double alpha_turb, double alpha_zero,
                                 double alpha, double beta, double mu, double phi,
                                 double nu, double g, double B, double u);
//****************************************


// this gets the peak biomass
double get_peak_biomass(double surface_elevation, double MHT,
		   double max_depth, double min_depth, double max_bmass, double temperatureincrease, int yr, double boriginal, double bfactor);

double get_peak_biomass_parab(double surface_elevation, double MHT,
		   double max_depth, double min_depth, double max_bmass);

// this gets the biomass growth and mortality for a given time of year
vector<double> biomass_and_mortality(double day, double peak_biomass, double min_biomass,
		      double day_peak_season, double peak_growth, double min_growth,
		      double dt);

vector<double> biomass_and_mortality2(double day, double peak_biomass, double min_biomass,
		      double day_peak_season, double peak_growth, double min_growth,
		      double dt,double phase_shift);

// the main model engine
double column_model(double RSLR, double kfactor, double bfactor, double tA, double conc_silt, double conc_fs,
				  double labile_frac, double refrac_frac,
				  double theta_root_efold, double effective_svel, double starting_depth_frac,
				  ofstream& data_out, ofstream& col_out);

// a utility function used to convert integers to strings
string itoa(int num);

// the function main is a wrapper function used to cycle through model runs
// with different parameter values. You simply compile with a start run number
// and end run number, the program then goes into directories ./data/column_run_i
// where i is the run number and the model extracts the parameter file and particle
// file from this directory and prints output data into this directory.

// to be added: the start and end run numbers can be
// changed to arguments to main() -> main(int argv, char* argc)
int main()
{
	time_t start,end;			// used for timing processes
	double dif;				// difference between machine times

	double kfactor=.25;				// MK- I added these two variables to set enhanced productivity and decay response to temperature
	double bfactor=.06;				// MK- I added these two variables to set enhanced productivity and decay response to temperature
	double SLR;
	double conc_fs;
	double conc_silt;
	double tot_conc;
	double labile_frac;
	double refrac_frac;
	double root_efold;
	double fine_silt_frac;
	double silt_frac;
	double gi_multiplier;
	double effective_svel;
	double tA;				// the tidal amplitude
	double peak_Bmass;		// this is returned from the model loop

	string s_temp;
	string num;
	string fname = "cSLR_tA1.tcout";

    vector<double> TC(10);
    TC[0] = 0.01;

    double kfactor_array[4] = {0,kfactor,kfactor,0};	// MK- To do multiple paramter runs
    double bfactor_array[4] = {0,bfactor,0,bfactor};	// MK- 1st is ambient, 2nd enhances both, 3rd enhances decay, 4th enhances productivity

    ofstream data_out;
    data_out.open(fname.c_str());

    string col_out_fname_prefix = "column_out_";
    string col_out_fname_suffix = ".cdata";
    string col_out_fname;

	double start_depth_frac;

	int i_start_frac=1;	//MK- I changed this from 8 to 4. Want the marsh to start low in the tide frame and build up.
	start_depth_frac = double(i_start_frac)*0.1;
	tA = 0.5;
	//tA = 0.1*double(i_tA)+0.4;		// this is the tidal amplitude
	tot_conc = .001;		// i_ssc is 2, so ssc = 0.01, or 10mg/L

	// these are the parameters for north inlet
	root_efold = 0.11;
	// effective_svel = 0.000037;			// UPDATE 9-sept-2011: tke model calcualtes effective settling directly
	labile_frac = 0.842;
	silt_frac = 1.0;
	conc_silt = tot_conc*silt_frac;
	conc_fs = 0;
	refrac_frac = 1-labile_frac;


	// shooting loop
	// the loop starts at a low rate of SLR (1 mm/yr)
	// each model run returns the peak biomass
	// if the biomass is 0, then the model has overshot
	// so the model increases the SLR by half the value it increased
	// previously and runs again#
	double old_SLR=0;
	double SLR_increase = 0.0001;

	//while (SLR_increase >= .0001)	 //MK- I commented this out
	for (int i = 0; i<=3; i++)
	{
		kfactor=kfactor_array[i];	// MK- Next 3 lines are to cycle through multiple parameters
		bfactor=bfactor_array[i];
		cout<< " Run #: " << i+1 << " kfactor= " << kfactor << " bfactor= " << bfactor << endl;

		//SLR = double(i)*0.001;
		SLR = old_SLR+SLR_increase;

		cout << endl << " fname: " << fname << endl
			 << "theta_gamma_roots: " << root_efold << endl
			 << "effective s vel: " << effective_svel << endl
			 << "labile frac: " << labile_frac << endl
			 << "SLR is: " << SLR << endl
			 << "tot_conc is: " << conc_silt << endl
			 << "tidal amplitude is: " << tA << endl;


		int SLR_int = int(SLR*10000);	// saving file with SLRR as the file number
		string SLR_str = itoa(SLR_int);
		col_out_fname = col_out_fname_prefix+SLR_str+col_out_fname_suffix;
		cout << "column fname is: " << col_out_fname << endl;
		ofstream col_out;
		col_out.open(col_out_fname.c_str());


		time (&start);				// get starting time
		peak_Bmass = column_model(SLR,kfactor, bfactor, tA,conc_silt,conc_fs,
				   labile_frac, refrac_frac,
				   root_efold, effective_svel, start_depth_frac,
				   data_out, col_out);
		time (&end);					// get ending time
		dif = difftime (end,start);
		cout << "\nruntime was: " << dif << " seconds\n";
		col_out.close();

		if (peak_Bmass > 0)			// we want SLR to go up when eq depth is approached, as long as plants remain. Stop everything when plants die.
		{
			old_SLR = SLR;
		}
		else
		{
			cout << "Plants died at SLR = " << SLR << " m/yr" << endl;
			//SLR_increase = SLR_increase/2;
			SLR_increase=0;
		}
	}



	data_out.close();
}


// this is the main model. It takes the int run number and loads
// model parameters from the filt ./data/column_run_i where i is
// the run number
// note, the implementation changes depending on whether you compile
// on a unix or windows machine
double column_model(double RSLR, double kfactor, double bfactor, double tA, double conc_silt, double conc_fs,
				  double labile_frac, double refrac_frac,
				  double theta_root_efold, double effective_svel, double start_depth_frac,
				  ofstream& data_out, ofstream& col_out)
{

	/// MK- I added this section to load temperature data from a file///////////////
	ifstream myfile;
	myfile.open("A2temp.txt", ios::in);
	vector<double> temperaturevector;
	double temperature;
	while (!myfile.eof())	// This makes a vector that contains one extra entry at end
		{
		myfile>>temperature;
		temperaturevector.push_back(temperature);
		}
	myfile.close();
	temperaturevector.pop_back(); // Removing that extra entry at end
	cout << "size of temperature vector: " << temperaturevector.size() << endl;

	///////MK- I added this section to load SL. Put all SL calculations in column_model function /////////////
	ifstream myfile2;
	myfile2.open("A2smoothsl.txt", ios::in);
	vector<double> slvector;
	double sealevel;
	while (!myfile2.eof())	// This makes a vector that contains one extra entry at end
		{
		myfile2>>sealevel;
		slvector.push_back(sealevel);
		}
	myfile2.close();
	slvector.pop_back(); // Removing that extra entry at end
	cout << "size of sea level vector: " << slvector.size() << endl;
        /////////////////
	/// 6. Figure out a good way to switch on and off productivity/decomp effects so that it uses the same basic code for each
	/// 7. Perhaps use main() to cycle through baseline, +productivity, +decomp, +both. Define bfactor and kfactor in main() and pass to needed functions. 			Setting bfactor to zero turns off enhanced productivity, setting kfactor to zero does the same.

	vector<double> TimeSeries_MHT;	//MK- I added these 7 lines for a new output file
	vector<double> TimeSeries_Biomass;
	vector<double> TimeSeries_Elevation;
	vector<double> TimeSeries_Accretion;
	vector<double> TimeSeries_Yr;
	vector<double> TimeSeries_Ref;
	vector<double> TimeSeries_Lab;
	vector<double> TimeSeries_Silt;
/////////////////////////////////////////////////////////



//***************  This section is new parameter values for Simon's TKE functions ////////////
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

        // peak B for fert 3280
        // alpha_T: 5.82x10^6
        // beta_T: 0.2476

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
//*******************************************************************************//

	double SS_acc_ratio_trigger = 0.0001;
	double cycles_per_month  = 1.0;
	double months_per_year = 12.0;
	int timesteps_per_year = int(double(months_per_year*cycles_per_month+0.5));
	double dt = 365.0/(cycles_per_month*months_per_year);			// time spacing
	double yr_time;

	vector<int> annual_band_layer_top;		// gives the location of the top
											// layer of an annual band
	int ab_pointer_sz;

	int end_year = 100;				// the final year	//MK- I changed this from 1500 to 10000

	int fixed_top;
	int fixed_bot;

	// some variables
	// load the particle info
	depo_particle_info dpi("dpart_20cd.dplist");
	int np_types = dpi.get_n_types();			// number of particle types
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
	// Simon- These five variables are not used with the new TKE calcs, or are previously declared
	//double flow_velocity = 0.01;		// flow velocity at column
	//double alpha_T = 1.0e6;		// parameter for trapping
	//double beta_T = 0.12;			// parameter for trapping
	//double epsilon = 2.08;		// parameter for trapping
	//double gamma = 0.718;			// parameter for trapping
	double t_ime;
	double tidal_periods_per_day;
	int number_of_particle_types;

	// some parameters for settling and trapping
	vector<double> particle_concentrations_0;
	vector<double> particle_diameters;
	vector<double> particle_settling_velocities;
	vector<double> mass_trap;

	// parameters for the biomass_component
	double temperatureincrease;		// MK- cumulative temperature increase: temp(t)-temp(o)
	double boriginal=400;			// MK- baseline biomass to calculate % change from
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
	double SLR;				// rate of sea level rise in m/yr
	double peak_Bmass;
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
							// match the indices in dpart.dlist
	vector<double> chi_carbon_types;
  							// the mass fractions of the carbon types
  							// the mass fractions should add to 1

	vector<double> bm_and_mort(2,0.0);
							// vector containing the aboveground biomass
							// and the death of belowground roots
	double BG_to_AG_ratio;
	double accretion_rate;
	double accretion_ratio;

	// These two lines are new from Simon's TKE code //
	vector<double> tmass;
    vector<double> smass;
	//

	// parameters for radioisotopes
	int Pb_type;				// index in the d_particle_info of Pb210
	double CRS_Pb_supply;		// this is in g/yr

	// the following algorthm loads the parameter file....................//
	ifstream paramfile;													//
	paramfile.open("N_inlet_final.param");								//
	string temp_string;													//
	double temp_double;													//
	int temp_int;														//
	paramfile >> temp_string >> Tidal_Period							//
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
			>> temp_string >> SLR										//
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
		paramfile >> temp_int;											//
		carbon_types.push_back(temp_int);								//
	}																	//
	paramfile >> temp_string;											//
	for (int i = 0; i< n_carbon_types; i++)
	{
		paramfile >> temp_double;
		chi_carbon_types.push_back(temp_double);
	}
	paramfile >> temp_string;
	//*********New from simon. This replaces particle_concentrations input loop
	for (int i = 0; i<number_of_particle_types; i++)
        {
                paramfile >> temp_double;
                particle_concentrations_0.push_back(temp_double);
                tmass.push_back(0.0);
                smass.push_back(0.0);
        }
	paramfile >> temp_string;
	//**************************
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



	// SMM 9-sept: calcualte the particle settling velocities explicitly
	effective_svel = calculate_w_s(particle_diameters[0]);


	// reset the sea level rise and concnetrations for this run
	SLR = RSLR;
	particle_concentrations_0[0] = conc_silt;
	particle_concentrations_0[1] = conc_fs;

	// set the particle settling velocity
	particle_settling_velocities[0] = effective_svel;

	cout << "You chave chose particles with diameter: " << particle_diameters[0]
	     << " their settling velocity is: " << effective_svel << endl;


	chi_carbon_types[0] = refrac_frac;
	chi_carbon_types[1] = labile_frac;

	Tidal_Amplitude = tA;
	gamma_bg_biomass = Tidal_Amplitude*theta_root_efold;


	// now we load a layer of particles
	// to do this, the constructor for a sediment layer is used that takes a vector of masses...
	// these particles are sand
	// first get the number of particle types (you can look at the "dpart.dplist" file to find out
	// the particles
	double initial_mass = 25.0;			// all masses are in kg per unit surface area (e.g., a
						// marsh column with aerial extent of 1m^2)
	vector<double> initial_masses;

	// now load the sand into the model
	for (int i = 0; i<np_types; i++)
	{
	if (i == 0)
		initial_masses.push_back(initial_mass);
	else
		initial_masses.push_back(0.0);
	}

	// now intialize the layer
	//cout <<"LINE 303 main initializing layer\n";
	sediment_layer sl(dpi,initial_masses);
	//cout <<"LINE 303 main initialized layer\n";
	//cout << "LINE 354 particle information:"<< endl;
	//sl.print_particle_info_to_screen();

	// lets initialize the stack. To create a stack object we give it the mean sea level
	// this mean sea level is above a datum...assumed to be the base of the sediment column
	sediment_stack sed_stack(Mean_Tide);

	// now add some layers to the sediment stack
	for (int i = 0; i< 100; i++)
	sed_stack.deposit_surface_sediment(sl);

	//cout << "LINE 365 printing layer 300: " << endl;
	//sed_stack.print_ind_particles_of_layer_to_screen(300);
	//cout << "LINE 367 printing layer 100: " << endl;
	//sed_stack.print_ind_particles_of_layer_to_screen(100);

	// compression
	// now calcualte the thickness and porosity of the layers
	sed_stack.calculate_overburden_and_thickness();

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

/*************************************************
**
** NOW FOR THE TIME LOOP
**
***********************************************************/

	t_ime = 0;						// initilaize the time


	// first, get the top layer of the initial pile of sediment
	annual_band_layer_top.push_back(sed_stack.get_n_layers()-1);
	cout << "LINE 359 starting n_layers: " << sed_stack.get_n_layers() << endl;

	// start a loop
	accretion_ratio = 10;
	int dead_biomass_counter = 0;
	int low_acc_ratio_counter = 0;
	int yr = 0;
	//while(yr < end_year && dead_biomass_counter < 10 && (accretion_ratio < .99 || accretion_ratio > 1.01 ))
	//while (yr<end_year && dead_biomass_counter<10)	//MK- I commented out these two lines?
	while (yr<end_year)
	{
		yr++;
		t_ime = yr*365;
		yr_time = t_ime/365;

		SLR= slvector[yr]-slvector[yr-1]; //MK- I added this line to cacluate the rate of sea level rise from the loaded file
		//cout << " Yr: "<< yr << " Sl: " << slvector[yr] << " Sl yr-1: " << slvector[yr-1] << " SLRR: " << SLR << endl;
		temperatureincrease=temperaturevector[yr]-temperaturevector[0]; //MK- I added this line too

		// reset the biomass
		Bmass = 0;

		// recalcualte the growth index. see notes in sediment_stack.cpp
		vector<double> growth_index = sed_stack.create_growth_index_inf(gamma_bg_biomass);
		//cout << " LINE 346 main size of growth index is: " << growth_index.size() << endl;

		marsh_surface_elevation = sed_stack.get_marsh_surface_elevation();
		marsh_surface_elevation_day0 = marsh_surface_elevation;
		//cout << "LINE 397 main marsh_surface is: " << marsh_surface_elevation << endl;

		// calcualte t./column.out
		//the peak biomass for the year
		peak_Bmass =  get_peak_biomass(marsh_surface_elevation, MHT,
		   max_depth, min_depth, max_bmass,temperatureincrease,yr,boriginal,bfactor);	//MK- I added temperature,yr,boriginal to this line

		if (yr==1) // MK- I added these 2 lines
		{boriginal=peak_Bmass;}

		//cout << "LINE 400 main peak_Bmass is: " << peak_Bmass
		//     << " and water depth is: " <<  MHT-marsh_surface_elevation << endl;

		min_bmass = peak_Bmass*theta_bmin;
		max_growth = nu_Gp*peak_Bmass;
		min_growth = nu_Gmin*peak_Bmass;
		//cout << "Line 416 peak bmass: " << peak_Bmass << " min_Bmass: " << min_bmass << endl
		//     << "max_growth: " << max_growth << " min_growth: "<< min_growth << endl;


		particle_masses = sed_stack.get_total_stack_mass_of_each_type(np_types);
		particle_masses_day0 = sed_stack.get_total_stack_mass_of_each_type(np_types);
		//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
		// now go day by day, calucalating the biomass, root growth and death
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

		ab_pointer_sz = annual_band_layer_top.size();

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


	}		// !! end annual loop

   cout << " ENDING THE LOOP, SLR: " << SLR << endl;

   // now print the data
   data_out << SLR << " " << " " << Tidal_Amplitude << " "
			<< start_depth_frac << " " << yr << " "
			<< conc_silt << " " << effective_svel << " "
			<< MHT << " " << MHT-marsh_surface_elevation << " "
			<< max_depth << " "
			<< marsh_surface_elevation-marsh_surface_elevation_day0
			<< " " <<peak_Bmass;

   // print the accretion rate for each particle type
   for (int i = 0; i< np_types; i++)
   {
	   data_out << " " << particle_masses[i]-particle_masses_day0[i];
   }
   data_out << endl;

	ab_pointer_sz = annual_band_layer_top.size();
	//cout << "n annual bands: " << ab_pointer_sz << endl;
	vector<double> pm_band;
	for (int i = 1; i<ab_pointer_sz; i++)
	{
		col_out << i << " " << sed_stack.get_layer_top_elevation(annual_band_layer_top[i]);
		pm_band = sed_stack.get_stack_slice_mass_of_each_type(np_types,
											annual_band_layer_top[i-1]+1,
									annual_band_layer_top[i]);

		for (int type = 0; type < np_types; type++)
		{
			col_out << " " << pm_band[type];
		}
		col_out << endl;



	}


   	// MK- New function to save time series data
	ofstream series_out;
	string fname2_prefix="series.";
	string fname2_suffix=".txt";
	string fname2_k=itoa(0);
	string fname2_b=itoa(0);
		if (kfactor!=0)
		{fname2_k=itoa(1);}
		if (bfactor!=0)
		{fname2_b=itoa(1);}
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


   return peak_Bmass;


}

/********************************************************
// this routine gives the total mass deposited on the surface of the marsh
// during one tidal period.
// particle concentrations are in kg/m^3
******************************************************/
vector<double> get_trapping_settling(double Tidal_Period, double Tidal_Amplitude,
                           double Mean_Tide, double marsh_surface_elevation,
                           vector<double> particle_concentrations_0,
                           vector<double> particle_diameters,
                           vector<double> particle_settling_velocities,
                           double flow_velocity,
                           double Biomass, double alpha_T, double beta_T,
                           double epsilon, double gamma_tr)
 {

  double Time_fractions = 10000; 	// the number of time steps
  					// in the simulation over one tidal
  					// cycle
  double dt = Tidal_Period/Time_fractions;
  double start_time;
  double water_depth;

  //cout << "LINE 727, depth_below_MHT: "
  //     << Mean_Tide + Tidal_Amplitude - marsh_surface_elevation << endl;

  int n_particle_types = particle_diameters.size();
  vector<double> trapping_mass(n_particle_types,0.0);
  vector<double> settling_mass(n_particle_types,0.0);
  vector<double> trapping_multiplier(n_particle_types);
  double ts_settling_mass;
  double ts_trapping_mass;

  //cout << "dt is: " << dt << endl;

  for (int i = 0; i<n_particle_types; i++)
   {
    trapping_multiplier[i] = 3600*alpha_T*pow(flow_velocity,gamma_tr+1)*
  				pow(particle_diameters[i],epsilon)*
  				particle_concentrations_0[i]*
  				pow(Biomass,beta_T)*dt;
  				// the 3600 is because the rest of the
                                // equation is in s but dt is in hours
    //cout << "multiplier is: " << trapping_multiplier[i] << endl;
   }

/****************************************************************
**
  if (marsh_surface_elevation > Mean_Tide+Tidal_Amplitude)
   {
    start_time = Tidal_Period;
   }
  else if (marsh_surface_elevation < Mean_Tide-Tidal_Amplitude)
   {
    start_time = 0;
   }
  else
   {
    start_time = Tidal_Period/(2*M_PI)*
                   asin((marsh_surface_elevation-Mean_Tide)/Tidal_Amplitude) +
                   Tidal_Period/4;
   }
**
****************************************************************/

  double t_ime = 0;
  double time_submerged = 0;
  while (t_ime < Tidal_Period)
   {
    t_ime += dt;
    water_depth = Tidal_Amplitude*sin(2*M_PI*(t_ime/Tidal_Period - 0.25))+
                  Mean_Tide - marsh_surface_elevation;
    if (water_depth < 0)
     water_depth = 0;
    else
     time_submerged += dt;

    //cout << "time is: " << t_ime << " and water depth is: " << water_depth << endl;
    for (int i = 0; i<n_particle_types; i++)
     {
      if (water_depth != 0)
       {
        ts_settling_mass = 3600*particle_concentrations_0[i]*dt*
                         particle_settling_velocities[i];
                         // the 3600 is becuase the settling velocity is in m/s
                         // and dt is in hours
       }
      else ts_settling_mass = 0;
      settling_mass[i] += ts_settling_mass;

      ts_trapping_mass = trapping_multiplier[i]*water_depth;
      trapping_mass[i] += ts_trapping_mass;
     }

   }
  //cout << "conc: " << particle_concentrations_0[0] << " w_s: " << particle_settling_velocities[0] << endl;
  //cout << "settling_mass: " << settling_mass[0] << " and trapping mass: " << trapping_mass[0] << endl;
  vector<double> total_mass(n_particle_types);
  //cout << "time submerged: " << time_submerged << " settling_mass: " << settling_mass[0] << endl;
  for (int i = 0; i<n_particle_types; i++)
   {
    total_mass[i] = trapping_mass[i]+settling_mass[i];
    //cout << "LINE 682 main: mass["<<i<<"]: " << total_mass[i] << " settling: "
    //     << settling_mass[i] << " trapping: " << trapping_mass[i] << endl;
   }


  return total_mass;
 }


/********************************************************
// this routine gives the total mass deposited on the surface of the marsh
// during one tidal period.
// particle concentrations are in kg/m^3
******************************************************/
vector<double> get_trapping_settling_mudd(double Tidal_Period, double Tidal_Amplitude,
                           double Mean_Tide, double marsh_surface_elevation,
                           vector<double> particle_concentrations_0,
                           vector<double> particle_diameters,
                           vector<double> particle_settling_velocities,
                           double flow_velocity,
                           double Biomass, double alpha_T, double beta_T,
                           double epsilon, double gamma_tr)
 {
	double dmin = 0.001;			// minimum depth for settling and trapping
									// necessary becasue in sediment continuity
									//
  	double Time_fractions = 2000; 	// the number of time steps
  					// in the simulation over one tidal
  					// cycle
  	double dt = Tidal_Period/Time_fractions;
  	double start_time;
  	double water_depth;
  	double dwater_depth_dt;			// derivative of the water depth
  	double flow_vel;

  	//cout << "LINE 727, depth_below_MHT: "
  	//     << Mean_Tide + Tidal_Amplitude - marsh_surface_elevation << endl;

  	int n_particle_types = particle_diameters.size();
  	vector<double> trapping_mass(n_particle_types,0.0);
  	vector<double> settling_mass(n_particle_types,0.0);
  	vector<double> trapping_multiplier(n_particle_types);
  	vector<double> particle_conc = particle_concentrations_0;
  	double ts_settling_mass;
  	double ts_trapping_mass;
  	vector<double> ts_conc_loss_settling(n_particle_types,0.0);
  	vector<double> ts_conc_loss_trapping(n_particle_types,0.0);
  	vector<double> starting_mass(n_particle_types,0.0);

  	//cout << "dt is: " << dt << endl;

  	for (int i = 0; i<n_particle_types; i++)
  	{
  	  	trapping_multiplier[i] = 3600*alpha_T*
  				pow(particle_diameters[i],epsilon)*
  				pow(Biomass,beta_T)*dt;
  				// the 3600 is because the rest of the
                                // equation is in s but dt is in hours
    	//cout << "multiplier is: " << trapping_multiplier[i] << endl;
   	}

  	double t_ime = 0;
  	double time_submerged = 0;
  	double water_depth_old;
  	double diff_dw_dt;
    water_depth = 0;
  	while (t_ime < Tidal_Period)
   	{
    	t_ime += dt;

    	//water_depth_old = water_depth;
    	water_depth = Tidal_Amplitude*sin(2*M_PI*(t_ime/Tidal_Period - 0.25))+
    	              Mean_Tide - marsh_surface_elevation;

    	//cout << "LINE 141 tims: " << t_ime << " TA: " << Tidal_Amplitude
    	//	<< " mse: " << marsh_surface_elevation << " MT: " << Mean_Tide << " wd: " << water_depth << endl;
    	//diff_dw_dt = (water_depth-water_depth_old)/dt;
        dwater_depth_dt =   2*M_PI*Tidal_Amplitude*cos(2*M_PI*(t_ime/Tidal_Period - 0.25))/Tidal_Period;

        //cout << "LINE 144 dwdt diff: " << diff_dw_dt << " and analytical: " << dwater_depth_dt << endl;

		// reset the concentration loss vectors
		for (int i = 0; i<n_particle_types; i++)
     	{
			ts_conc_loss_settling[i] = 0;
			ts_conc_loss_trapping[i] = 0;
		}

		if (water_depth >dmin)
      	{
			// if the tide is coming in, water has a starting concntration of C_0
			if (dwater_depth_dt >0)
			{
				for (int i = 0; i<n_particle_types; i++)
     			{
					starting_mass[i] = particle_concentrations_0[i]*water_depth;

					//cout << "LINE 463, part conc["<<i<<"]: " << particle_conc[i]
					//     << " and wd is: " << water_depth<< endl;
					ts_settling_mass = 3600*particle_concentrations_0[i]*dt*
					                         particle_settling_velocities[i];
					                         // the 3600 is becuase the settling velocity is in m/s
                         					 // and dt is in hours
                    ts_trapping_mass = trapping_multiplier[i]*particle_concentrations_0[i]
                    					*pow(flow_velocity,gamma_tr+1)
  										*water_depth;

                    if (starting_mass[i] < ts_trapping_mass)
                    {
						ts_trapping_mass = starting_mass[i];
						ts_settling_mass = 0;
					}
					else if ( starting_mass[i] -ts_trapping_mass < ts_settling_mass)
					{
						ts_settling_mass = starting_mass[i] -ts_trapping_mass;
					}

					trapping_mass[i] += ts_trapping_mass;
                    settling_mass[i] += ts_settling_mass;


					//cout << "LINE 485, ts_sett_mass["<<i<<"]: " << ts_settling_mass << endl;
					//cout << "LINE 486, ts_trap_mass["<<i<<"]: " << ts_trapping_mass << endl;
                    //cout << "LINE 487 ts: " <<t_ime << " part_conc["<<i<<"]: " << particle_conc[i] <<endl;
				}
			}
			// if the tide is going out, the concentration can get depleted
			else
			{
				for (int i = 0; i<n_particle_types; i++)
     			{
					starting_mass[i] = particle_conc[i]*water_depth;

					//cout << "LINE 497, part conc["<<i<<"]: " << particle_conc[i]
					//     << " and wd is: " <<water_depth
					//     << " and starting mass: " << starting_mass[i]  <<  endl;
					ts_settling_mass = 3600*particle_conc[i]*dt*
					                         particle_settling_velocities[i];
					                         // the 3600 is becuase the settling velocity is in m/s
                         					 // and dt is in hours
                    ts_trapping_mass = trapping_multiplier[i]*particle_conc[i]
                    					*pow(flow_velocity,gamma_tr+1)
  										*water_depth;

                    if (starting_mass[i] < ts_trapping_mass)
                    {
						ts_trapping_mass = starting_mass[i];
						ts_settling_mass = 0;
					}
					else if ( starting_mass[i] -ts_trapping_mass < ts_settling_mass)
					{
						ts_settling_mass = starting_mass[i] -ts_trapping_mass;
					}

					trapping_mass[i] += ts_trapping_mass;
                    settling_mass[i] += ts_settling_mass;


					//cout << "LINE 519, ts_sett_mass["<<i<<"]: " << ts_settling_mass << endl;
					//cout << "LINE 520, ts_trap_mass["<<i<<"]: " << ts_trapping_mass << endl;

                    particle_conc[i] = particle_conc[i]
                    		-dt*particle_conc[i]*dwater_depth_dt/water_depth
                    		-(ts_settling_mass-ts_trapping_mass)/water_depth;
                    if (particle_conc[i] < 0)
                     particle_conc[i] = 0;

                    //cout << "LINE 528 ts: " <<t_ime << " part_conc["<<i<<"]: " << particle_conc[i] <<endl;

				}
			}

		}
		//cout << "time is: " << t_ime << " and water depth is: " << water_depth << endl;



	}
  	//cout << "conc: " << particle_concentrations_0[0] << " w_s: " << particle_settling_velocities[0] << endl;
  	//cout << "settling_mass: " << settling_mass[0] << " and trapping mass: " << trapping_mass[0] << endl;
  	vector<double> total_mass(n_particle_types);
  	//cout << "time submerged: " << time_submerged << " settling_mass: " << settling_mass[0] << endl;
  	for (int i = 0; i<n_particle_types; i++)
  	{
    	total_mass[i] = trapping_mass[i]+settling_mass[i];
    	//cout << "LINE 682 main: mass["<<i<<"]: " << total_mass[i] << " settling: "
    	//     << settling_mass[i] << " trapping: " << trapping_mass[i] << endl;
   	}


  	return total_mass;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=




/********************************************************
// this returns the peak biomass is g/m^2
******************************************************/
double get_peak_biomass(double surface_elevation, double MHT,
		   double max_depth, double min_depth, double max_bmass, double temperatureincrease, int yr, double boriginal, double bfactor)
{
  double depth_range = max_depth-min_depth;
  double water_depth = MHT-surface_elevation;
  double B_ps;

  //double bfactor=.06; // MK- I added this line. Set to 0.06 for standard effect, Set to zero to turn off temperature effect.


  //cout << "LINE 688 MHT: " << MHT << " z: " << surface_elevation
  //     << " d: " << water_depth << " max_MB: " << max_bmass << endl;
  //cout << "LINE 690 multiplier: " << (water_depth - min_depth)/(depth_range) << endl;
  if (water_depth > max_depth)
   B_ps = 0;
  else if (water_depth < min_depth)
   B_ps = 0;
  else
    {
    B_ps = max_bmass*(water_depth - min_depth)/(depth_range);

    cout << "year: " << yr << " B before temp: "<< B_ps << " Temp inc: " << temperatureincrease<<"       "; //MK- I added these 4 lines
    B_ps = B_ps+(temperatureincrease*bfactor*B_ps);	// MK- temperature changes peak biomass by same % regardless of ambient biomass
    //B_ps = B_ps + (temperatureincrease*bfactor*boriginal); //MK- temperature changes biomass by same absolute amount (a function of initial biomass)
    //cout<<" B after temp: "<< B_ps <<endl; //MK
    }

  return B_ps;
 }

/********************************************************
// this returns the peak biomass is g/m^2
******************************************************/
double get_peak_biomass_parab(double surface_elevation, double MHT,
		   double max_depth, double min_depth, double max_bmass)
{
  double depth_range = max_depth-min_depth;
  double water_depth = MHT-surface_elevation;
  double B_ps;

  //cout << "LINE 688 MHT: " << MHT << " z: " << surface_elevation
  //     << " d: " << water_depth << " max_MB: " << max_bmass << endl;
  //cout << "LINE 690 multiplier: " << (water_depth - min_depth)/(depth_range) << endl;
  if (water_depth > max_depth)
   B_ps = 0;
  else if (water_depth < min_depth)
   B_ps = 0;
  else
   B_ps = max_bmass*(water_depth-min_depth)*(max_depth-water_depth);

  //cout << "Line 698 min_depth: " << min_depth << " B_ps: " << B_ps << endl;
  return B_ps;
 }

string itoa(int num)
{
    stringstream converter;
    converter << num;
    return converter.str();
}


/********************************************************
// this function gives the growth
// this is the change in biomass (aboveground) during the
// interval dt for derivation, see notes_on_growth.nb
// returns a vector<double> where
// element[0] = biomass	(in mass)
// element[1] = mortality (in mass)
******************************************************/
vector<double> biomass_and_mortality(double day, double peak_biomass, double min_biomass,
		      double day_peak_season, double peak_growth, double min_growth,
		      double dt)
 {
  //cout << "line 724 peak_biomass: " << peak_biomass << endl;

  vector<double> bm_and_mort(2);
  //if (peak_biomass < min_biomass)
  double biomass = 0.5*(min_biomass+peak_biomass+(peak_biomass-min_biomass)*
                        cos( 2*M_PI*(day-day_peak_season)/365) );
  double mortality = (1/(4*M_PI))*(2*dt*(min_growth+peak_growth)*M_PI+
                      2*(min_biomass-peak_biomass)*M_PI*
                      (cos(2*M_PI*(day-day_peak_season)/365) -
                       cos(2*M_PI*(dt-day+day_peak_season)/365)) -
                      365*(min_growth-peak_growth)*
                      (sin(2*M_PI*(day-day_peak_season)/365) +
                       sin(2*M_PI*(dt-day+day_peak_season)/365)));

   //cout << "LINE 838 tidal_colum2 day: " << day << " , bmass: " << biomass
   //     << " and mort: " << mortality << endl;
   bm_and_mort[0] = biomass;
   bm_and_mort[1] = mortality;

   if (mortality < 0)
    cout << "LINE 754 mortality is neagitve! " << bm_and_mort[1] << endl;



   return bm_and_mort;
 }

/********************************************************
// this function gives the growth
// this is the change in biomass (aboveground) during the
// interval dt for derivation, see notes_on_growth.nb
// returns a vector<double> where
// element[0] = biomass	(in mass)
// element[1] = mortality (in mass)
******************************************************/
vector<double> biomass_and_mortality2(double day, double peak_biomass, double min_biomass,
		      double day_peak_season, double peak_growth, double min_growth,
		      double dt,double phase_shift)
 {
  //cout << "line 724 peak_biomass: " << peak_biomass << endl;

  vector<double> bm_and_mort(2);
  //if (peak_biomass < min_biomass)
  double biomass = 0.5*(min_biomass+peak_biomass+(peak_biomass-min_biomass)*
                        cos( 2*M_PI*(day-day_peak_season)/365) );
  double mortality = (1/(4*M_PI))*(2*dt*(min_growth+peak_growth)*M_PI+
                      2*(min_biomass-peak_biomass)*M_PI*
                      (cos(2*M_PI*(day-day_peak_season)/365) -
                       cos(2*M_PI*(dt-day+day_peak_season)/365)) -
                      365*(min_growth-peak_growth)*
                      (sin(2*M_PI*(day-day_peak_season+phase_shift)/365) +
                       sin(2*M_PI*(dt-day+day_peak_season-phase_shift)/365)));

   //cout << "LINE 838 tidal_colum2 day: " << day << " , bmass: " << biomass
   //     << " and mort: " << mortality << endl;
   bm_and_mort[0] = biomass;
   bm_and_mort[1] = mortality;

   if (mortality < 0)
    cout << "LINE 754 mortality is negative! " << bm_and_mort[1] << endl;



   return bm_and_mort;
 }



//***********Simon- Functions from here to end are new, part of the TKE code
double get_slope(double aa, double bb, double alpha_turb, double alpha_zero,
                                 double alpha, double beta, double mu, double phi,
                                 double nu, double g, double B, double u)
{
        double S;

        // for a given veolicty, what is the slope?
        double term1 = aa*pow(B,beta)*u*u*alpha;
        double term2 = 0.25*pow(B,2*beta+phi)*bb*M_PI*u*u*alpha*alpha*mu;
        double term3 = pow(B,beta-phi)*u*alpha*alpha_zero*nu/mu;
        double term4 = (1-0.25*alpha*mu*M_PI*pow(B,beta+phi))*g;

        S = -(term1+term2+term3)/term4;
        return S;
}

double get_u(double aa, double bb, double alpha_turb, double alpha_zero,
                                 double alpha, double beta, double mu, double phi,
                                 double nu, double g, double B, double S)
{
        double u;

        // now get the velocity for a given slope
        double term5 = -4*pow(B,beta)*alpha*alpha_zero*nu;
        double term6 = -4*pow(B,phi)*g*S*mu+pow(B,beta+2*phi)*g*M_PI*S*alpha*mu*mu;
        double term7 = -4*aa*pow(B,beta+phi)*alpha*mu-pow(B,2*beta+2*phi)*bb*M_PI*alpha*alpha*mu*mu;
        double term8 = 16*pow(B,2*beta)*alpha*alpha*alpha_zero*alpha_zero*nu*nu;
        double term9 = 2*(4*aa*pow(B,beta+phi)*alpha*mu+pow(B,2*beta+2*phi)*bb*M_PI*alpha*alpha*mu*mu);

        u = (term5+sqrt(-4*term6*term7+term8))/term9;
        return u;
}

double get_k(double aa, double bb, double alpha_turb, double alpha_zero,
                                 double alpha, double beta, double mu, double phi,
                                 double nu, double g, double B, double u)
{
        double k;
        double term1 = alpha_zero*nu/(u*mu*pow(B,phi));
        double term2 = aa+bb*alpha*mu*M_PI*0.25*pow(B,beta+phi);
        double term3 = alpha*mu*pow(B,beta+phi);

        k = alpha_turb*alpha_turb*u*u*pow(2*(term1+term2)*term3,2/3);
        return k;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this routine gives the total mass deposited on the surface of the marsh
// during one tidal period.
// particle concentrations are in kg/m^3
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<double> get_trap_TKE_eff_sett(double Tidal_Period, double Tidal_Amplitude,
                           double Mean_Tide, double marsh_surface_elevation,
                           vector<double> particle_concentrations_0,
                           vector<double> particle_diameters,
                           vector<double> particle_settling_velocities,
                           double peak_flow_velocity,
                           double Biomass, double alpha_T, double beta_T,
                           double epsilon, double gamma_tr,
                           double aa, double bb, double alpha_turb, double alpha_zero,
                                                         double alpha, double beta, double mu, double phi,
                                                         double nu, double g,
                           vector<double>& smass, vector<double>& tmass)
 {
        int counter = 0;
        double dmin = 0.002;                                // minimum depth for settling and trapping
                                                                        // necessary becasue in sediment continuity
                                                                        //
          double Time_fractions = 15000;         // the number of time steps
                                          // in the simulation over one tidal
                                          // cycle
          double dt = Tidal_Period/Time_fractions;
          double start_time;
          double water_depth;
          double dwater_depth_dt;                        // derivative of the water depth
          double flow_vel;
          double TKE;                                                // turbulent kinetic energy
          double Slope;                                        // surface water slope

          double von_karman = 0.4;                // von karman constant
          double rho_w = 1000;                        // density of water
          double TKE_coefficient = 0.2;        // coeffieicent relating TKE to
                                                                          // shear stress
          double w_up;                                        // upward velocity due to turbulence

          //cout << "LINE 727, depth_below_MHT: "
          //     << Mean_Tide + Tidal_Amplitude - marsh_surface_elevation << endl;

          int n_particle_types = particle_diameters.size();
          vector<double> trapping_mass(n_particle_types,0.0);
          vector<double> settling_mass(n_particle_types,0.0);
          vector<double> settling_mass_no_reduc(n_particle_types,0.0);
          vector<double> trapping_multiplier(n_particle_types);
          vector<double> particle_conc(n_particle_types,0.0);
          double ts_settling_mass;
          double ts_settling_mass_no_reduc;
          double ts_trapping_mass;
          vector<double> ts_conc_loss_settling(n_particle_types,0.0);
          vector<double> ts_conc_loss_trapping(n_particle_types,0.0);
          vector<double> starting_mass(n_particle_types,0.0);
          vector<double> mass_present(n_particle_types,0.0);

          //cout << "dt is: " << dt << endl;

          for (int i = 0; i<n_particle_types; i++)
          {
                    trapping_multiplier[i] = 3600*alpha_T*
                                  pow(particle_diameters[i],epsilon)*
                                  pow(Biomass,beta_T)*dt;
                                  // the 3600 is because the rest of the
                                // equation is in s but dt is in hours
            //cout << "multiplier is: " << trapping_multiplier[i] << endl;
            particle_conc[i] = particle_concentrations_0[i];
           }

        double water_depth_amp = 2*Tidal_Amplitude*M_PI/Tidal_Period;
          double t_ime = 0;
          double time_submerged = 0;
          double water_depth_old;
          double diff_dw_dt;
          double diff_flow_vel;
          double mass_in;
          double mass_out;
    water_depth = 0;

    int ws_reduc_counter = 0;
    double ws_mean_reduc = 0;

          while (t_ime < Tidal_Period)
           {
            t_ime += dt;

            water_depth_old = water_depth;
            water_depth = Tidal_Amplitude*sin(2*M_PI*(t_ime/Tidal_Period - 0.25))+
                          Mean_Tide - marsh_surface_elevation;
        dwater_depth_dt =   2*M_PI*Tidal_Amplitude*cos(2*M_PI*(t_ime/Tidal_Period - 0.25))/Tidal_Period;
                flow_vel = peak_flow_velocity*fabs(dwater_depth_dt/water_depth_amp);

                TKE = get_k(aa, bb, alpha_turb, alpha_zero,
                                 alpha, beta, mu, phi, nu, g, Biomass, flow_vel);
                //Slope = get_slope(aa, bb, alpha_turb, alpha_zero,
                //                 alpha, beta, mu, phi, nu, g, Biomass, flow_vel);



                w_up = von_karman*sqrt(TKE_coefficient*TKE/rho_w);
                //cout << "TKE is: " << TKE << " and w_up is: " << w_up << endl;

                //ws_mean_reduc += (particle_settling_velocities[0]-w_up)/
                                                        particle_settling_velocities[0];
                //ws_reduc_counter++;

                //cout << "LINE 352, flow velocity is: " << flow_vel << endl;

                // reset the concentration loss vectors
                for (int i = 0; i<n_particle_types; i++)
             {
                        ts_conc_loss_settling[i] = 0;
                        ts_conc_loss_trapping[i] = 0;
                }

                if (water_depth >dmin)
              {
                        // if the tide is coming in, water has a starting concntration of C_0
                        if (dwater_depth_dt >0)
                        {
                                for (int i = 0; i<n_particle_types; i++)
                             {
                                        time_submerged+=dt;
                                        starting_mass[i] = particle_concentrations_0[i]*water_depth;

                                        //cout << "LINE 463, part conc["<<i<<"]: " << particle_conc[i]
                                        //     << " and wd is: " << water_depth<< endl;
                                        ts_settling_mass = 3600*particle_concentrations_0[i]*dt*
                                                                 (particle_settling_velocities[i]-w_up);
                                        if (ts_settling_mass < 0)
                                        {
                                                ts_settling_mass = 0;
                                        }
                                        ts_settling_mass_no_reduc = 3600*particle_concentrations_0[i]*dt*
                                                                 (particle_settling_velocities[i]);
                                                                 // the 3600 is becuase the settling velocity is in m/s
                                                                  // and dt is in hours
                    ts_trapping_mass = trapping_multiplier[i]*particle_concentrations_0[i]
                                                            *pow(flow_vel,gamma_tr+1)
                                                                                  *water_depth;

                    if (starting_mass[i] < ts_trapping_mass)
                    {
                                                ts_trapping_mass = starting_mass[i];
                                                ts_settling_mass = 0;
                                        }
                                        else if ( starting_mass[i] -ts_trapping_mass < ts_settling_mass)
                                        {
                                                ts_settling_mass = starting_mass[i] -ts_trapping_mass;
                                        }

                                        settling_mass_no_reduc[i] += ts_settling_mass_no_reduc;
                                        trapping_mass[i] += ts_trapping_mass;
                    settling_mass[i] += ts_settling_mass;


                                        //cout << "LINE 485, ts_sett_mass["<<i<<"]: " << ts_settling_mass << endl;
                                        //cout << "LINE 486, ts_trap_mass["<<i<<"]: " << ts_trapping_mass << endl;
                    //cout << "LINE 487 ts: " <<t_ime << " part_conc["<<i<<"]: " << particle_conc[i] <<endl;
                                }
                        }
                        // if the tide is going out, the concentration can get depleted



                        else
                        {
                                for (int i = 0; i<n_particle_types; i++)
                             {
                                        starting_mass[i] = particle_conc[i]*water_depth;

                                        //if (counter == 0)
                                        //{
                                        //        cout << "LINE 497, part conc["<<i<<"]: " << particle_conc[i]
                                        //             << " and wd is: " <<water_depth
                                        //             << " and starting mass: " << starting_mass[i]  <<  endl;
                                        //}
                                        ts_settling_mass = 3600*particle_conc[i]*dt*
                                                                 (particle_settling_velocities[i]-w_up);
                                                                 // the 3600 is becuase the settling velocity is in m/s
                                                                  // and dt is in hours
                    if (ts_settling_mass < 0)
                                        {
                                                ts_settling_mass = 0;
                                        }
                    ts_trapping_mass = trapping_multiplier[i]*particle_conc[i]
                                                            *pow(flow_vel,gamma_tr+1)
                                                                                  *water_depth;

                    if (starting_mass[i] < ts_trapping_mass)
                    {
                                                ts_trapping_mass = starting_mass[i];
                                                ts_settling_mass = 0;
                                        }
                                        else if ( starting_mass[i] -ts_trapping_mass < ts_settling_mass)
                                        {
                                                ts_settling_mass = starting_mass[i] -ts_trapping_mass;
                                        }

                                        trapping_mass[i] += ts_trapping_mass;
                    settling_mass[i] += ts_settling_mass;


                                        //cout << "LINE 519, ts_sett_mass["<<i<<"]: " << ts_settling_mass << endl;
                                        //cout << "LINE 520, ts_trap_mass["<<i<<"]: " << ts_trapping_mass << endl;

                    particle_conc[i] = particle_conc[i]
                                    -dt*particle_conc[i]*dwater_depth_dt/water_depth
                                    -(ts_settling_mass+ts_trapping_mass)/water_depth;
                    if (particle_conc[i] < 0)
                     particle_conc[i] = 0;

                    //if (particle_conc[i] > 0)
                    //        cout << "LINE 528 ts: " <<t_ime << " part_conc["<<i<<"]: " << particle_conc[i] <<endl;

                           if (counter == 0)
                            counter ++;

                                }
                        }


                }
                //cout << "time is: " << t_ime << " and water depth is: " << water_depth << endl;




        }

        //cout << "time submerged is: " << time_submerged << endl;
        //cout << "w_mean_reduc = " << ws_mean_reduc/double(ws_reduc_counter) << endl;
        //cout << "settling: " << settling_mass[0] << " settling no reduc: " << settling_mass_no_reduc[0] << endl;

          //cout << "conc: " << particle_concentrations_0[0] << " w_s: " << particle_settling_velocities[0] << endl;
          //cout << "settling_mass: " << settling_mass[0] << " and trapping mass: " << trapping_mass[0] << endl;
          vector<double> total_mass(n_particle_types);
          //cout << "time submerged: " << time_submerged << " settling_mass: " << settling_mass[0] << endl;
          for (int i = 0; i<n_particle_types; i++)
          {
            total_mass[i] = trapping_mass[i]+settling_mass[i];
            //cout << "LINE 682 main: mass["<<i<<"]: " << total_mass[i] << " settling: "
            //     << settling_mass[i] << " trapping: " << trapping_mass[i] << endl;
           }

          tmass = trapping_mass;
          smass = settling_mass;
          //cout << "smass: " << smass[0] << " and tmass: " << tmass[0] << endl;
          return total_mass;


}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


// this function calculates the settling velocity
double calculate_w_s(double d_p)
{
        double w_s;
        double A = 38.0;
        double F = 3.55;
        double m = 1.12;
        double g = 9.80;
        double rho_s = 2600;
        double rho_w = 1000;
        double s = rho_s/rho_w;
        double nu = pow(10,-6);
        double term1,term2,term3,term4,term5,term6;

        term1 = 0.25*pow(A/F,2.0/m);
        term2 = 4*d_p*d_p*d_p*g*(s-1)/(3*F*nu*nu);
        term3 = pow(term2,1/m);
        term4 = sqrt(term1+term3);
        term5 = 0.5*pow(A/F,1/m);
        term6 = pow((term4-term5),m);
        w_s = nu*term6/d_p;

        return w_s;
}



