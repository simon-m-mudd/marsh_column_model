//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDStatsTools
// Land Surface Dynamics StatsTools
//
// A collection of statistical routines for use with the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//
// Developed by:
//  Simon M. Mudd
//  Martin D. Hurst
//  David T. Milodowski
//  Stuart W.D. Grieve
//  Declan A. Valters
//  Fiona Clubb
//
// Copyright (C) 2013 Simon M. Mudd 2013
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
// either version 2 of the License, or (at your option) any later version.
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
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//-----------------------------------------------------------------
//DOCUMENTATION URL: http://www.geos.ed.ac.uk/~s0675405/LSD_Docs/
//-----------------------------------------------------------------

#include <vector>
#include <map>
#include <math.h>
#include "TNT/tnt.h"
using namespace std;
using namespace TNT;



// Initiates a stack of sediment made exclusively of sand. 
// note that it is assumed the depo_particle is initiated with a sand fraction in the 
// first element of particle sizes. 
// the column will actually be made of whatever particle size is in the first element
sediment_stack initiate_column(depo_particle_info dpi,double intitial_mass, double mean_tide);

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
		   double max_depth, double min_depth, double max_bmass, double temperatureincrease, int yr, double bfactor);

double get_peak_biomass_parab(double surface_elevation, double MHT,
		   double max_depth, double min_depth, double max_bmass);

// this gets the biomass growth and mortality for a given time of year
vector<double> biomass_and_mortality(double day, double peak_biomass, double min_biomass,
		      double day_peak_season, double peak_growth, double min_growth,
		      double dt);

// similar to above but includes a phase shift. 
vector<double> biomass_and_mortality2(double day, double peak_biomass, double min_biomass,
		      double day_peak_season, double peak_growth, double min_growth,
		      double dt,double phase_shift);

#ifndef MarshUtil_H
#define MarshUtil_H




#endif
