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
#include <fstream>
#include <string>
#include <algorithm>
#include <functional> // For string splitter
#include <iostream>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <ctime>
#include <map>
#include <list>
#include <math.h>
#include "TNT/tnt.h"
#include "TNT/jama_lu.h"
#include "LSDStatsTools.hpp"
#include "sediment_stack.hpp"
#include "sediment_layer.hpp"
#include "depo_particle.hpp"

using namespace std;
using namespace TNT;
using namespace JAMA;

#ifndef MarshUtil_CPP
#define MarshUtil_CPP



sediment_stack initiate_column(depo_particle_info dpi, double intitial_mass, double mean_tide)
{
	// now we load a layer of particles
	// to do this, the constructor for a sediment layer is used that takes a vector of masses...
	// these particles are sand
	// first get the number of particle types (you can look at the "dpart.dplist" file to find out
	// the particles
	// all masses are in kg per unit surface area (e.g., a
	// marsh column with aerial extent of 1m^2)
  vector<double> initial_masses;

	int np_types = dpi.get_n_types();			// number of particle types

	// now load the sand into the model
	for (int i = 0; i<np_types; i++)
	{
	  if (i == 0)
		{
      initial_masses.push_back(intitial_mass);
    }
	  else
    {
		  initial_masses.push_back(0.0);
    }
	}

	// now intialize the layer
	//cout <<"LINE 303 main initializing layer\n";
	sediment_layer sl(dpi,initial_masses);
	//cout <<"LINE 303 main initialized layer\n";
	//cout << "LINE 354 particle information:"<< endl;
	//sl.print_particle_info_to_screen();

	// lets initialise the stack. To create a stack object we give it the mean sea level
	// this mean sea level is above a datum...assumed to be the base of the sediment column
	sediment_stack sed_stack(mean_tide);

	// now add some layers to the sediment stack
  cout <<"Adding 100 layers of " << intitial_mass << " kg each." << endl;
	for (int i = 0; i< 100; i++)
  {
	  sed_stack.deposit_surface_sediment(sl);
  }

	//cout << "LINE 365 printing layer 300: " << endl;
	//sed_stack.print_ind_particles_of_layer_to_screen(300);
	//cout << "LINE 367 printing layer 100: " << endl;
	//sed_stack.print_ind_particles_of_layer_to_screen(100);

	// compression
	// now calcualte the thickness and porosity of the layers
  cout << "I now need to compress your intitial column. " << endl;
	sed_stack.calculate_overburden_and_thickness();

  cout << "Your column is initiated!" << endl;

  return sed_stack;
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
		   double max_depth, double min_depth, double max_bmass, double temperatureincrease, int yr, double bfactor)
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
   {
    cout << "LINE 754 mortality is neagitve! " << bm_and_mort[1] << endl;
   }



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
                                                        // in metres
    double Time_fractions = 15000;       // the number of time steps
                                         // in the simulation over one tidal
                                         // cycle
    double dt = Tidal_Period/Time_fractions;
    double start_time;
    double water_depth;
    double dwater_depth_dt;              // derivative of the water depth
    double flow_vel;
    double TKE;                          // turbulent kinetic energy
    double Slope;                        // surface water slope

    double von_karman = 0.4;             // von karman constant
    double rho_w = 1000;                 // density of water
    double TKE_coefficient = 0.2;        // coeffieicent relating TKE to
                                         // shear stress
    double w_up;                         // upward velocity due to turbulence

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
            {
              particle_conc[i] = 0;
            }

            //if (particle_conc[i] > 0)
            //        cout << "LINE 528 ts: " <<t_ime << " part_conc["<<i<<"]: " << particle_conc[i] <<endl;

            if (counter == 0)
            {
              counter ++;
            }

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


#endif
