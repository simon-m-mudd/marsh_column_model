//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// MuddColTrapping
// The trapping object from the MuddCol model
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for manipulating
//  and analysing raster data, with a particular focus on topography
//
// Developed by:
//  Simon M. Mudd
// Copyright (C) 2015 Simon M. Mudd 2015
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


#ifndef MuddColTrapping_CPP
#define MuddColTrapping_CPP

#include <string>
#include <cmath>
#include "MuddColTrapping.hpp"
using namespace std;

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void MuddColTrapping::create()
{
  alpha = 0.55;
  beta = 0.40;
  mu = 0.00066;
  phi = 0.29;
  kappa = 0.224;
  nu = 10e-6;                
  gamma_tr = 0.718;        
  epsilon = 2.08;
  g = 9.80;

  
  double alpha_T = alpha*kappa*pow(mu,gamma_tr-epsilon)/pow(nu,gamma_tr);
  double beta_T = phi*(gamma_tr-epsilon)+beta;
  double aa = 0.46;             
  double bb = 3.8;                
  double alpha_turb = 0.9;
  double alpha_zero = 11;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function calculates the settling velocity
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double MuddColTrapping::calculate_w_s(double d_p)
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
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


#endif