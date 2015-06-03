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


#ifndef MuddColTrapping_H
#define MuddColTrapping_H

#include <string>
using namespace std;

///@brief Main analysis object to interface with other LSD objects.
class MuddColTrapping
{
  public:
    /// @brief The default constructor, just sets parameter values
    /// @author SMM
    /// @date 03/06/2015
    MuddColTrapping()    { create(); }

    /// @brief This function calculates settling velocity
    /// @param d_p the particle diameter in units (??? m??)
    /// @author SMM
    /// @date 03/06/2015
    double calculate_w_s(double d_p);
  protected:
  


  
    double alpha;
    double beta;
    double mu;
    double phi;
    double kappa;
    
    /// kinematic viscosity of water a 20C
    double nu;
    
    /// exponent from Nepf (1999)
    double gamma_tr;
    
    /// exponent from Nepf (1999)
    double epsilon;
    
    /// Gravitational acceleration
    double g;
    
    double alpha_T;
    double beta_T;
    
    /// coefficient in eq of alpha_1 as fxn of solid fraction in Tanino and Nepf eq 13. 
    double aa;
    
    /// coefficient in eq of alpha_1 as fxn of solid fraction in Tanino and Nepf eq 13. 
    double bb;
    
    /// turbulence coefficient Nepf (1999)
    double alpha_turb;
    
    /// from tanino and nepf (2008)
    double alpha_zero;


  private:
    void create();


};

#endif
