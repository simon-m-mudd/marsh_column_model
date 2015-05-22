//depo_particle.cpp
// implementation of the dParticle object
#include <fstream>
#include <iostream>
#include <math.h>
#include "depo_particle.hpp"
using namespace std;


/**************************************************************************\
**
**  depo_particle::depo_particle  Constructor for particles.
**			   The intial age of the particle is zero
**			   The location is assigned
\**************************************************************************/
void depo_particle::create()
 {
  //cout << "LINE 18 depo_particle.cpp empty set " << endl;
  Type = 0;
  Age = 0.0;
  OSLage = -9999.0;
  Mass = 0.0;
  //cout << "Type: " << Type << endl;
  //cout << "yo " << endl;
  //cout << "Age: " << Age << endl;
  //cout << "OSL_age: " << OSLage << endl;
  //cout << "Mass: " << Mass << endl;
 }


void depo_particle::create( int StartType )
 {
  Type = StartType;
  Age = 0;
  OSLage = -9999;
  Mass = 0;
 }

void depo_particle::create( int StartType, double StartMass )
 {
  Type = StartType;
  Age = 0;
  OSLage = -9999;
  Mass = StartMass;
 }


void depo_particle::create( int StartType, double StartAge, double StartOSLage, double StartMass )
 {
  Type = StartType;
  Age = StartAge;
  OSLage = StartOSLage;
  Mass = StartMass;
 }

/**************************************************************************\
**
**  depo_particle::operators
** operators for the depo_particle object
**
\**************************************************************************/
depo_particle& depo_particle::operator=(const depo_particle& rhs)
 {
  if (&rhs != this)
   {
    create(rhs.getType(),rhs.getAge(),rhs.getOSLage(), rhs.getMass());
   }
  return *this;
 }

std::ostream& operator<<(std::ostream& os, const depo_particle& tP)
 {
  os <<   tP.getType() << " " << tP.getAge() << " " << tP.getOSLage() << " " << tP.getMass();
  return os;
 }


/**************************************************************************\
**
**  depo_particle::accessor functions:
** depo_particle::incrementAge  increment the particles age
** depo_particle::changeMass(double deltaMass this adds mass (if delataMass is positive)
**						     or removes mass (if deltaMass is negative)
** depo_particle::updateMass(double newMass) updates the mass with newMass
**
\**************************************************************************/
void depo_particle::incrementAge(double dt)
 {
  Age += dt;
  if (OSLage > 0)
   OSLage += dt;
 }

void depo_particle::changeMass(double deltaMass)
 { Mass += deltaMass; }

void depo_particle::updateMass(double newMass)
 { Mass = newMass; }

/**************************************************************************\
**
**  depo_particle::printing tools
**
\**************************************************************************/
void depo_particle::dp_properties_to_screen()
 {
  cout << "Type: " << Type << endl;
  cout << "Age: " << Age << endl;
  cout << "OSL_age: " << OSLage << endl;
  cout << "Mass: " << Mass << endl;
 }


/**************************************************************************\
**
**  depo_particle::decay
** depo_particle::linear_decay_second_order_accurate
** this function calculates the decay over a time dt to second order accuracy.
**
** depo_particle::linear_decay_multiplier
**  this updates mass with an externally clculated multiplier (could be
**  second order of higher accurate). The purpose of this function is
**  to save computational time
**
** depo_particle::linear_decay_integration_order_accurate
**  perfact accuracy, but higher computational costs
\**************************************************************************/
//
// this calucaltion is just the O(1) and O(2) terms from the taylor expantion of the
// exponential function
void depo_particle::linear_decay_second_order_accurate(double kdt, double kdt_squared_over_2)
 {
  // second order
  Mass = Mass*(1-kdt+kdt_squared_over_2);
 }

void depo_particle::linear_decay_multiplier(double multiplier)
 {
  // second order, but mutiplier is calculated outside
  Mass = Mass*multiplier;
 }

void depo_particle::linear_decay_integration_order_accurate(double k, double dt)
 {
  // do this second order
  Mass = Mass*exp(-k*dt);
 }


