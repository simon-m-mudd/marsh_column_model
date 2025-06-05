//depo_Particle.h
// header file for the discrete particle object
// this version of the dicrete particle keeps track
// of two integers, the time and the position

#include<iostream>
#include<vector>
using namespace std;


#ifndef Depo_Particle_H
#define Depo_Particle_H

class depo_particle
{
 public:
    depo_particle()					{ create(); }
    depo_particle( int StartType)			{ create(StartType); }
    depo_particle( int StartType, double StartMass)	{ create(StartType, StartMass); }

    // accessor functions
    int    getType() const			{ return Type; }
    double getAge() const			{ return Age; }
    double getOSLage() const			{ return OSLage; }
    double getMass() const			{ return Mass; }

    // this function changes the mass of the particle. Used for growth and death of roots.
    // decay of radioisotopes and carbon should use the decay functions
    void changeMass(double deltaMass);
    void updateMass(double newMass);



    depo_particle(const depo_particle& tP)
    	{ create(tP.getType(),tP.getAge(),tP.getOSLage(), tP.getMass() ); }
    depo_particle(depo_particle& tP)
    	{ create(tP.getType(),tP.getAge(),tP.getOSLage(), tP.getMass() ); }

    depo_particle& operator=(const depo_particle& tP);

    std::ostream& operator<< (std::ostream&);

    void incrementAge(double dt);					// adds to both the OSL and age
    void linear_decay_second_order_accurate(double kdt, double kdt_squared_over_2);
    									// calculate decay with a second order
    									// accurate approximation
    void linear_decay_integration_order_accurate(double k, double dt);	// calculate decay exactly
    void linear_decay_multiplier(double multiplier);			// calculate decay with an externally
    									// calculated decay multiplyer (typically the
    									// 2nd order accurate multiplier)
    void dp_properties_to_screen();					// prints depo particle properties to screen


 private:
    int Type;					// sediment type
    double Age;					// age of particle
    double OSLage;				// OSL age of particle
    double Mass;				// mass of particle


    void create();
    void create(int);
    void create(int, double);
    void create(int, double, double, double);

};

#endif
