//sediment_stack.cpp
// implementation of the dParticle object
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include "sediment_stack.hpp"
#include "sediment_layer.hpp"
#include "depo_particle.hpp"
#include "depo_particle_info.hpp"
using namespace std;


/**************************************************************************\
**
**  sediment_stack::sediment_stack:  Constructor for the stack.
**
\**************************************************************************/
void sediment_stack::create()
 {
  cout << "you need to give more information to the sediment stack!\n";
 }

void sediment_stack::create(double MSL)
 {
  basement_elevation = 0;
  mean_sea_level_elevation = MSL;
  cout << "you are starting with an empty stack\n";
 }


/**************************************************************************\
**
**  sediment_stack::deposit_surface sediment
**			   this function deposits sediment on the marsh surface
\**************************************************************************/
void sediment_stack::deposit_surface_sediment(depo_particle_info DPI, vector<double> dpart_masses)
 {
  sediment_layer temp_sl(DPI,dpart_masses);
  layers.push_back(temp_sl);
  total_layer_mass.push_back(temp_sl.get_total_layer_solid_mass());
  layer_top_elevations.push_back(-9999);
  layer_thickness.push_back(-9999);
  layer_porosity.push_back(-9999);
 }

void sediment_stack::deposit_surface_sediment(sediment_layer sl)
 {
  layers.push_back(sl);
  total_layer_mass.push_back(sl.get_total_layer_solid_mass());
  layer_top_elevations.push_back(-9999);
  layer_thickness.push_back(-9999);
  layer_porosity.push_back(-9999);
 }

double sediment_stack::get_layer_thickness(int lay)
{
	return layer_thickness[lay];
}
double sediment_stack::get_layer_top_elevation(int lay)
{
	return layer_top_elevations[lay];
}
/**************************************************************************\
**
**  sediment_stack::create_growth_index:  each year, roots grow and die
** 			root growth is a funciton of the depth below
** 			the surface. This is only calcualted once a year,
**			however, so sediment accumulating after january
**			has no roots until the following year when the roots reset.
**			this index sets a multiplier.
**			Then, for each layer, the root mass is caluclated by
**			multiplying the layer multiplier*dt*the surface growth rate
**			The death rates are also caluclated in the way, exept for
**			the additional chi_organic multiplier which determines
**			the partitioning of dead material into the various carbon pools
**
\**************************************************************************/
vector<double> sediment_stack::create_growth_index(double gamma)
 {
 // first get the depths:
  int n_layers = layers.size();
  double layer_midpoint_depth;
  double top;
  double bottom;
  vector<double> growth_index;

  // the exp multplier is e^2/(e^2-1) and is derived from
  // the relationship between total biomass or mortality and the
  // mortality or biomass at the surface. (see notes_on_growth.nb
  double exp_multiplier = 1.15652;

  double sum = 0;
  for (int i = 0; i< n_layers; i++)
   {
    // get the top elevation
    top = layer_top_elevations[n_layers-1]-layer_top_elevations[i];

    // get the bottom elevation
    bottom = top + layer_thickness[i];

    // if the layer is below 2*gamma, the growth index is zero
    // note there will be a slight error introduced into these calucaltions because the
    // total biomass or growth is assumed to lie between a depth of zero and 2*gamma, but
    // the boundary of the lowest layer will not neccesarily be at 2*gamma
    // there will also be and lower than normal biomass growth if the seidment is less than
    // 2*gamma thick. This may be altered in the future with additional logic.
    if (top > 2*gamma)
     growth_index.push_back(0);
    else
     growth_index.push_back(exp_multiplier*(exp(-top/gamma)-exp(-bottom/gamma)));

    //cout << "layer: " << i << " growth_index: " << growth_index[i] << endl;

    sum+=growth_index[i];

  }

  // check to see how close the index is to unity
  //cout << "LINE 110 sediment_stack.cpp: sum of growth index is: " << sum << endl;

  return growth_index;


 }

vector<double> sediment_stack::create_growth_index_deep(double gamma, double efolding_multiplier,
                                                        double exp_multiplier)
 {
  // the exp multplier is e^n/(e^n-1) where n
  // is teh efolding_multiplier and is derived from
  // the relationship between total biomass or mortality and the
  // mortality or biomass at the surface. (see notes_on_growth.nb

 // first get the depths:
  int n_layers = layers.size();
  double layer_midpoint_depth;
  double top;
  double bottom;
  vector<double> growth_index;


  double sum = 0;
  for (int i = 0; i< n_layers; i++)
   {
    // get the top elevation
    top = layer_top_elevations[n_layers-1]-layer_top_elevations[i];

    // get the bottom elevation
    bottom = top + layer_thickness[i];

    // if the layer is below 2*gamma, the growth index is zero
    // note there will be a slight error introduced into these calucaltions because the
    // total biomass or growth is assumed to lie between a depth of zero and 2*gamma, but
    // the boundary of the lowest layer will not neccesarily be at 3*gamma
    // there will also be and lower than normal biomass growth if the seidment is less than
    // 2*gamma thick. This may be altered in the future with additional logic.
    if (top > efolding_multiplier*gamma)
     growth_index.push_back(0);
    else
     growth_index.push_back(exp_multiplier*(exp(-top/gamma)-exp(-bottom/gamma)));

    //cout << "layer: " << i << " growth_index: " << growth_index[i] << endl;

    sum+=growth_index[i];

  }

  // check to see how close the index is to unity
  //cout << "LINE 110 sediment_stack.cpp: sum of growth index is: " << sum << endl;

  return growth_index;


 }

vector<double> sediment_stack::create_growth_index_inf(double gamma)
 {
  // the exp multplier is e^n/(e^n-1) where n
  // is teh efolding_multiplier and is derived from
  // the relationship between total biomass or mortality and the
  // mortality or biomass at the surface. (see notes_on_growth.nb

 // first get the depths:
  int n_layers = layers.size();
  double layer_midpoint_depth;
  double top;
  double bottom;
  vector<double> growth_index;


  double sum = 0;
  for (int i = 0; i< n_layers; i++)
   {
    // get the top elevation
    top = layer_top_elevations[n_layers-1]-layer_top_elevations[i];

    // get the bottom elevation
    bottom = top + layer_thickness[i];

    // if the layer is below 2*gamma, the growth index is zero
    // note there will be a slight error introduced into these calucaltions because the
    // total biomass or growth is assumed to lie between a depth of zero and 2*gamma, but
    // the boundary of the lowest layer will not neccesarily be at 3*gamma
    // there will also be and lower than normal biomass growth if the seidment is less than
    // 2*gamma thick. This may be altered in the future with additional logic.
    growth_index.push_back((exp(-top/gamma)-exp(-bottom/gamma)));

    //cout << "layer: " << i << " growth_index: " << growth_index[i] << endl;

    sum+=growth_index[i];

  }

  // check to see how close the index is to unity
  //cout << "LINE 110 sediment_stack.cpp: sum of growth index is: " << sum << endl;

  return growth_index;


 }

/**************************************************************************\
**
**  sediment_stack::root_growth_and_death
**			   this calcualtes both the growth and the death of
*			   belowground biomass
\**************************************************************************/
void sediment_stack::root_growth_and_death(vector<double> growth_index, int root_type,
                        vector<int> carbon_types, vector<double> chi_org,
                        double BG_Bmass, double tot_mort_mass)
 {
  //cout << "LINE 128 sediemnt_stack.cpp, printing particles at layer 10\n";
  //print_ind_particles_of_layer_to_screen(10);

  int n_carbon_types = carbon_types.size();
  double layer_root_mass;
  vector<double> layer_carbon_ingrowth(n_carbon_types);
  									// change in carbon mass
  //cout <<"LINE 135 sed_stack.cpp root type: " << root_type << endl;
  //for (int i = 0; i< n_carbon_types ; i++)
  // cout << "ctype["<<i<<"]: " << carbon_types[i] << endl;

  // the root mass depends on the total aboveground biomass and theta_bg, which is the ratio
  // between belowground and aboveground biomass
  // first we see the number of nodes in the growth index:
  int n_g_layers = growth_index.size();
  //cout << "LINE 137 sed_stack.cpp n_g_layers: "<<n_g_layers << endl;
  //cout << BG_Bmass*growth_index[150] << endl;
  for (int i = 0; i< n_g_layers; i++)
   {
    layer_root_mass = BG_Bmass*growth_index[i];
    //cout << "Line 147 sed_stack.cpp layer["<<i<<"] root mass: " << layer_root_mass << endl;

    //cout << "LINE 149 sediemnt_stack.cpp, printing particles at layer " << i << endl;
    //print_ind_particles_of_layer_to_screen(i);

    for (int j = 0; j<n_carbon_types; j++)
     {
	  layer_carbon_ingrowth[j] = tot_mort_mass*chi_org[j]*growth_index[i];
	  //cout << " LINE 156 sed_stack.cpp: c" << j <<": "<< layer_carbon_ingrowth[j];
     }
    //cout << endl;

    //depo_particle dp;
    //cout << "LINE 153 sed_stack.cpp initiating depo particle: " << endl;
    //dp.dp_properties_to_screen();
    //cout << endl;

    //cout << "LINE 149 sed_stack.cpp editing layers, there are " << layers.size() << " layers " << endl;
    layers[i].root_growth_and_death(root_type, carbon_types,
                                    layer_root_mass, layer_carbon_ingrowth);
    //cout << "LINE 152 sed_stack.cpp done with layer " << endl;
   }

 }


/**************************************************************************\
**
**  sediment_stack::get_marsh_surface_elevation
**
\**************************************************************************/
double sediment_stack::get_marsh_surface_elevation()
 {
  int n_layers = layers.size();
  return layer_top_elevations[n_layers-1];
 }

/**************************************************************************\
**
**  sediment_stack::get_n_layers
**
\**************************************************************************/
int sediment_stack::get_n_layers()
 {
  int n_layers = layers.size();
  return n_layers;
 }

/**************************************************************************\
**
**  sediment_stack::get_total_stack_mass_of_each_type()
**
\**************************************************************************/
vector<double> sediment_stack::get_total_stack_mass_of_each_type(int n_p_types)
 {

  vector<double> stack_masses(n_p_types,0.0);
  int n_layers = layers.size();

  //cout << "LINE 198 sed_stack.cpp n_layers: " << n_layers
  //     << " np_types: " << n_p_types << endl;
  for (int i = 0; i< n_layers; i++)
   {
	vector<double> part_masses = layers[i].get_mass_of_particle_types();
    for (int j = 0; j< n_p_types; j++)
	 stack_masses[j]+= part_masses[j];
   }

  return stack_masses;
 }




/**************************************************************************
**
**  sediment_stack::get_stack_slice_mass_of_each_type()
**
\**************************************************************************/
vector<double> sediment_stack::get_stack_slice_mass_of_each_type(int n_p_types,
						double slice_top_depth, double slice_bottom_depth)
 {

  vector<double> stack_masses(n_p_types,0.0);
  int n_layers = layers.size();

  double top_elev = layer_top_elevations[n_layers-1];
  double layer_top_depth;
  double layer_bottom_depth;

  int counting_start_trigger = 0;
  int counting_end_trigger = 1;

  double true_slice_bottom;
  double true_slice_top;

  //cout << "LINE 198 sed_stack.cpp n_layers: " << n_layers
  //     << " np_types: " << n_p_types << endl;
  for (int i = 0; i< n_layers; i++)
   {
	// get the top and bottom depths of the layer
	layer_top_depth = top_elev-layer_top_elevations[i];
	if (i == 0)
	 layer_bottom_depth = top_elev;
	else
	 layer_bottom_depth = top_elev-layer_top_elevations[i-1];

	//cout << "Layer_top_depth: " << layer_top_depth
	//     << " and bottom_depth: " << layer_bottom_depth << endl;



	// when counting start trigger == 0, the bottom of the slice
	// has yet to be triggered.
	if (counting_start_trigger == 0)
	{
		// if
		if (slice_bottom_depth <= layer_bottom_depth &&
		    slice_bottom_depth >= layer_top_depth)
		{
			counting_start_trigger = 1;
			true_slice_bottom = layer_bottom_depth;
		}
	}
	if (counting_end_trigger == 1)
	{
		if (slice_top_depth <= layer_bottom_depth &&
		    slice_top_depth >= layer_top_depth)
		{
			counting_end_trigger = 0;
			true_slice_top = layer_bottom_depth;
		}
	}



    if (counting_start_trigger == 1 && counting_end_trigger == 1)
    {
		//cout << "Adding it up!" << endl;
		vector<double> part_masses = layers[i].get_mass_of_particle_types();
    	for (int j = 0; j< n_p_types; j++)
	    {
			stack_masses[j]+= part_masses[j];
		}
	}
   }

  //cout << "slice_top: " << slice_top_depth << endl
  //     << "true_slice_top: " << true_slice_top << endl
  //     << "slice_bottom: " << slice_bottom_depth << endl
  //     << "true_slice_bottom: " << true_slice_bottom << endl
  //     << "slice_thick: " << slice_bottom_depth-slice_top_depth << endl
  //     << "true_thick: " << true_slice_bottom-true_slice_top << endl;

  //for (int i = 0; i<n_p_types; i++)
  // cout << "type: " << i << " mass: " << stack_masses[i] << endl;

  return stack_masses;
 }

 vector<double> sediment_stack::get_stack_slice_mass_of_each_type(int n_p_types,
 										int bottom_layer, int top_layer)
  {
   int n_layers = layers.size();
   vector<double> particle_masses;
   particle_masses = layers[0].get_mass_of_particle_types();
   int np = particle_masses.size();
   vector<double> band_particle_masses(np,0.0);


   //cout << "bot lay: " << bottom_layer << " and top lay: " << top_layer << endl;
   if (bottom_layer > 0 && top_layer >0)
   {
	   if (bottom_layer < top_layer && top_layer < n_layers)
	   {
		   for (int lay = bottom_layer; lay<=top_layer; lay++)
		   {
			   // get the particle masses for the layer
               particle_masses = layers[lay].get_mass_of_particle_types();

               for (int ptype = 0; ptype<np; ptype++)
			   {
			    	band_particle_masses[ptype]+= particle_masses[ptype];
			   }
		   }
	   }
   }
   return band_particle_masses;
 }

/**************************************************************************\
**
**  sediment_stack::organic_decay:  this decays the non-refractory carbon
** 					in the sediment column
**
\**************************************************************************/
void sediment_stack::decay_and_mass_loss(double dt, double temperatureincrease, double kfactor) // MK- I added temperatureincrease and kfactor
 {
  // first get the depths:
  int n_layers = layers.size();
  double layer_midpoint_depth;

  //cout << "\n\n starting layers\n";
  for (int i = 0; i< n_layers; i++)
   {

    layer_midpoint_depth = layer_top_elevations[n_layers-1]-
                           layer_top_elevations[i]+layer_thickness[i]*0.5;

    layers[i].decay_and_mass_loss(dt, layer_midpoint_depth, temperatureincrease, kfactor); // MK- I added temperatureincrease and kfactor
    total_layer_mass[i] = layers[i].get_total_layer_solid_mass();
  }
 }


/**************************************************************************\
**
**  sediment_stack::update_thickness:  recalcualte
**			the layer thicknesses based on new overburden
**
\**************************************************************************/
void sediment_stack::calculate_overburden_and_thickness()
 {
  int n_layers = layers.size();
  vector<double> overburden(n_layers);
  overburden[n_layers-1] = 0;			// the overburden at the surface is zero
  double bw = 0;

  // prgressivly add overburden for each layer
  if (n_layers > 1)
   {
    for (int i = n_layers-2; i>=0; i--)
     {
      bw += layers[i+1].calculate_buoyant_weight();
      overburden[i] = bw;
     }
   }

  vector<double> compaction_vec(2);		// this is the vector returned
  						// by the compaction function in sediment_layer
  						// element [0] is porosity and element [1] is layer thickness
  						// (this assumes the column is 1m^2)

  // calculate thicknesses and porosities
  double lte = 0;
  for (int i = 0; i<n_layers; i++)
   {
    compaction_vec = layers[i].instant_compaction_compression_index(overburden[i]);
    layer_porosity[i] = compaction_vec[0];
    layer_thickness[i] = compaction_vec[1];
    lte+=layer_thickness[i];
    layer_top_elevations[i] = lte;

    //cout << "layer number: " << i << " thickness: " << layer_thickness[i] << endl
    //     << " porosity: " << layer_porosity[i] << " top elevation: " << layer_top_elevations[i] << endl;
   }

 }


/**************************************************************************\
**
**  sediment_stack::print...:  printing functions
**
\**************************************************************************/
void sediment_stack::print_particle_masses_to_screen(double t_ime)
 {
  int n_layers = layers.size();
  cout << endl << "time is: " <<t_ime << endl;
  vector<double> masses;
  for (int i = 0; i<n_layers; i++)
   {
    cout << "layer " << i << " ";
    masses = layers[i].get_mass_of_particle_types();
    int nparts = masses.size();
    for (int j = 0; j<nparts; j++)
     cout << " " << masses[j];
    cout << endl;
   }
  cout << endl;
 }


void sediment_stack::print_layer_thickness_to_screen(double t_ime)
 {
  int n_layers = layers.size();
  cout << t_ime;
  for (int i = 0; i<n_layers; i++)
   cout << " " << layer_thickness[i];
  cout << endl;
 }

void sediment_stack::print_layer_mass_to_screen(double t_ime)
 {
  int n_layers = layers.size();
  cout << t_ime;
  for (int i = 0; i<n_layers; i++)
   cout << " " << total_layer_mass[i];
  cout << endl;
 }

void sediment_stack::print_layer_porosity_to_screen(double t_ime)
 {
  int n_layers = layers.size();
  cout << t_ime;
  for (int i = 0; i<n_layers; i++)
   cout << " " << layer_porosity[i];
  cout << endl;
 }

void sediment_stack::print_layer_top_elevations_to_screen(double t_ime)
 {
  int n_layers = layers.size();
  cout << t_ime;
  for (int i = 0; i<n_layers; i++)
   cout << " " << layer_top_elevations[i];
  cout << endl;
 }

void sediment_stack::print_ind_particles_of_layer_to_screen(int layer_number)
 {
  int sz = layers.size();
  if (layer_number <0 || layer_number > sz-1)
   cout << "LINE 345 sedment_stack.cpp you have called a layer that doesn't exist\n";
  else
   layers[layer_number].print_particle_info_to_screen();
 }

/**************************************************************************\
**
**  sediment_stack::print...:  printing functions
**
\**************************************************************************/
void sediment_stack::print_layer_thickness(double t_ime, ofstream& fout)
 {
  int n_layers = layers.size();
  fout << t_ime << " " << n_layers;
  for (int i = 0; i<n_layers; i++)
   fout << " " << layer_thickness[i];
  fout << endl;
 }

void sediment_stack::print_layer_mass(double t_ime, ofstream& fout)
 {
  int n_layers = layers.size();
  fout << t_ime << " " << n_layers;
  for (int i = 0; i<n_layers; i++)
   fout << " " << total_layer_mass[i];
  fout << endl;
 }

void sediment_stack::print_layer_porosity(double t_ime, ofstream& fout)
 {
  int n_layers = layers.size();
  fout << t_ime << " " << n_layers;
  for (int i = 0; i<n_layers; i++)
   fout << " " << layer_porosity[i];
  fout << endl;
 }

void sediment_stack::print_layer_top_elevations(double t_ime, ofstream& fout)
 {
  int n_layers = layers.size();
  fout << t_ime << " " << n_layers;
  for (int i = 0; i<n_layers; i++)
   fout << " " << layer_top_elevations[i];
  fout << endl;
 }

void sediment_stack::print_layer_particle_proportions(double t_ime, ofstream& fout)
 {
  int n_layers = layers.size();
  vector<double> particle_masses;
  particle_masses = layers[0].get_mass_of_particle_types();
  int np = particle_masses.size();

  vector< vector<double> > prop_vecs(np);



  for (int i = 0; i<n_layers; i++)
   {
    // get the particle masses for the layer
    particle_masses = layers[i].get_mass_of_particle_types();

    // put these masses in the prop_vecs vectors
    for (int ptype = 0; ptype<np; ptype++)
     prop_vecs[ptype].push_back(particle_masses[ptype]);
   }

  vector<double> p_masses;
  // now print a line for each particle type
  for (int ptype = 0; ptype<np; ptype++)
   {
    p_masses = prop_vecs[ptype];
    fout << t_ime << " " << n_layers;
    for (int i = 0; i<n_layers; i++)
     fout << " " << p_masses[i]/total_layer_mass[i];
    fout << endl;
   }

 }

void sediment_stack::print_layer_particle_masses(double t_ime, ofstream& fout)
 {
  int n_layers = layers.size();
  vector<double> particle_masses;
  particle_masses = layers[0].get_mass_of_particle_types();
  int np = particle_masses.size();

  vector< vector<double> > prop_vecs(np);

  for (int i = 0; i<n_layers; i++)
   {
    // get the particle masses for the layer
    particle_masses = layers[i].get_mass_of_particle_types();

    // put these masses in the prop_vecs vectors
    for (int ptype = 0; ptype<np; ptype++)
     prop_vecs[ptype].push_back(particle_masses[ptype]);
   }

  vector<double> p_masses;
  // now print a line for each particle type
  for (int ptype = 0; ptype<np; ptype++)
   {
    p_masses = prop_vecs[ptype];
    fout << t_ime << " " << n_layers;
    for (int i = 0; i<n_layers; i++)
     fout << " " << p_masses[i];
    fout << endl;
   }

 }


void sediment_stack::print_single_layer_particle_masses(double t_ime, ofstream& fout, int layer_number)
 {
  int n_layers = layers.size();
  vector<double> particle_masses;

  int np = particle_masses.size();


  if (layer_number < n_layers)
  {
	  particle_masses = layers[layer_number].get_mass_of_particle_types();
	  fout << t_ime;
	  for (int ptype = 0; ptype<np; ptype++)
      {
		  fout << " " << particle_masses[ptype];
	  }

      fout << endl;
  }


 }

 void sediment_stack::print_band_layer_particle_masses(double t_ime, ofstream& fout,
 										int bottom_layer, int top_layer)
  {
   int n_layers = layers.size();
   vector<double> particle_masses;
   particle_masses = layers[0].get_mass_of_particle_types();
   int np = particle_masses.size();


   //cout << "bot lay: " << bottom_layer << " and top lay: " << top_layer << endl;
   if (bottom_layer > 0 && top_layer >0)
   {
	   if (bottom_layer < top_layer && top_layer < n_layers)
	   {
		   for (int lay = bottom_layer; lay<=top_layer; lay++)
		   {
			   // get the particle masses for the layer
               particle_masses = layers[lay].get_mass_of_particle_types();

               fout << t_ime << " " << lay;
               for (int ptype = 0; ptype<np; ptype++)
			   {
			    	fout << " " << particle_masses[ptype];
			   }
			   fout << endl;
		   }
	   }
   }
 }

void sediment_stack::print_stack_stats(double t_ime, ofstream& fout)

{
	  int n_layers = layers.size();
	  vector<double> particle_masses;
	  particle_masses = layers[0].get_mass_of_particle_types();
	  int np = particle_masses.size();
	  int yr = 1;

	  for (int i = 0; i<n_layers; i++)
	  {
		  if (i>100)
		  {
			//cout << "i-100: " << i-100 << " (i-100)%24: " << (i-100)%24 << endl;
		  	if ( (i-100)%24 == 0)
		  	{
				  yr++;
		 	 }
		  }

		  particle_masses = layers[i].get_mass_of_particle_types();
		  fout << t_ime << " " << yr << " " << i << " " << layer_top_elevations[i] << " "
		       << layer_porosity[i] << " " << total_layer_mass[i];
	      for (int ptype = 0; ptype<np; ptype++)
		  {
		  		fout << " " << particle_masses[ptype];
		  }
		  fout << endl;
	  }
 }
