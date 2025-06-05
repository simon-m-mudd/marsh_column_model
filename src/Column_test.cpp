//unix directory routing
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include "./class_headers/depo_particle.hpp"
#include "./class_headers/depo_particle_info.hpp"
#include "./class_headers/sediment_layer.hpp"
#include "./class_headers/sediment_stack.hpp"
#include "./class_headers/sediment_stack.hpp"
#include "./class_headers/MuddColTrapping.hpp"
using namespace std;

int main (int nNumberofArgs,char *argv[])
{
  // some paramters
  //Test for correct input arguments
  if (nNumberofArgs!=3)
  {
    cout << "=========================================================" << endl;
    cout << "|| Welcome to the Marsh Column Model!  ||" << endl;
    cout << "=========================================================" << endl;
    cout << "This program requires two inputs: " << endl;
    cout << "* First the path to the parameter files." << endl;
    cout << "   The path must have a slash at the end." << endl;
    cout << "  (Either \\ or / depending on your operating system.)" << endl;
    cout << "* Second the prefix of the parameter files." << endl;
    cout << "=========================================================" << endl;
    exit(EXIT_SUCCESS);
  }
  
  string path_name = argv[1];
  string param_name_prefix = argv[2];
  
  cout << "The path name is: " << path_name << endl 
       << " and the paramter prefix is: " << param_name_prefix << endl;
       
  // now initiate a MuddColTrapping object
  MuddColTrapping MCT;
  
  //calculate a settling velocity
  double d_p = 0.0001;    // particle diamter in m
  double w_s = MCT.calculate_w_s(d_p);
  cout << "The settling velocity is: " << w_s << endl;
}
