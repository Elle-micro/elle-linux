#ifndef _E_phreeqc_h
#define _E_phreeqc_h

#define dimX 50
#define dimY 50

#include <boost/algorithm/string.hpp>
#include <boost/array.hpp>
#include <boost/lexical_cast.hpp>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
// header of the Iphreeqc class
#include <IPhreeqc.hpp>
// include lattice particle class
#include "particle.h"

#include <omp.h>

using namespace std;
using std::vector;

class PhreeqC {
  /***********************************************************************
   * Now we build the class which can access the Iphreec functions 	   *
   **********************************************************************/
public:
  PhreeqC();
  IPhreeqc iphreeqc;
  IPhreeqc *iphreeqc_ptr;
  ~PhreeqC(){};

  PhreeqC(const PhreeqC &rhs) {}
  PhreeqC &operator=(const PhreeqC &rhs) {};

  // for convinience we define two vector types.
  typedef vector<std::pair<const char *, double>>
      FLUID; // vector to hold concentrations of a fluid
  typedef vector<std::pair<const char *, boost::array<double, (dimY * dimX)>>>
      FLUID_array; // this is the container that holds all the matrices
  FLUID Conc_inf;  // the infiltrating fluid
  FLUID Conc_por;  // the pore fluid
  FLUID_array Conc_matrices; // the pore fluid matrices for teh transport
  string units; // this holds the units (haas to be consitnt in the input files)
  const char *p_database; // this guy holds the name of the databse to use
  const char *pore_fluid; // this is the name of the file containing the pore
                          // fluid comosition
  const char *inf_fluid; // this is the name of the file containing the
                         // infiltrating fluid comosition

  int nb_threads;
  vector<double> Fluid;
  vector<double>::iterator F;
  vector<VAR> fluid;

  void Test();

  // basic functions-------------------------------------------------------
  void Version();
  void ERROR();

  int loader(const char *);          // load the databse
  void SetPoreFluid(const char *f1); // check and set the file of teh pore fluid
  void SetInfFluid(
      const char *f2); // check and set the file of the infiltrating fluid
  void SetDatabase(const char *database); // check and set the database
  double ReturnValue(
      VAR vs,
      string name); // return the value of 'name' from the phreeqc output vector

  // these are some helper functions to read values
  double FindConcentration(const char *type);
  FLUID GetFluid(int x, int y);
  void GetConcentration(const char *type, double (&Conc)[dimX][dimY]);
  void SetConcentration(const char *type, double Conc[dimX][dimY]);

  // test function to read phreeqc outputs
  void Selected(VAR v);
  // Initialisation--------------------------------------------------------
  void InitializeFluidVector(); // this generates the fluid vector for the
                                // simulations and the list of concentrations in
                                // the domain
  void AdjustUnits(vector<std::pair<const char *, double>> &fluid);
  void UpdateVector(VAR v, vector<std::pair<const char *, double>> &vec);
  // Main------------------------------------------------------------------
  void Equilibrate(vector<std::pair<const char *, double>> &fluid,
                   string phase);
  double SI(FLUID fluid, string phase, bool write_file);
  void React(FLUID fluid, double k, double n, double A0, double V, double m0,
             double time_step, double tot_time, string mineral, string formula,
             bool write_file);
  double dissCalcite(vector<double> Fluid, double radius, double porosity,
                     int NB, double &AdCalcite);
};
#endif
