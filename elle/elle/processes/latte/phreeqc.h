#ifndef _E_phreeqc_h
#define _E_phreeqc_h

#define dimX 100
#define dimY 100

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
#include <vector>
// header of the IPhreeqc class (prefer system, fallback to local stub)
#if __has_include(<IPhreeqc.hpp>)
#include <IPhreeqc.hpp>
#else
#include "iphreeqc_stub/IPhreeqc.hpp"
#endif
// include lattice particle class
#include "particle.h"

// OpenMP is optional; include only if available
#if __has_include(<omp.h>)
#include <omp.h>
#endif

using namespace std;
using std::vector;

class PhreeqC {
  /***********************************************************************
   * Now we build the class which can access the Iphreec functions 	   *
   **********************************************************************/
public:
  PhreeqC();
  IPhreeqc iphreeqc;
  ~PhreeqC() {};

  PhreeqC(const PhreeqC &rhs) = default;
  PhreeqC &operator=(const PhreeqC &rhs) = default;

  // for convinience we define two vector types.
  typedef vector<std::pair<const char *, double>>
      FLUID; // vector to hold concentrations of a fluid
  typedef vector<std::pair<const char *, boost::array<double, (dimY * dimX)>>>
      FLUID_array; // this is the container that holds all the matrices
  FLUID Conc_inf;  // the infiltrating fluid
  FLUID Conc_por, Conc_por_b; // the pore fluid(s)
  FLUID_array Conc_matrices;  // the pore fluid matrices for the transport
  string units; // this holds the units (haas to be consitnt in the input files)
  const char *p_database; // this guy holds the name of the databse to use
  const char *pore_fluid; // this is the name of the file containing the pore
                          // fluid comosition
  const char *inf_fluid;  // this is the name of the file containing the
                          // infiltrating fluid comosition

  int time_step;
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
  double TermalDiffusivity(double thermal_conductivity, double heat_capacity,
                           double density);
  double Diffusivity(double D, double Ea, double T);

  // these are some helper functions to read values
  double FindConcentration(const char *type, FLUID fluid);
  FLUID GetFluid(int x, int y, FLUID_array fluid_c);
  void GetConcentration(const char *type, double (&Conc)[dimX][dimY]);
  void SetConcentration(const char *type, double Conc[dimX][dimY]);
  void SetFluidMatrices(vector<std::pair<const char *, double>> &fluid,
                        int boundary);

  void AddPoreFluidLayer(vector<std::pair<const char *, double>> fluid,
                         int y_min, int y_max);

  // test function to read phreeqc outputs
  void Selected(VAR v);
  // Initialisation--------------------------------------------------------
  void InitializeFluidVector(); // this generates the fluid vector for the
                                // simulations and the list of concentrations in
                                // the domain
  void AdjustUnits(vector<std::pair<const char *, double>> &fluid);
  void ChargeBalance(vector<std::pair<const char *, double>> &fluid);
  void UpdateVector(VAR v, vector<std::pair<const char *, double>> &vec);
  // Main------------------------------------------------------------------
  void Equilibrate(vector<std::pair<const char *, double>> &fluid,
                   string phase);
  void OUTPUT(int time_step, vector<std::pair<const char *, double>> c_fluid,
              vector<std::pair<std::string, double>> &local_result,
              bool write_file);
  double SI(int time_step, FLUID fluid, string phase, bool write_file);
  void React(FLUID fluid, double k, double n, double A0, double V, double m0,
             double time_step, double tot_time, string mineral, string formula,
             bool write_file);
  double dissCalcite(vector<double> Fluid, double radius, double porosity,
                     int NB, double &AdCalcite);
};
#endif
