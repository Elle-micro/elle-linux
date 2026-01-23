#ifndef _Fluid_Lattice_H
#define _Fluid_Lattice_H

// change the resolution of the elle file her:
#define dimX 100 // num particles along x
#define dimY 100 // num particles along y
#define dimM 100 // matrix stacker

#define PI 4 * atan(1)

#include <cstdlib>
#include <math.h>
#include <time.h>

#include "particle.h"
#include "phreeqc.h"          // header for the phreeqc functions
#include <boost/progress.hpp> //a boost libar that shows progress in terminal

class Fluid_Lattice {

  friend class PhreeqC;

public:
  int dim;

  int pl_size;
  double rho, rho_s, width, delta_x, delta_y;

  int box; // counters & box: the corresponding four repulsion boxes for each
           // pressure node

  double x_diff, y_diff; // difference b/w x & y positions of particles and the
                         // pressure nodes resp.

  /******************************************************
   * Fluid lattice parameters
   ******************************************************/

  double area_node; // area of pressure node
  double
      rad; // (smooth function) radius of circle centred at respective pressure
           // node and edged at the centres of four surrounding nodes.

  double rho_fluid;
  double mu; // visocsity for calculation of alpha & beta variables of pressure
             // lattice
  double compressibility; // compressibility for calculation of alpha & beta
                          // variables of pressure lattice
  double compressibility_c;
  double compressibility_T;

  double alpha[dimX][dimY]; // aplha (calculated) variable for each pressure
                            // node of pressure lattice

  double a_inv[dimX][dimY]; // inverse matrix for comparing with a[dim][dim] &
                            // b[dim][dim] matrices

  double al[dimX][dimY][dimM];
  double bl[dimX][dimY][dimM];
  double a_invl[dimX][dimY][dimM];

  double interim_Pf[dimX][dimY]; // corrosponding intermediate pressure value of
                                 // each pressure mode during multiplication of
                                 // two pressure matrices

  double Pf[dimX][dimY];
  double oldPf[dimX][dimY];

  double Con[dimX][dimY];
  double conc[dimX][dimY];
  double Metal_A[dimX][dimY];
  double Temperature[dimX][dimY];
  double Con_S6[dimX][dimY];
  double Con_Zn[dimX][dimY];
  double Con_Ca[dimX][dimY];
  double Con_Mg[dimX][dimY];
  double Si[dimX][dimY];

  double density_factor[dimX][dimY];

  double ddyc[dimX][dimY];
  double ddxc[dimX][dimY];
  double dyc[dimX][dimY];
  double dxc[dimX][dimY];
  double vel_darcyx[dimX][dimY];
  double vel_darcyy[dimX][dimY];
  double vel_darcyx_av[dimX][dimY];
  double vel_darcyy_av[dimX][dimY];

  double solid_fraction[dimX][dimY];
  double kappa[dimX][dimY];
  double fracture[dimX][dimY];

  double react[dimX][dimY]; // for PhreqC

  double par_area;
  double par_rad;

  double smooth_func; // tent function or smooth function

  double vel_x[dimX][dimY]; // velocity field (average velocity) along
                            // x-direction associated with each pressure node
  double vel_y[dimX][dimY]; // velocity field along y-direction

  double source_a[dimX][dimY];
  double source_b[dimX][dimY];
  double source_c[dimX][dimY];
  double source[dimX][dimY];

  double time_a; // time for calculation of alpha & beta variables of pressure
                 // lattice

  int scale_constt; // total dimension in meter of real system to model down

  //*******************************************

  Fluid_Lattice(
      int lattice_size); // constructor, function with pointer argument, points
                         // to particle class pointer, points to particle class
  Particle *help; // help pointer, points to particle class

  PhreeqC iphreeqc_obj;  // object of pgreeqc
  PhreeqC *iphreeqc_ptr; // pointer to phreeqc

  ~Fluid_Lattice(); // destructor

  int alpha_set(int adjust); // function for the construction of alpha and beta
                             // variables of pressure lattice

  double alpha_set_time();
  void alpha_set_fluid(double nonlin_time);

  void matrices_set_1(); // functions for the construction of a[dim][dim] &
                         // b[dim][dim] matrices of pressure lattice
  void matrices_set_2();

  // function to take transpose of pressure matrix during ADI

  void transpose_trans_P();

  void invert(); // function to take inversion of a[dim][dim] pressure matrix

  void copy_matrix();
  void multiplication();
  void multiplication_inv();
  void test_pressure_with_source(int byo);
  void Pass_Back_Gradients(Particle **list, int av_conc);
  void Background_Pressure(double pressure, double concentration);
  void Set_Scale(double scale, double time_factor);

  void Set_Concentration(int x, int y, double concentration);
  void Set_hydrostatic_gradient(int depth);
  void Read_Porosity_Movement_Matrix(double fracture_effect, Particle **list,
                                     int boundary);
  void Calculate_Permeability_Matrix(double Coseny_grain_size);

  void Fluid_Input(double pressure, int x, int y);

  void Background_Concentration(double concentration);
  void Get_Fluid_Velocity();
  void SolveAdvectionLax_Wenddroff(double time, double space, int boundary);
  void SolveAdvectionLax_Wenddroff_Phreeqc(double time, double space,
                                           int boundary);
  void SolveAdvectionLax_Wenddroff_Metal(double time, double space,
                                         int boundary);
  void SolveAdvectionLax_Wenddroff_Temperature(double time, double space,
                                               int boundary);
  void SolveDiffusionImplicit(double dif_const, double time, double space);
  void SolveAdvectionFromm(double space, double time);
  void alpha_set_dif(double dif_constant, double time_dif);
  void multiplication_dif();
  void transpose_trans_dif();
  void multiplication_inv_dif();
  void SolveDiffusionImplicit_Phreeqc(double dif_const, double time,
                                      double space);
  void alpha_set_Phreeqc(double dif_constant, double time_dif);
  void multiplication_Phreeqc();
  void transpose_trans_Phreeqc();
  void multiplication_inv_Phreeqc();
  void SolveMetalDiffusionImplicit(double dif_const, double time, double space);
  void alpha_set_met(double dif_constant, double time_dif);
  void multiplication_met();
  void transpose_trans_met();
  void multiplication_inv_met();
  void SolveTemperatureDiffusionImplicit(double dif_const, double time,
                                         double space);
  void alpha_set_temp(double dif_constant, double time_dif);
  void multiplication_temp();
  void transpose_trans_temp();
  void multiplication_inv_temp();

  void SolveadvDisExplicitTransient(double fluid_v, double dif_const);
  void ChangeConcentration(int box_x, int box_y, double change);

  void Set_Boundaries(double UpperBoundaryPressure,
                      double LowerBoundaryPressure, int double_sides,
                      double conc, double additional);

  /**************************************
   * Phreeqc Hydrogeochemical simulations
   **************************************/
  void SetFluidPhreeqc(const char *pore_fluid, const char *inf_fluid,
                       const char *d_base);
  void TestPhreeqc(double timestep, double space);
};
#endif
