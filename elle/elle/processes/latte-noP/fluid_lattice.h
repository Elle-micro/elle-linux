#ifndef _Fluid_Lattice_H
#define _Fluid_Lattice_H

#define PI 4 * atan(1)

#include <cstdlib>
#include <math.h>
#include <time.h>

#include "particle.h"

class Fluid_Lattice {

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

  double alpha[200][200];
  double beta[200][200]; // aplha (calculated) variable for each pressure node
                         // of pressure lattice

  double a_inv[200][200]; // inverse matrix for comparing with a[dim][dim] &
                          // b[dim][dim] matrices

  double al[200][200][200];
  double bl[200][200][200];
  double a_invl[200][200][200];

  double interim_Pf[200][200]; // corrosponding intermediate pressure value of
                               // each pressure mode during multiplication of
                               // two pressure matrices

  double Pf[200][200];
  double oldPf[200][200];

  double Mineral_Diffusion[200][200];
  double Fluid_Boundary[200][200];
  int Fluid[200][200];
  int Is_Reactive_Grain[200][200];

  double Mg_conc[200][200];
  double Fe_conc[200][200];
  double Mg_conc_b[200][200];
  double Fe_conc_b[200][200];

  double HeatSheet[200][100];
  double IntHeatSheet[200][100];

  double Con[200][200];
  double Metal_A[200][200];
  double Temperature[200][200];

  double density_factor[200][200];

  double ddyc[200][200];
  double ddxc[200][200];
  double dyc[200][200];
  double dxc[200][200];
  double vel_darcyx[200][200];
  double vel_darcyy[200][200];
  double vel_darcyx_av[200][200];
  double vel_darcyy_av[200][200];

  double solid_fraction[200][200];
  double kappa[200][200];
  double fracture[200][200];

  double par_area;
  double par_rad;

  double smooth_func; // tent function or smooth function

  double vel_x[200][200]; // velocity field (average velocity) along x-direction
                          // associated with each pressure node
  double vel_y[200][200]; // velocity field along y-direction

  double source_a[200][200];
  double source_b[200][200];
  double source_c[200][200];
  double source[200][200];

  double time_a; // time for calculation of alpha & beta variables of pressure
                 // lattice

  double scale_constt; // total dimension in meter of real system to model down

  //*******************************************

  Fluid_Lattice(
      int lattice_size); // constructor, function with pointer argument, points
                         // to particle class pointer, points to particle class
  Particle *help; // help pointer, points to particle class

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

  // exchange functions new Daniel and Simon 2021
  void Init_Mineral(double constant);
  void Set_Mineral_Fluid(int gr, int min, Particle **list);
  void Pass_Back_Concentration(Particle **list);
  void SolveDiffusionImplicitSolid(double time, double space);
  void alpha_set_solid(double time_dif);
  void alpha_set_dif_por(double time_dif);
  void beta_set_dif_por(double time_dif);
  void matrices_set_1_fluid();
  void matrices_set_2_fluid();
  void multiplication_solid_Mg();
  void multiplication_inv_solid_Mg();
  void multiplication_solid_Fe();
  void multiplication_inv_solid_Fe();
  void SolveDiffusionImplicitFluid(double time, double space);
  void multiplication_fluid_Mg();
  void multiplication_inv_fluid_Mg();
  void multiplication_fluid_Fe();
  void multiplication_inv_fluid_Fe();
  void SolveDiffusionExplicitFluid(double time, double space);
  void Exchange(double time, double space, Particle **list);
  void Set_Boundary_Fluid(double background, double bound, Particle **list);

  void SetHeatSheet(double pos, double pos2);
  void SolveHeatSheet(double space, double time, double dif);
  void Pass_Back_SheetT(Particle **list, double temp, double critAge);
  void Read_TSheet(Particle **list, int boundary);
  void SetHeatSheetbox(double x, double x2, double y, double y2, int Temp);
  void InitHeatSheet(int Temp);
};
#endif
