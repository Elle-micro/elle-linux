#ifndef _Pressure_Lattice_H
#define _Pressure_Lattice_H

#define PI 4 * atan(1)

#include <cstdlib>
#include <math.h>
#include <time.h>

#include "particle.h"

class Pressure_Lattice {

private:
  const static int dim =
      103; // 53; //103;			// dimension of pressure lattice
           //  any function with extension _a will function with dimension 50*50
           //  functions with extension _b will work with 52*52 dimensions
           //  functions with extension _bb will work with 102*102 dimensions
           //  functions with extension _c will work with 51*51 dimensions
           //  functions with extension _cc will work with 101*101 dimensions
           //  functions with extension _d will work with 53*53 dimensions
           //  functions with extension _dd will work with 103*103 dimensions

  // static const int Low = 10000, High = 50000;     // Lower & higher bounds
  // for random increment of specific Pressure nodes in Pascals
  int i, j, k, box, v; // counters & box: the corresponding four repulsion boxes
                       // for each pressure node

  float x_wall,
      y_wall; // distance b/w side walls & upper-lower walls of pressure lattice
  // float delta_x, delta_y;			// pressure lattice constants in
  // x & y directions
  float x_diff, y_diff; // difference b/w x & y positions of particles and the
                        // pressure nodes resp.

  /******************************************************
   * Particle parameters
   * *****************************************************/
  float rad_particle; // particle radius for calculation of kappa in pressure
                      // lattice
  float rad_particle_b;
  // float volume_particle;

  float area_particle;   // area of single particle
  float area_particle_b; // area of single particle w.r.t given particle's
                         // radius for specific defined lattice zones

  // float mass_particle;				// grain particle mass
  float rho_rock; // mass density of the material that makes up the grain
                  // particles
  float g_accl;   // gravity
  float z_depth;  // vertical distance of interconnectivity of the fluid
  float delta_z_depth;
  double CozenyGrainSize;

  /******************************************************
   * Fluid lattice parameters
   ******************************************************/
  int i_lower, i_upper; // lower and upper limits of rows(y-axis) calling by
                        // experiment.cc
  int j_lower, j_upper; // lower and upper limits of columns(x-axis) calling by
                        // experiment.cc

  float area_node; // area of pressure node
  float
      rad; // (smooth function) radius of circle centred at respective pressure
           // node and edged at the centres of four surrounding nodes.
  float area_circle; // area of circle with radius 'rad'

  float rho_fluid;
  double mu; // visocsity for calculation of alpha & beta variables of pressure
             // lattice
  double compressibility; // compressibility for calculation of alpha & beta
                          // variables of pressure lattice
  double pressure_val;    // pressure (background) parameter of pressure lattice

  double source[dim]
               [dim]; // value of pressure source term in the giverning equation
  double source_a[dim][dim];
  double source_b[dim][dim];

  double PArea[dim][dim];
  float PRadius[dim][dim];

  double alpha[dim][dim]; // aplha (calculated) variable for each pressure node
                          // of pressure lattice
  double beta[dim][dim]; // beta (calculated) variable for each pressure node of
                         // pressure lattice

  double a[dim][dim]; // matrices built up on alpha & beta variables of pressure
                      // lattice
  double b[dim][dim];
  double a_inv[dim][dim]; // inverse matrix for comparing with a[dim][dim] &
                          // b[dim][dim] matrices

  double al[dim][dim][dim];
  double bl[dim][dim][dim];
  double a_invl[dim][dim][dim];

  double sourced_pressure[dim][dim];
  double interim_pressure_dist[dim]
                              [dim]; // corrosponding intermediate pressure
                                     // value of each pressure mode during
                                     // multiplication of two pressure matrices

  double dxc[dim][dim];
  double ddxc[dim][dim];
  double ddyc[dim][dim];
  double dyc[dim][dim];

  double press_background[dim][dim];

  /******************************************************
   * variables, to give pressure to each particle of
   * particle class from pressure lattice through repulsion boxes
   ******************************************************/
  int particle_in_box; // total number of particles in each of four repulsion
                       // boxes associated with pressure nodes
  int Particles;
  int particle[dim][dim];

  double par_area;
  double par_rad;

  double smooth_func;     // tent function or smooth function
  double rho_s[dim][dim]; // solid fraction associated to each pressure node
  double rho_s_circle[dim][dim];
  double rho[dim][dim]; // number density associated to each pressure node

  double vel_x[dim][dim]; // velocity field (average velocity) along x-direction
                          // associated with each pressure node
  double vel_y[dim][dim]; // velocity field along y-direction

  double F_P_x; // pressure force on each particle of particle lattice from
                // pressure lattice
  double F_P_y;

  double delta_P_x;
  double delta_P_y;

  double smooth_i_low;
  double smooth_i_up;
  double smooth_j_back;
  double smooth_j_fwd;

  double smooth_1st;
  double smooth_2nd;
  double smooth_3rd;
  double smooth_4th;

  double delta_P_x_1st;
  double delta_P_x_2nd;
  double delta_P_x_3rd;
  double delta_P_x_4th;

  double delta_P_y_1st;
  double delta_P_y_2nd;
  double delta_P_y_3rd;
  double delta_P_y_4th;

  double F_P_x_1st;
  double F_P_x_2nd;
  double F_P_x_3rd;
  double F_P_x_4th;

  double F_P_y_1st;
  double F_P_y_2nd;
  double F_P_y_3rd;
  double F_P_y_4th;

  double delta_P_i_low;
  double delta_P_i_up;
  double delta_P_j_back;
  double delta_P_j_fwd;

  double delta_P_j_low;
  double delta_P_j_up;
  double delta_P_i_back;
  double delta_P_i_fwd;

  double delta_P_ij_low;
  double delta_P_ij_up;
  double delta_P_ij_back;
  double delta_P_ij_fwd;

  double F_P_i_low;
  double F_P_i_up;
  double F_P_j_back;
  double F_P_j_fwd;

  double F_P_j_low;
  double F_P_j_up;
  double F_P_i_back;
  double F_P_i_fwd;

  double F_P_ij_low;
  double F_P_ij_up;
  double F_P_ij_back;
  double F_P_ij_fwd;

  double dl;
  double dl_i_low;
  double dl_i_up;
  double dl_j_back;
  double dl_j_fwd;

  double deta;
  double sin_deta;
  double cos_deta;

  double deta_i_low;
  double sin_deta_i_low;
  double cos_deta_i_low;

  double deta_i_up;
  double sin_deta_i_up;
  double cos_deta_i_up;

  double deta_j_back;
  double cos_deta_j_back;
  double sin_deta_j_back;

  double deta_j_fwd;
  double cos_deta_j_fwd;
  double sin_deta_j_fwd;

  double F_P_ijx;
  double F_P_ijy;

  double F_P_ix_low;
  double F_P_iy_low;

  double F_P_ix_up;
  double F_P_iy_up;

  double F_P_jx_back;
  double F_P_jy_back;

  double F_P_jx_fwd;
  double F_P_jy_fwd;

  double F_P_ij_x;
  double F_P_ij_y;

  int seal_y_limit;

public:
  double time_a; // time for calculation of alpha & beta variables of pressure
                 // lattice

  int scale_constt; // total dimension in meter of real system to model down

  static const int min_size = 0;
  static const int max_size = dim;

  double pressure_dist[dim][dim]; // pressure value of each pressure node in
                                  // pressure lattice

  float Pressure_x[dim][dim]; // x & y positions of pressure nodes of pressure
                              // lattice
  float Pressure_y[dim][dim];

  float delta_x, delta_y;
  float press_pos_y, press_pos_x;
  double porosity[dim][dim];
  double kappa[dim][dim];
  double displment;

  bool Seal_hor;
  bool ModelHalve;

  float LowSeal_y;
  float UppSeal_y;

  float ModHalf_y;

  // variables for Advection Diffusion

  double concent[dim][dim];

  double vel_darcyx[dim][dim];
  double vel_darcyy[dim][dim];

  double v_darcyx_1st, v_darcyx_2nd, v_darcyx_3rd, v_darcyx_4th, v_darcyy_1st,
      v_darcyy_2nd, v_darcyy_3rd, v_darcyy_4th;

  double conc_1, conc_2, conc_3, conc_4;

  //*******************************************

  Pressure_Lattice(Particle *first); // constructor, function with pointer
                                     // argument, points to particle class
  // Particle *f;
  // // pointer, points to particle class
  Particle *help; // help pointer, points to particle class

  ~Pressure_Lattice(); // destructor

  /********************************
   * function for lattice of dimension
   * 50*50. where rep_Box start from
   * node(0,0), in a way that this continuum
   * node covers four rep_Box
   * ******************************/
  void a_pressure_node_pos_n_area(); // function to locate the x & y positions
                                     // of pressure nodes in pressure lattice

  void Change_Time(float change);

  /********************************
   * function for lattice of dimension
   * 52*52. where rep_Box start from
   * node(1,1) and not from boundary in
   * a way that this node(1,1) covers
   * four rep_Box
   * ******************************/
  void b_pressure_node_pos_n_area();

  /********************************
   * function for lattice of dimension
   * 51*51. where rep_Box start from
   * node(0,0), in a way that this node
   * covers only first rep_Box
   * ******************************/
  void c_pressure_node_pos_n_area();

  /********************************
   * function for lattice of dimension
   * 51*51. where rep_Box start from
   * node(0,0), in a way that this node
   * covers only first rep_Box
   * ******************************/
  void cc_pressure_node_pos_n_area();

  /********************************
   * function for lattice of dimension
   * 53*53. where rep_Box start from
   * node(1,1), in a way that this node
   * covers only first rep_Box
   * ******************************/
  void d_pressure_node_pos_n_area(int resolution);

  void Make_Melt(double change_mu, double change_compressibility);

  void pressure_background(double init_pressure_val, int z_distance);
  void boundary_fixed(int boundry_sourc, double source_val_boundry,
                      int z_distance);
  void boundary_set_hydrostatic(int box_size, int depth, double increase);
  void boundary_set_hydrostatic_gradient(int box_size, int depth);

  void alpha_beta_set(double time); // function for the construction of alpha
                                    // and beta variables of pressure lattice

  void matrices_set_1(); // functions for the construction of a[dim][dim] &
                         // b[dim][dim] matrices of pressure lattice
  void matrices_set_2();

  void matrix_a_1st();
  void matrix_a_2_step();
  void multiplication_x();
  void multiplication_y();

  void
  transpose_P(); // function to take transpose of pressure matrix during ADI
  void transpose_trans_P();

  void invert(); // function to take inversion of a[dim][dim] pressure matrix

  void multiplication_1st();
  void multiplication_inv();

  void dan_matrices_set_1();
  void dan_multiplication_1st();
  void dan_invert();
  void dan_multiplication_inv();
  void dan_matrices_set_2();

  void press_nodes(int press_y, int press_x);

  void
  numberDensity_porosity_kappa(); // function for particle density, porosity
                                  // and permeability by Kozney-carman equation

  void source_pressure();
  void pressure_with_source(double time);
  void background_press();
  void ohne_background_press();
  void mit_background_press();

  void
  a_solid_frac_n_vel(Particle **list); // consider lattice edges as boundary
                                       // (dim = 50)+ square smooth function
  void cir_solid_frac_n_vel(
      Particle **list); // consider lattice edges as boundary (dim = 50) +
                        // circular smooth function

  void b_solid_frac_n_vel(
      Particle **list); // with extra rows and columns for boundary (dim = 52) +
                        // square smooth function
  void b_dasmal_solid_frac_n_vel(Particle **list);
  void bb_dasmal_solid_frac_n_vel(Particle **list);

  void
  c_solid_frac_n_vel(Particle **list); // consider lattice edges as boundary
                                       // (dim = 51)+ square smooth function
  void c_dasmal_solid_frac_n_vel(
      Particle **list); // consider lattice edges as boundary (dim = 51)+ square
                        // smooth function + high resolution dasmal.elle file
  void cc_dasmal_solid_frac_n_vel(
      Particle *
          *list); // consider lattice edges as boundary (dim = 101)+ square
                  // smooth function + high resolution dasmal.elle file

  void
  d_solid_frac_n_vel(Particle **list); // consider lattice edges as boundary
                                       // (dim = 53)+ square smooth function
  void
  cir_d_solid_frac_n_vel(Particle **list); // consider lattice edges as boundary
                                           // (dim = 53)+ square smooth function
  void d_dasmal_solid_frac_n_vel(Particle **list);
  void dd_dasmal_solid_frac_n_vel(Particle **list);

  void e_solid_frac_n_vel(Particle **list); // way of solid fraction calculation
                                            // is different than above functions

  void solid_frac_n_vel(
      Particle **list); // smoothing function w.r.t circle edged at the centre
                        // of neighbouring continuous nodes divide by the area
                        // of cicle, excluding the boxes out of boundaries
  void solid_frac_circled(
      Particle *
          *list); //// smoothing function w.r.t circle edged at the centre of
                  ///neighbouring continuous nodes divide by the area of boxes
                  ///(included), excluding the boxes out of boundaries

  void pressure_particles(Particle **list);
  void a_pressure_particles(Particle **list);
  void circled_pressure_particles(Particle **list);

  void b_pressure_particles(Particle **list);
  void b_dasmal_pressure_particles(Particle **list);
  void bb_dasmal_pressure_particles(Particle **list);

  void c_pressure_particles(Particle **list);
  void c_dasmal_pressure_particles(Particle **list);
  void cc_dasmal_pressure_particles(Particle **list);

  void d_pressure_particles(Particle **list);
  void cir_d_pressure_particles(Particle **list);
  void d_dasmal_pressure_particles(Particle **list);
  void dd_dasmal_pressure_particles(Particle **list);

  void pressure_incre_source(int press_y, int press_x, double additional_press);
  void Pressure_incre_twoXsources(int y_press, int x_press_diff,
                                  double press_additional);
  void Pressure_incre_twoYsources(int x_press, int y_press_diff,
                                  double press_additional);
  void Pressure_source_n_sink(int press_y_a, int press_x_a, int press_y_b,
                              int press_x_b, double additional_press,
                              double reduced_press);

  void pressure_point_sink(int press_y, int press_x, double reduced_press);
  void pressure_sink(float y_press_min, float y_press_max, float x_press_min,
                     float x_press_max, double reducing_press);
  void pressure_sink_cont(float y_press_min, float y_press_max,
                          float x_press_min, float x_press_max,
                          double reducing_press);

  void rand_source_PressIncrement(int press_y, int press_x,
                                  double pressure_val);

  void source_hor_PressIncrement(int press_y, int press_min_x, int press_max_x,
                                 double pressure_val);
  void source_ver_PressIncrement(int press_x, int press_min_y, int press_max_y,
                                 double pressure_val);
  void source_area_PressIncrement(int press_min_y, int press_max_y,
                                  int press_min_x, int press_max_x,
                                  double pressure_val);

  void test_pressure_with_source(double time);
  void test_c_pressure_particles(Particle **list);
  void test_c_dasmal_pressure_particles(Particle **list);
  void test_cc_dasmal_pressure_particles(Particle **list);

  void grad_pressure(float y_limit);

  void d_1x4_Boxes(Particle **list, int position, int l, int m);
  void d_2x4_Boxes(Particle **list, int position, int l, int m);
  void d_4x4_Boxes(Particle **list, int position, int l, int m);

  void d_LowLeft_4a_Boxes(int pos_count, int q, int r);
  void d_LowRight_4b_Boxes(int pos_count, int q, int r);
  void d_UppLeft_4c_Boxes(int pos_count, int q, int r);
  void d_UppRight_4d_Boxes(int pos_count, int q, int r);

  void d_a_Boxes(Particle **list, int n, int o);
  void d_b_Boxes(Particle **list, int n, int o);
  void d_c_Boxes(Particle **list, int n, int o);
  void d_d_Boxes(Particle **list, int n, int o);

  // advection diffusion boxes

  void conc_1x4_Boxes(Particle **list, int position, int l, int m);
  void conc_2x4_Boxes(Particle **list, int position, int l, int m);
  void conc_4x4_Boxes(Particle **list, int position, int l, int m);

  void conc_LowLeft_4a_Boxes(int pos_count, int q, int r);
  void conc_LowRight_4b_Boxes(int pos_count, int q, int r);
  void conc_UppLeft_4c_Boxes(int pos_count, int q, int r);
  void conc_UppRight_4d_Boxes(int pos_count, int q, int r);

  void conc_a_Boxes(Particle **list, int n, int o);
  void conc_b_Boxes(Particle **list, int n, int o);
  void conc_c_Boxes(Particle **list, int n, int o);
  void conc_d_Boxes(Particle **list, int n, int o);

  /***************************
   * advection diffusion
   *
   * **************************/

  void setconcent(double concentration);
  void concentration_particles(Particle **list);
  void get_concentration(Particle **list);
  void setDarcyVelocity(Particle **list);
  void SolveadvDisExplicitTransient(double fluid_v, double dif_const);
  void ChangeConcentration(int box_x, int box_y, double change);
};
#endif
