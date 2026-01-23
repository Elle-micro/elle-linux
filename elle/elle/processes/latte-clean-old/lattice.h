/******************************************************
 * Lattice Spring Code 2.0
 *
 * Class Lattice in lattice.h
 *
 * Basic header for Fracture Class Lattice
 *
 * Daniel and Jochen oct. 2001
 * Daniel and Jochen Feb 2002 to Feb 2003
 * Oslo, Mainz
 *
 * Daniel Dec. 2003
 *
 * We thank Anders Malthe-Srenssen for his enormous help
 * and introduction to these codes
 *
 * new 2004/2005
 ******************************************************/
#ifndef _E_lattice_h
#define _E_lattice_h

#include "fluid_lattice.h" // header for Pressure Lattice class
#include "heat_lattice.h"  // header for the heat diffusion class
#include "particle.h" // include the particle class so we can have pointers of type particle
#include "pressure_lattice.h" // header for Pressure Lattice class

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>

using namespace std;
using std::vector;

/******************************************************
 * the lattice class is the highest class at moment
 * is called from the main function (well, in elle the mike_elle.cc function
 * it will build the lattice in its constructor including
 * the particle list. This lists can be
 * reached from lattice which has Runpointers and pointers
 * pointing to the beginning and end of list. In order
 * to reach it have to refer to a refParticle.  The list is
 * circular, each member points to the next in the list and
 * last member points to the first one, and also backwards.
 * the lattice class can also make a full relaxation and
 * can deform i.e. apply boundary conditions
 * Should also define the walls in this class, not done yet
 *
 * introduced several deformation times in order to preload
 * system, daniel March 4 2003
 *******************************************************/

class Lattice {
public:
  //------------------------------------------------------------------------------
  // first define the variables for the class that are important for all
  // the functions used for the class
  //------------------------------------------------------------------------------

  float dimx;    // lattice x dimension, 1.0 as default
  float dimy;    // lattice y dimension, 1.0 as default
  float xyratio; // xy ratio depending on lattice type

  double default_Young, default_Poisson;

  int nbBreak;       // counts broken springs in Lattice after current time step
  int local_nbBreak; // local counter for fracture plot
  int internal_break; // counter for dumping fracture pictures

  int nbBreak_last; // counts broken springs in Lattice before current time step
  int nbBreak_current; // counts broken springs in current time step

  int LowSeal_bnd_nb, Seal_bnd_nb,
      UppSeal_bnd_nb; // counter for zonally broken bond nr. in experiment for
                      // hor_seal
  int nbLowSeal_current, nbSeal_current,
      nbUppSeal_current; // counts zonally broken springs in current time step
                         // in experiment for hor_seal
  int nbLowSeal_last, nbSeal_last,
      nbUppSeal_last; // // counts zonally broken springs in Lattice before
                      // current time step in experiment for hor_seal

  int ModLowr_bnd_nb, ModUppr_bnd_nb; // counter for bronken bonds in upper and
                                      // lower parts of Model.
  int nbModLowr_current, ModUppr_current; // in current time step
  int nbModLowr_last, nbModUppr_last;     // before current time step

  int def_time;   // time for fractures from experiment
  int def_time_p; // time counter for steps of progressive deformation (pure
                  // shear)
  int def_time_u; // time counter for steps of progressive def. (uniaxial
                  // compression)
  int fract_Plot_Numb; // after how many fracture do I plot
  float relaxthres;    // relaxation threshold
  int relaxmax;
  int layer[10000]; // numbers of grains or Flynns in a layer

  int Ulist[100000]; // list of unused Unodes for ridges
  int Ulistcount;

  int particle_age;

  bool set_max_pict; // set a maximum (default false)
  int max_pict;      // maximum number of pictures dumped
  int num_pict;      // current number of pictures dumped (counter)

  float Pi; // number pi

  int visc_rel; // flag for viscous relax routines in relaxation

  float pressure_scale; // scale for pressure (at moment 1 = 10 GPa)
  float pascal_scale;   // scales MPa to Pascal

  double area_box; // area of whole box for pure shear

  double average_break;
  double average_springs;
  double average_springf;

  int particlex; // number of particles in x
  int particley; // number of particles in y

  float boundary_strength;
  float boundary_constant;

  //--------------------------------------------------------------------------------
  // walls of the lattice
  //--------------------------------------------------------------------------------

  bool walls;
  float wall_rep;
  float upper_wall_pos;
  float lower_wall_pos;
  float right_wall_pos;
  float left_wall_pos;

  int grain_counter; // initial grain number

  float polyg_area;  // area for viscous step
  float polyg_area8; // area for the octagon (BoxRad())

  bool adjust; // for debugging only

  bool debugflag, heat;

  bool shear;

  int debug_nb;

  float alpha;

  //-----------------------------------------------------------------------------
  // and define some pointers
  //-----------------------------------------------------------------------------

  Particle *firstParticle;  // pointer to particle list begin
  Particle *lastParticle;   // pointer to last particle in list
  Particle *runParticle;    // run-pointer for list
  Particle *preRunParticle; // 2. run-pointer

  //---------------------------------------------------------------------------
  // these arrays have to have the same dimensions !
  //---------------------------------------------------------------------------

  Particle *repBox[2500000]; // The neighbour box for Repulsion

  int node_Box[2500000]; // neighbour box for nodes

  float grain_young[100000]; // saves youngs moduli of grains

  //-----------------------------------------------------------------
  //  maximum Elle nodes
  //-----------------------------------------------------------------

  int next_inBox[1000000]; // vector for nodes to save node number for node box

  bool nodes_Check[1000000]; // node vector for node check if connected to
                             // lattice

  //-----------------------------------------------------------------------------
  // and define a reference particle. this will be the first particle in the
  // list note that this has to change if that particle is deleted for whatever
  // reason
  //-----------------------------------------------------------------------------

  Particle *driftlist[200]; // The list for drifting

  Particle refParticle; // the reference particle

  int numParticles; // number of particle

  int highest_grain; // highest grain number
  int olivine_grain_nb, spinel_grain_nb;

  bool visc_flag, transition;
  float surf_pos;
  float volume;
  bool sheet;

  bool grav;
  bool grav_Press;

  bool wrapping;
  float right_lattice_wall_pos;
  float left_lattice_wall_pos;

  double cal_time, press_cal_time;

  void SetRightLatticeWallPos();

  //-----------------------------------------------------------------
  //  variables of Pressure_lattice
  //-----------------------------------------------------------------
  Pressure_Lattice
      *P; // pointer to call variables and functions of Pressure_Lattice class
  Fluid_Lattice *Fluid;
  double
      diff; // error check to get steady state pressure in Calculate_Pressure()
  // double err_last;
  // double err_new;
  int i, j; // counter for Pressure_lattice functions in Lattice.cc

  float sourc_Ymin, sourc_Ymax, sourc_Xmin, sourc_Xmax; // seal edges
  int P_x_mean; // Gaussian mean of pressure distribution.

  // functions for CSV files

  bool doneWithHeaders;
  void setMapForCSV(int i);
  void writeCSVfile(std::map<std::string, double> attribMap);
  ofstream myExpCsv;
  int experimentTime, saveCSVInterval;
  void getTimeForCsv(int expTime, int saveInt);

  //------------------------------------------------------------------------------
  // now define the functions
  //------------------------------------------------------------------------------

  Lattice();    // Constructor
  ~Lattice(){}; // Destructor

  double Standard_Breaking_Probability(Particle *prtcl, int nb);
  bool Standard_Spring_Failure();
  void Set_TimeFrac(int time);
  void Activate_Lattice();
  void MakeLattice(int type);   // build lattice
  void InitiateGrainParticle(); // find grains,boundaries,nodes etc.
  void MakeRepulsionBox();      // build repulsion box
  void SetBoundaries();         // boundary conditions
  void FindUnodeElle();         // connect to Unodes
  bool FullRelax();             // simple full relaxation(overrelaxation)
  void Relaxation();            // Relaxation
  int HighestGrain();           // find highest grain in Elle
  void UpdateElle();            // talk a bit to Elle
  void UpdateElle(bool movenode);
  int TimeStep;

  //------------------------------------------------------------------
  // functions for external use called from mike.elle
  //------------------------------------------------------------------

  void SetSpringConstants();
  void SetWeightForce();

  // set wall boundaries

  void SetWallBoundaries(int both, float constant);

  // set a sinusoidal variation of youngs moduli

  void SetSinAnisotropy(float nb, float ampl);

  void Set_Absolute_Box(double size, int gr);

  // reread the elle file

  void GetNewElleStructure();
  void GetNewElleStructureTwo();

  void ApplyGravity(float height);

  void Gravity(float depth, float density);
  void Gravity_no_Move(float depth, float density);
  void Adjust_Gravity();
  void Gravity_noload_nofluid();

  // security stop for pict dump

  void Set_Max_Pict(int max);

  void RelaxSheet(float t_step);

  void ActivateSheet(float visc, float young);

  void ReleaseBoundaryParticlesY(int fracture, int release);

  void ChangeYoung(float change);

  // changes grain constant and tensile strength

  void WeakenGrain(int nb, float constant, float visc, float break_strength);

  void WeakenGrainAttribute(double nb, float constant, float visc,
                            float break_strength);

  // changes grainboundary constant and tensile strength

  void MakeGrainBoundaries(float constant, float break_strength);

  void SetViscous(float max);

  void Only_Extension(bool extension);

  // function to adjust modulus of grainboundary springs

  void CalculateGravityOverburden(float depth, float height);

  void AdjustConstantGrainBoundaries();

  // function to adjust constants of particles so that springs have mean value
  // of neighbours

  void AdjustParticleConstants();

  // plot after certain amount of fractures or not

  void SetFracturePlot(int numbfrac, int local);

  void MarkMohoParticles(double y);

  void DamageModel(float change, float stress, int damage);

  void OutputTopParticles(double time);
  void OutputMohoParticles(double time);

  float DeformLatticeExtend(float move, int plot);

  float DeformLatticeXMid(float move, int plot, float neck_height, float left,
                          float right);

  float DeformLatticeXMidTwo(float move, float ymov, int plot,
                             float left_margin_width, float right_margin_width,
                             float left_rift_width, float right_rift_width);

  float DeformLatticeXMidTwo_b(float move, float ymov, int plot);
  // Boundary conditions deform from upper side, do an average (uniaxial
  // compression)

  float DeformLatticedriftmid(float move, int plot, float change, float leftup,
                              float rightup, float leftdown, float rightdown);

  float DeformLatticedriftmidnew(float move, int plot, float change,
                                 float leftup, float rightup, float leftdown,
                                 float rightdown);

  float DeformLatticedriftmidnew2(float move, int plot, float change,
                                  float leftup, float rightup, float leftdown,
                                  float rightdown);
  float DeformLatticedriftmidnew3(float move, int plot, float change,
                                  float leftup, float rightup, float leftdown,
                                  float rightdown);
  float DeformLatticedriftmidnew4(float move, int time, float changeone,
                                  float leftup, float changetwo, float lefttwo,
                                  float changethree, float leftthree,
                                  float leftdown, float grad);

  float DeformLatticedriftmidnew5(float move, int time, float changeone,
                                  float leftup, float changetwo, float lefttwo,
                                  float changethree, float leftthree,
                                  float leftdown, float grad);

  float DeformLatticedrift6(float move, int time, float changeone, float leftup,
                            float changetwo, float lefttwo);
  float InitDriftList();

  float DeformLattice(float move, int plot);

  // deform volume constant with average

  float DeformLatticePureShear(float move, int plot);

  // deform volume not constant

  float DeformLatticePureShearAreaChange(float move, float area_change,
                                         int plot);

  // Deform Lattice without an average

  float DeformLatticeNoAverage(float move, int plot);

  // Deform Lattice from upper and lower boundaries, no average

  float DeformLatticeNoAverage2side(float move, int plot);

  void DeformLatticeExtension(float move, int plot, float dx, float shear);

  float Save_Particle_Velocity();

  // Deform Lattice from upper and lower boundaries, do an average

  float DeformLatticeNewAverage2side(float move, int plot);

  // do a viscous retardation step

  void ViscousRelax(int dump, float timestep);

  // functions for the viscous retardation step

  void SetSpringRatio();

  double VirtualNeighbour(int spring, int xy);

  void Adjust_Radii(); // adjust radii so that area is kept constant

  void Copy_NeigP_List(); // help function

  void BoxRad();

  // change the relaxation threshold

  void ChangeRelaxThreshold(float change);

  /*********************************************************
   * CHANGES AND DISTRIBUTIONS
   *********************************************************/

  // shrink or grow one particle

  void ChangeParticle(int number, float shrink, int plot);

  // shrink particle that is stressed most (highest mean stress)

  void ChangeParticleStress(float shrink, int plot);

  // shrink or expand a whole grain

  void ShrinkGrain(int nb, float shrink, int plot);

  void SetCosViscSheetMod(float factor, float factorb);

  // shrink the whole box

  void ShrinkBox(float shrink, int plot, int yes);

  // define layer using Elle grains

  void WeakenHorizontalLayer(double y_min, double y_max, float constant,
                             float break_strength);

  void WeakenHorizontalBox(double y_min, double y_max, double x_min,
                           double x_max, float constant, float vis,
                           float break_strength, float rad_frac);

  // horizontal layer using only particles, independent of elle grains
  void WeakenTiltedParticleLayerabsolute(double y_min, double y_max,
                                         double shift, float constant,
                                         float vis, float break_strength,
                                         float rad_frac);
  void WeakenHorizontalParticleLayer(double y_min, double y_max, double x_min,
                                     double x_max, float constant, float vis,
                                     float break_strength, float rad_frac);
  void WeakenHorizontalParticleLayerX(double y_min, double y_max, double x_min,
                                      double x_max, float constant, float vis,
                                      float break_strength);
  void PoreDepnttBrkStr_Seal(float y_min, float y_max, double parRad_modl);

  // a whole row of even horizontal layers

  void MakeHorizontalLayers(float young_a, float viscous_a, float break_a,
                            float thick_a, float young_b, float viscous_b,
                            float break_b, float thick_b, int av);

  // make a tilted layer starting at ymin to ymax

  void WeakenTiltedParticleLayer(double y_min, double y_max, double shift,
                                 float constant, float vis,
                                 float break_strength, float rad_frac);

  // make a tilted layer starting at xmin to xmax

  void WeakenTiltedParticleLayerX(double x_min, double x_max, double shift,
                                  float constant, float vis);

  // change constants of all particles

  void WeakenAll(float constant, float viscosity, float break_strength);

  // set a horizontal anisotropy using small particles rows of 1 - 3 particles

  void SetAnisotropy(float break_strength, float youngs_mod, float viscosity,
                     float ratio, float length);

  // set a horizontal anisotropy, little mica flakes with 3 to 7 particle length
  // and varying youngs moduli

  void SetAnisotropyRandom(float max, float boundary);
  void SetDistributionPorosity(float dis);

  // Change the Young modulus of grains by a factor

  void Change_Young(float young_factor);

  // sets a Gaussian distribution on spring constants of whole grains

  void SetGaussianSpringDistribution(double g_mean, double g_sigma);

  // set a Gaussian distribution on breaking stengths of all springs

  void SetGaussianStrengthDistribution(double g_mean, double g_sigma);

  // set a Gaussian distribution on spring constant of all particles

  void SetGaussianYoungDistribution(double g_mean, double g_sigma);

  // sets pseudorandom distribution of spring constants of grains and of
  // breaking strength of all springs

  void SetPhase(float str_mean, float str_size, float br_mean, float br_size);

  void SetVariationBreakingThreshold(float var);
  void Dumplayer(double strain);
  /*********************************************************
   *  STATISTICS DATADUMPS
   *********************************************************/

  void DumpStatisticPorosity(double y_box_min, double y_box_max,
                             double x_box_min, double x_box_max, double strain);

  void DumpStatisticGravityStress(int particle_row);

  // write out the x and y coordinates plus strain of particles on the surface

  void DumpStatisticSurface(double strain);

  // write out x and y coordinates plus strain of particles on surface in a
  // region

  void DumpTimeStatisticSurface(double strain, double xmin, double xmax);

  // function to dump statistics of a single grain

  void DumpStatisticStressGrain(double strain, int nb);

  // function to dump statistics of two grains

  void DumpStatisticStressTwoGrains(double strain, int nb1, int nb2);

  // function to dump statistics of a predefined box

  void DumpStatisticStressBox(double y_box_min, double y_box_max,
                              double x_box_min, double x_box_max,
                              double strain);

  void Pressure_force_range_text(double y_box_min, double y_box_max,
                                 double x_box_min, double x_box_max);

  void DumpStat_Press_stress_pore_perm(float x_box_min, float x_box_max,
                                       float y_box_min, float y_box_max);
  void DumpStat_TwoSourc_Pres_stres_por_perm_BrkBnd(float y_ModHalf);

  void DumpStat_PressMDMxMnStress_LHRes_c(int press_node_y, int press_node_x);
  void DumpStat_PressMDMxMnStress_HRes_c(int press_node_y, int press_node_x);
  void DumpStat_PressMDMxMnStress_LHRes_d(int press_node_y, int press_node_x);
  void DumpStat_PressMDMxMnStress_HRes_d(int press_node_y, int press_node_x);

  void DumpStat_PressMDStress_LHRes_c_d_verSeal(float min_press_limit_x,
                                                float max_press_limit_x,
                                                int y_nod, int x_nod);
  void DumpStat_PressMDStress_HRes_c_d_verSeal(float min_press_limit_x,
                                               float max_press_limit_x,
                                               int y_nod, int x_nod);
  void DumpStat_PressMDStress_LHRes_c_d_horSeal(float min_press_limit_y,
                                                float max_press_limit_y,
                                                int y_nod, int x_nod);
  void DumpStat_PressMDStress_HRes_c_d_horSeal(float min_press_limit_y,
                                               float max_press_limit_y,
                                               int y_nod, int x_nod);

  void
  DumpStat_PressMDStress_LHRes_c_d_horSeal_randSource(float min_press_limit_y,
                                                      float max_press_limit_y);
  void DumpStat_PressMDStress_LHRes_c_d_randSource();

  void DumpStat_Boxes_Stress(int press_node_y, int press_node_x);
  void DumpStat_Boxes_Stress_LHRes_c_d(int press_node_y, int press_node_x);
  void DumpStat_UnderSealUpper_Stress_LHRes_c_d(float sealy_min,
                                                float sealy_max);
  void DumpStat_Boxes_Stress_HRes_c_d(int press_node_y, int press_node_x);

  void DumpStat_Boxes_Press_brkBond(int press_node_y, int press_node_x);
  void DumpStat_Boxes_Press_brkBond_LHRes_c_d(int press_node_y,
                                              int press_node_x);
  void DumpStat_UnderSealUpper_Press_brkBond_LHRes_c_d(float sealy_min,
                                                       float sealy_max);
  void DumpStat_Boxes_Press_brkBond_HRes_c_d(int press_node_y,
                                             int press_node_x);

  void DumpStat_Boxes_Pore_Perm(int press_node_y, int press_node_x);
  void DumpStat_Boxes_Pore_Perm_LHRes_c_d(int press_node_y, int press_node_x);
  void DumpStat_Boxes_Pore_Perm_HRes_c_d(int press_node_y, int press_node_x);

  void DumpStat_UnderSealUpper_Pore_Perm_LHRes_c_d(float sealy_min,
                                                   float sealy_max);

  // additional functions for the viscoelastic part, hidden

  float Alen(int i, int j);

  void SetSpringAngle();

  void UpdateNeighList();

  void EraseNeighList();

  void UpdateSpringLength();

  void Tilt();

  float DeformLatticeSimpleShear(float move, int plot);

  float DeformLatticeX(float move, int plot);

  void WeakenVerticalParticleLayer(double x_min, double x_max, double y_min,
                                   double y_max, float constant, float vis,
                                   float break_strength);

  void SetRealDensityYoung(double x_max, double x_min, double y_max,
                           double y_min, double den, double you);

  void InitRealDensityYoung(double den, double you);

  void InitInternalYoung(double young_standard);

  double SpreadingStrain(double m, double box_size, double years);

  void SetViscosity(double x_max, double x_min, double y_max, double y_min,
                    double visc);

  void InitViscosity(double visc);

  void SetBreakingStrength(double x_max, double x_min, double y_max,
                           double y_min, double factor);

  void UnsetNobreakUpperRow();

  void SetBreakingStrengthGradient(double x_max, double x_min, double y_max,
                                   double y_min, double c);

  void SetWallBoundariesX(int release, int fracture, int walls, float constant);

  /*****************************************************************
   *
   * ------------functions for Pressure_Lattice--------------------
   *
   * ***************************************************************/

  void Pressure_Insert_Random_Node(int res, double press_increment);

  void Initialize_Pressure_Lattice(double Press_Background, int scale);

  void Pressure_init_RandSource(int rand_min_y, int rand_max_y, int rand_min_x,
                                int rand_max_x, double pressure_min,
                                double pressure_max);
  void Pressure_init_RandSource_hor_lay(int hor_lay, int rand_min_x,
                                        int rand_max_x, double pressure_min,
                                        double pressure_max);

  void Pressure_RandSource_fix(int rand_min_y, int rand_max_y, int rand_min_x,
                               int rand_max_x, double press_increment);
  void Pressure_RandSource(int rand_min_y, int rand_max_y, int rand_min_x,
                           int rand_max_x, double pressure_min,
                           double pressure_max);

  void Pressure_RandFixSource_hor_lay(int lay_y, int rand_min_x, int rand_max_x,
                                      double pressfix);

  void Pressure_RandSource_hor_lay(int lay_y, int rand_min_x, int rand_max_x,
                                   double pressure_min, double pressure_max);
  void Pressure_RandSource_ver_lay(int lay_x, int rand_min_y, int rand_max_y,
                                   double pressure_min, double pressure_max);

  void Pressure_xGauss_RandSource_area(int y_rand_min, int y_rand_max,
                                       int x_rand_min, int x_rand_max,
                                       double press_increment, float P_sigma);
  void Pressure_yGauss_RandSource_area(int y_rand_min, int y_rand_max,
                                       int x_rand_min, int x_rand_max,
                                       double press_increment, float P_y_mean,
                                       float P_sigma);
  void Pressure_2dCircularGauss_RandSource_area(
      float y_rand_min, float y_rand_max, float x_rand_min, float x_rand_max,
      double press_increment, float P_y_mean, float P_sigma);
  void Pressure_2dEllipticGauss_RandSource_area(
      float y_rand_min, float y_rand_max, float x_rand_min, float x_rand_max,
      double press_increment, float P_y_mean, float P_y_sigma, float P_x_sigma);

  void Pressure_init_source(
      int press_node_y, int press_node_x,
      double press_increment); // initial value of pressure source
  void Pressure_init_twoXsources(
      int y_press_node, int x_press_node_diff,
      double increment_press); // initial values of two pressure point sources
                               // at same y but different x positions
  void Pressure_init_twoYsources(int x_press_node, int y_press_node_diff,
                                 double increment_press);
  void Pressure_init_source_n_sink(int press_node_y_a, int press_node_x_a,
                                   int press_node_y_b, int press_node_x_b,
                                   double press_decrement,
                                   double press_increment);

  void pressure_for_run(float seal_y_min);

  void Pressure_init_sink(int press_node_y, int press_node_x,
                          double press_decrement);
  void Pressure_constt_SinkArea(float y_box_min, float y_box_max,
                                float x_box_min, float x_box_max,
                                double press_decrement);
  void Pressure_SinkArea_cont(float y_box_min, float y_box_max, float x_box_min,
                              float x_box_max, double press_decrement);

  void Pressure_BuildUp(int press_node_y, int press_node_x,
                        double press_increment);
  void Pressure_BuildUp_twoXsources(int y_press_node, int x_press_node_diff,
                                    double increment_press);
  void Pressure_BuildUp_twoYsources(int x_press_node, int y_press_node_diff,
                                    double increment_press);
  void Pressure_build_Up_n_Down(int press_node_y_a, int press_node_x_a,
                                int press_node_y_b, int press_node_x_b,
                                double press_increment, double press_decrement);

  void Pressure_BuildDown_point(int press_node_y, int press_node_x,
                                double press_decrement);
  void Pressure_BuildDown(float y_box_min, float y_box_max, float x_box_min,
                          float x_box_max, double press_decrement);
  void Pressure_cont_BuildDown(float y_box_min, float y_box_max,
                               float x_box_min, float x_box_max,
                               double press_decrement);

  void Calculate_Pressure_1st_time(int stable, double step_time);
  void Calculate_Pressure();

  void Particle_Fluid_parameters();

  void Pressure_Insert_Random_Nodey(double ymin, double ymax, int res,
                                    double press_increment);

  void Pressure_init_Boundary_cond(int LboundrySourc, double boundry_val,
                                   int vertical_z_depth);
  void Pressure_Boundary_cond(int LboundrySourc, double boundry_val,
                              int vertical_z_depth);
  void Pressure_initialize_Hydrostatic(int box_size, int depth,
                                       double add_pressure);
  void Pressure_initialize_Hydrostatic_gradient(int box_size, int depth);

  void Pressure_HorSource(int lay_y, int min_x, int max_x,
                          double press_increment);
  void Pressure_VerSource(int lay_x, int min_y, int max_y,
                          double press_increment);
  void Pressure_AreaSource(int min_y, int max_y, int min_x, int max_x,
                           double press_increment);

  void Make_Seal(float seal_y_min, float seal_y_max, float seal_x_min,
                 float seal_x_max, float rad_mean, float rad_sigma);
  void Make_Seal_Slot(float seal_y_min, float seal_y_max, float seal_x_min,
                      float seal_x_max, float rad_frac);

  void WeakenHorizontal_Slot(double y_min, double y_max, double x_min,
                             double x_max, float constant, float vis,
                             float break_strength);
  void Fluid_Insert_Random_Nodexy(double ymin, double ymax, double xmin,
                                  double xmax, double press_increment);

  void Zonally_brkBond_Cal(); // function to count nr. of zonally broken bonds.
  void HalfMod_brkBond_Cal(); // function to count breaking bonds in two halves
                              // of model.

  void Real_Density_Young(float seal_y_min, float seal_y_max, double dens_seal,
                          double dens, double yong_seal, double yong, int scal);
  void Press_Initial_Gravity(int depth, int height);

  void DumpStat_EdgeSeal_Pore_Perm_LHRes_c_d(float sealy_min, float sealy_max);
  void DumpStat_EdgeSeal_Press_brkBond_LHRes_c_d(float sealy_min,
                                                 float sealy_max);
  void DumpStat_EdgeSeal_Stress_LHRes_c_d(float sealy_min, float sealy_max);
  void Melt_Parameters(double viscosity, double compressibility);
  void Healing(double heal_dist, float prob, int time_heal, float springfactor,
               float springstrength);
  void HealingTMol(double heal_dist, float prob, int time_heal,
                   float springfactor, float springstrength);
  void Fracture_Reaction(float fluid_threshold, float prob, float max_young);

  void Loose_Springs(bool boundary);
  void Shift_Particles(float variation);
  void Connect_New_Lattice(bool changeRad);
  void Adjust_Springs_Particles();

  /******************
   * Advection Diffusion
   *
   * ***************************/

  void Calculate_Concentration(double vel, double dif);
  void Initialize_Concentration_Lattice(double Press_Background);
  void SetDarcyforAdv();
  void Interpolation_Concentration();
  void Get_Concentration_Directly();
  void Set_Calcite();
  void Add_Concentration(int nodex, int nodey, double concentration);
  void Make_Dolomite(double reaction_time, int step, double shrink,
                     double volume, double scale);
  void Make_Metal(double reaction_time, int step, double shrink, double volume,
                  double scale);
  void Reaction_Arrhenius(double reaction_time, double scale);
  void Permeable_Grain_Boundary(double factor);
  void Permeable_Grain(double factor, int grain);

  void Initialize_Fluid_Lattice(double pressure, double scale, double conc,
                                double time_factor);
  void Calculate_Fluid_Pressure(double upper_pressure, double lower_pressure,
                                int iteration, int sides, int average,
                                double additional);
  double Calculate_Fluid_Pressure_timeadjust(double upper_pressure,
                                             double lower_pressure,
                                             int iteration, int sides,
                                             int average, double additional);
  void Input_Fluid_Lattice(double pressure, int x, int y);
  void Fluid_Parameters(double fracture_Perm, int boundary, double grain_size);
  void Hydrostatic_Fluid(int depth);
  void Fluid_Insert_Random_Node(double ymin, double ymax,
                                double press_increment);
  void Calculate_Con(double vel, double dif);
  void Calculate_Advection(double time, double space, int boundary, int method);
  void Calculate_Advection_Metal(double time, double space, int boundary);
  void Calculate_Advection_Temperature(double time, double space, int boundary);
  void Calculate_Diffusion(double constant, double time, double space);
  void Calculate_Diffusion_Metal(double constant, double time, double space);
  void Calculate_Diffusion_Temperature(double constant, double time,
                                       double space);
  void Dump_concentration_growth(double strain);

  // Exchange Reactions

  void Initiate_Mineral(double background);
  void Set_Mineral(int Grain, int Mineral);
  void Calculate_Diffusion_Solid(double time, double space);
  void Calculate_Diffusion_Fluid(double time, double space);
  void Calculate_Diffusion_Fluid_Explicit(double time, double space);
  void Exchange_Solid_Fluid(double constant, double time, double space);
  void Set_Fluid(double back, double boundary);

  void InitHeat(int T, double critAge);
  void Initialize_Heat(double xpos, double xpos2, double temp, double critAge);
  void Initialize_Heat_Box(double xpos, double xpos2, double ypos, double ypos2,
                           int T, double critAge);
  void solve_Heat(double space, double time, double constant, double temp,
                  double critAge, int tfac);
  void Adjust_break();
  void Loose_Unodes(double xmin, double xmax);
  void Grow_Particle(double posx, double posy, int boxpos, double sheetYoung,
                     double sheetVisc, int age, int T);
  void Find_Empty_Box(double age, int temp);
  void Cooling(double ref);
  void CoolingDirect(double ref, double refspring, double visc);
  void CoolingDirectExp(double ref, double refspring, double visc);

  void Dump_Young(int start, int end);

  void Get_Aperture();
  void ShrinkTemperature(int shrink);
  void ViscousStep(double time);
  void ChangePoisson(double poisson);
  void ChangeRelaxMax(int max);
};

#endif
