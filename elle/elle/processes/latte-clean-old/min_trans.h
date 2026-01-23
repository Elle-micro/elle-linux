/******************************************************
 * Spring Code Mike 1.1
 *
 * Class Phase_Lattice in phase_lattice.h
 *
 * Basic header for Phase-transition Class Phase_Lattice
 *
 * Daniel and Jochen, March 2003
 * Mainz
 *
 * Daniel Dec. 2003
 ******************************************************/
#ifndef _E_min_trans_lattice_h
#define _E_min_trans_lattice_h

#include "heat_lattice.h"
#include "min_trans.h"
#include "phase_base.h"
#include "phase_lattice.h"

#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>

#include <vector>

using namespace std;

class Min_Trans_Lattice : public Phase_Lattice {
public:
  // till-variablen

  double actual_time, now, time_interval, old_time;
  bool timeflag,         // multipurpose
      reaction_finished, // flag for heat-conduction and time-management in
                         // gbm()
      logicalflag, nucleation_occured;
  double mole_per_particle, highest_prob, lowest_prob;
  double pressure, org_wall_pos_y, org_wall_pos_x;

  int beware_grain1, beware_grain2, beware_grain3, count_steps, timestep;

  double static_activation_energy;

  enum { OL_SPIN, SPIN_OL };

  bool set_act_energy,
      set_pressure_barrier; // whether the activation-energy and the
                            // pressure-barrier for the transformation will be
                            // given by the user
  double act_energy,
      pressure_barrier; // activation_energy and pressure barrier as set by user

  bool no_heat_flag;

  Heat_Lattice heat_distribution; // the boundary-condition

  Particle *nucleus, *help_pointer;

  bool include_elastic_energy;
  bool include_surface_energy;
  bool use_curvature;
  bool use_dynamic_activation_energy;
  bool include_negative_rates_into_elastic_relax;
  double allowed_error_factor;
  bool numeric_compensation;
  bool protect_flickering_particles;
  int secondary_timestep;
  bool plot_grains;
  bool use_homogeneous_spin;
  bool check_percolation;

  bool system_percolation;

  bool reset_pressure_scale;

  bool keep_pressure;

  int plot_grain_nb;

  double finite_strain;

  double overpressure;

  vector<double> grain_young;

  Min_Trans_Lattice();    // Constructor
  ~Min_Trans_Lattice(){}; // Destructor

  void Activate_MinTrans();

  // ------------------------------------------------------------------
  // usr functions for the run, reactions
  // ------------------------------------------------------------------

  void Reaction_One();

  void Reaction_Two();

  void Start_Reactions(double pres_bar);

  void Energy_Calculation(Particle *partikel);

  void Provisoric_Mineral_Transformation();

  void Mineral_Transformation();
  void Mineral_Transformation_No_Spin_Distri();

  void Make_Phase_Boundaries();

  void Heat_Flow(double);

  void Find_Local_Maxima();

  void Read_Heat_Lattice();

  void Nucleation();

  void Gbm();

  double Pressure();

  double Latent_Heat_Release(Particle *, int);

  void Set_Heat_Lattice(Particle *);

  double Adjust_Molare_Volume();

  double Molare_Volume(int mineral);

  void Invert_Provisoric_Mineral_Transformation();

  void Save_Heat_Lattice();

  void Invert_Heat_Lattice();

  void Change_Grain_Mineral(int, int, int);

  void Change_Timestep();

  void Exchange_Probabilities();

  void Call_Reaction_Two();

  double Activation_Energy(Particle *);

  double Growth_Rate_For_Nucleation();

  double Nucleation_Rate();

  double Normal_Stress(Particle *fixed_neig, int);

  void Change_Young();

  void Restore_Young();

  double Excess_Pressure();

  double Undercooling();

  double Energy_Barrier();

  double Driving_Force();

  void Check_Nucleus_Stability();

  void Update_Surf_Neighbour_List();

  double CalcTime();

  double CalcArea(double);

  void DumpDiffStressAndTime(Particle *, double, double);

  void SetGaussianYoungDistribution_2(double, double);

  bool Change_Particle();

  bool Change_Particle_Random(int nb_grains);

  bool Change_Particle_Distri();

  bool Change_Particle_Stick();

  void SetHeatFlowEnabled(int);

  void SetHeatFlowParameters(double, double, double, double);

  void HeatGrain(int, int);

  void SetHeatLatticeHeatFlowExample();

  void SetReactions(double, double);

  void InitProvisoricElasticProperties();

  void SetProvisoricElasticProperties();

  double clausius_clapeyron_threshold(double); //

  void SetTimestep(int);

  int Determine_grain_number(Particle *);

  void Insert_Hardbodies(vector<vector<double>> hardbodies);

  double growth_rate(double lowest_time, Particle *prtcl);

  void Provisoric_Mineral_Transformation_No_Spin_Distri();

  void CalcPressureScale(double bulkmod, double required_mean_stress);

  double Average_Temperature();

  double Calc_Curvature_Surface_Energy(Particle *prtcl, bool);

  double Calc_Particle_Surface_Energy(Particle *prtcl, Particle *lneig,
                                      Particle *rneig, double &surfE1,
                                      double &surfE2, double &surfE3,
                                      bool provisoric_transformation);

  bool Is_Same_Surface(Particle *prtcl1, Particle *prtcl2);

  bool Is_Same_Surface_Grain(Particle *prtcl1, Particle *prtcl2);

  void Plot_Grains();

  void Init_Grain_Young();

  bool Transform_System_Horizontal();

  bool Transform_System_Diagonal();

  void Check_Percolation();

  bool Transform_Grains(int, int, int);

  bool Transform_System_Vertical();

  bool Insert_Heterogeneous_Nuclei(double press_bar, int nb_nuclei,
                                   double radius, bool boundary_nucleation);

  void ResetPressureScale(double value, double t);

  void KeepPressure(double value, double t);

  void SetOverPressure(bool, double);

  void SetMolePerParticle();

  double OlivinePressure();

  void Check_Percolation_All();
};

#endif
