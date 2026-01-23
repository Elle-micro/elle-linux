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
#include "min_trans_lattice.h"
#include "phase_base.h"

class Min_Trans_Lattice : public Phase_Base {
public:
  // till-variablen
  double actual_time, now, time_interval, old_time;
  bool timeflag,         // multipurpose
      reaction_finished, // flag for heat-conduction and time-management in
                         // gbm()
      logicalflag, nucleation_occured;
  float mole_per_particle, //
      initial_diameter, highest_prob, lowest_prob;
  float pressure;

  int beware_grain1, beware_grain2, beware_grain3, count_steps;

  Heat_Lattice heat_distribution;

  Min_Trans_Lattice();    // Constructor
  ~Min_Trans_Lattice(){}; // Destructor

  // ------------------------------------------------------------------
  // usr functions for the run, reactions
  // ------------------------------------------------------------------

  void Reaction_One();

  void Reaction_Two();

  void Start_Reactions();

  void Prob_Calculation();

  void Energy_Calculation(Particle *partikel);

  void Provisoric_Mineral_Transformation();

  void Mineral_Transformation();

  void Make_Phase_Boundaries();

  void Heat_Flow();

  void Find_Local_Maxima();

  void Read_Heat_Lattice();

  void Nucleation();

  void Gbm();

  float Pressure();

  float Latent_Heat_Release();

  void Set_Heat_Lattice();

  float Adjust_Molare_Volume();

  void Invert_Provisoric_Mineral_Transformation();

  void Save_Heat_Lattice();

  void Invert_Heat_Lattice();

  float Rate_Constant_Nucleation();

  void Change_Grain_Mineral(int, int, int);

  void Change_Timestep();

  void Exchange_Probabilities();

  void Compare_Particles();

  void Check_Reaction();

  void Call_Reaction_Two();

  float Activation_Energy();

  float Growth_Rate();

  float Growth_Rate_For_Nucleation();

  double Nucleation_Rate();

  float Normal_Stress(Particle *fixed_neig);

  void Change_Young();

  void Restore_Young();

  float Normal_Stress_2();

  float Neig_Energy();

  float Excess_Pressure();

  float Undercooling();

  float Energy_Barrier();

  float Driving_Force();

  void Exchange_Probabilities_Nucleation();

  void Monte_Carlo();

  void Check_Nucleus_Stability();
};

#endif
