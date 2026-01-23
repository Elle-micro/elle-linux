/******************************************************
 * Spring Code Mike 2.0
 *
 * Functions for min_trans_lattice class in min_trans_lattice.cc
 *
 * Expansion of Spring Code for Reactions
 *
 *
 *
 *
 * Daniel Koehn and Jochen Arnold March 2003
 * Mainz
 *
 * Daniel Dec. 2003
 *
 * We thank Anders Malthe-Srenssen for his enormous help
 *
 * Daniel Koehn and Till Sachau 2004/5
 ******************************************************/

// srand(std::time(0));   Problem on the mac ! (line 201)

// runParticle->temperature = round(runParticle->temperature); Problem on Mac
// line 2074

// ------------------------------------
// system headers
// ------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <list>
#include <vector>

// ------------------------------------
// elle headers
// ------------------------------------

#include "lattice.h" // include lattice class
#include "min_trans.h"
#include "phase_base.h"
// inherits lattice class
#include "unodes.h" // include unode funct. plus undodesP.h
// for c++
#include "attrib.h"
#include "attribarray.h"
#include "attribute.h"
#include "check.h"
#include "convert.h"
#include "display.h"
#include "error.h"
#include "file.h"
#include "general.h"
#include "interface.h"
#include "nodes.h"
#include "nodesP.h"
#include "polygon.h"
#include "runopts.h"
#include "tripoly.h"
#include "update.h"

// CONSTRUCTOR
/*******************************************************
 * Constructor for the phase_lattice class
 * Does not do much at the moment, defines some variables
 * Phase lattice inherits the lattice class that does
 * most of the intialization work.
 *
 *
 * Daniel March 2003
 ********************************************************/

// ---------------------------------------------------------------
// Constructor of Min_Trans_Lattice class
// ---------------------------------------------------------------
Min_Trans_Lattice::Min_Trans_Lattice() {
  transition = true;

  // create the file for the rate
  ofstream rate;
  rate.open("Rate.txt", ios::out);
  rate << "Type ;\t Timestep ;\t percent transformed ;\t total area ;\t "
          "lowest_time ;\t nb spin partikel ;\t total spinarea ;\t new area "
          "(total) ;\t new_area (total) / lowest_time ;\t average rate from "
          "new area ;\t percentage of new area ;\t new_area / (the_size * "
          "lowest_time) ;\t percentage transformed area (partikel) ;\t "
          "partikel_area / (the_size * lowest_time) ;\t "
          "average_normal_stress_on_transformed_particle ;\t average "
          "dir_rate\n\n";
  rate.close();

  include_negative_rates_into_elastic_relax = false;

  // set as default: use the elastic energy
  include_elastic_energy = true;

  // default: use surface energy
  include_surface_energy = true;
  // whether particle by particle or rather 'average' over particle + 2
  // neighbours
  use_curvature = true;

  use_dynamic_activation_energy = true;

  numeric_compensation = false;

  // just some unrealistically large number
  beware_grain1 = beware_grain2 = beware_grain3 = 1000000;

  // the default, if the activation energy isn't dynamic or set otherwise
  static_activation_energy = 350000.0;

  // the default for the error which will be allowed by the program in order to
  // circumvent numerical imprecision (1 percent)
  allowed_error_factor = 0.01;

  // is increased by one every time gbm starts over
  secondary_timestep = 0;

  // not sure if this works, anyway
  protect_flickering_particles = false;

  // produce a file with grain-informations by default
  plot_grains = true;

  // default: let spinel grains have the same youngs modulus as their nucleus
  use_homogeneous_spin = true;

  // default: make a plot-grain file after every single transformation
  plot_grain_nb = 1;

  // default: don't plot the statisticstressbox by setting the finite strain,
  // passed to the function, to zero
  finite_strain = 0.0;

  check_percolation = false;

  system_percolation = false;

  reset_pressure_scale = false;

  keep_pressure = false;

  overpressure = 0.0;
}

void Min_Trans_Lattice::Init_Grain_Young() {

  int i, j;

  grain_young.clear();

  for (i = 0; i <= highest_grain; i++) {
    grain_young.push_back(0.0);
  }

  runParticle = &refParticle;
  for (i = 0; i < numParticles; i++) {
    grain_young[runParticle->grain] = runParticle->young;
    runParticle = runParticle->nextP;
  }
}

void Min_Trans_Lattice::Activate_MinTrans() {
  int i;
  // --------------------------------------------------------------------
  // Set some more starting values for particles
  // --------------------------------------------------------------------

  Activate_Lattice();

  for (i = 0; i < numParticles; i++) {
    if (runParticle->mineral == 1) {
      runParticle->surf_free_E = 0.6;

      // set the initial temperature and the molare volume
      runParticle->temperature = 1000.0;

      runParticle->mV = Adjust_Molare_Volume(); // check this !

      runParticle->previous_prob = 0;
      runParticle->spinel_content = 0;
      runParticle->previous_spinel_content = 0;

      runParticle->real_radius = the_size * runParticle->radius;

      runParticle->nucleus = false;

      runParticle->area = pow(double(runParticle->real_radius), 2.0) * 3.1415;

      runParticle->no_reaction = false;

      if (runParticle->xpos > 0.95 || runParticle->xpos < 0.05 ||
          runParticle->ypos > 0.95 || runParticle->ypos < 0.05)
        runParticle->no_nucleation = true;

      if (wrapping) {
        if (!runParticle->is_top_lattice_boundary &&
            !runParticle->is_bottom_lattice_boundary)
          runParticle->is_lattice_boundary = false;
      }
    }
    runParticle = runParticle->nextP;
  }

  // Mineral_Transformation();
  initial_diameter = 2 * runParticle->real_radius;

  // fuer die zeit-minimierung
  logicalflag = true;
  timeflag = false;

  // set initial parameters for heat-flow
  runParticle = &refParticle;       // start of list
  runParticle = runParticle->prevP; // go back one
  org_wall_pos_y = runParticle->ypos;

  runParticle = runParticle->prevP;
  org_wall_pos_x = runParticle->xpos;
  for (i = 0; i < particlex; i++) {
    runParticle = runParticle->prevP;
  }
  if (runParticle->xpos > org_wall_pos_x) {
    org_wall_pos_x = runParticle->xpos;
  }

  heat_distribution.heat_flow_enabled =
      false; // in standard settings heat-flow is turned off
  heat_distribution.numParticles = numParticles;

  old_time = time;

  no_heat_flag = false;

  // activation-energy set by user? default == no
  set_act_energy = false;
  act_energy = 400000.0;

  // pressure barrier ste by user ? if not->
  set_pressure_barrier = false;
  pressure_barrier = 12e9;

  // contains the youngs moduli of the grains, to assign it more easily to
  // transformed particles
  Init_Grain_Young();
}

// concept stolen from Update_Fluid_Neighbour_List()
void Min_Trans_Lattice::Update_Surf_Neighbour_List() {
  int i, j, count;
  double slope, slope_a, average_x, average_y;

  for (i = 0; i < numParticles; i++) {
    runParticle->fluid_P = 0.0;
    runParticle->rightNeighbour = NULL;
    runParticle->leftNeighbour = NULL;

    if (runParticle->is_boundary && !runParticle->is_lattice_boundary) {
      average_x = 0.0;
      average_y = 0.0;
      count = 0;
      for (j = 0; j < 6; j++) {
        if (runParticle->neigP[j]) {
          if (runParticle->neigP[j]->grain == runParticle->grain) {
            runParticle->rightNeighbour = runParticle->neigP[j];
            break;
          }
        }
      }
      for (j = 0; j < 6; j++) {
        if (runParticle->neigP[j]) {
          if (runParticle->neigP[j]->grain == runParticle->grain) {
            if (runParticle->neigP[j]->nb != runParticle->rightNeighbour->nb)
              runParticle->leftNeighbour = runParticle->neigP[j];
          }
        }
      }
      if (!runParticle->rightNeighbour) {
        // cout <<" right not set" << endl;   // if particle is totally alone
      } else if (!runParticle->leftNeighbour) {
        // cout << "left not set" << endl;    // can happen at edges, happens
        // quite often.
        for (j = 0; j < 6; j++) {
          if (runParticle->neigP[j]) {
            if (runParticle->neigP[j]->grain == runParticle->grain) {
              average_x = average_x + runParticle->neigP[j]->xpos;
              average_y = average_y + runParticle->neigP[j]->ypos;
              count = count + 1;
            }
          }
        }
        average_x = average_x / count;
        average_y = average_y / count;

        if (runParticle->xpos < runParticle->rightNeighbour->xpos) {
          if (average_y < runParticle->ypos) {
            runParticle->leftNeighbour = runParticle->rightNeighbour;
            runParticle->rightNeighbour = NULL;
          }
        } else if (runParticle->xpos > runParticle->rightNeighbour->xpos) {
          if (average_y > runParticle->ypos) {
            runParticle->leftNeighbour = runParticle->rightNeighbour;
            runParticle->rightNeighbour = NULL;
          }
        }
      }

      else if (runParticle->leftNeighbour->nb ==
               runParticle->rightNeighbour->nb) {
        // cout <<"left and right same !"<< endl;  //should be ok;
        ;
      } else {
        count = 0;
        average_x = 0;
        average_y = 0;
        for (j = 0; j < 8; j++) {
          if (runParticle->neigP[j]) {
            if (runParticle->neigP[j]->grain == runParticle->grain) {
              average_x = average_x + runParticle->neigP[j]->xpos;
              average_y = average_y + runParticle->neigP[j]->ypos;
              count = count + 1;
            }
          }
        }
        average_x = average_x / count;
        average_y = average_y / count;

        if (runParticle->rightNeighbour->xpos <
            (runParticle->leftNeighbour->xpos - (runParticle->radius))) {
          if (runParticle->rightNeighbour->ypos <
              (runParticle->leftNeighbour->ypos - (runParticle->radius))) {
            if (average_x < runParticle->xpos) {
              runParticle->neig = runParticle->leftNeighbour;
              runParticle->leftNeighbour = runParticle->rightNeighbour;
              runParticle->rightNeighbour = runParticle->neig;
            }
          } else if (runParticle->rightNeighbour->ypos >
                     (runParticle->leftNeighbour->ypos +
                      (runParticle->radius))) {
            if (runParticle->xpos < average_x) {
              runParticle->neig = runParticle->leftNeighbour;
              runParticle->leftNeighbour = runParticle->rightNeighbour;
              runParticle->rightNeighbour = runParticle->neig;
            }
          } else {
            if (runParticle->ypos < average_y) {
              // cout << "in" << endl;
              runParticle->neig = runParticle->leftNeighbour;
              runParticle->leftNeighbour = runParticle->rightNeighbour;
              runParticle->rightNeighbour = runParticle->neig;
            }
          }

        } else if (runParticle->rightNeighbour->xpos >
                   (runParticle->leftNeighbour->xpos + (runParticle->radius))) {
          if (runParticle->rightNeighbour->ypos <
              (runParticle->leftNeighbour->ypos - (runParticle->radius))) {
            if (average_x < runParticle->xpos) {

              runParticle->neig = runParticle->leftNeighbour;
              runParticle->leftNeighbour = runParticle->rightNeighbour;
              runParticle->rightNeighbour = runParticle->neig;
            }
          } else if (runParticle->rightNeighbour->ypos >
                     (runParticle->leftNeighbour->ypos +
                      (runParticle->radius))) {
            if (average_x > runParticle->xpos) {

              runParticle->neig = runParticle->leftNeighbour;
              runParticle->leftNeighbour = runParticle->rightNeighbour;
              runParticle->rightNeighbour = runParticle->neig;
            }
          } else // y the same
          {
            if (average_y < runParticle->ypos) {
              // cout << " y same " << endl;
              runParticle->neig = runParticle->leftNeighbour;
              runParticle->leftNeighbour = runParticle->rightNeighbour;
              runParticle->rightNeighbour = runParticle->neig;
            }
          }

        } else // x the same
        {
          // cout << "in x same" << endl;
          if (runParticle->rightNeighbour->ypos >
              runParticle->leftNeighbour->ypos) {
            if (average_x > runParticle->xpos) {
              runParticle->neig = runParticle->leftNeighbour;
              runParticle->leftNeighbour = runParticle->rightNeighbour;
              runParticle->rightNeighbour = runParticle->neig;
            }
          } else {
            if (average_x < runParticle->xpos) {
              runParticle->neig = runParticle->leftNeighbour;
              runParticle->leftNeighbour = runParticle->rightNeighbour;
              runParticle->rightNeighbour = runParticle->neig;
            }
          }
        }
      }
    }

    runParticle = runParticle->nextP;
  }
}

// "main" function for phase-transitions. is called in experiment.cc,
// coordinates the rest.
void Min_Trans_Lattice::Start_Reactions(double pres_bar) {

  int i, // for loop
      j;
  double interim_time, // stores time during gbm
      pressure, factor, av_temp;

  double wall_pos, wall_pos_x;

  cout << "mole_per_particle: " << mole_per_particle << endl;
  cout << "numParticles: " << numParticles << " " << endl;
  cout << "rate_constant: " << rate_constant << endl;

  pressure = Pressure();

  cout << "Pressure: " << pressure << endl;

  // secondary_timestep = 0;

  /************************************
   * preset variables:
   ************************************/

  runParticle = &refParticle; // for particle loops

  pois = 0.3333333; // poisson number for the triangular lattice

  Make_Phase_Boundaries();

  actual_time = time;

  for (i = 0; i < numParticles; i++) {
    runParticle->mV = Adjust_Molare_Volume();

    runParticle->real_radius = runParticle->radius * the_size;
    initial_diameter = runParticle->real_radius * 2;
    runParticle->area =
        runParticle->real_radius * runParticle->real_radius * Pi;

    runParticle->transformation_timestep = 0;

    runParticle = runParticle->nextP;
  }

  mole_per_particle =
      pow(refParticle.real_radius, 3.0) * 1.33333333 * 3.14 / refParticle.mV;
  SetMolePerParticle();

  time_interval = actual_time;

  nucleus = NULL;

  Relaxation();

  if (pressure > pres_bar) // start if barrier is oversteped
  {

    for (i = 0; i < numParticles; i++) {
      runParticle->spinel_content = 0.0;
      runParticle->rate = 0.0;

      runParticle->prob = 0.0;

      runParticle->neig3 = NULL;
      runParticle->right_now_transformed = false;
      runParticle->previously_transformed = false;

      runParticle = runParticle->nextP;
    }

    cout << "time interval:" << time_interval << endl;

    do {
      highest_prob = 0.0;
      srand(std::time(0));

      for (i = 0; i < numParticles; i++) {
        runParticle->phase_change_required = false;
        runParticle->previous_sigman = 0;

        runParticle->rate = 0.0;

        runParticle->prob = 0.0;
        runParticle->previous_deltaG = 0.0;
        runParticle->neig3 = NULL;

        runParticle->mV = Adjust_Molare_Volume();

        runParticle = runParticle->nextP;
      }

      Gbm();

      for (i = 0; i < numParticles; i++) {

        if (!runParticle->no_reaction)
          runParticle->phase = runParticle->mineral;

        runParticle = runParticle->nextP;
      }

      UpdateElle();

      av_temp = Average_Temperature();

      cout << "Average temperature: " << av_temp << endl;

      // plot
      if (plot_grains)
        Plot_Grains();

      // stress
      runParticle = &refParticle;       // start of list
      runParticle = runParticle->prevP; // go back one
      wall_pos =
          runParticle->ypos + runParticle->radius; // define wall position

      runParticle = &refParticle;       // start of list
      runParticle = runParticle->prevP; // go back one
      wall_pos_x = runParticle->xpos;   // define wall position
      for (i = 0; i < particlex; i++)
        runParticle = runParticle->prevP;
      if (runParticle->xpos > wall_pos)
        wall_pos_x = runParticle->xpos;

      if (finite_strain != 0.0)
        DumpStatisticStressBox(0.1, wall_pos - 0.1, 0.1, wall_pos_x - 0.1,
                               finite_strain);

      // percolation
      if (check_percolation)
        Check_Percolation();

      // keep pressure as it is, using the scaling parameter
      // the first value will be added to the clausius-clapeyron-threshold,
      // which will be calculated by the temperature overgiven in the second
      // value
      if (keep_pressure)
        KeepPressure(overpressure, av_temp);

    } while (highest_prob > 1.0 &&
             Pressure() > clausius_clapeyron_threshold(av_temp));
  }
}

double Min_Trans_Lattice::Average_Temperature() {

  int i, j, counter = 0;
  double temperature = 0.0;

  runParticle = &refParticle;

  for (i = 0; i < numParticles; i++) {

    temperature += runParticle->temperature;
    counter++;

    runParticle = runParticle->nextP;
  }

  return (temperature / double(counter));
}

double Min_Trans_Lattice::clausius_clapeyron_threshold(double temperature) {

  double critical_pressure;
  double P_0 = 1.08e10;
  double slope = 2436708.86;

  critical_pressure = P_0 + slope * temperature;

  return critical_pressure;
}

// if nb_nuclei == 0, the complete list of possible nuclei will be transformed
bool Min_Trans_Lattice::Insert_Heterogeneous_Nuclei(double press_bar,
                                                    int nb_nuclei,
                                                    double radius,
                                                    bool boundary_nucleation) {

  int k, i, j, index;
  double x_dist, y_dist;
  bool lattice;
  Particle *ref_particle;
  vector<Particle *> list;
  vector<int> todelete;

  mole_per_particle =
      pow(refParticle.real_radius, 3.0) * 1.33333333 * 3.14 / refParticle.mV;
  SetMolePerParticle();
  cout << "mole per particle: " << mole_per_particle << endl;

  runParticle = &refParticle;

  for (i = 0; i < numParticles; i++) {

    if (!runParticle->no_nucleation)
      runParticle->local_max = true;
    else
      runParticle->local_max = false;

    runParticle->prob = 0.0;

    runParticle = runParticle->nextP;
  }

  for (i = 0; i < numParticles; i++) {

    runParticle->average_pressure =
        -(runParticle->sxx + runParticle->syy) / 2.0;

    if (runParticle->no_nucleation)
      runParticle->average_pressure = 0.0;

    runParticle = runParticle->nextP;
  }

  for (i = 0; i < numParticles; i++) {

    if (!runParticle->no_nucleation) {

      ref_particle = &refParticle;

      for (j = 0; j < numParticles; j++) {

        if (!runParticle->is_lattice_boundary &&
            !runParticle->is_left_lattice_boundary &&
            !runParticle->is_left_lattice_boundary) {

          x_dist = runParticle->xpos - ref_particle->xpos;
          y_dist = runParticle->ypos - ref_particle->ypos;

          if (sqrt(y_dist * y_dist + x_dist * x_dist) < radius) {

            if (runParticle->average_pressure <
                ref_particle->average_pressure) {

              runParticle->local_max = false;
              break;
            }
          }
        }

        ref_particle = ref_particle->nextP;
      }
    }

    runParticle = runParticle->nextP;
  }

  list.clear();

  for (i = 0; i < numParticles; i++) {

    if (runParticle->local_max) {

      // dont want lattice-boundaries
      if (!runParticle->is_lattice_boundary) {

        lattice = false;

        for (j = 0; j < 6; j++) {
          if (runParticle->neigP[j]) {
            if (runParticle->neigP[j]->is_lattice_boundary) {
              lattice = true;
              break;
            }
          }
        }
        if (!lattice)
          list.push_back(runParticle);
      }
    }

    runParticle = runParticle->nextP;
  }

  cout << "max nb of particles: " << list.size() << endl;

  // sort list!

  for (k = 0; k < list.size(); k++) {
    index = k;
    ref_particle = list[k];

    for (j = 0; j < list.size() - k; j++) {
      if ((-(list[k + j]->sxx + list[k + j]->syy) / 2.0) >
          (-(ref_particle->sxx + ref_particle->syy) / 2.0)) {
        ref_particle = list[k + j];
        index = k + j;
      }
    }
    list[index] = list[k];
    list[k] = ref_particle;
  }

  cout << "max nb of particles: " << list.size() << endl;

  k = 0;
  for (i = 0; i < list.size(); i++) {
    if ((list[i]->average_pressure * pressure_scale * pascal_scale) <
        press_bar) {
      k = i;
      break;
    }
  }

  cout << "max nb of particles: " << list.size() << endl;

  for (i = list.size() - k; i < list.size(); i++) {
    list.pop_back();
  }

  if (nb_nuclei == 0)
    nb_nuclei = list.size() - 1;

  if (nb_nuclei < list.size()) {

    cout << "max nb of particles: " << list.size() << endl;

    if (boundary_nucleation) {
      for (i = 0; i < nb_nuclei; i++) {
        if (list[i]->is_boundary) {
          list[i]->prob = 2.0;
          list[i]->spin_young = runParticle->young * 1.4;
          list[i]->newly_overgrown = runParticle->area;
          list[i]->temperature += Latent_Heat_Release(runParticle, OL_SPIN);
          list[i]->no_reaction = true;
        }
      }
    } else {
      for (i = 0; i < nb_nuclei; i++) {
        list[i]->prob = 2.0;
        list[i]->spin_young = runParticle->young * 1.4;
        list[i]->newly_overgrown = runParticle->area;
        list[i]->temperature += Latent_Heat_Release(runParticle, OL_SPIN);
        list[i]->no_reaction = true;
      }
    }

    UpdateElle();

    no_heat_flag = true;
    Mineral_Transformation_No_Spin_Distri();
    no_heat_flag = false;

    // adjust changed springs on grain boundaries
    AdjustConstantGrainBoundaries();
    Make_Phase_Boundaries();
    // and relax
    Relaxation();
    // Exchange_Probabilities();

    return true;

  } else {

    return false;
  }
}

bool Min_Trans_Lattice::Transform_System_Vertical() {

  int i;

  mole_per_particle =
      pow(refParticle.real_radius, 3.0) * 1.33333333 * 3.14 / refParticle.mV;
  SetMolePerParticle();

  runParticle = &refParticle;

  for (i = 0; i < numParticles; i++) {

    if (runParticle->xpos < 0.5) {

      if (!runParticle->is_lattice_boundary) {

        runParticle->prob = 2.0;
        runParticle->spin_young = runParticle->young * 1.4;
        runParticle->newly_overgrown = runParticle->area;
        runParticle->temperature += Latent_Heat_Release(runParticle, OL_SPIN);
      }
    }

    runParticle = runParticle->nextP;
  }

  no_heat_flag = true;
  // Mineral_Transformation();
  Mineral_Transformation_No_Spin_Distri();
  no_heat_flag = false;

  // adjust changed springs on grain boundaries
  AdjustConstantGrainBoundaries();
  Make_Phase_Boundaries();
  // and relax
  Relaxation();
  // Exchange_Probabilities()

  UpdateElle();

  return (true);
}

bool Min_Trans_Lattice::Transform_Grains(int g1, int g2, int g3) {

  int i;

  mole_per_particle =
      pow(refParticle.real_radius, 3.0) * 1.33333333 * 3.14 / refParticle.mV;
  SetMolePerParticle();

  runParticle = &refParticle;

  for (i = 0; i < numParticles; i++) {

    if (runParticle->grain == g1 || runParticle->grain == g2 ||
        runParticle->grain == g3) {

      if (!runParticle->is_lattice_boundary) {

        runParticle->prob = 2.0;
        runParticle->spin_young = runParticle->young * 1.4;
        runParticle->newly_overgrown = runParticle->area;
        runParticle->temperature += Latent_Heat_Release(runParticle, OL_SPIN);
      }
    }

    runParticle = runParticle->nextP;
  }

  no_heat_flag = true;
  // Mineral_Transformation();
  Mineral_Transformation_No_Spin_Distri();
  no_heat_flag = false;

  // adjust changed springs on grain boundaries
  AdjustConstantGrainBoundaries();
  Make_Phase_Boundaries();
  // and relax
  Relaxation();
  // Exchange_Probabilities()

  UpdateElle();

  return (true);
}

bool Min_Trans_Lattice::Transform_System_Diagonal() {

  int i;

  mole_per_particle =
      pow(refParticle.real_radius, 3.0) * 1.33333333 * 3.14 / refParticle.mV;
  SetMolePerParticle();

  runParticle = &refParticle;

  for (i = 0; i < numParticles; i++) {

    if (runParticle->ypos > runParticle->xpos) {

      if (!runParticle->is_lattice_boundary) {

        runParticle->prob = 2.0;
        runParticle->spin_young = runParticle->young * 1.4;
        runParticle->newly_overgrown = runParticle->area;
        runParticle->temperature += Latent_Heat_Release(runParticle, OL_SPIN);
      }
    }

    runParticle = runParticle->nextP;
  }

  no_heat_flag = true;
  // Mineral_Transformation();
  Mineral_Transformation_No_Spin_Distri();
  no_heat_flag = false;

  // adjust changed springs on grain boundaries
  AdjustConstantGrainBoundaries();
  Make_Phase_Boundaries();
  // and relax
  Relaxation();
  // Exchange_Probabilities()

  UpdateElle();

  return (true);
}

bool Min_Trans_Lattice::Transform_System_Horizontal() {

  int i;

  mole_per_particle =
      pow(refParticle.real_radius, 3.0) * 1.33333333 * 3.14 / refParticle.mV;
  SetMolePerParticle();

  runParticle = &refParticle;

  for (i = 0; i < numParticles; i++) {

    if (runParticle->ypos < 0.5) {

      if (!runParticle->is_lattice_boundary) {

        runParticle->prob = 2.0;
        runParticle->spin_young = runParticle->young * 1.4;
        runParticle->newly_overgrown = runParticle->area;
        runParticle->temperature += Latent_Heat_Release(runParticle, OL_SPIN);
      }
    }

    runParticle = runParticle->nextP;
  }

  no_heat_flag = true;
  // Mineral_Transformation();
  Mineral_Transformation_No_Spin_Distri();
  no_heat_flag = false;

  // adjust changed springs on grain boundaries
  AdjustConstantGrainBoundaries();
  Make_Phase_Boundaries();
  // and relax
  Relaxation();
  // Exchange_Probabilities()

  UpdateElle();

  return (true);
}

bool Min_Trans_Lattice::Change_Particle() {
  int i, j, one, two, row, column1, column2;

  column1 = int(sqrt(numParticles));
  row = int(column1 / 2);

  column1 = int(column1 / 3);
  column2 = column1 * 2;

  one = (row - 1) * int(sqrt(numParticles)) + column1;
  two = (row + 1) * int(sqrt(numParticles)) + column2;

  mole_per_particle =
      pow(refParticle.real_radius, 3.0) * 1.33333333 * 3.14 / refParticle.mV;
  SetMolePerParticle();

  for (i = 0; i < numParticles; i++) {
    if (runParticle->nb ==
        23500 /*|| runParticle->nb == 23501 || runParticle->nb == 5057*/) {
      if (runParticle->mineral == 1) {
        runParticle->prob = 2.0;
        runParticle->spin_young = runParticle->young * 1.4;
        runParticle->newly_overgrown = runParticle->area;
        runParticle->temperature += Latent_Heat_Release(runParticle, OL_SPIN);

        //		        for (j=0;j<6;j=j++){
        //		            runParticle->neigP[j]->prob = 2.0;
        //		            runParticle->neigP[j]->spin_young =
        //runParticle->young*1.4; 		            runParticle->neigP[j]->newly_overgrown =
        //runParticle->area; 		            runParticle->temperature +=
        //Latent_Heat_Release(runParticle, OL_SPIN);
        //		        }
      }
    }

    runParticle = runParticle->nextP;
  }

  no_heat_flag = true;
  // Mineral_Transformation();
  Mineral_Transformation_No_Spin_Distri();
  no_heat_flag = false;

  // adjust changed springs on grain boundaries
  AdjustConstantGrainBoundaries();
  Make_Phase_Boundaries();
  // and relax
  Relaxation();
  // Exchange_Probabilities()

  UpdateElle();

  return (true);
}

// insert "nuclei", from which the grainboundary migration will start
bool Min_Trans_Lattice::Change_Particle_Stick() {
  int i, j, one, two, row, column1, column2;

  column1 = int(sqrt(numParticles));
  row = int(column1 / 2);

  column1 = int(column1 / 3);
  column2 = column1 * 2;

  one = (row - 1) * int(sqrt(numParticles)) + column1;
  two = (row + 1) * int(sqrt(numParticles)) + column2;

  mole_per_particle =
      pow(refParticle.real_radius, 3.0) * 1.33333333 * 3.14 / refParticle.mV;
  SetMolePerParticle();

  for (i = 0; i < numParticles; i++) {
    if (runParticle->nb == 1325 || runParticle->nb == 1324 ||
        runParticle->nb == 1326 || runParticle->nb == 1323 ||
        runParticle->nb == 1327) {
      if (runParticle->mineral == 1) {
        runParticle->prob = 2.0;
        runParticle->spin_young = runParticle->young * 1.4;
        runParticle->newly_overgrown = runParticle->area;
        runParticle->temperature += Latent_Heat_Release(runParticle, OL_SPIN);
      }
    }

    runParticle = runParticle->nextP;
  }

  no_heat_flag = true;
  // Mineral_Transformation();
  Mineral_Transformation_No_Spin_Distri();
  no_heat_flag = false;

  // adjust changed springs on grain boundaries
  AdjustConstantGrainBoundaries();
  Make_Phase_Boundaries();
  // and relax
  Relaxation();
  // Exchange_Probabilities()

  UpdateElle();

  return (true);
}

bool Min_Trans_Lattice::Change_Particle_Distri() {
  int i, one, two, row, column1, column2, j;

  column1 = int(sqrt(numParticles));
  row = int(column1 / 2);

  column1 = int(column1 / 3);
  column2 = column1 * 2;

  one = (row - 1) * int(sqrt(numParticles)) + column1;
  two = (row + 1) * int(sqrt(numParticles)) + column2;

  mole_per_particle =
      pow(refParticle.real_radius, 3.0) * 1.33333333 * 3.14 / refParticle.mV;
  SetMolePerParticle();

  for (i = 0; i < numParticles; i++) {
    if (runParticle->nb == 5050 || runParticle->nb == 4545 ||
        runParticle->nb == 4555) {
      if (runParticle->mineral == 1) {
        runParticle->prob = 2.0;
        runParticle->spin_young = runParticle->young * 1.4;
        runParticle->newly_overgrown = runParticle->area;
      }
    }
    runParticle = runParticle->nextP;
  }

  no_heat_flag = true;
  // Mineral_Transformation();
  Mineral_Transformation_No_Spin_Distri();
  no_heat_flag = false;

  // adjust changed springs on grain boundaries
  AdjustConstantGrainBoundaries();
  Make_Phase_Boundaries();
  // and relax
  Relaxation();
  // Exchange_Probabilities();
  UpdateElle();

  return (true);
}

// insert "nuclei", from which the grainboundary migration will start
bool Min_Trans_Lattice::Change_Particle_Random(int nb_grains) {
  int i, j, one, two, row, column1, column2, rd_nb, k;

  column1 = int(sqrt(numParticles));
  row = int(column1 / 2);

  column1 = int(column1 / 3);
  column2 = column1 * 2;

  one = (row - 1) * int(sqrt(numParticles)) + column1;
  two = (row + 1) * int(sqrt(numParticles)) + column2;

  srand(std::time(0));

  mole_per_particle =
      pow(refParticle.real_radius, 3.0) * 1.33333333 * 3.14 / refParticle.mV;
  SetMolePerParticle();

  for (j = 0; j < nb_grains; j++) {

    rd_nb = rand() % numParticles;

    for (i = 0; i < numParticles; i++) {
      if (runParticle->nb == rd_nb) {
        if (runParticle->xpos > 0.05 && runParticle->xpos < 0.95 &&
            runParticle->ypos > 0.01 && runParticle->ypos < 0.8) {

          if (runParticle->mineral == 1 && !runParticle->is_lattice_boundary &&
              !runParticle->no_nucleation) {
            runParticle->prob = 2.0;
            runParticle->spin_young = runParticle->young * 1.4;
            runParticle->newly_overgrown = runParticle->area;
            runParticle->no_reaction = true;
          }
          //		           for (k=0;k<6;k++){
          //		               runParticle->neigP[k]->prob = 2.0;
          //		               runParticle->neigP[k]->spin_young =
          //runParticle->young*1.4; 		               runParticle->neigP[k]->newly_overgrown =
          //runParticle->area; 		               runParticle->temperature +=
          //Latent_Heat_Release(runParticle, OL_SPIN);
          //		           }
        }
      }
      runParticle = runParticle->nextP;
    }
  }

  no_heat_flag = true;
  // Mineral_Transformation();
  Mineral_Transformation_No_Spin_Distri();
  no_heat_flag = false;

  // adjust changed springs on grain boundaries
  AdjustConstantGrainBoundaries();
  Make_Phase_Boundaries();
  // and relax
  Relaxation();
  // Exchange_Probabilities();
  UpdateElle();

  return (true);
}

// calculates the energy difference for boundary-particles due to a phase-change
void Min_Trans_Lattice::Reaction_Two() {
  int i, j;

  // neue variablen im wesentlichen fr rate-law / normal stress
  double merker1, // multipurpose
      merker2,    // dito
      numerator, denominator, molareVolume, delta_Ed,
      delta_Ev, // delta elastic energy
      sigman,
      deltav, // volumendifferenz, molare
      inverse_mean_stress,
      time2, // fr zeitzwischenspeicherung
      deltaV;

  double surfE1, surfE2,
      surfE3; // for the surface-energy. surfeE2/3 are mainly used as dummies

  // variablen fr boundary-steigung
  bool neigh_found; // fr abfrage in prob-calculation

  double rate_constant, rate, rate2, exponent, N, k, R;

  double delta_Ev_pos = 0.0, delta_Ev_neg = 0.0, delta_Ed_pos = 0.0,
         delta_Ed_neg = 0.0, delta_surf_pos = 0.0, delta_surf_neg = 0.0,
         deltaG_pos = 0.0, deltaG_neg = 0.0;
  // if (runParticle->neig3) runParticle->time = runParticle->neig3->time;

  double eEl_dif_pos, eEl_dif_neg, surfE_dif_pos, surfE_dif_neg;

  N = 6e23;
  k = 1.5 * 1e-23;
  R = N * k;

  Particle *helper;

  runParticle = &refParticle;

  /******************************
   * surface energy difference
   ******************************/

  if (include_surface_energy) {

    runParticle = &refParticle;

    // first calculate the current state
    for (i = 0; i < numParticles; i++) {
      if (runParticle->is_phase_boundary) {
        if (!runParticle->is_lattice_boundary) {
          if (use_curvature)
            Calc_Curvature_Surface_Energy(runParticle, false);
          else
            Calc_Particle_Surface_Energy(runParticle, NULL, NULL, surfE1,
                                         surfE2, surfE3, false);
        }

        runParticle->surfE = surfE1;
      } else
        runParticle->surfE = 0.0;

      runParticle = runParticle->nextP;
    }

    runParticle = &refParticle;

    // now calculate the provisoric state
    for (i = 0; i < numParticles; i++) {
      if (runParticle->is_phase_boundary) {
        if (!runParticle->is_lattice_boundary) {
          if (use_curvature)
            Calc_Curvature_Surface_Energy(runParticle, true);
          else
            Calc_Particle_Surface_Energy(runParticle, NULL, NULL, surfE1,
                                         surfE2, surfE3, true);
        }

        runParticle->delta_surfE = surfE1 - runParticle->surfE;
      } else
        runParticle->delta_surfE = 0.0;

      runParticle = runParticle->nextP;
    }
  }

  /******************************
   * normal stresses
   ******************************/

  for (i = 0; i < numParticles; i++) {
    if (runParticle->phase_change_required) {

      // Calc_Curvature_New_Surf_Energy(runParticle);

      for (j = 0; j < 6; j++) {
        runParticle->dir_rate[j] = 0.0;

        if (runParticle->neigP[j]) {
          if (runParticle->mineral != runParticle->neigP[j]->mineral) {
            runParticle->normal_stress[j] =
                Normal_Stress(runParticle->neigP[j], j);
          } else
            runParticle->normal_stress[j] = 0.0;
        } else
          runParticle->normal_stress[j] = 0.0;

        if (isnan(runParticle->normal_stress[j]))
          cout << "isnogud" << endl;
      }
    }

    runParticle = runParticle->nextP;
  }

  /******************************
   * calculate driving force
   ******************************/

  for (i = 0; i < numParticles; i++) {
    if (runParticle->phase_change_required) {

      for (j = 0; j < 6; j++) {
        if (runParticle->neigP[j]) {
          if (runParticle->normal_stress[j] != 0.0) {

            sigman =
                0.5 * (-runParticle->sxx - runParticle->syy) +
                0.5 * (-runParticle->sxx + runParticle->syy) *
                    cos(2.0 * runParticle->boundary_angle[j]) -
                runParticle->sxy * sin(2.0 * runParticle->boundary_angle[j]);

            sigman *= pressure_scale * pascal_scale;

            runParticle->normal_stress[j] =
                clausius_clapeyron_threshold(runParticle->temperature) -
                runParticle->normal_stress[j];
          }
        }
        if (isnan(runParticle->normal_stress[j]))
          cout << "isnogud" << endl;
      }
    }

    runParticle = runParticle->nextP;
  }

  /******************************
   * initial and secondary elastic energies
   ******************************/

  if (include_elastic_energy) {

    for (i = 0; i < numParticles; i++) {
      Energy_Calculation(runParticle);

      runParticle->previous_eEl = runParticle->eEl;
      runParticle->previous_eEl_dilatation = runParticle->eEl_dilatation;
      runParticle->previous_eEl_distortion = runParticle->eEl_distortion;

      runParticle->previous_mV = runParticle->mV;

      runParticle = runParticle->nextP;
    }

    Provisoric_Mineral_Transformation_No_Spin_Distri();

    // save current values, avoid a second relaxation!
    for (i = 0; i < numParticles; i++) {
      // instead of relaxation 1
      runParticle->prev_xpos = runParticle->xpos;
      runParticle->prev_ypos = runParticle->ypos;
      runParticle->prev_sxx = runParticle->sxx;
      // instead of relaxation 2
      runParticle->prev_syy = runParticle->syy;
      runParticle->prev_sxy = runParticle->sxy;

      runParticle = runParticle->nextP;
    }

    Relaxation();

    for (i = 0; i < numParticles; i++) {
      Energy_Calculation(runParticle);

      runParticle = runParticle->nextP;
    }
  }

  for (i = 0; i < numParticles; i++) {
    if (runParticle->phase_change_required) {

      if (include_elastic_energy) {

        deltaV = runParticle->mV - runParticle->previous_mV;

        runParticle->eEl_dif = (runParticle->eEl - runParticle->previous_eEl);

        delta_Ev =
            runParticle->eEl_dilatation - runParticle->previous_eEl_dilatation;
        delta_Ed =
            runParticle->eEl_distortion - runParticle->previous_eEl_distortion;

        // energy stored in the spherical deformation
        if (delta_Ev != 0.0) {
          if (delta_Ev > 0.0) {
            delta_Ev_pos = delta_Ev;
            delta_Ev_neg = 0.0;
          } else {
            delta_Ev_neg = delta_Ev;
            delta_Ev_pos = 0.0;
          }
        } else {
          delta_Ev_neg = 0.0;
          delta_Ev_pos = 0.0;
        }

        // energy of distortion
        if (delta_Ed != 0.0) {
          if (delta_Ed > 0.0) {
            delta_Ed_pos = delta_Ev;
            delta_Ed_neg = 0.0;
          } else {
            delta_Ed_neg = delta_Ev;
            delta_Ed_pos = 0.0;
          }
        } else {
          delta_Ed_neg = 0.0;
          delta_Ed_pos = 0.0;
        }

        // this should be the sum of the two energies above
        if (runParticle->eEl_dif != 0.0) {
          if (runParticle->eEl_dif > 0.0) {
            delta_Ed_pos = runParticle->eEl_dif;
            delta_Ed_neg = 0.0;
          } else {
            delta_Ed_pos = 0.0;
            delta_Ed_neg = runParticle->eEl_dif;
          }
        } else {

          delta_Ed_pos = 0.0;
          delta_Ed_neg = 0.0;
        }

      } else {

        delta_Ed_pos = 0.0;
        delta_Ed_neg = 0.0;

        if (runParticle->mineral == 1)
          deltaV = Molare_Volume(2) - runParticle->mV;
        else
          deltaV = Molare_Volume(1) - runParticle->mV;
      }

      if (include_surface_energy) {

        if (runParticle->delta_surfE != 0.0) {

          if (runParticle->delta_surfE > 0.0) {
            surfE_dif_pos = runParticle->delta_surfE;
            surfE_dif_neg = 0.0;
          } else {
            surfE_dif_pos = 0.0;
            surfE_dif_neg = runParticle->delta_surfE;
          }

        } else {

          surfE_dif_pos = surfE_dif_neg = 0.0;
        }

      } else {

        surfE_dif_pos = surfE_dif_neg = 0.0;
      }

      for (j = 0; j < 6; j++) {
        deltaG_pos = -(runParticle->normal_stress[j]) * deltaV * 2;
        deltaG_neg = 0.0;

        if (deltaG_pos < 0.0) {
          deltaG_neg = deltaG_pos;
          deltaG_pos = 0.0;
        }

        if (runParticle->normal_stress[j] != 0.0) {
          runParticle->dir_rate[j] =
              /*162000.0*/ 1.58e12 * runParticle->temperature *
              (exp(-(Activation_Energy(runParticle) / N + deltaG_neg / N +
                     surfE_dif_neg / N + delta_Ed_neg / N + delta_Ev_neg / N) /
                   (k * runParticle->temperature)) -
               exp(-(Activation_Energy(runParticle) / N - deltaG_pos / N -
                     surfE_dif_pos / N + delta_Ed_pos / N + delta_Ev_pos / N) /
                   (k * runParticle->temperature)));
        }
      }
    }

    runParticle = runParticle->nextP;
  }

  /******************************
   * redo the provisoric mineral transformation that was used to calculated
   *delta E_d/v
   ******************************/

  if (include_elastic_energy) {

    Invert_Provisoric_Mineral_Transformation();

    // instead of a second relaxation, use the saved values
    for (i = 0; i < numParticles; i++) {

      runParticle->xpos = runParticle->prev_xpos;
      runParticle->ypos = runParticle->prev_ypos;
      runParticle->sxx = runParticle->prev_sxx;
      runParticle->syy = runParticle->prev_syy;
      runParticle->sxy = runParticle->prev_sxy;

      runParticle->phase_change_required = false;

      runParticle = runParticle->nextP;
    }
  }
}

double Min_Trans_Lattice::Normal_Stress(Particle *fixed_neig, int neignb) {

  double merker2, merker1, mean_stress, r, sigman, grain_angle, boundary_angle,
      phi, dd_neigh,
      dd_neigh_max, // zur Berechnung der Entfernung zweier Nachbarpartikeln
      neigh1_x,
      neigh1_y, // speichert pos des vorgger-neighbours
      neigh2_x, neigh2_y, boundary_slope, average_normal_stress;
  double x1, y1, norm;
  Particle *neig1, *neig2;
  bool case2 = true;
  int i, j,
      numNeig; // number of neighbours of other phase
  bool done, ok;
  vector<Particle *> list, list1, list2;

  list.clear();
  list1.clear();
  list2.clear();

  list.push_back(fixed_neig);

  //	for (i=0; i<6; i++)
  //	{
  //		if (fixed_neig->neigP[i])
  //		{
  //			if (fixed_neig->neigP[i]->is_phase_boundary &&
  //fixed_neig->neigP[i]->mineral != runParticle->mineral)
  //			{
  //				for (j=0; j<6; j++)
  //				{
  //					if (fixed_neig->neigP[i]->neigP[j] ==
  //runParticle)
  //					{
  //						list1.push_back(fixed_neig->neigP[i]);
  //					}
  //				}
  //			}
  //		}
  //	}

  //	for (i=0; i<6; i++)
  //	{
  //		if (fixed_neig->neigP[i])
  //		{
  //			if (fixed_neig->neigP[i]->is_phase_boundary &&
  //fixed_neig->neigP[i]->mineral == runParticle->mineral)
  //			{
  //				for (j=0; j<6; j++)
  //				{
  //					if (fixed_neig->neigP[i]->neigP[j] ==
  //runParticle)
  //					{
  //						list2.push_back(fixed_neig->neigP[i]);
  //					}
  //				}
  //			}
  //		}
  //	}

  //	if (list1.size() > 0 || list2.size() > 0) {
  //		if (list2.size() > list1.size()) {
  //			for (i=0; i<list2.size(); i++) {
  //				list.push_back(list2[i]);
  //			}
  //		} else {
  //			for (i=0; i<list1.size(); i++) {
  //				list.push_back(list1[i]);
  //			}
  //		}
  //	}

  for (i = 0; i < 6; i++) {
    if (fixed_neig->neigP[i]) {
      if (fixed_neig->neigP[i]->is_phase_boundary) {
        for (j = 0; j < 6; j++) {
          if (fixed_neig->neigP[i]->neigP[j] == runParticle) {
            if (fixed_neig->neigP[i]->mineral != runParticle->mineral)
              list1.push_back(fixed_neig->neigP[i]);
            // else
            // list2.push_back(fixed_neig->neigP[i]);
          }
        }
      }
    }
  }

  //	if (list1.size() > 1) {
  for (i = 0; i < list1.size(); i++) {
    list.push_back(list1[i]);
  }
  //	} else if (list2.size() > 1) {
  //		for (i=0; i<list2.size(); i++) {
  //			list.push_back(list2[i]);
  //		}
  //	} else {
  //		for (i=0; i<list2.size(); i++) {
  //			list.push_back(list2[i]);
  //		}
  //	}

  average_normal_stress = 0.0;
  j = 0;
  boundary_angle = 0;
  x1 = y1 = 0;
  for (i = 0; i < list.size(); i++) {

    j++;

    y1 += list[i]->ypos;
    x1 += list[i]->xpos;
  }

  x1 /= list.size();
  y1 /= list.size();

  x1 -= runParticle->xpos;
  y1 -= runParticle->ypos;

  //	x1 = list[0]->xpos - runParticle->xpos;
  //	y1 = list[0]->ypos - runParticle->ypos;

  norm = sqrt(x1 * x1 + y1 * y1);
  boundary_angle = acos(x1 / norm);

  if (y1 < 0)
    boundary_angle = (2.0 * Pi - boundary_angle);

  average_normal_stress =
      -runParticle->sxx * cos(boundary_angle) * cos(boundary_angle) -
      runParticle->syy * sin(boundary_angle) * sin(boundary_angle) -
      runParticle->sxy * sin(2.0 * boundary_angle);

  // average_normal_stress /= double(j);

  average_normal_stress *= pressure_scale * pascal_scale;

  return (average_normal_stress);
}

// function calculates normal stress on grain boundaries, where it uses
// neighbouring gb-particles to estimate the true grain-boundary direction
// double Min_Trans_Lattice::Normal_Stress(Particle * fixed_neig,int neignb)
//{
//  double      merker2,
//  merker1,
//  mean_stress,
//  r,
//  sigman,
//  grain_angle,
//  boundary_angle,
//  phi,
//  dd_neigh,
//  dd_neigh_max, // zur Berechnung der Entfernung zweier Nachbarpartikeln
//  neigh1_x,
//  neigh1_y,        // speichert pos des vorg�ger-neighbours
//  neigh2_x,
//  neigh2_y,
//  boundary_slope;

//  Particle   *neig1,
//  *neig2;

//  int i,
//  j,
//  numNeig; //number of neighbours of other phase

//  bool done,ok;

//  //-------------------------------------------------------------------------------------------
//  //find out about number of neig-particles of other phase, to distinguish
//  types of boundary
//  //-------------------------------------------------------------------------------------------
//  numNeig = 0;
//  neig1 = NULL;
//  neig2 = NULL;
//  done = false;

//  for (i=0; i<6; i++)
//    {
//      if (fixed_neig->neigP[i])
//        {
//          if (fixed_neig->neigP[i]->is_phase_boundary &&
//          fixed_neig->neigP[i]->grain == runParticle->grain)
//            {
//              if (fixed_neig->neigP[i] != runParticle)
//                {
//                  for (j=0; j<6; j++)
//                    {
//                      if (fixed_neig->neigP[i]->neigP[j] == runParticle)
//                        {
//                          if (!done)
//                            {
//                              neig1 = fixed_neig->neigP[i];
//                              done = true;
//                            }
//                          else
//                            {
//                              neig2 = fixed_neig->neigP[i];
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }

//  if (!neig2)
//    {
//      done = false;
//      for (i=0; i<6; i++)
//        {
//          if (fixed_neig->neigP[i])
//            {
//              if (fixed_neig->neigP[i]->is_phase_boundary &&
//              fixed_neig->neigP[i]->grain != runParticle->grain)
//                {
//                  for (j=0; j<6; j++)
//                    {
//                      if (fixed_neig->neigP[i]->neigP[j] == runParticle)
//                        {
//                          if (!done)
//                            {
//                              neig1 = fixed_neig->neigP[i];
//                              done = true;
//                            }
//                          else
//                            {
//                              neig2 = fixed_neig->neigP[i];
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }

//  if (!neig2)
//    {
//      done = false;
//      //find connected neighbours
//      for (i=0; i<6; i++)
//        {
//          if (fixed_neig->neigP[i])
//            {
//              if (fixed_neig->neigP[i] != runParticle)
//                {
//                  if (fixed_neig->neigP[i]->is_phase_boundary &&
//                  fixed_neig->neigP[i]->grain == runParticle->grain)
//                    {
//                      neig1 = fixed_neig->neigP[i];
//                      done = true;
//                      neig2 = runParticle;
//                    }
//                }
//            }
//        }

//      if (!done)
//        {
//          for (i=0; i<6; i++)
//            {
//              if (fixed_neig->neigP[i])
//                {
//                  if (fixed_neig->neigP[i]->is_phase_boundary &&
//                  fixed_neig->neigP[i]->grain == fixed_neig->grain)
//                    {
//                      neig1 = fixed_neig->neigP[i];
//                      done = true;
//                      neig2 = fixed_neig;
//                    }
//                }
//            }
//        }
//    }

//  if (neig2)
//    {
//      neigh1_x = neig1->xpos;
//      neigh2_x = neig2->xpos;
//      neigh1_y = neig1->ypos;
//      neigh2_y = neig2->ypos;

//      //-------------------------------------------------------------------------------------------
//      // slope-determination
//      //-------------------------------------------------------------------------------------------

//      //vermeidet 0 im nenner:
//      if (neigh1_x == neigh2_x)
//        {
//          neigh2_x = neigh1_x + 0.0000000001;
//        }

//      boundary_slope =(neigh1_y - neigh2_y) / (neigh1_x - neigh2_x);

//      // winkelabweichung von boundary auf horizontale
//      boundary_angle = atan(boundary_slope);
//      boundary_angle += Pi/2;

//      runParticle->boundary_angle[neignb] = boundary_angle;

//      sigman = 0.5*(-runParticle->sxx - runParticle->syy) +
//      0.5*(-runParticle->sxx + runParticle->syy)*cos(2.0*boundary_angle) -
//      runParticle->sxy*sin(2.0*boundary_angle);

//      sigman *= pressure_scale * pascal_scale;

//      if (isnan(sigman))
//        cout << "isnogud" << endl;

//      return(sigman);
//    }
//  else
//    {
//      return (NULL); //dummy-value if no result was obtained
//    }
//}

double Min_Trans_Lattice::Calc_Curvature_Surface_Energy(
    Particle *prtcl, bool provisoric_transformation) {
  int i, j;
  double E, surfE1, surfE2, surfE3;
  Particle *lneig = NULL, *rneig = NULL;
  bool phase_boundary[6] = {
      false, false, false, false, false, false,
  };

  if (provisoric_transformation) {

    if (prtcl->mineral == 1)
      prtcl->mineral = 2;
    else
      prtcl->mineral = 1;

    for (i = 0; i < 6; i++) {
      if (prtcl->neigP[i]) {
        for (j = 0; j < 6; j++) {
          phase_boundary[i] = prtcl->neigP[i]->is_phase_boundary;
          if (prtcl->neigP[i]->neigP[j]) {
            if (prtcl->neigP[i]->neigP[j]->mineral !=
                prtcl->neigP[i]->mineral) {
              prtcl->neigP[i]->is_phase_boundary = true;
              break;
            }
          }
        }
      }
    }
  }

  for (i = 0; i < 6; i++) {
    if (prtcl->neigP[i]) {
      if (prtcl->neigP[i]->is_phase_boundary &&
          prtcl->neigP[i]->mineral == prtcl->mineral) {
        if (!lneig)
          lneig = prtcl->neigP[i];
        else
          rneig = prtcl->neigP[i];
      }
    }
  }

  // if it's a 'stick'
  if (!rneig)
    rneig = lneig;

  // control
  if (rneig) {
    if (rneig->is_lattice_boundary)
      rneig = NULL;
  }

  if (lneig) {
    if (lneig->is_lattice_boundary)
      lneig = NULL;
  }

  Calc_Particle_Surface_Energy(prtcl, lneig, rneig, surfE1, surfE2, surfE3,
                               false);

  if (provisoric_transformation) {

    if (prtcl->mineral == 1)
      prtcl->mineral = 2;
    else
      prtcl->mineral = 1;

    for (i = 0; i < 6; i++) {
      if (prtcl->neigP[i]) {
        prtcl->neigP[i]->is_phase_boundary = phase_boundary[i];
      }
    }
  }

  // formula see daniels paper, 2007, about stylolithe teeth
  E = (surfE1 + (surfE2 + surfE3) / 2.0) / 2.0;
  // E = E / 9.0;
  E *= 1.728572;

  return (E);
}

double Min_Trans_Lattice::Calc_Particle_Surface_Energy(
    Particle *prtcl, Particle *lneig, Particle *rneig, double &surfE1,
    double &surfE2, double &surfE3, bool provisoric_transformation) {

  int i, j;
  int mark1, mark2;
  int phase_marker, counter;
  Particle *vec[3] = {prtcl, lneig, rneig};
  double surfE;
  bool phase_boundary[6] = {
      false, false, false, false, false, false,
  };

  if (provisoric_transformation) {

    if (prtcl->mineral == 1)
      prtcl->mineral = 2;
    else
      prtcl->mineral = 1;

    for (i = 0; i < 6; i++) {

      if (prtcl->neigP[i]) {
        phase_boundary[i] = prtcl->neigP[i]->is_phase_boundary;
        for (j = 0; j < 6; j++) {
          if (prtcl->neigP[i]->neigP[j]) {
            if (prtcl->neigP[i]->neigP[j]->mineral !=
                prtcl->neigP[i]->mineral) {
              prtcl->neigP[i]->is_phase_boundary = true;
              break;
            }
          }
        }
      }
    }
  }

  surfE1 = surfE2 = surfE3 = 0.0;

  for (i = 0; i < 3; i++) {

    if (vec[i]) {

      // find out if there are NONE open springs: this demands a differnet
      // treatment
      counter = 0;
      for (j = 0; j < 6; j++) {
        if (vec[i]->neigP[j]) {
          if (vec[i]->neigP[j]->mineral != prtcl->mineral) {
            counter++;
          }
        }
      }
      if (counter != 6 && counter != 0) {
        counter = 0;

        phase_marker = vec[i]->neigP[0]->mineral;

        if (phase_marker == prtcl->mineral) {
          mark1 = 1;
          while (vec[i]->neigP[mark1]->mineral == prtcl->mineral) {
            mark1++;
            if (mark1 == 6)
              mark1 = 0;
          }
        } else {
          mark1 = 5;
          while (vec[i]->neigP[mark1]->mineral != prtcl->mineral) {
            mark1--;
            if (mark1 == -1)
              mark1 = 5;
          }
          mark1++;
          if (mark1 == 6)
            mark1 = 0;
        }

        mark2 = mark1;
        counter = 0;
        do {

          if (vec[i]->neigP[mark2]->mineral != prtcl->mineral) {

            counter++;

          } else if (counter != 0) {

            surfE = vec[i]->surf_free_E / (vec[i]->real_radius);
            surfE *= ((counter - 2.0) / 4.0);
            surfE *= vec[i]->mV;

            switch (i) {
            case 0:
              surfE1 += surfE;
              break;
            case 1:
              surfE2 += surfE;
              break;
            case 2:
              surfE3 += surfE;
              break;
            }

            counter = 0;
          }

          mark2++;
          if (mark2 == 6)
            mark2 = 0;

        } while (mark2 != mark1);
      } else if (counter == 6) {

        surfE = vec[i]->surf_free_E / (vec[i]->real_radius);
        surfE *= ((counter - 2.0) / 4.0);
        surfE *= vec[i]->mV;

        switch (i) {
        case 0:
          surfE1 = surfE;
          break;
        case 1:
          surfE2 = surfE;
          break;
        case 2:
          surfE3 = surfE;
          break;
        }
      } else if (counter == 0) {

        switch (i) {
        case 0:
          surfE1 = 0.0;
          break;
        case 1:
          surfE2 = 0.0;
          break;
        case 2:
          surfE3 = 0.0;
          break;
        }
      }
    }
  }

  if (provisoric_transformation) {

    if (prtcl->mineral == 1)
      prtcl->mineral = 2;
    else
      prtcl->mineral = 1;

    for (i = 0; i < 6; i++) {
      if (prtcl->neigP[i]) {
        prtcl->neigP[i]->is_phase_boundary = phase_boundary[i];
      }
    }
  }
}

// note: the particles given to this function should be neighbours already. the
// function doesn't check this!
bool Min_Trans_Lattice::Is_Same_Surface(Particle *prtcl1, Particle *prtcl2) {

  int i, j;
  Particle *common_neig;

  if (prtcl1->mineral != prtcl2->mineral)
    return (false);

  // the particles are on the same surface (given they are NEIGBORS!), if there
  // is one particle of the other phase which is connected to both of them
  for (i = 0; i < 6; i++) {
    if (prtcl1->neigP[i]) {
      if (prtcl1->neigP[i]->mineral != prtcl1->mineral) {
        for (j = 0; j < 6; j++) {
          if (prtcl1->neigP[i]->neigP[j]) {
            if (prtcl1->neigP[i]->neigP[j] == prtcl2) {
              return (true);
            }
          }
        }
      }
    }
  }

  return (false);
}

// does basic energy-calculaton, i.e. elastic energy,
// surface energy
void Min_Trans_Lattice::Energy_Calculation(Particle *partikel) {
  int j, ii,
      jj,           // counters for loops
      spring_count, // counter of open springs
      box,
      pos;      // dummies for relaxation box position
  double lame1, // first lame constant
      lame2,    // second lame constant
      uxx,      // x component of infinitesimal strain tensor
      uyy,      // y component of infinitesimal strain tensor
      uxy,      // xy component of infinitesimal strain tensor
      dx,
      dy, // delta x, delta y
      alen,
      dd; // equilibrium length, normal length
  double fn, fx,
      fy, // normal force, force in x, force in y
      uxi,
      uyi;          // specific components of strain tensor
  double rate,      // rate for reaction
      surfE,        // surface energy
      rep_constant, // repulsion constant
      eEl_distortion, eEl_dilatation, bulk;
  double sum_Energy;
  double bound_stress;

  Particle *neig; // particle object pointer

  double surface_right[10];
  double surface_left[10];
  double mean, dxx, dyy, dxy;
  Particle *help;

  // for the surface-energy
  Update_Surf_Neighbour_List();

  partikel->open_springs = 0;
  partikel->bound_stress = 0.0;
  partikel->fluid_neig_count = 0;
  spring_count = 0;

  partikel->done = true;

  for (j = 0; j < 6; j++) {
    if (partikel->neigP[j]) {
      partikel->neigP[j]->done = true; // set flag for repulsion

      spring_count = spring_count + 1; // count number of active connections

      if (partikel->grain == partikel->neigP[j]->grain) {

        if (partikel->neigP[j]->is_boundary) {

          if (!partikel->neigP[j]->is_lattice_boundary) {

            partikel->fluid_neigP[partikel->fluid_neig_count] =
                partikel->neigP[j];
            partikel->fluid_neig_count = partikel->fluid_neig_count + 1;
          }
        }
      }
      // else
      {
        //
        if (partikel->neigP[j]->isHole)
          partikel->open_springs = partikel->open_springs + 1;
        if (!partikel->neigP[j]->isHole) {
          if (partikel->is_boundary) {
            if (partikel->neigP[j]->grain != partikel->grain) {

              partikel->open_springs += 1;

              dx = partikel->neigP[j]->xpos - partikel->xpos;
              dy = partikel->neigP[j]->ypos - partikel->ypos;

              dd = sqrt((dx * dx) + (dy * dy));
              alen = partikel->radius + partikel->neigP[j]->radius;

              // get unit length in x and y

              if (dd != 0.0) {
                uxi = dx / dd;
                uyi = dy / dd;
              } else {
                cout << " zero divide in Reaction2 " << endl;
                uxi = 0.0;
                uyi = 0.0;
              }

              fn = partikel->springf[j] * (dd - alen); // normal force
              fx = fn * uxi;                           // force x
              fy = fn * uyi;                           // force y

              fn = fx * uxi + fy * uyi; // force in

              fn = fn * 2.0 * partikel->radius; // scaling

              fn = fn / ((partikel->radius * partikel->radius * Pi) / 6.0);
            }
          }
        }
      }
    }
  }

  if (spring_count < 6) {
    partikel->open_springs = partikel->open_springs + 6 - spring_count;

    for (ii = 0; ii < 3; ii++) {
      for (jj = 0; jj < 3; jj++) {

        if (jj == 0)
          box = particlex * (-1); // lower row of
        // boxes
        if (jj == 1)
          box = 0; // middle row of boxes
        if (jj == 2)
          box = particlex; // upper row of
        // boxes

        pos = partikel->box_pos + (ii - 1) + box; // define
        // box
        // position
        //

        if (pos > 0) {
          neig = 0;
          neig = repBox[pos]; // get particle at
          // that position

          while (neig) {

            if (!neig->done) {  // if not already done
              if (neig->xpos) { // just in case

                dx = neig->xpos - partikel->xpos;
                dy = neig->ypos - partikel->ypos;

                dd = sqrt((dx * dx) + (dy * dy));

                alen = partikel->radius + neig->radius;

                if (dd != 0.0) {
                  uxi = dx / dd;
                  uyi = dy / dd;
                } else {
                  cout << " zero divide in Relax " << endl;
                  uxi = 0.0;
                  uyi = 0.0;
                }

                rep_constant = (partikel->young + neig->young) / 2.0;

                // scale back to spring
                // constant

                rep_constant = rep_constant * sqrt(3.0) / 2.0;

                fn = rep_constant *
                     (dd -
                      alen); // normal"force per length" (pressure times length)
                fx = fn * uxi; // " force x"
                fy = fn * uyi; // " force y"

                if (fn < 0.0) // if compressive
                {

                  fn = fx * uxi + fy * uyi; // force in direction of spring
                  fn = fn * 2.0 * partikel->radius; // scaling to particle size
                  fn = fn / ((partikel->radius * partikel->radius * Pi) / 6.0);

                  // -----------------------------------------------------------
                  // and add that to the
                  // particles boundary
                  // stress
                  // -----------------------------------------------------------

                  partikel->bound_stress = partikel->bound_stress + fn;
                }
              }
            }
            if (neig->next_inBox) {    // if there is another one in box
              neig = neig->next_inBox; // move pointer
            } else {
              break;
            }
          }
        }
      }
    }
  }
  // ---------------------------------------------------
  // now turn off the flags
  // ---------------------------------------------------

  partikel->done = false;

  for (j = 0; j < 6; j++) {
    if (partikel->neigP[j]) {
      partikel->neigP[j]->done = false;
    }
  }

  lame1 = partikel->young * pois / ((1.0 - 2.0 * pois) * (1.0 + pois));
  lame2 = partikel->young / (2.0 + 2.0 * pois);

  // -----------------------------------------------------------------------
  // getting strain tensor
  // -----------------------------------------------------------------------

  uxy = ((1 + pois) / partikel->young) * partikel->sxy;
  uxx = (1 / partikel->young) * (partikel->sxx - (partikel->syy * pois));
  uyy = (1 / partikel->young) * (partikel->syy - (partikel->sxx * pois));

  // -----------------------------------------------------------------------
  // calculate elastic energy for particle
  // -----------------------------------------------------------------------

  partikel->eEl = (0.5 * lame1 * (uxx + uyy) * (uxx + uyy) +
                   (lame2 * ((uxx * uxx) + (uyy * uyy) + 2 * (uxy * uxy))));
  // scaling:
  partikel->eEl *= partikel->mV * pressure_scale * pascal_scale;

  mean = -(partikel->sxx + partikel->syy) / 2.0;
  dxx = (-partikel->sxx) - mean;
  dxy = -partikel->sxy;
  dyy = (-partikel->syy) - mean;

  bulk = partikel->young / (3 * (1 - 2 * pois));

  eEl_distortion = (dxx * dxx + dyy * dyy + 2.0 * dxy * dxy) / (4.0 * lame2);
  eEl_dilatation = 0.5 * (mean * mean / bulk);

  partikel->eEl_distortion =
      eEl_distortion * partikel->mV * pressure_scale * pascal_scale;
  partikel->eEl_dilatation =
      eEl_dilatation * partikel->mV * pressure_scale * pascal_scale;
}

// called by reaction_two(), to restore the original spinel-distribution
// without having to relax again. saves some time, costs some memory.
void Min_Trans_Lattice::Invert_Provisoric_Mineral_Transformation() {
  int i, j;

  // firstly, restore the old heat-distribution
  // Invert_Heat_Lattice();

  // secondly, restore old particle-volumes, radii and young-moduli
  for (i = 0; i < numParticles; i++) {

    runParticle->mineral = runParticle->previous_mineral;
    runParticle->mV = runParticle->previous_mV; // restore the molecular volume
    runParticle->young = runParticle->previous_young; // restore young's modulus

    runParticle->viscosity = runParticle->prev_viscosity;

    runParticle->original_springf = runParticle->previous_original_springf;
    for (j = 0; j < 6; j++) // and of course the spring force
    {
      runParticle->springf[j] = runParticle->previous_springf[j];
      runParticle->rad[j] = runParticle->previous_rad[j];
    }

    runParticle->grain = runParticle->previous_grain;

    runParticle->radius = runParticle->previous_radius; // and change radius
    runParticle->real_radius = runParticle->previous_real_radius;
    runParticle->area = runParticle->previous_area;

    // cout << "temperature (reverse): " << runParticle->temperature << endl;

    runParticle->phase_change_required = false;

    runParticle = runParticle->nextP;
  }

  Make_Phase_Boundaries();

  AdjustParticleConstants();
}

void Min_Trans_Lattice::InitProvisoricElasticProperties() {

  runParticle = &refParticle;

  for (int i = 0; i < numParticles; i++) {

    runParticle->spin_young = runParticle->olivine_young = runParticle->young;

    for (int j = 0; j < 6; j++) {

      if (runParticle->neigP[j])
        runParticle->spin_springf[j] = runParticle->olivine_springf[j] =
            runParticle->springf[j];
      else
        runParticle->spin_springf[j] = runParticle->olivine_springf[j] = 0.0;
    }

    runParticle = runParticle->nextP;
  }
}

void Min_Trans_Lattice::Mineral_Transformation_No_Spin_Distri() {
  int i, j, counter;

  bool new_grain; // flag for new grain or not
  double factor;

  // add latent heat release (even if there is no complete transformation)
  if (heat_distribution.heat_flow_enabled && !no_heat_flag) {

    runParticle = &refParticle;

    for (i = 0; i < numParticles; i++) {

      if (runParticle->newly_overgrown != 0.0) {
        switch (runParticle->mineral) {
        case 1:
          runParticle->temperature += Latent_Heat_Release(runParticle, OL_SPIN);
          Set_Heat_Lattice(runParticle);
          // cout << runParticle->temperature << endl;
          break;
        case 2:
          runParticle->temperature += Latent_Heat_Release(runParticle, SPIN_OL);
          Set_Heat_Lattice(runParticle);
          break;
        }
      }
      runParticle = runParticle->nextP;
    }
  }

  // set new grain / new mineral
  for (i = 0; i < numParticles; i++) {
    if (runParticle->prob >= 1.0) {

      cout << "particle-number: " << runParticle->nb << endl;
      cout << "prob: " << runParticle->prob << endl;

      runParticle->overgrown_area = 0.0;
      // fuer die zeit-anpassung nach erster umwandlung

      switch (runParticle->mineral) {
      case 1:
        new_grain = true;         // first we assume its a new grain of spinell
        runParticle->mineral = 2; // change mineral index

        runParticle->viscosity *= 1.2;

        // runParticle->young = runParticle->spin_young; // change repulsion
        // constant

        if (!use_homogeneous_spin)
          runParticle->young = runParticle->young * 1.4;

        for (j = 0; j < 8; j++) {
          if (runParticle->neigP[j]) { // if spring is active

            // -----------------------------------------------------------
            // now check what the neighbours are, if they are
            // of a
            // different mineral spring is a boundary, if not
            // this is
            // not a new grain of spinell in olivine.
            // -----------------------------------------------------------

            if (runParticle->neigP[j]->mineral == 1) {
              runParticle->neigP[j]->is_boundary = true; // real boundary
            } else {
              new_grain = false;
              runParticle->grain =
                  runParticle->neigP[j]
                      ->grain; // give it the name of the other grain of spinell
              runParticle->neigP[j]->is_phase_boundary = true;
            }
          }
        }

        if (new_grain) {                          // if this is a new grain
          runParticle->grain = highest_grain + 1; // new number for grain
          highest_grain = highest_grain + 1;      // and adjust grain counter
          spinel_grain_nb += 1;
          grain_young.push_back(runParticle->young * 1.4);
        } else if (runParticle->neig_of_two_spinel_grains) {
          runParticle->grain = Determine_grain_number(runParticle);
        }

        break;

      case 2:
        new_grain = true;         // first we assume its a new grain of olivine
        runParticle->mineral = 1; // change mineral index0

        runParticle->viscosity /= 1.2;

        // runParticle->young = runParticle->olivine_young;  // change repulsion
        // constant

        for (j = 0; j < 8; j++) {
          if (runParticle->neigP[j]) // if spring is active
          {

            // -----------------------------------------------------------
            // now check what the neighbours are, if they are
            // of a
            // different mineral spring is a boundary, if not
            // this is
            // not a new grain of spinell in olivine.
            // -----------------------------------------------------------

            runParticle->young = runParticle->olivine_young;

            if (runParticle->neigP[j]->mineral == 2) {
              runParticle->neigP[j]->is_boundary = true; // real boundary
              runParticle->neigP[j]->is_phase_boundary = true;
            } else {
              new_grain = false;
              runParticle->grain =
                  runParticle->neigP[j]
                      ->grain; // give it the name of the other grain of spinell
            }
          }
        }

        if (new_grain) {                          // if this is a new grain
          runParticle->grain = highest_grain + 1; // new number for grain
          highest_grain = highest_grain + 1;      // and adjust grain counter
          grain_young.push_back(runParticle->olivine_young);
        }

        break;
      }

      // zur zeit-minimierung (s.a. Set_Mineral_Parameters, Time-change
      // funktion)
      timeflag = true;

      for (j = 0; j < 4; j++)
        runParticle->l[j] = 0.0;
    }

    // weil die uhren jetzt auch fuer addition zurckgedreht werden wollen:
    // runParticle->previous_prob = 0.0;
    // runParticle->prob = 0.0;

    runParticle = runParticle->nextP;
  }

  // assign youngs modulus
  if (use_homogeneous_spin) {
    for (i = 0; i < numParticles; i++) {

      if (runParticle->mineral == 2)
        runParticle->young = grain_young[runParticle->grain];

      runParticle = runParticle->nextP;
    }
  }

  counter = 0;

  for (i = 0; i < numParticles; i++) {
    if (runParticle->grain != beware_grain1 &&
        runParticle->grain != beware_grain2 &&
        runParticle->grain != beware_grain3) {

      runParticle->mV = Adjust_Molare_Volume();
      runParticle->previous_real_radius = runParticle->real_radius;
      runParticle->real_radius =
          pow(((3.0 * runParticle->mole_per_particle * runParticle->mV) /
               (4.0 * 3.1415)),
              (1.0 / 3.0)); // and change radius
      factor = (runParticle->real_radius / runParticle->previous_real_radius);

      runParticle->radius *= factor;
      for (j = 0; j < 6; j++) {
        if (runParticle->neigP[j])
          runParticle->rad[j] *= factor;
      }

      runParticle->area = pow(double(runParticle->real_radius), 2.0) * 3.1415;
    }

    runParticle = runParticle->nextP;
  }

  for (i = 0; i < numParticles; i++) {

    runParticle->is_boundary = false;
    for (j = 0; j < 6; j++) {
      if (runParticle->neigP[j]) {
        if (runParticle->neigP[j]->grain != runParticle->grain) {
          runParticle->is_boundary = true;
        }
      }
    }

    runParticle = runParticle->nextP;
  }

  AdjustParticleConstants();
}

int Min_Trans_Lattice::Determine_grain_number(Particle *prtcl) {

  // bestimme, aus welcher Richtung her der Partikel am stärksten gewachsen ist

  int i, buffer, counter;
  double dir[6];
  bool tmp[6] = {
      false, false, false, false, false, false,
  };

  dir[0] = prtcl->l[0];
  dir[1] = prtcl->l[0] * cos(Pi / 3) + prtcl->l[1] * cos(Pi / 6);
  dir[2] = prtcl->l[1] * cos(Pi / 6) + prtcl->l[2] * cos(Pi / 3);
  dir[3] = prtcl->l[2];
  dir[4] = prtcl->l[2] * cos(Pi / 3) + prtcl->l[3] * cos(Pi / 6);
  dir[5] = prtcl->l[3] * cos(Pi / 6) + prtcl->l[0] * cos(Pi / 3);

  // if the first try doesnt work, we have a while loop
  buffer = 0;
  while (prtcl->neigP[buffer]->mineral != 2) {
    buffer = 0;
    for (i = 1; i < 6; i++) {
      if (!tmp[i]) {
        if (dir[i] > dir[buffer]) {
          buffer = i;
        }
      }
    }
    tmp[buffer] = true;
    counter++;
    if (counter > 5)
      break;
  }

  // cout << "fatal error in Determine_grain_number" << endl;

  return (prtcl->neigP[buffer]->grain);
}

// function sets the "is_phase_boundary"-flag of
// all particles.
void Min_Trans_Lattice::Make_Phase_Boundaries() {
  int i, j, ii, jj;
  Particle *neig;
  int box, pos; // dummies for relaxation box position tensor

  for (i = 0; i < numParticles; i++) { // and loop through particles
    for (j = 0; j < 8; j++) {          // loop through neighbours of particle
      if (runParticle->neigP[j]) {     // if there is a neighbour
        if ((runParticle->neigP[j]->mineral == 1 &&
             runParticle->mineral == 2) ||
            (runParticle->neigP[j]->mineral == 2 &&
             runParticle->mineral == 1)) { // if I am a boundary
          runParticle->is_phase_boundary = true;
          break;
        } else {
          runParticle->is_phase_boundary = false;
        }
      }
    }
    runParticle = runParticle->nextP; // go on looping
  }
}

// function returns mean stress over all particles

double Min_Trans_Lattice::Pressure() {
  int i, p_counter;
  double smax, smin, mean_stress, sxx, sxy, syy;
  // set values

  // first find out about stress: mean_stress
  p_counter = 0;
  mean_stress = 0.0;
  runParticle = &refParticle;            // start
  for (i = 0; i < numParticles; i++) {   // look which particles are in the box
    if (runParticle->ypos > 0.1) {       // larger ymin
      if (runParticle->ypos < 0.9) {     // smaller ymax
        if (runParticle->xpos > 0.1) {   // larger xmin
          if (runParticle->xpos < 0.9) { // smaller xmax
            p_counter = p_counter + 1;   // count particles
            // set stress tensor
            sxx = runParticle->sxx;
            syy = runParticle->syy;
            sxy = runParticle->sxy;
            // eigenvalues
            smax = ((sxx + syy) / 2.0) +
                   sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) + sxy * sxy);
            smin = ((sxx + syy) / 2.0) -
                   sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) + sxy * sxy);
            mean_stress += (smax + smin) / 2.0;
          }
        }
      }
    }
    runParticle = runParticle->nextP;
  }
  // and divide by number of particles
  mean_stress = mean_stress / p_counter;

  return (-mean_stress * pascal_scale * pressure_scale);
}

double Min_Trans_Lattice::OlivinePressure() {
  int i, p_counter;
  double smax, smin, mean_stress, sxx, sxy, syy;
  // set values

  // first find out about stress: mean_stress
  p_counter = 0;
  mean_stress = 0.0;
  runParticle = &refParticle;          // start
  for (i = 0; i < numParticles; i++) { // look which particles are in the box

    if (runParticle->mineral == 1 && !runParticle->is_top_lattice_boundary &&
        !runParticle->is_bottom_lattice_boundary &&
        !runParticle->is_left_lattice_boundary &&
        !runParticle->is_right_lattice_boundary) {

      p_counter = p_counter + 1; // count particles
      // set stress tensor
      sxx = runParticle->sxx;
      syy = runParticle->syy;
      sxy = runParticle->sxy;
      // eigenvalues
      smax = ((sxx + syy) / 2.0) +
             sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) + sxy * sxy);
      smin = ((sxx + syy) / 2.0) -
             sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) + sxy * sxy);
      mean_stress += (smax + smin) / 2.0;
    }
    runParticle = runParticle->nextP;
  }
  // and divide by number of particles
  mean_stress = mean_stress / p_counter;

  return (-mean_stress * pascal_scale * pressure_scale);
}

// calculates the time it takes for the fastest particle to transform
// called in start_reactions()
void Min_Trans_Lattice::Gbm() {

  int i, j, k, l;

  double rate, time2, prob, lowest_time, lowest_time2, local_time, v, t,
      time_interval2, dd;

  double random_number;

  Particle *neig;

  vector<Particle *> neig_list;

  neig = NULL;

  lowest_time = time_interval; // zeit nach nukleiierung, wenn keine
                               // nukleiierung, dann actual_time
  highest_prob = 0.0;

  // cout << "grainboundary migration, time = " << actual_time << ",
  // highest_prob = " << highest_prob << endl;

  cout << "gbm 1" << endl;
  cout << "pressure: " << Pressure() << endl;

  secondary_timestep++;

  if (time_interval > 0.0) {
    for (i = 0; i < numParticles; i++) {
      runParticle->neig3 = NULL;
      runParticle->prob = 0.0;
      runParticle->spinel_content = 0.0;

      runParticle->delta_eEl_distortion = 0.0;
      runParticle->delta_eEl_dilatation = 0.0;

      runParticle->newly_overgrown = 0.0;

      runParticle->transform = false;

      for (j = 0; j < 6; j++) {
        runParticle->dir_rate[j] = 0.0;
      }

      for (j = 0; j < 4; j++) {
        runParticle->save_rate[j] = 0.0;
      }

      if (protect_flickering_particles) {
        if (runParticle->transformation_timestep != 0) {
          if (runParticle->transformation_timestep == secondary_timestep - 1) {
            runParticle->transformed_in_last_timestep = true;
          } else {
            runParticle->transformed_in_last_timestep = false;
          }
        } else {
          runParticle->transformed_in_last_timestep = false;
        }
      }

      runParticle = runParticle->nextP;
    }

    SetProvisoricElasticProperties();

    neig_list.clear();

    cout << "time_interval: " << time_interval
         << " refParticle->time: " << refParticle.time << endl;

    // transform boundary-particles, calculate delta G
    Call_Reaction_Two();

    // initialize random-number generator:
    srand(std::time(0));

    // set a standard

    lowest_time = time_interval;

    for (i = 0; i < numParticles; i++) {
      if (runParticle->is_phase_boundary && !runParticle->is_lattice_boundary) {
        if (runParticle->grain != beware_grain1 &&
            runParticle->grain != beware_grain2 &&
            runParticle->grain != beware_grain3 && !runParticle->no_reaction) {
          t = CalcTime();

          if (t < lowest_time) {
            lowest_time = t;
            neig = runParticle;

            if (t == 0.0)
              break;
          }
        }
      }
      runParticle = runParticle->nextP;
    }

    //	neig_list.push_back (neig);

    //  for (i=0; i<numParticles; i++)
    //    {
    //      if (runParticle->is_phase_boundary &&
    //      !runParticle->is_lattice_boundary)
    //        {
    //          if (runParticle->grain != beware_grain1 && runParticle->grain !=
    //          beware_grain2 && runParticle->grain != beware_grain3 &&
    //          !runParticle->no_reaction)
    //            {
    //              t = CalcTime();

    //              if (t == lowest_time)
    //                {
    //                  neig_list.push_back (runParticle);
    //                }

    //            }
    //        }
    //      runParticle = runParticle->nextP;
    //    }

    if (lowest_time == 0.0)
      cout << "warning: lowest time == 0.0" << endl;

    cout << "lowest time: " << lowest_time << endl;

    if (lowest_time > 0.0 && neig) {

      neig->prob = 2.0;
      highest_prob = 2.0;

      for (i = 0; i < numParticles; i++) {
        if (runParticle->is_phase_boundary &&
            !runParticle->is_lattice_boundary) {
          CalcArea(lowest_time);
        }
        runParticle = runParticle->nextP;
      }

      for (i = 0; i < numParticles; i++) {
        if (runParticle->is_phase_boundary &&
            !runParticle->is_lattice_boundary) {
          if (runParticle->transform && runParticle->grain != beware_grain1 &&
              runParticle->grain != beware_grain2 &&
              runParticle->grain != beware_grain3 &&
              !runParticle->no_reaction) {
            neig_list.push_back(runParticle);

            if (runParticle->transformed_in_last_timestep) {
              neig_list.pop_back();
              runParticle->no_reaction = true;
            }
          }
        }
        runParticle = runParticle->nextP;
      }

      if (lowest_time < time_interval && neig_list.size() > 0) {
        for (i = 0; i < neig_list.size(); i++) {
          neig_list[i]->prob = 2.0;
          highest_prob = 2.0;
          runParticle->transformation_timestep = secondary_timestep;
        }

        DumpDiffStressAndTime(neig_list[0], lowest_time, neig_list[0]->eEl_dif);
      } else
        highest_prob = 0.0;

      // time setting
      time_interval -= lowest_time;

      if (time_interval < 0.0) {
        for (i = 0; i < numParticles; i++) {
          runParticle->prob = 0.0;
          runParticle = runParticle->nextP;
        }
        time_interval = 0.0;
        highest_prob = 0.0;
      }

      cout << time_interval << endl;

      for (i = 0; i < numParticles; i++) {
        if (runParticle->prob > 1.0) {
          k = 0;
          for (j = 0; j < 6; j++) {
            if (runParticle->neigP[j]) {
              if (runParticle->neigP[j]->mineral == runParticle->mineral) {
                k++;
              } else
                l = j;
            }
          }
          if (k == 5) {
            if (l != 0 && l != 3) {

              if (l == 2 || l == 4)
                runParticle->neigP[3]->prob = 2.0;
              else
                runParticle->neigP[0]->prob = 2.0;

              cout << "did it" << endl;
            }
          }
        }
        runParticle = runParticle->nextP;
      }

      if (neig_list.size() > 0) {
        Mineral_Transformation_No_Spin_Distri();
        for (i = 0; i < neig_list.size(); i++) {
          cout << "already transformed particles in this loop: " << i + 1
               << endl;
          for (j = 0; j < 4; j++) {
            neig->l[j] = 0.0;
          }
        }
      }
    } else if (lowest_time == 0.0) {
      neig->prob = 2.0;
      highest_prob = 2.0;
      time_interval -= lowest_time;
      Mineral_Transformation_No_Spin_Distri();
    } else
      cout << "t < 0" << endl;

    // adjust changed springs on grain boundaries
    AdjustConstantGrainBoundaries();
    Make_Phase_Boundaries();

    // save the rates to a file
    if (neig_list.size() > 0)
      growth_rate(lowest_time, neig_list[0]); // direction of the transformation

    if (heat_distribution.heat_flow_enabled)
      Heat_Flow(lowest_time);

    // and relax
    Relaxation();
  } else {
    highest_prob = 0.0;
  }

  /*if (highest_prob != 0)
      ViscousRelax(1, lowest_time);
  else
      ViscousRelax(1, time_interval);*/

  // Exchange_Probabilities();
  UpdateElle();
}

// calculates the energy for different grain-growth scenarios
// (only spinel grows vs only olvine growths)
// called in gbm()
void Min_Trans_Lattice::Call_Reaction_Two() {
  int i, j, actual_grain_nb;

  runParticle = &refParticle;

  // preliminary settings
  for (i = 0; i < numParticles; i++) {

    runParticle->phase_change_required = false;
    runParticle->neig_of_two_spinel_grains = false;

    // dummy number, to large to be realistic
    actual_grain_nb = 100000000;

    // tag olivine-particles which have TWO spinel-grains as neighbours
    if (!runParticle->is_lattice_boundary) {
      if (runParticle->is_phase_boundary) {
        if (runParticle->mineral == 1) {
          for (j = 0; j < 6; j++) {
            if (runParticle->neigP[j]) {
              if (runParticle->neigP[j]->mineral == 2) {
                actual_grain_nb = runParticle->neigP[j]->grain;
                break;
              }
            }
          }
          for (j = 0; j < 6; j++) {
            if (runParticle->neigP[j]) {
              if (runParticle->neigP[j]->mineral == 2) {
                if (actual_grain_nb != runParticle->neigP[j]->grain) {
                  runParticle->neig_of_two_spinel_grains = true;
                }
              }
            }
          }
        }
      }
    }

    runParticle = runParticle->nextP;
  }

  if (include_elastic_energy) {

    for (i = 0; i < numParticles; i++) {
      if (!runParticle->is_lattice_boundary) {
        if (runParticle->is_phase_boundary) {
          if (runParticle->mineral == 1) {
            runParticle->phase_change_required = true;
          }
        }
      }
      runParticle = runParticle->nextP;
    }

    for (i = 0; i < numParticles; i++) {

      if (runParticle->no_reaction)
        runParticle->phase_change_required = false;

      runParticle = runParticle->nextP;
    }

    Reaction_Two();

    if (include_negative_rates_into_elastic_relax) {

      for (i = 0; i < numParticles; i++) {

        runParticle->phase_change_required = false;

        runParticle = runParticle->nextP;
      }

      for (i = 0; i < numParticles; i++) {
        if (!runParticle->is_lattice_boundary) {
          if (runParticle->is_phase_boundary) {
            if (runParticle->mineral == 2) {

              runParticle->phase_change_required = true;
            }
          }
        }
        runParticle = runParticle->nextP;
      }

      for (i = 0; i < numParticles; i++) {

        if (runParticle->no_reaction)
          runParticle->phase_change_required = false;

        runParticle = runParticle->nextP;
      }

      Reaction_Two();
    }

  }
  //	if (include_elastic_energy) {
  //	  // from now on change particles spinel-grain by spinel-grain
  //	  for ( actual_grain_nb = olivine_grain_nb+1; actual_grain_nb <
  //highest_grain+1; actual_grain_nb++ ) {
  //
  //	  	cout << "probing grain nb " << actual_grain_nb << endl;
  //
  //	  	for (i=0; i<numParticles; i++) {
  //	  		if (!runParticle->is_lattice_boundary) {
  //	  			if (runParticle->is_phase_boundary) {
  //	//  				if
  //(!runParticle->neig_of_two_spinel_grains){ 		  				for (j=0; j<6; j++) { 		  					if
  //(runParticle->neigP[j]->grain == actual_grain_nb) {
  //		  						runParticle->phase_change_required
  //= true;
  //		  					}
  //	  					}
  //	//  				}
  //	  			}
  //	  		}
  //	  		runParticle = runParticle->nextP;
  //	  	}
  //
  //	  	for (i=0; i<numParticles; i++) {
  //
  //	  		if (runParticle->no_reaction)
  //		  		runParticle->phase_change_required = false;
  //
  //	  		runParticle = runParticle->nextP;
  //	  	}
  //
  //	  	Reaction_Two();
  //
  //	  	for (i=0; i<numParticles; i++) {
  //
  //	  		runParticle->phase_change_required = false;
  //
  //	  		runParticle = runParticle->nextP;
  //	  	}
  //	  }
  //  }
  else {
    for (i = 0; i < numParticles; i++) {

      if (runParticle->mineral == 1 && runParticle->is_phase_boundary &&
          !runParticle->is_lattice_boundary) {
        runParticle->phase_change_required = true;
      } else {
        runParticle->phase_change_required = false;
      }

      if (runParticle->no_reaction)
        runParticle->phase_change_required = false;

      runParticle = runParticle->nextP;
    }

    Reaction_Two();
  }

  // clean up
  for (i = 0; i < numParticles; i++) {

    runParticle->phase_change_required = false;

    runParticle = runParticle->nextP;
  }
}

double Min_Trans_Lattice::CalcArea(double lowest_time) {
  int i, k;
  double rate[4];             // indicate directions of rate in the rectangle
  double l, l1, l2, ln1, ln2; // the length of a basic square size
  double x1, x2, F0,
      F1; // transformed area variables (F0=partikel-area, F1=remaining area
          // after timestep, x1,x2=new side-lengthes

  F0 = runParticle->real_radius * runParticle->real_radius * Pi;
  l = sqrt(F0);

  l1 = l - runParticle->l[0] - runParticle->l[2];
  l2 = l - runParticle->l[1] - runParticle->l[3];

  if (l1 < 0.0)
    l1 = 0.0;

  if (l2 < 0.0)
    l2 = 0.0;

  for (i = 0; i < 4; i++)
    rate[i] = 0.0;

  runParticle->l[0] += lowest_time * runParticle->save_rate[0];
  runParticle->l[1] += lowest_time * runParticle->save_rate[1];
  runParticle->l[2] += lowest_time * runParticle->save_rate[2];
  runParticle->l[3] += lowest_time * runParticle->save_rate[3];

  for (i = 0; i < 4; i++) {
    if (runParticle->l[i] < 0.0) {
      runParticle->l[i] = 0.0;
    }

    if (runParticle->l[i] > l) {
      runParticle->l[i] = l;
      runParticle->transform = true;
    }
  }

  // for output: goes into grate2 in ::growth_rate
  ln1 = l - runParticle->l[0] - runParticle->l[2];
  ln2 = l - runParticle->l[1] - runParticle->l[3];

  runParticle->newly_overgrown = l1 * l2 - ln1 * ln2;

  runParticle->overgrown_area = F0 - ln1 * ln2;

  // just in case
  if (runParticle->overgrown_area < 0.0)
    runParticle->overgrown_area = 0.0;

  if (runParticle->overgrown_area >= F0) {
    runParticle->overgrown_area = F0;
    runParticle->transform = true;
  }

  if (numeric_compensation)
    if (F0 - runParticle->overgrown_area < F0 * allowed_error_factor)
      runParticle->transform = true;
}

double Min_Trans_Lattice::CalcTime() {

  int i, j, buffer1, buffer2;
  double rate[4];   // indicate directions of rate in the rectangle
  double l1, l2, l; // the length of a basic square size
  double a;         // the transformed area
  double lowest_time = time_interval, t = time_interval;
  int k = 0;
  double time1, time2, max_time = 0.0;
  double cur_angle, dx1, dy1, dx2, dy2, length, scalarproduct;

  l1 = l2 = l = sqrt(runParticle->real_radius * runParticle->real_radius * Pi);

  l1 -= (runParticle->l[0] + runParticle->l[2]);
  l2 -= (runParticle->l[1] + runParticle->l[3]);

  if (l1 < 0.0)
    l1 = 0.0;

  if (l2 < 0.0)
    l2 = 0.0;

  for (i = 0; i < 4; i++)
    rate[i] = 0.0;

  // convert hexagonaly directed rates to rates relating to a square
  // rate[i] += cos(cur_angle) * runParticle->dir_rate[j];

  for (i = 0; i < 4; i++) {

    switch (i) {
    case 0:
      dx1 = 1;
      dy1 = 0;
      break;
    case 1:
      dx1 = 0;
      dy1 = 1;
      break;
    case 2:
      dx1 = -1;
      dy1 = 0;
      break;
    case 3:
      dx1 = 0;
      dy1 = -1;
      break;
    }

    for (j = 0; j < 6; j++) {

      if (runParticle->neigP[j]) {

        dx2 = runParticle->neigP[j]->xpos - runParticle->xpos;
        dy2 = runParticle->neigP[j]->ypos - runParticle->ypos;

        length = sqrt(dx2 * dx2 + dy2 * dy2);

        scalarproduct = dx1 * dx2 + dy1 * dy2;

        cur_angle = acos(scalarproduct / length);

        if (cur_angle < (Pi / 2.0))
          rate[i] += cos(cur_angle) * runParticle->dir_rate[j];
      }
    }
  }

  if (!(rate[0] <= 0.0 && rate[2] <= 0.0)) {

    if (rate[0] + rate[2] > 0) {

      lowest_time = l1 / (rate[0] + rate[2]);

      if (rate[0] < 0)
        if (lowest_time * (-rate[0]) > runParticle->l[0])
          lowest_time = (l1 + runParticle->l[0]) / rate[2];

        else if (rate[2] < 0)
          if (lowest_time * (-rate[2]) > runParticle->l[2])
            lowest_time = (l1 + runParticle->l[2]) / rate[0];
    }
  }

  if (!(rate[0] <= 0.0 && rate[2] <= 0.0)) {
    if (rate[1] + rate[3] > 0) {

      t = l2 / (rate[1] + rate[3]);

      if (rate[1] < 0)
        if (t * (-rate[1]) > runParticle->l[1])
          t = (l2 + runParticle->l[1]) / rate[3];

        else if (rate[3] < 0)
          if (t * (-rate[3]) > runParticle->l[3])
            t = (l2 + runParticle->l[3]) / rate[1];
    }
  }

  /************************/

  if (t < lowest_time)
    lowest_time = t;

  if (t == 0)
    cout << "bad" << endl;

  for (i = 0; i < 4; i++)
    runParticle->save_rate[i] = rate[i];

  if (lowest_time < 0)
    cout << "doof" << endl;

  return lowest_time;
}

// this function calculates the time, that it takes to overgrow the
// area of a particle, depending on the number of neighbours (the
// contact area) and the rates into this direction
// idea: calculate the time for either two sides to meet at a certain
// point, choose the lowest calculated time
//  called in gbm()
// double
// Min_Trans_Lattice::CalcTime()
//{
//   int i,j, buffer1, buffer2;
//   double rate[4];         //indicate directions of rate in the rectangle
//   double l1,l2,l;               //the length of a basic square size
//   double a;               //the transformed area
//   double lowest_time, t;
//   int k = 0;
//   double time1, time2;

//  l1 = l2 = l = sqrt(runParticle->real_radius*runParticle->real_radius*Pi);

//	if (l1 > l || l2 > l)
//		cout << "WTF" << endl;

//  l1 -= (runParticle->l[0] + runParticle->l[2]);
//  l2 -= (runParticle->l[1] + runParticle->l[3]);

//  for (i=0; i<4; i++)
//    rate[i] = 0.0;

//  //for (i=0; i<6; i++)
//  //  runParticle->dir_rate[i] *= -1;

//  //convert hexagonaly directed rates to rates relating to a square

//	rate[0] += (runParticle->dir_rate[0]);

//	rate[0] += cos(Pi/3) * runParticle->dir_rate[1];
//	rate[1] += cos(Pi/6) * runParticle->dir_rate[1];

//	rate[1] += cos(Pi/6) * runParticle->dir_rate[2];
//	rate[2] += cos(Pi/3) * runParticle->dir_rate[2];

//	rate[2] += runParticle->dir_rate[3];

//	rate[2] += cos(Pi/3) * runParticle->dir_rate[4];
//	rate[3] += cos(Pi/6) * runParticle->dir_rate[4];

//	rate[3] += cos(Pi/6) * runParticle->dir_rate[5];
//	rate[0] += cos(Pi/3) * runParticle->dir_rate[5];

//
////  if (runParticle->neig_of_two_spinel_grains) {
////  	cout << "debug" << endl;
////  }

//  //now calculate the time of growth to a common meeting point
//  if ( rate[0] > 0.0 && rate[2] > 0.0 ) {
//
//  	lowest_time = l1 / (rate[0]+rate[2]);
//
//// 	if ( (lowest_time * rate[0] > l) || (lowest_time * rate[2] > l) )
//// 		cout << "! to large" << endl;
//
//  }
//  else if ( rate[0] > 0.0 && rate[2] == 0.0 ) {
//
//	lowest_time = l1 / rate[0];
//
////	  if ( (lowest_time * rate[0] > l) || (lowest_time * rate[2] > l) )
////  	cout << "! to large" << endl;
//
//  }
//  else if ( rate[0] == 0.0 && rate[2] > 0.0 ) {
//
//  	lowest_time = l1 / rate[2];
//
////  	if ( (lowest_time * rate[0] > l) || (lowest_time * rate[2] > l) )
////  		cout << "! to large" << endl;
//
//  }
//
//  else if ( rate[0] > 0.0 && rate[2] < 0.0 ) {
//
//  	time1 = (l - runParticle->l[0]) / rate[0];
//
//  	time2 = - runParticle->l[2] / rate[2];
//
//  	if ( time1 > time2 )
//  		lowest_time = time1;
//  	else
//  		lowest_time = l1 / (rate[0]+rate[2]);
//
////  	if ( (lowest_time * rate[0] > l) || (lowest_time * rate[2] > l) )
////  		cout << "! to large" << endl;
//
//  }
//
//  else if ( rate[0] < 0.0 && rate[2] > 0.0 ) {
//
//  	time1 = (l - runParticle->l[2]) / rate[2];
//
//  	time2 = - runParticle->l[0] / rate[0];
//
//  	if ( time1 > time2 )
//  		lowest_time = time1;
//  	else
//  		lowest_time = l1 / (rate[0]+rate[2]);
//
//// 	if ( (lowest_time * rate[0] > l) || (lowest_time * rate[2] > l) )
////  		cout << "! to large" << endl;
//
//  }
//
//  else if ( rate[0] < 0.0 && rate[2] == 0.0 ) {
//  	lowest_time = 10.0 * time_interval;
//  }
//
//  else if ( rate[0] < 0.0 && rate[2] < 0.0 ) {
//  	lowest_time = 10.0 * time_interval;
//  }

//  else if ( rate[0] == 0.0 && rate[2] == 0.0 ) {
//  	lowest_time = 10.0 * time_interval;
//  }
//
//  else if ( rate[0] == 0.0 && rate[2] < 0.0 ) {
//  	lowest_time = 10.0 * time_interval;
//  }
//
//	/************************/
//
//  if ( rate[1] > 0.0 && rate[3] > 0.0 ) {
//
//  	t = l2 / (rate[1]+rate[3]);
//
////  	if ( (t * rate[1] > l) || (t * rate[3] > l) )
////  		cout << "! to large" << endl;
//
//  }
//  else if ( rate[1] > 0.0 && rate[3] == 0.0 ) {
//
//	t = l2 / rate[1];
//
////	if ( (t * rate[1] > l) || (t * rate[3] > l) )
////		cout << "! to large" << endl;
//
//  }
//  else if ( rate[1] == 0.0 && rate[3] > 0.0 ) {
//
//  	t = l2 / rate[3];
//
////  	if ( (t * rate[1] > l) || (t * rate[3] > l) )
////  		cout << "! to large" << endl;
//
//  }
//
//  else if ( rate[1] > 0.0 && rate[3] < 0.0 ) {
//
//  	time1 = (l - runParticle->l[1]) / rate[1];
//
//  	time2 = - runParticle->l[3] / rate[3];
//
//  	if ( time1 > time2 )
//  		t = time1;
//  	else
//  		t = l1 / (rate[1]+rate[3]);
//
////  	if ( (t * rate[1] > l) || (t * rate[3] > l) )
////  		cout << "! to large" << endl;
//
//  }
//
//  else if ( rate[1] < 0.0 && rate[3] > 0.0 ) {
//
//  	time1 = (l - runParticle->l[3]) / rate[3];
//
//  	time2 = - runParticle->l[1] / rate[1];
//
//  	if ( time1 > time2 )
//  		t = time1;
//  	else
//  		t = l1 / (rate[1]+rate[3]);
//
//// 	if ( (t * rate[1] > l) || (t * rate[3] > l) )
//// 		cout << "! to large" << endl;
//
//  }
//
//  else if ( rate[1] < 0.0 && rate[3] == 0.0 ) {
//  	t = 10.0 * time_interval;
//  }
//
//  else if ( rate[1] < 0.0 && rate[3] < 0.0 ) {
//  	t = 10.0 * time_interval;
//  }

//  else if ( rate[1] == 0.0 && rate[3] == 0.0 ) {
//  	t = 10.0 * time_interval;
//  }
//
//  else if ( rate[1] == 0.0 && rate[3] < 0.0 ) {
//  	t = 10.0 * time_interval;
//  }

//	/************************/
//
//  for (i=0; i<4; i++)
//  	runParticle->save_rate[i] = rate[i];
//
//  if ( t < 0 || lowest_time < 0 )
//  	cout << "doof" << endl;

//  if (t < lowest_time)
//    return t;
//  else
//  	return lowest_time;

//}

//// calculate the area, that will be overgrown during a timestep,
//// which is defined in the CalcTime()-function
//// called in gbm()
// double
// Min_Trans_Lattice::CalcArea(double lowest_time)
//{
//   int i;
//   double rate[4];        //indicate directions of rate in the rectangle
//   double l,l1,l2;        //the length of a basic square size
//   double x1, x2, F0, F1; //transformed area variables (F0=partikel-area,
//   F1=remaining area after timestep, x1,x2=new side-lengthes

//  F0 = runParticle->real_radius*runParticle->real_radius*Pi;
//  l = sqrt(F0);

//  l1 = l-runParticle->l[0]-runParticle->l[2];
//  l2 = l-runParticle->l[1]-runParticle->l[3];

//  for (i=0; i<4; i++)
//    rate[i] = 0.0;

//  //pass hexagonaly directed rates over to rates in a square
//  for (i=0; i<6; i++)
//    {
//      if (runParticle->dir_rate[i] < 0.0 && runParticle->mineral == 2)
//        runParticle->dir_rate[i] = 0.0;

//      switch (i)
//        {
//        case 0:
//          rate[0] += runParticle->dir_rate[0];
//          break;

//        case 1:
//          rate[0] += sin(Pi/6) * runParticle->dir_rate[1];
//          rate[1] += cos(Pi/6) * runParticle->dir_rate[1];
//          break;

//        case 2:
//          rate[1] += cos(Pi/6) * runParticle->dir_rate[2];
//          rate[2] += sin(Pi/6) * runParticle->dir_rate[2];
//          break;

//        case 3:
//          rate[2] += runParticle->dir_rate[3];
//          break;

//        case 4:
//          rate[2] += sin(Pi/6) * runParticle->dir_rate[4];
//          rate[3] += cos(Pi/6) * runParticle->dir_rate[4];
//          break;

//        case 5:
//          rate[3] += cos(Pi/6) * runParticle->dir_rate[5];
//          rate[0] += sin(Pi/6) * runParticle->dir_rate[5];
//          break;
//        }

//    }

//  runParticle->l[0] += lowest_time*rate[0];
//  runParticle->l[1] += lowest_time*rate[1];
//  runParticle->l[2] += lowest_time*rate[2];
//  runParticle->l[3] += lowest_time*rate[3];

//}

////this function calculates the time, that it takes to overgrow the
////area of a particle, depending on the number of neighbours (the
////contact area) and the rates into this direction
////idea: calculate the time for either two sides to meet at a certain
////point, choose the lowest calculated time
//// called in gbm()
// double
// Min_Trans_Lattice::CalcTime()
//{
//   int i,j, buffer1, buffer2;
//   double rate[4];         //indicate directions of rate in the rectangle
//   double l1,l2;               //the length of a basic square size
//   double a;               //the transformed area
//   double lowest_time, t;

//  l1 = l2 = sqrt(runParticle->real_radius*runParticle->real_radius*Pi);

//  l1 -= (runParticle->l[0] + runParticle->l[2]);
//  l2 -= (runParticle->l[1] + runParticle->l[3]);

//  for (i=0; i<4; i++)
//    rate[i] = 0.0;

//  //for (i=0; i<6; i++)
//    //  runParticle->dir_rate[i] *= -1;

//  //pass hexagonaly directed rates over to rates in a square
//  for (i=0; i<6; i++)
//    {

//      //if (runParticle->dir_rate[i] < 0.0 && runParticle->mineral == 2)
//	  //runParticle->dir_rate[i] = 0.0;

//      if (isnan(runParticle->dir_rate[i]))
//        cout << "isnogud" << endl;

//      switch (i)
//        {
//        case 0:
//          rate[0] += runParticle->dir_rate[0];
//          break;

//        case 1:
//          rate[0] += sin(Pi/6) * runParticle->dir_rate[1];
//          rate[1] += cos(Pi/6) * runParticle->dir_rate[1];
//          break;

//        case 2:
//          rate[1] += cos(Pi/6) * runParticle->dir_rate[2];
//          rate[2] += sin(Pi/6) * runParticle->dir_rate[2];
//          break;

//        case 3:
//          rate[2] += runParticle->dir_rate[3];
//          break;

//        case 4:
//          rate[2] += sin(Pi/6) * runParticle->dir_rate[4];
//          rate[3] += cos(Pi/6) * runParticle->dir_rate[4];
//          break;

//        case 5:
//          rate[3] += cos(Pi/6) * runParticle->dir_rate[5];
//          rate[0] += sin(Pi/6) * runParticle->dir_rate[5];
//          break;
//        }

//    }

//  //now calculate the time of growth to a common meeting point
//  lowest_time = l1 / (rate[0]+rate[2]);
//  if (lowest_time < 0)
//    lowest_time = 10*time_interval;

//  t = l2 / (rate[1]+rate[3]);
//  if (t<0)
//    t = 10*time_interval;

//  if (t < lowest_time)
//    lowest_time = t;

//  //cout << lowest_time << endl;

//  return (lowest_time);
//}

// calculates the latent-heat realease after a phase-change. called in
// fctn "mineral_transformation"

// Akaogie, 1989, does not take into account the effect of pressure on the heat
// capacity

// returns the NEW TEMPERATURE (NOT the Delta T)
double Min_Trans_Lattice::Latent_Heat_Release(Particle *prtcl, int type) {
  int i, j;

  double C_OL, C_SPIN, T, P, V0_OL, V0_SPIN, DeltaC, ALPHA_OL, ALPHA_SPIN, V_OL,
      V_SPIN;
  double DeltaH0, DeltaH_TP, DeltaV_T, DeltaH_T, DeltaT, area, factor1, factor2;
  double initial_c, new_c;

  T = prtcl->temperature;
  P = ((prtcl->sxx + prtcl->syy) / 2.0) * pressure_scale * pascal_scale;

  if (runParticle->newly_overgrown != 0.0) {

    C_OL = 217.8 - 1408 * pow(T, -0.5) - 5.264e8 * pow(T, -3.0);
    C_SPIN = 216.4 - 1472 * pow(T, -0.5) - 5.579e8 * pow(T, -3.0);
    DeltaC = C_SPIN - C_OL;

    V0_OL = 0.00004367;
    V0_SPIN = 0.00004051;
    ALPHA_OL = 0.00003052 + 8.504e-9 * T - 0.5824 * pow(T, -2.0);
    ALPHA_SPIN = 0.00002711 + 6.885e-9 * T - 0.5767 * pow(T, -2.0);
    V_OL = V0_OL * exp(ALPHA_OL * (T - 298.0));
    V_SPIN = V0_SPIN * exp(ALPHA_SPIN * (T - 298.0));
    DeltaV_T = V_SPIN - V_OL;

    area = runParticle->real_radius * runParticle->real_radius * Pi;
    factor1 =
        ((runParticle->overgrown_area - runParticle->newly_overgrown) / area);
    factor2 = ((runParticle->overgrown_area) / area);

    switch (type) {

    case OL_SPIN:

      initial_c = C_SPIN * factor1 + C_OL * (1.0 - factor1);
      new_c = C_SPIN * factor2 + C_OL * (1.0 - factor2);

      DeltaH0 = -29970;

      DeltaH_T = DeltaH0 + DeltaC * (T - 975.0);

      DeltaH_TP = DeltaH_T + DeltaV_T * (P - 1);

      DeltaH_TP = DeltaH_TP * (prtcl->newly_overgrown / area);

      DeltaT = ((T * initial_c + fabs(DeltaH_TP)) / new_c) - T;

      DeltaT = fabs(DeltaT);

      if (prtcl->newly_overgrown < 0)
        DeltaT *= -1.0;

      break;

    case SPIN_OL:

      initial_c = C_OL * factor1 + C_SPIN * (1.0 - factor1);
      new_c = C_OL * factor2 + C_SPIN * (1.0 - factor2);

      DeltaH0 = -29970;

      DeltaH_T = DeltaH0 + DeltaC * (T - 975.0);

      DeltaH_TP = DeltaH_T + DeltaV_T * (P - 1);

      DeltaH_TP = DeltaH_TP * (prtcl->newly_overgrown / area);

      DeltaT = ((T * initial_c - fabs(DeltaH_TP)) / new_c) - T;

      DeltaT = fabs(DeltaT);

      if (prtcl->newly_overgrown > 0)
        DeltaT *= -1.0;

      break;
    }
  } else
    DeltaT = 0.0;

  return (DeltaT);
}

// called in reaction_two(). changes boundary-particles
void Min_Trans_Lattice::Provisoric_Mineral_Transformation_No_Spin_Distri() {

  int i, j, k;
  double factor;

  runParticle = &refParticle;

  for (i = 0; i < numParticles; i++) {

    runParticle->previous_mineral = runParticle->mineral;
    runParticle->previous_mV = runParticle->mV; // save the molecular volume
    runParticle->previous_young = runParticle->young; // save young
    runParticle->previous_grain = runParticle->grain;
    runParticle->previous_radius = runParticle->radius;
    runParticle->previous_viscosity = runParticle->viscosity;
    runParticle->previous_real_radius = runParticle->real_radius;
    runParticle->previous_area = runParticle->area;

    for (j = 0; j < 6; j++) {
      runParticle->previous_springf[j] = runParticle->springf[j];
      runParticle->previous_rad[j] = runParticle->rad[j];
    }

    runParticle = runParticle->nextP;
  }

  runParticle = &refParticle;

  // set new grain / new mineral
  for (i = 0; i < numParticles; i++) {
    if (runParticle->phase_change_required) {

      switch (runParticle->mineral) {
      case 1:

        runParticle->mineral = 2; // change mineral index

        runParticle->viscosity *= 1.2;

        if (!use_homogeneous_spin)
          runParticle->young = runParticle->young * 1.4;
        else
          runParticle->young =
              runParticle->spin_young; // change repulsion constant

        break;

      case 2:

        runParticle->mineral = 1; // change mineral index0

        runParticle->viscosity /= 1.2;

        runParticle->young =
            runParticle->olivine_young; // change repulsion constant

        break;
      }
    }

    runParticle = runParticle->nextP;
  }

  for (i = 0; i < numParticles; i++) {
    if (runParticle->grain != beware_grain1 &&
        runParticle->grain != beware_grain2 &&
        runParticle->grain != beware_grain3) {
      runParticle->mV = Adjust_Molare_Volume();
      runParticle->previous_real_radius = runParticle->real_radius;
      runParticle->real_radius =
          pow(((3.0 * runParticle->mole_per_particle * runParticle->mV) /
               (4.0 * 3.1415)),
              (1.0 / 3.0)); // and change radius
      factor = (runParticle->real_radius / runParticle->previous_real_radius);

      runParticle->radius *= factor;
      for (j = 0; j < 6; j++) {
        if (runParticle->neigP[j])
          runParticle->rad[j] *= factor;
      }

      runParticle->area = pow(double(runParticle->real_radius), 2.0) * 3.1415;
    }

    runParticle = runParticle->nextP;
  }

  AdjustParticleConstants();
}

// set heat distribution in class heat_lattice by
// passing x/y-position of particles + their temperatures
// performed only in the beginning for all particles and
// if a particle actually changed the phase.
// called in heat_flow
void Min_Trans_Lattice::Set_Heat_Lattice(Particle *prtl) {
  double delta_x, delta_y, Time_Step_Heat;

  Time_Step_Heat = 5.0; // dummy-value

  runParticle = &refParticle;
  runParticle = runParticle->prevP;
  delta_y = runParticle->ypos / org_wall_pos_y;

  runParticle = runParticle->prevP;
  delta_x = runParticle->xpos / org_wall_pos_x;

  heat_distribution.initial_diameter = initial_diameter;

  heat_distribution.set_Parameters(delta_y, delta_x, Time_Step_Heat);

  runParticle = prtl;

  heat_distribution.set_cell_temperature(runParticle->xpos, runParticle->ypos,
                                         runParticle->temperature);
}

// read the actual temperature-distribution, apply it
// to the particles. called in "heat-flow()" function
void Min_Trans_Lattice::Read_Heat_Lattice() {
  int i;

  for (i = 0; i < numParticles; i++) {
    runParticle->temperature = heat_distribution.set_particle_temperature(
        runParticle->xpos, runParticle->ypos);
    // cout << "temperature of particle: " << runParticle->temperature << endl;
    runParticle = runParticle->nextP;
  }
}

// "main"-function for the heat flow
// is called in the experiment class or
// - if turned on - in the gbm-routine
void Min_Trans_Lattice::Heat_Flow(double Time_Step_Heat) {
  int p_counter, i, col;
  double delta_t, wall_pos, delta_x, delta_y, density, specific_heat, lambda,
      highest_temp;

  if (Time_Step_Heat > 0.0) {

    runParticle = &refParticle;
    runParticle = runParticle->prevP;
    delta_y = runParticle->ypos / org_wall_pos_y;

    runParticle = runParticle->prevP;
    delta_x = runParticle->xpos / org_wall_pos_x;

    heat_distribution.initial_diameter = initial_diameter;

    heat_distribution.set_Parameters(delta_y, delta_x, Time_Step_Heat);

    // perform heat flow, read and set particle-properties
    // Set_Heat_Lattice();		//temperaturen der partikel uebergeben
    heat_distribution.heat_flow(); // ausfhren der waermeleitungs-routine
    Read_Heat_Lattice(); // rueckuebergabe der temperaturen an partikel
  }
}

// A function to adjust the molare volume, using the
// temperature. It does NOT include the compressibility.
// Called after a phase-change of a particle.
double Min_Trans_Lattice::Adjust_Molare_Volume() {

  double alpha, mean_stress, mV, t, delta, gamma;

  mean_stress = ((-runParticle->sxx) + (-runParticle->syy)) / 2;
  mean_stress *= pascal_scale * pressure_scale;

  switch (runParticle->mineral) {
  case 1:
    t = runParticle->temperature;
    // after temperature
    // runParticle->temperature = 1000.0;
    alpha = 0.00003052 * t + 0.000000004 * pow(double(t), 2.0) +
            0.5824 * pow(double(t), -1.0);
    t = 298;
    alpha -= 0.00003052 * t + 0.000000004 * pow(double(t), 2.0) +
             0.5824 * pow(double(t), -1.0);

    mV = 0.00004367 * exp(alpha);

    // after pressure
    // adjustment occurs using the bulk-modulus
    // t = runParticle->temperature;  //formel nach °C
    // mV -= (1.0/(135000000000.0 + 16000000.0 * (t))) * mV *
    // fabs(mean_stress);//s.Fei et al.

    break;

  case 2:
    t = runParticle->temperature;
    // after temperature
    // runParticle->temperature = 1000.0;

    alpha = 0.00002711 * t + 0.00000000344 * pow(double(t), 2.0) +
            0.5767 * pow(double(t), -1.0);
    t = 298;
    alpha -= 0.00002711 * t + 0.00000000344 * pow(double(t), 2.0) +
             0.5767 * pow(double(t), -1.0);

    mV = 0.00003953 * exp(alpha);

    // after pressure
    // adjustment occurs using the bulk-modulus
    // t = runParticle->temperature - 27;//formel nach Raumtemperatur in °K
    // mV -= (1.0/(174000000000.0 - 27000000.0 * t)) * mV *
    // fabs(mean_stress);//s.Fei et al.

    break;

  case 3:
    // pyroxen
    alpha = 32.2e-6;
    mV = 3.1292e-5 * (1.0 + alpha * (runParticle->temperature - 293.0));
    break;

  case 4:
    // pyrop
    alpha = 2.36e-5;
    delta = 6.88;

    alpha = alpha *
            (1.0 / (1 - delta * alpha * (runParticle->temperature - 293.0)));

    mV = 1.13163e-4 * (1.0 + alpha * (runParticle->temperature - 293.0));

    break;
  }

  mV = fabs(mV);
  return (mV);
}

// A function to adjust the molare volume, using the
// temperature. It does NOT include the compressibility.
// Called after a phase-change of a particle.
double Min_Trans_Lattice::Molare_Volume(int mineral) {

  double alpha, mean_stress, mV, t;

  mean_stress = ((-runParticle->sxx) + (-runParticle->syy)) / 2;
  mean_stress *= pascal_scale * pressure_scale;

  switch (mineral) {
  case 1:
    t = runParticle->temperature;
    // after temperature
    // runParticle->temperature = 1000.0;
    alpha = 0.00003052 * t + 0.000000004 * pow(double(t), 2.0) +
            0.5824 * pow(double(t), -1.0);
    t = 298;
    alpha -= 0.00003052 * t + 0.000000004 * pow(double(t), 2.0) +
             0.5824 * pow(double(t), -1.0);

    mV = 0.00004367 * exp(alpha);

    // after pressure
    // adjustment occurs using the bulk-modulus
    // t = runParticle->temperature;  //formel nach °C
    // mV -= (1.0/(135000000000.0 + 16000000.0 * (t))) * mV *
    // fabs(mean_stress);//s.Fei et al.

    break;

  case 2:
    t = runParticle->temperature;
    // after temperature
    // runParticle->temperature = 1000.0;

    alpha = 0.00002711 * t + 0.00000000344 * pow(double(t), 2.0) +
            0.5767 * pow(double(t), -1.0);
    t = 298;
    alpha -= 0.00002711 * t + 0.00000000344 * pow(double(t), 2.0) +
             0.5767 * pow(double(t), -1.0);

    mV = 0.00003953 * exp(alpha);

    // after pressure
    // adjustment occurs using the bulk-modulus
    // t = runParticle->temperature - 27;//formel nach Raumtemperatur in °K
    // mV -= (1.0/(174000000000.0 - 27000000.0 * t)) * mV *
    // fabs(mean_stress);//s.Fei et al.

    break;
  }

  mV = fabs(mV);
  return (mV);
}

// The activation energy used for the calculation of the
// rate of reaction. Currently it is fixed to a certain value,
// commenting this out will calculate the activation energy depending
// including the activation-volume.
// The activation energy used for the calculation of the
// rate of reaction. Currently it is fixed to a certain value,
// commenting this out will calculate the activation energy depending
// including the activation-volume.
double Min_Trans_Lattice::Activation_Energy(Particle *prtcl) {
  double activation_energy, ActV, alpha, Tm, deltaH1, deltaH2, mean_stress;

  // default for set_act_energy == false
  if (!use_dynamic_activation_energy) {

    // set in experiment.cc, default in constructor: 350000
    activation_energy = static_activation_energy;

  } else {

    mean_stress = -(prtcl->sxx + prtcl->syy) / 2.0;

    mean_stress *= pressure_scale;

    ActV = 1.19e4;
    ActV = ActV * (1.0 / (1.0 + mean_stress / (1e9 * 8.0)));

    alpha = 159.0;
    Tm = 2163.0;

    deltaH1 = alpha * Tm;

    deltaH2 = (mean_stress / 1e9) * ActV;

    activation_energy = deltaH2 + deltaH1;
  }

  return (activation_energy);
}

// presets an ectivation-energy and the pressure-barrier where gbm starts
// for olivine->spinel!
void Min_Trans_Lattice::SetReactions(double energy, double barrier) {
  act_energy = energy;
  set_act_energy = true;

  pressure_barrier = barrier;
  set_pressure_barrier = true;

  heat_distribution.SetHeatFlowParameters(4.2, 1005.0, 0.000001, 1000.0);
}

// changes the phase of complete grain as basis for grain-growth
// works with up to 3 grains
void Min_Trans_Lattice::Change_Grain_Mineral(int grain1, int grain2,
                                             int grain3) {
  // aufruf in main():
  // uses grain-number to identify, changes phase to spinel

  int i;
  double factor;

  beware_grain1 = grain1;
  beware_grain2 = grain2;
  beware_grain3 = grain3;

  //**********************************************
  // define grain:
  //**********************************************

  for (i = 0; i < numParticles; i++) {
    if (runParticle->grain == beware_grain1 ||
        runParticle->grain == beware_grain2 ||
        runParticle->grain == beware_grain3) {
      runParticle->mineral = 2;
      runParticle->young *= 1.4;

      runParticle->mV = Adjust_Molare_Volume();
      runParticle->previous_real_radius = runParticle->radius;
      runParticle->real_radius =
          pow(((3.0 * runParticle->mole_per_particle * runParticle->mV) /
               (4.0 * 3.1415)),
              (1.0 / 3.0)); // and change radius
      factor = (runParticle->real_radius / runParticle->previous_real_radius);
      runParticle->radius *= factor;

      runParticle->area = pow(double(runParticle->real_radius), 2.0) * 3.1415;
    }

    runParticle = runParticle->nextP;
  }
}

/*****************************************************
 * Convenience-function, to approach the relevant pressure
 * with large deformation- and time-steps and to continue
 * with smaller steps.
 *****************************************************/
void Min_Trans_Lattice::Change_Timestep() {
  int i;

  if (timeflag && logicalflag) {
    time /= 20.0;
    logicalflag = false;
  }
}

void Min_Trans_Lattice::Exchange_Probabilities() {

  // used for the addition of probabilities to change phase

  int i;

  for (i = 0; i < numParticles; i++) {
    if (runParticle->prob < 1.0) {
      runParticle->previous_spinel_content +=
          pow(double(runParticle->rate * time_interval), 3);
      // cout << "prev. prob: " << runParticle->previous_prob << endl;
    } else {
      runParticle->previous_spinel_content = runParticle->spinel_content;
    }

    runParticle = runParticle->nextP;
  }
}

void Min_Trans_Lattice::Change_Young() {
  if (runParticle->mineral == 1) {
    runParticle->young *= 1.4;
  } else {
    runParticle->young *= 1.0 / 1.4;
  }
}

void Min_Trans_Lattice::Restore_Young() {
  if (runParticle->mineral == 1) {
    runParticle->young *= 1.0 / 1.4;
  } else {
    runParticle->young *= 1.4;
  }
}

/**********************************************************************************************************
******************************* Functions for Nucleation, not in use
***************************************
**********************************************************************************************************/

// function to nucleate. still experimental
void Min_Trans_Lattice::Nucleation() {
  int i, j;
  double rate, I, Y, v, lowest_time, buffer, prob;

  // voreinstellungen, vielleicht sind nicht alle davon nötig...
  // damit der partikel beim nucleation_list_pointer nicht mitgezählt wird
  for (i = 0; i < numParticles; i++) {
    runParticle->prob = 0;
    runParticle = runParticle->nextP;
  }

  runParticle = &refParticle;

  lowest_time = actual_time; // zeit nach wärmeleitung

  nucleation_occured = false;

  help_pointer = NULL;

  for (i = 0; i < numParticles; i++) {
    if (!runParticle->is_lattice_boundary && runParticle->mineral == 1) {
      Y = Growth_Rate_For_Nucleation();
      I = Nucleation_Rate();

      v = pow(double(actual_time), 4.0) * pow(double(Y), 3.0) * I * Pi / 3.0;
      runParticle->rate = pow(double(Y), 3.0) * I;

      // standardisiert für ein mol, daher kein faktor
      runParticle->spinel_content = v;
    }
    runParticle = runParticle->nextP;
  }

  // sortiere nukleiierte partikel entsprechend der hoechsten volumina
  Find_Local_Maxima();

  for (i = 0; i < numParticles; i++) {
    if (runParticle->mineral == 1 && !runParticle->is_lattice_boundary) {
      if (runParticle->prob >= 1.0 && runParticle != nucleus) {
        buffer = pow((runParticle->mole_per_particle * runParticle->mV) /
                         (runParticle->rate * Pi) * 3.0,
                     0.25);

        if (buffer < actual_time) {
          nucleation_occured = true;
          runParticle->time = buffer;
        } else
          runParticle->time = actual_time;
      } else
        runParticle->time = actual_time;
    } else
      runParticle->time = actual_time;

    runParticle = runParticle->nextP;
  }

  if (!nucleation_occured)
    help_pointer = NULL;

  if (nucleation_occured) {
    if (!nucleus) {
      time = actual_time;

      for (i = 0; i < numParticles; i++) {
        if (runParticle->time < time && runParticle->time > 0.0 &&
            runParticle->prob > 1.0) {
          nucleus = runParticle;
          time = runParticle->time;
        }
        runParticle = runParticle->nextP;
      }
    }

    time = actual_time;

    for (i = 0; i < numParticles; i++) {
      if (runParticle->time < time && runParticle != nucleus &&
          runParticle->time > 0.0 && runParticle->prob > 1.0) {
        help_pointer = runParticle;
        time = runParticle->time;
      }
      runParticle = runParticle->nextP;
    }

    if (help_pointer) {
      time_interval = help_pointer->time - nucleus->time;
      actual_time -= time_interval;
    } else {
      time_interval = 0;
      actual_time -= nucleus->time;
    }
  }

  if (help_pointer) {
    for (i = 0; i < numParticles; i++) {
      if (runParticle != nucleus) {
        runParticle->prob = 0.0;
      }
      runParticle = runParticle->nextP;
    }

    // nucleus = help_pointer;

    nucleus = NULL;

    // Mineral_Transformation();
    Mineral_Transformation_No_Spin_Distri();

    AdjustConstantGrainBoundaries();

    Make_Phase_Boundaries();

    Relaxation();
  }

  UpdateElle();
}

// calculates energy barrier for nucleation

double Min_Trans_Lattice::Energy_Barrier() {
  double prefactor, surface_energy, driving_force, molare_volume,
      energy_barrier, mean_stress;

  mean_stress = ((-runParticle->sxx) + (-runParticle->syy)) / 2;
  mean_stress *= pascal_scale * pressure_scale;

  if (runParticle->mineral == 1) {
    molare_volume = 0.00004367;
  } else if (runParticle->mineral == 2) {
    molare_volume = 0.000039;
  }

  prefactor = (16.0 * Pi) / 3.0;
  surface_energy = 0.6;
  driving_force = Driving_Force();

  surface_energy = pow(double(surface_energy), 3.0);
  molare_volume = pow(double(molare_volume), 2.0);
  driving_force = pow(double(driving_force), 2.0);

  energy_barrier =
      prefactor * ((surface_energy * molare_volume) / driving_force);

  return (energy_barrier);
}

double Min_Trans_Lattice::Excess_Pressure() {

  double delta_P, // differenz to the equilibrium pressure for the given
                  // temperature (akaogie, daessler)
      P, // at phase boundary with temperature T
      T, mean_stress;

  double clapeyron_constant, enthalpy_differenz, volume_differenz;

  enthalpy_differenz = 29970.0;
  volume_differenz = 0.0000316; // volumen ist KORRIGIERT fr J/(pa*mol)!!!!!
                                // vgl. stosch, akaogie table 3
  T = runParticle->temperature;

  P = 14400000000.0 + (enthalpy_differenz / volume_differenz) *
                          log(T / 1473); // P(1) und T(1) = 14.4 GPa, 1473 K

  mean_stress = ((-runParticle->sxx) + (-runParticle->syy)) / 2;
  mean_stress *= pascal_scale * pressure_scale;

  delta_P = P - mean_stress;

  // cout << "deltaP: " << delta_P << endl;
  // cout << "mean stress: " << mean_stress << endl;

  return (delta_P);
}

double Min_Trans_Lattice::Undercooling() {

  double delta_T, // differenz to the equilibrium temperature for the given
                  // pressure (akaogie, daessler)
      P, // at phase boundary with temperature T
      T, mean_stress;

  double clapeyron_constant, enthalpy_differenz, volume_differenz;

  enthalpy_differenz = 29970.0;
  volume_differenz = 0.00000316; // volumen ist KORRIGIERT fr J/(pa*mol)!!!!!
                                 // vgl. stosch, akaogie table 3
  // TODO: volumenkorrektur auf beta-phasen umwandeln!!!

  mean_stress = ((-runParticle->sxx) + (-runParticle->syy)) / 2;
  mean_stress *= pascal_scale * pressure_scale;

  T = 1473 * exp((mean_stress - 14400000000.0) *
                 (volume_differenz / enthalpy_differenz));

  delta_T = T - runParticle->temperature;

  // cout << "deltaT: " << delta_T << endl;

  return (delta_T);
}

double Min_Trans_Lattice::Driving_Force() {
  double driving_force, entropy_change, mean_stress, mean_stress_change,
      volume_change, excess_pressure, undercooling;

  excess_pressure = Excess_Pressure();
  undercooling = Undercooling();

  // cout << "excess P: " << excess_pressure << " Delta T: " << undercooling <<
  // endl;

  mean_stress = ((-runParticle->sxx) + (-runParticle->syy)) / 2;
  mean_stress *= pascal_scale * pressure_scale;

  if (runParticle->mineral == 1) {
    volume_change = 0.00000316; // hier: v in m^3
    entropy_change = 7.7;
  } else if (runParticle->mineral == 2) {
    volume_change = 0.00000316; // hier: v in m^3
    entropy_change = 7.7;
  }

  // runParticle->mV * mole_per_particledriving_force = excess_pressure *
  // volume_change - undercooling * entropy_change;
  driving_force =
      mean_stress * volume_change - runParticle->temperature * (entropy_change);

  // cout <<"driving force: " << driving_force << endl;

  return (driving_force);
}

double Min_Trans_Lattice::Growth_Rate_For_Nucleation() {
  double rate, mean_stress, kelvin;

  kelvin = runParticle->temperature;

  mean_stress = mean_stress = ((-runParticle->sxx) + (-runParticle->syy)) / 2;
  mean_stress *= pascal_scale * pressure_scale;

  runParticle->delta_pV = Driving_Force();

  rate = 162000 * kelvin *
         exp(-(Activation_Energy(runParticle) / (gas_const * kelvin)));
  rate *= (1 - exp(-(runParticle->delta_pV) / (gas_const * kelvin)));

  return (rate);
}

double Min_Trans_Lattice::Nucleation_Rate() {
  // eEL UND SURFe ÜERARBEITEN: EINHEITEN
  double k, rate, deltaG_hom, kelvin;

  double I_zero;

  int i;

  kelvin = runParticle->temperature;

  I_zero = 1.0e+40;

  k = 1.38 * pow(10.0, -23); // boltzmann-constant

  deltaG_hom = Energy_Barrier(); // geeicht auf spinel-nukleiierung in olivin!!
  // cout << "Energy Barrier: " << deltaG_hom << endl;

  rate = I_zero * kelvin * exp(-(450000 / (gas_const * kelvin)));
  // cout << "rate: " << rate << endl;
  // cout << "deltaG_hom / (k * kelvin): " << deltaG_hom / (k * kelvin) << endl;

  rate *= exp(-(deltaG_hom / (k * kelvin)));
  // cout << "rate: " << rate << endl;

  return (rate);
}

void Min_Trans_Lattice::Find_Local_Maxima() {
  int i, j;

  for (i = 0; i < numParticles; i++) { // loop through particles and find out
                                       // local maxima among remaining particles
    if (runParticle->spinel_content > 0.0) {
      for (j = 0; j < 8; j++) {
        if (runParticle->neigP[j] &&
            runParticle->neigP[j]->spinel_content > 0.0) {
          if (runParticle->neigP[j]->spinel_content >
              runParticle->spinel_content) {
            runParticle->local_maximum = false;
            break;
          } else {
            runParticle->local_maximum = true;
          }
        }
      }
    }
    runParticle = runParticle->nextP;
  }

  for (i = 0; i < numParticles; i++) {
    if (!runParticle->local_maximum) {
      runParticle->spinel_content = runParticle->previous_spinel_content;
      runParticle->prob = 0.0;
    } else
      runParticle->prob = 2.0;

    runParticle = runParticle->nextP;
  }
}

/**********************************************************************************************************
************************************ Additional transition-stuff
*******************************************
**********************************************************************************************************/

// dumps diff-stress and time for the last overgrown particle
void Min_Trans_Lattice::DumpDiffStressAndTime(Particle *particle, double time,
                                              double eEl_diff) {
  FILE *stat;
  double sxx, syy, sxy, differential, smax, smin;

  sxx = particle->sxx;
  sxy = particle->sxy;
  syy = particle->syy;

  smax = ((sxx + syy) / 2.0) +
         sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) + sxy * sxy);
  smin = ((sxx + syy) / 2.0) -
         sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) + sxy * sxy);
  differential = smax - smin;

  stat = fopen("Transition", "a");

  fprintf(stat, " differential stress:");
  fprintf(stat, " %f", differential);
  fprintf(stat, " time:");
  fprintf(stat, " %f\n", time);
  fprintf(stat, " elast. energy diff.:");
  fprintf(stat, " %f\n", eEl_diff);

  fclose(stat); // close file
}

// gaussian distribution for the youngs-moduli of grains, additional function
void Min_Trans_Lattice::SetGaussianYoungDistribution_2(double g_mean,
                                                       double g_sigma) {

  double prob;     // probability from gauss
  double k_spring; // pict spring
  int i, j;        // counters
  double ran_nb;   // rand number

  // grain_counter counter how many grain there were intially
  // times two is for security, not all grains may be there
  // this can certainly be made nicer

  runParticle = &refParticle;

  // initialize random-number generator:
  srand(std::time(0));

  for (i = 0; i < numParticles; i++) // loop through the particles now
  {
    do {
      k_spring = rand() /
                 (double)RAND_MAX; // pic a pseudorandom double between 0 and 1

      //----------------------------------------------------------------
      // now convert this distribution to a distribution from
      // zero minus 8 times sigma to plus 8 times sigma
      //-----------------------------------------------------------------

      k_spring = (k_spring - 0.5) * 8.0 * (8.0 * g_sigma);

      //--------------------------------------------------------------
      // and shift it to mean plus minus 2 times sigma
      //--------------------------------------------------------------

      k_spring = k_spring + g_mean;

      //---------------------------------------------------------------
      // now apply gauss function to determine a probability for this
      // reaction constant to occur
      //---------------------------------------------------------------

      prob = (1 / (g_sigma *
                   sqrt(2.0 * 3.1415927))); // part one (right part) of function

      prob = prob * (exp(-0.5 * (((k_spring)-g_mean) / g_sigma) *
                         (((k_spring)-g_mean) / g_sigma))); // rest

      //--------------------------------------------------------------
      // now adjust probablity to run from 0 to 1.0
      //--------------------------------------------------------------

      prob = prob * sqrt(2.0 * 3.1415927) * g_sigma;

      //---------------------------------------------------------------
      // pic the second pseudorandom number
      //---------------------------------------------------------------

      ran_nb = rand() / (double)RAND_MAX;

      //--------------------------------------------------------------
      // if the number picted is smaller or the same as the probability
      // from the gauss function accept the rate constant, if not go on
      // looping.
      //-------------------------------------------------------------------

      if (ran_nb <= prob) // if smaller
      {
        if (k_spring <= 0.0) // if not 0.0 or smaller (dont want that)
          ;
        else {
          // set the constant in the list !
          runParticle->rate_factor = k_spring;

          break;
        }
      }
    } while (1); // loop until break

    runParticle = runParticle->nextP;
  }

  for (i = 0; i < numParticles; i++) {
    runParticle->young *= runParticle->rate_factor;

    for (j = 0; j < 6; j++) {
      if (runParticle->springf[j])
        runParticle->springf[j] *= runParticle->rate_factor;
    }

    runParticle = runParticle->nextP;
  }
}

// convieniance function for heat flow-example
void Min_Trans_Lattice::HeatGrain(int T, int nb) {
  int i;

  runParticle = &refParticle;

  for (i = 0; i < numParticles; i++) {
    if (runParticle->grain == nb) {
      runParticle->temperature = T;
      Set_Heat_Lattice(runParticle);
    }

    runParticle = runParticle->nextP;
  }
}

// convieniance function for heat flow-example
void Min_Trans_Lattice::SetHeatLatticeHeatFlowExample() {
  double delta_x, delta_y, Time_Step_Heat;
  int i;

  Time_Step_Heat = 0.5;

  runParticle = &refParticle;
  runParticle = runParticle->prevP;
  delta_y = runParticle->ypos / org_wall_pos_y;

  runParticle = runParticle->prevP;
  delta_x = runParticle->xpos / org_wall_pos_x;

  heat_distribution.initial_diameter = initial_diameter;

  heat_distribution.set_Parameters(delta_y, delta_x, Time_Step_Heat);

  HeatGrain(1500, 3);

  HeatGrain(1500, 5);

  runParticle = &refParticle;

  for (i = 0; i < numParticles; i++) {
    heat_distribution.set_cell_temperature(runParticle->xpos, runParticle->ypos,
                                           runParticle->temperature);
    //       if (runParticle->temperature == 1500)
    //       	cout << "1500" << endl;

    runParticle = runParticle->nextP;
  }
}

// here, the
void Min_Trans_Lattice::SetProvisoricElasticProperties() {

  int i, j, rand_nb, counter;
  double young_buffer;
  vector<Particle *> spin_neig;

  runParticle = &refParticle;

  srand(std::time(0));

  // set back to basic mode, olivine can stay as it is
  for (i = 0; i < numParticles; i++) {

    runParticle->spin_young = runParticle->young;

    runParticle = runParticle->nextP;
  }

  // just to be on the safe side:
  Make_Phase_Boundaries();

  for (i = 0; i < numParticles; i++) {

    if (runParticle->is_phase_boundary) {
      if (runParticle->mineral == 1) {

        spin_neig.clear();

        for (j = 0; j < 6; j++) {
          if (runParticle->neigP[j]) {
            if (runParticle->neigP[j]->mineral == 2) {

              spin_neig.push_back(runParticle->neigP[j]);
            }
          }
        }

        // take the list from above and find a random-entry // assign to
        // olivine-particle this is now only for the provisoric mode...

        // this assigns a random spinel-neighbour-youngs modulus to the spinel
        // young

        /*
        rand_nb = rand() % spin_neig.size();
        runParticle->spin_young = spin_neig[rand_nb]->young;
        */

        // this assigns the average spinel-neighbour-youngs modulus to the
        // spinel young

        counter = 0;
        runParticle->spin_young = 0.0;
        for (j = 0; j < spin_neig.size(); j++) {
          counter++;
          runParticle->spin_young += spin_neig[j]->young;
        }
        runParticle->spin_young /= double(counter);
      }
    }
    runParticle = runParticle->nextP;
  }
}

void Min_Trans_Lattice::SetTimestep(int i) { timestep = i; }

void Min_Trans_Lattice::Insert_Hardbodies(vector<vector<double>> hardbodies) {

  runParticle = &refParticle;

  for (int i = 0; i < numParticles; i++) {
    for (int j = 0; j < hardbodies.size(); j++) {

      if (runParticle->grain == (int)hardbodies[j][0]) {

        runParticle->no_reaction = true;
        runParticle->no_nucleation = true;
        runParticle->young *= hardbodies[j][1];

        for (int k = 0; k < 6; k++) {
          runParticle->springf[k] *= hardbodies[j][1];
        }

        runParticle->viscosity *= hardbodies[j][2];

        runParticle->phase = int(hardbodies[j][3]);
        runParticle->mineral = int(hardbodies[j][3]);

        runParticle->mV = hardbodies[j][4];

        Adjust_Molare_Volume();

        break;
      }
    }
    runParticle = runParticle->nextP;
  }

  AdjustParticleConstants();
}

double Min_Trans_Lattice::growth_rate(double lowest_time, Particle *prtcl) {

  ofstream rate;
  int i, j, counter1, counter2;
  int spin_gr_nb = 0;
  int divisor = 0;
  double rad = 0, dir1, dir2;
  double average_normal_stress = 0.0;

  rate.open("Rate.txt", ios::out | ios::app);

  // stupid variable names, countercheck
  double total_area = 0.0;
  double particle_area = sqrt(pow(refParticle.real_radius, 2.0) * Pi);

  double grate = 0.0;
  double average_rate = 0.0;
  double new_area = 0.0;

  runParticle = &refParticle;

  for (i = 0; i < numParticles; i++) {

    if (runParticle->mineral == 2) {

      spin_gr_nb++;
    }

    if (runParticle->newly_overgrown != 0.0) {

      new_area += runParticle->newly_overgrown;
      average_rate += sqrt(runParticle->newly_overgrown);
      divisor++;
    }

    runParticle = runParticle->nextP;
  }

  if (divisor != 0 && lowest_time != 0.0) {
    average_rate /= divisor;
    average_rate /= lowest_time;
  }

  // get the average dir_rate
  counter1 = 0;
  for (i = 0; i < numParticles; i++) {

    counter2 = 0;
    dir2 = 0.0;

    for (j = 0; j < 6; j++) {
      if (runParticle->dir_rate[j] != 0.0) {
        counter2++;
        dir2 += runParticle->dir_rate[j];
      }
    }

    if (counter2 != 0) {
      counter1++;
      dir1 += dir2 / double(counter2);
    }

    runParticle = runParticle->nextP;
  }

  if (counter1 != 0)
    dir1 /= double(counter1);
  else
    dir1 = 0.0;

  rad = pow(((3.0 * prtcl->mole_per_particle * prtcl->mV) / (4.0 * 3.1415)),
            (1.0 / 3.0));
  total_area = spin_gr_nb * Pi * rad * rad;
  new_area = Pi * rad * rad;

  if (prtcl) {
    if (prtcl->mineral == 1)
      rate << "spin->ol ;\t";
    else
      rate << "ol->spin ;\t";

    cout << "new radius is: " << prtcl->radius << endl;
  }

  divisor = 0;
  for (i = 0; i < 6; i++) {
    if (prtcl->normal_stress[i] != 0.0) {
      divisor++;
      average_normal_stress += prtcl->normal_stress[i];
    }
  }

  if (divisor != 0)
    average_normal_stress /= divisor;

  // still to do: areas / area-rates in percent of the system-size

  rate << timestep << " ;\t " << spin_gr_nb / numParticles << " ;\t "
       << the_size << " ;\t " << lowest_time << " ;\t " << spin_gr_nb << " ;\t "
       << total_area << " ;\t " << new_area << " ;\t " << new_area / lowest_time
       << " ;\t " << average_rate << " ;\t ";
  rate << new_area / the_size << " ;\t " << new_area / (the_size * lowest_time)
       << " ;\t " << (Pi * rad * rad) / the_size << " ;\t "
       << (Pi * rad * rad) / (the_size * lowest_time) << " ;\t "
       << average_normal_stress << " ;\t " << dir1 << endl;

  rate.close();
}

void Min_Trans_Lattice::CalcPressureScale(double bulkmod,
                                          double required_mean_stress) {

  double lambda, my;
  double poisson = 1.0 / 3.0;
  double dilatation, young = 0.0;
  double K, mean_stress;
  int i;

  for (i = 0; i < numParticles; i++) {
    young += runParticle->young;

    runParticle = runParticle->nextP;
  }

  young /= double(numParticles);

  lambda = -(poisson * young) / (2.0 * poisson * poisson + poisson - 1.0);
  my = young / (2.0 * poisson + 2.0);

  K = lambda + (2.0 / 3.0) * my;

  // reset the second scaling factor
  pascal_scale = 1.0;

  // define main scaling factor
  // we assume cubical dilatation (better: compression)
  dilatation = required_mean_stress / bulkmod;

  mean_stress = K * dilatation;

  pressure_scale = required_mean_stress / mean_stress;

  cout << "pressure-scale: " << pressure_scale << endl;
}

void Min_Trans_Lattice::Plot_Grains() {

  int i, j, k, m, grain_nb;
  ofstream grainfile;
  char filename[15] = "Grains";
  Particle *mark;
  bool cond;
  double xpos, ypos;
  bool allgrains[200];
  int counter, last_spring;

  if (!(secondary_timestep % plot_grain_nb)) {

    for (i = 0; i < 200; i++)
      allgrains[i] = false;

    for (i = 0; i < numParticles; i++) {
      if (runParticle->mineral == 2)
        allgrains[runParticle->grain] = true;

      runParticle->done = false;
      runParticle->marker = false;

      runParticle = runParticle->nextP;
    }

    sprintf(filename, "Grains.%d", secondary_timestep);

    grainfile.open(filename, ios::out | ios::trunc);

    for (i = 0; i < numParticles; i++) {

      if (!runParticle->is_lattice_boundary && runParticle->is_boundary &&
          !runParticle->done && allgrains[runParticle->grain] &&
          runParticle->mineral == 2) {

        counter = 0;
        for (j = 0; j < 6; j++) {
          if (runParticle->neigP[j]->grain == runParticle->grain &&
              runParticle->neigP[j]->mineral == runParticle->mineral)
            counter++;
        }

        grain_nb = runParticle->grain;

        if (counter > 0) {

          // preliminary settings for the first particle
          for (j = 0; j < 6; j++) {
            if (runParticle->neigP[j] == mark) {
              last_spring = j;
              break;
            }
          }

          mark = runParticle;
          mark->done = true;

          j = (j >= 5) ? 0 : j + 1;
          while (mark->neigP[j]->grain != grain_nb &&
                 mark->neigP[j]->mineral != mark->mineral)
            j = (j >= 5) ? 0 : j + 1;

          counter = last_spring = j;
          counter = (counter == 0) ? 5 : j - 1;
          while (mark->neigP[counter]->grain != grain_nb &&
                 mark->neigP[counter]->mineral != mark->mineral)
            counter = (counter == 0) ? 5 : counter - 1;

          j = (counter == 5) ? 0 : counter + 1;
          while (mark->neigP[j]->grain != grain_nb &&
                 mark->neigP[j]->mineral != mark->mineral) {
            xpos = (mark->xpos + mark->neigP[j]->xpos) / 2.0;
            ypos = (mark->ypos + mark->neigP[j]->ypos) / 2.0;

            grainfile << mark->grain << " " << xpos << " " << ypos << endl;

            j = (j >= 5) ? 0 : j + 1;
          }

          counter = last_spring;
          last_spring = mark->neig_spring[counter];
          mark = mark->neigP[counter];

          cond = true;
          while (cond) {

            j = (last_spring >= 5) ? 0 : last_spring + 1;
            while (mark->neigP[j]->grain != grain_nb &&
                   mark->neigP[j]->mineral != mark->mineral) {
              xpos = (mark->xpos + mark->neigP[j]->xpos) / 2.0;
              ypos = (mark->ypos + mark->neigP[j]->ypos) / 2.0;

              grainfile << mark->grain << " " << xpos << " " << ypos << endl;

              j = (j >= 5) ? 0 : j + 1;
            }

            counter = j;
            last_spring = mark->neig_spring[counter];

            mark = mark->neigP[j];

            if (mark->done)
              cond = false;
          }

        } else {
          for (j = 0; j < 6; j++) {
            xpos = (runParticle->xpos + runParticle->neigP[j]->xpos) / 2.0;
            ypos = (runParticle->ypos + runParticle->neigP[j]->ypos) / 2.0;

            grainfile << mark->grain << " " << xpos << " " << ypos << endl;
          }
        }

        mark = &refParticle;
        for (j = 0; j < numParticles; j++) {
          if (mark->grain == grain_nb)
            mark->done = true;

          mark = mark->nextP;
        }
      }

      mark = runParticle;
      runParticle = runParticle->nextP;
    }

    grainfile.close();
  }
}

bool Min_Trans_Lattice::Is_Same_Surface_Grain(Particle *prtcl1,
                                              Particle *prtcl2) {

  int i, j;
  Particle *common_neig;

  if (prtcl1->grain != prtcl2->grain)
    return (false);

  // the particles are on the same surface (given they are NEIGBORS!), if there
  // is one particle of the other phase which is connected to both of them
  for (i = 0; i < 6; i++) {
    if (prtcl1->neigP[i]) {
      if (prtcl1->neigP[i]->grain != prtcl1->grain) {
        for (j = 0; j < 6; j++) {
          if (prtcl1->neigP[i]->neigP[j]) {
            if (prtcl1->neigP[i]->neigP[j] == prtcl2) {
              return (true);
            }
          }
        }
      }
    }
  }

  return (false);
}

// called in Min_Trans_Lattice::Start_Reactions, after GBM() got called
// essentially this is the same routine as the Plot_Grain function, just not
// checking for Grains but for clusters
void Min_Trans_Lattice::Check_Percolation() {

  int i, j, k, spinel_particles = 0;
  vector<Particle *> clusterlist;
  list<int> clustersize;
  bool percolated = false;
  bool top_boundary, bottom_boundary, left_boundary, right_boundary;
  int counter, nb_particles;
  list<int>::iterator it;
  double tmp;

  if (!system_percolation) {

    for (i = 0; i < numParticles; i++) {

      runParticle->done = false;
      runParticle->cluster_nb = 0;

      runParticle = runParticle->nextP;
    }

    runParticle = &refParticle;

    spinel_particles = 0;
    clustersize.clear();

    for (i = 0; i < numParticles; i++) {

      if (runParticle->mineral == 2 && !runParticle->done) {

        clusterlist.clear();

        clusterlist.push_back(runParticle);

        runParticle->done = true;

        top_boundary = bottom_boundary = left_boundary = right_boundary = false;

        for (j = 0; j < clusterlist.size(); j++) {

          for (k = 0; k < 6; k++) {
            if (clusterlist[j]->neigP[k]) {
              if (clusterlist[j]->neigP[k]->mineral == 2 &&
                  !clusterlist[j]->neigP[k]->done) {
                clusterlist.push_back(clusterlist[j]->neigP[k]);
                clusterlist[j]->neigP[k]->done = true;
              }

              if (clusterlist[j]->neigP[k]->is_top_lattice_boundary)
                top_boundary = true;
              if (clusterlist[j]->neigP[k]->is_bottom_lattice_boundary)
                bottom_boundary = true;
              if (clusterlist[j]->neigP[k]->is_left_lattice_boundary)
                left_boundary = true;
              if (clusterlist[j]->neigP[k]->is_right_lattice_boundary)
                right_boundary = true;
            }
          }

          spinel_particles++;
        }

        // append the size of the cluster to the clustersize-vector
        clustersize.push_back(clusterlist.size());

        if (top_boundary && bottom_boundary && left_boundary && right_boundary)
          percolated = true;
      }

      runParticle = runParticle->nextP;
    }

    if (percolated) {
      system_percolation = true;
      cout << "Percolated" << endl;
    }

    if (percolated) {

      ofstream percolationfile;
      percolationfile.open("Spinel_Percolation", ios::out | ios::trunc);

      tmp = spinel_particles / numParticles;

      // total nb of spinel particles in the moment of percolation and the
      // spinel percentage of the total system
      percolationfile << spinel_particles << " \t" << numParticles << endl;

      // sort list
      while (clustersize.size() > 0) {

        counter = 0;
        nb_particles = *clustersize.begin();

        for (it = clustersize.begin(); it != clustersize.end(); ++it) {
          if (*it == nb_particles)
            counter++;
        }

        // number and size of clusters
        percolationfile << nb_particles << " \t" << counter << endl;

        clustersize.remove(nb_particles);
      }

      percolationfile.close();
    }
  }

  // now check the percolation over all minerals but olivine
  Check_Percolation_All();
}

// the other percolation-routine checks, if the spinel forms a throughgoing
// network. this one, instead, checks if all-but-olivine forms a throughgoing
// network
void Min_Trans_Lattice::Check_Percolation_All() {

  int i, j, k, spinel_particles = 0;
  vector<Particle *> clusterlist;
  list<int> clustersize;
  bool percolated = false;
  bool top_boundary, bottom_boundary, left_boundary, right_boundary;
  int counter, nb_particles;
  list<int>::iterator it;
  double tmp;

  if (!system_percolation) {

    for (i = 0; i < numParticles; i++) {

      runParticle->done = false;
      runParticle->cluster_nb = 0;

      runParticle = runParticle->nextP;
    }

    runParticle = &refParticle;

    spinel_particles = 0;
    clustersize.clear();

    for (i = 0; i < numParticles; i++) {

      if (runParticle->mineral != 1 && !runParticle->done) {

        clusterlist.clear();

        clusterlist.push_back(runParticle);

        runParticle->done = true;

        top_boundary = bottom_boundary = left_boundary = right_boundary = false;

        for (j = 0; j < clusterlist.size(); j++) {

          for (k = 0; k < 6; k++) {
            if (clusterlist[j]->neigP[k]) {
              if (clusterlist[j]->neigP[k]->mineral != 1 &&
                  !clusterlist[j]->neigP[k]->done) {
                clusterlist.push_back(clusterlist[j]->neigP[k]);
                clusterlist[j]->neigP[k]->done = true;
              }

              if (clusterlist[j]->neigP[k]->is_top_lattice_boundary)
                top_boundary = true;
              if (clusterlist[j]->neigP[k]->is_bottom_lattice_boundary)
                bottom_boundary = true;
              if (clusterlist[j]->neigP[k]->is_left_lattice_boundary)
                left_boundary = true;
              if (clusterlist[j]->neigP[k]->is_right_lattice_boundary)
                right_boundary = true;
            }
          }

          spinel_particles++;
        }

        // append the size of the cluster to the clustersize-vector
        clustersize.push_back(clusterlist.size());

        if (top_boundary && bottom_boundary && left_boundary && right_boundary)
          percolated = true;
      }

      runParticle = runParticle->nextP;
    }

    if (percolated) {
      system_percolation = true;
      cout << "Percolated" << endl;
    }

    if (percolated) {

      ofstream percolationfile;
      percolationfile.open("Total_Percolation", ios::out | ios::trunc);

      tmp = spinel_particles / numParticles;

      // total nb of spinel particles in the moment of percolation and the
      // spinel percentage of the total system
      percolationfile << spinel_particles << " \t" << numParticles << endl;

      // sort list
      while (clustersize.size() > 0) {

        counter = 0;
        nb_particles = *clustersize.begin();

        for (it = clustersize.begin(); it != clustersize.end(); ++it) {
          if (*it == nb_particles)
            counter++;
        }

        // number and size of clusters
        percolationfile << nb_particles << " \t" << counter << endl;

        clustersize.remove(nb_particles);
      }

      percolationfile.close();
    }
  }
}

void Min_Trans_Lattice::ResetPressureScale(double value, double t) {

  double mean_stress = 0.0;
  double required_mean_stress = clausius_clapeyron_threshold(t) + value;
  int i, counter = 0;

  for (i = 0; i < numParticles; i++) {
    if (runParticle->mineral == 1 && !runParticle->is_lattice_boundary) {
      mean_stress = mean_stress + (runParticle->sxx + runParticle->syy) / 2.0;
      counter++;
    }

    runParticle = runParticle->nextP;
  }

  mean_stress /= double(counter);

  // reset the second scaling factor
  pascal_scale = 1.0;
  pressure_scale = -(required_mean_stress) / mean_stress;

  cout << "pressure-scale: " << pressure_scale << endl;
}

void Min_Trans_Lattice::KeepPressure(double value, double t) {

  double mean_stress = 0.0;
  double required_mean_stress = clausius_clapeyron_threshold(t) + value;
  int i, counter = 0;

  for (i = 0; i < numParticles; i++) {
    if (runParticle->mineral == 1 && !runParticle->is_lattice_boundary) {
      mean_stress = mean_stress - (runParticle->sxx + runParticle->syy) / 2.0;
      counter++;
    }

    runParticle = runParticle->nextP;
  }

  mean_stress /= double(counter);

  if (mean_stress * pressure_scale < required_mean_stress)
    ResetPressureScale(value, t);
}

void Min_Trans_Lattice::SetOverPressure(bool yes, double op) {

  keep_pressure = yes;
  overpressure = op;
}

void Min_Trans_Lattice::SetMolePerParticle() {

  int i;

  for (i = 0; i < numParticles; i++) {

    runParticle->mole_per_particle = pow(runParticle->real_radius, 3.0) *
                                     1.33333333 * 3.14 / runParticle->mV;

    runParticle = runParticle->nextP;
  }
}
