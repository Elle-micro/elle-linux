
/******************************************************
 * Experiment Class for Elle/Latte
 *
 * works now mainly with the lattice spring code
 *
 * Daniel and Till 2005
 *
 * Latte Version 2.0
 * koehn_uni-mainz.de
 *
 * Latte Version 3.0
 * Erlangen 2021
 *
 * Clean Version 4.0
 * Erlangen 2025
 *
 ******************************************************/

// ------------------------------------
// system headers
// ------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <locale.h>
#include <vector>

// ------------------------------------
// elle headers
// ------------------------------------

#include "experiment.h"
#include "unodes.h" // include unode funct. plus undodesP.h
// for c++
#include "attrib.h" // enth�lt struct coord (hier: xy)
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

// have to define these in the new Elle version

using std::cout;
using std::endl;
using std::vector;

// CONSTRUCTOR
/*******************************************************
 * Dont really construct anything here at the moment
 *
 * The experiment class is directly called from the
 * Elle main function
 ********************************************************/

// ---------------------------------------------------------------
// Constructor of Experiment class
// ---------------------------------------------------------------

Experiment::Experiment()

// -----------------------------------------------------------------------------
//  first there is just some talking to the world plus
// -----------------------------------------------------------------------------

{
  cout << "Oh, my experiment starts !" << endl;
  cout << "Latte Version 4.0 2025" << endl;
  setlocale(LC_ALL, "en_US");
  experiment_time = 0;
}

double time_step;
double nonlinear_time;
double size;

/*************************************************************
 * Now we start with the initialization function
 * This function initialized some Elle basics
 * reads input file if one is there,
 * opens the interface
 * and starts the local initialization functions for the
 * desired process
 *
 *************************************************************/

void Experiment::Init() {
  //-------------------------------------------------
  // local variables
  //-------------------------------------------------

  char *infile; // input
  int err = 0;  // pass errors
  int i, j;     // counter
  float set;
  float visc;
  UserData udata; // Elle Structure for data from input
  int process;    // variable for the data (which process)
  int strength, value1, value2, value3, value4;

  //*-----------------------------------------------
  //* clear the data structures
  //*-----------------------------------------------

  ElleReinit();

  //*------------------------------------------------
  //* read the data
  //*-----------------------------------------------

  infile = ElleFile(); // input file specified by -i
  file = infile;
  ElleUserData(
      udata); // reads in data from the initial call of the program behind -u

  process = (int)udata[0]; // get this data from Elle, first number behind -u in
                           // Terminal during call

  // here we can define more udata values that are read directly into the script
  // from the interface this means for these values we do not have to recompile
  // which we normally have to do for a change in values. Here the numbers can
  // be changed during the call from the terminal directly remember that the
  // very first number after -u is the process below

  value1 = (int)udata[1]; // second number in Terminal
  value2 = (int)udata[2]; // third number in Terminal
  value3 = (int)udata[3]; // fourth number in Terminal
  value4 = (int)udata[4]; // fifth number in Terminal

  //----------------------------------------------------
  // processes are defined by simple integers
  //----------------------------------------------------
  wrapping = false; // not important right now

  if (strlen(infile) > 0) {
    //--------------------------------------------
    // only go in here if an input file is there
    //--------------------------------------------

    if (err = ElleReadData(infile))
      OnError(infile, err);

    //--------------------------------------------------
    // the switch function is used for manual input.
    // The processes can all be called from the
    // new interface.
    //
    // Without interface specify a process using the
    // -u command after calling the experiment.
    // ./my_experiment -u 1
    // -u 1 fracturing with aperture
    // -u 2 ridges with set geometry
    // -u 3 ridges with selforganized nucleation
    // -u 4 setup for grain fracturing
    //--------------------------------------------------

    switch (process) {

      /*******************************************************************************
       *
       *	PROCESSES Initialization
       *
       * so you can set up as many cases as you want for your
       * different ideas. That has the advantage that we dont have
       * to delete old ones when we want to make new stuff
       *
       * This will vary depending on the use of Latte, there should be
       * several versions for different processes, otherwise its too chaotic
       * this version has the fracturing and aperture (1) (Mathur et al., 2013)
       * the drifting for defined MORs and for the random ones (2) and (3)
       *(Hafermaas) and the new versions for fracturing grains in bands (Mathur
       *and Karimi)
       ********************************************************************************/

    case 0:
      // SetUpFromFile();
      break;

      /*******************************************************************************
       *	FRACTURING AND APERTURE (Mathur et al., 2023)
       *
       * This code simulates fracturing of a rock with or without grains and
       * measures the aperture of fractures during pure shear
       *
       * Deformation can include an initial compaction of the rock by uniaxial
       * compaction where normally nothing fractures. This simulates different
       * reservoir depth before the rock is deformed
       *
       ********************************************************************************/

    case 1: // fracturing

      cout << "Fracturing Aperture Mathur et al., 2023" << endl;
      cout << "Lattice Version 4.0, 2025" << endl;

      Activate_Lattice(); // construct the lattice

      // function from Giulia Fedrizzi to save Csv files for Latte
      getTimeForCsv(experiment_time, 1);

      // function in lattice.cc that is used to set a Variation on the breaking
      // threshold of springs.
      SetVariationBreakingThreshold(0.6);

      // function to change properties of all particles (Young,
      // Viscosity,breaking_strength)
      WeakenAll(0.4, 1.0, 1.5);

      // change the relaxation threshold (<1 = more precise)
      ChangeRelaxThreshold(0.5);

      // MakeGrainBoundaries(1.0, 0.8);		// Make grain boundaries weaker

      // Allowing only tensile failure, not shear
      // typically bonds can break in shear and tension following
      // Sachau+Koehn2013

      Only_Extension(true);

      // SetFracturePlot(50,1); 	// plot fractures after 50 bonds broken

      break;

      /*******************************************************************************
       *	DEVEPMENT OF RIDGES AND TRANSFORMS (Hafermaas and Koehn)
       *
       * This code simulates the growth of mid ocean ridges, transform faults
       *and microplates. This is a relatively complex code where particles can
       *grow and the model is coupled to a temperature field (finite difference
       *grid)
       *
       *
       * A healing function reconnects particles as a function of temperature
       * and a weakening function weakens elasticity, breaking strength and
       *viscosity as a function of temperature. Both thermo-mechanical couplings
       *are applied using an Arhenius like law (Hafermaas et al.).
       *
       * This version uses initially set T anomalies and a predefined offset
       * geometry for two plates that move apart from each other
       *
       ********************************************************************************/

    case 2: // ridges

      cout << "Mid Ocean Ridges and Transforms" << endl;
      cout << "Lattice Version 4.0, 2022-2025" << endl;

      Activate_Lattice(); // construct the lattice

      // function from Giulia Fedrizzi to save Csv files for Latte
      getTimeForCsv(experiment_time, 1);

      ChangePoisson(0.333333333);

      SetVariationBreakingThreshold(0.5); // Variation on breaking strings

      // This initializes an object fluid lattice in fluid_lattice.cc, normally
      // used for fluid pressure diffusion, but here for temperature diffusion

      Initialize_Fluid_Lattice(1000000, 1000, 0,
                               400.0); // (pressure, scale, conc, timefactor)

      // this allows the upper and lower row to also break,
      // normally they are not allowed to break

      ReleaseBoundaryParticlesY(1, 0); // 1st: break? 2nd: elastic wall

      // change properties of all particles (properties may be reset later,
      // careful)
      WeakenAll(
          0.5, 1 * 1e20,
          40); // WeakenAll(0.2,1,2); auf 0.1, (young,viscosity,break strength)

      // background temperature 200 degrees
      // critical age means new hot particles stay hot for x model stages
      InitHeat(200, 5); // 200 (temperature, critical age)

      // This is used for the self organized list, only active if called in
      // deformation function in the Run part
      InitDriftList();

      // activate the second layer that is used as analogue for the mantle
      // its a visco-elastic sheet pinning particles in the third dimension to
      // visco-elastic points

      ActivateSheet(1 * 1e20, 2 * 1e9); // das hier

      // intialize some hot spots, x1, x2, y1, y2, Temperature, critical age
      // x and y as boxes in heat lattice (here 100 x 100
      // careful, this depends on the size of the heat lattice

      Initialize_Heat_Box(29, 31, 75, 100, 1200, 1);
      Initialize_Heat_Box(44, 46, 50, 75, 1200, 1);
      Initialize_Heat_Box(29, 31, 25, 50, 1200, 1);
      Initialize_Heat_Box(44, 46, 0, 25, 1200, 1);

      // Thermomechanical weakening that adjusts elasticity, breaking strength
      // and viscosity to temperature using an Arhenius equation. Numbers are
      // prefactors for the equation first number for breaking strength, then
      // elasticity and viscosity

      CoolingDirectExp(0.3, 1.0, 1.0);

      // function that makes left and right side of box unbreakable, no action
      // happens here ymin ymax xmin xmax elasticity, viscosity, breaking
      // strength, porosity changes

      WeakenHorizontalBox(0, 1, 0.0, 0.25, 1.0, 1.0, 100, 1);
      WeakenHorizontalBox(0, 1, 0.75, 1.0, 1.0, 1.0, 100, 1);

      // Unodes on the left and right hand side are taken out in a list that is
      // used for the new particles that are inserted. This is why during the
      // simulation you will see them disappearing on the left and right hand
      // side. First number end of left side, second one beginning of right side
      // (start at 0, and at 1.0)

      Loose_Unodes(0.2, 0.8);

      // make relaxation a bit easier
      ChangeRelaxThreshold(5.0);

      break;

      /*******************************************************************************
       *	DEVEPMENT OF RIDGES AND TRANSFORMS (Hafermaas and Koehn)
       *
       * This code simulates the growth of mid ocean ridges, transform faults
       *and microplates. This is a relatively complex code where particles can
       *grow and the model is coupled to a temperature field (finite difference
       *grid)
       *
       *
       * A healing function reconnects particles as a function of temperature
       * and a weakening function weakens elasticity, breaking strength and
       *viscosity as a function of temperature. Both thermo-mechanical couplings
       *are applied using an Arhenius like law (Hafermaas et al.).
       *
       * This version uses the self organizing ridges with an intially vertical
       * seed of high T for better nucleation
       *
       ********************************************************************************/
    case 3: // ridges

      cout << "Mid Ocean Ridges and Transforms, self organizing part" << endl;
      cout << "Lattice Version4.0, 2022-2025" << endl;

      Activate_Lattice(); // construct the lattice

      // function from Giulia Fedrizzi to save Csv files for Latte
      getTimeForCsv(experiment_time, 1);

      SetVariationBreakingThreshold(0.5); // Variation on breaking strings (<1),

      // This initializes an object fluid lattice in fluid_lattice.cc, normally
      // used for fluid pressure diffusion, but here for temperature diffusion
      Initialize_Fluid_Lattice(1000000, 1000, 0,
                               400.0); // before (1000000, 1000, 0, 4000.0)
                                       // (pressure, scale, conc, timefactor)

      // this allows the upper and lower row to also break,
      // normally they are not allowed to break
      ReleaseBoundaryParticlesY(1, 0); // 1st: break? 2nd: elastic wall

      // change properties of all particles (properties may be reset later,
      // careful)
      WeakenAll(0.5, 1 * 1e20, 40); // WeakenAll(0.2,1,2); auf 0.1

      // background temperature 200 degrees
      // critical age means new hot particles stay hot for x model stages
      InitHeat(200, 5); // 200 (temperature, critical age)

      // This is used for the self organized list, only active if called in
      // deformation function in the Run part
      InitDriftList();

      // activate the second layer that is used as analogue for the mantle
      // its a visco-elastic sheet pinning particles in the third dimension to
      // visco-elastic points
      ActivateSheet(1 * 1e20, 2 * 1e9); // das hier

      // set a heat anomalie in the centre as a seed
      Initialize_Heat_Box(48, 52, 5, 95, 1200, 10);

      // Thermomechanical weakening that adjusts elasticity, breaking strength
      // and viscosity to temperature using an Arhenius equation. Numbers are
      // prefactors for the equation first number for breaking strength, then
      // elasticity and viscosity
      CoolingDirectExp(0.03, 0.50, 1.0);

      // function that makes left and right side of box unbreakable, no action
      // happens here ymin ymax xmin xmax elasticity, viscosity, breaking
      // strength, porosity changes
      WeakenHorizontalBox(0, 1, 0.0, 0.25, 1.0, 1.0, 100, 1);
      WeakenHorizontalBox(0, 1, 0.75, 1.0, 1.0, 1.0, 100, 1);

      // Unodes on the left and right hand side are taken out in a list that is
      // used for the new particles that are inserted. This is why during the
      // simulation you will see them disappearing on the left and right hand
      // side. First number end of left side, second one beginning of right side
      // (start at 0, and at 1.0)
      Loose_Unodes(0.2, 0.8);

      // make relaxation easier
      ChangeRelaxThreshold(5.0);

      break;

    case 4: // fracturing
      cout << "Fracturing Deformation Bands" << endl;
      cout << "Lattice Version 4.0 2025" << endl;

      //------------------------------------------------
      // default breaking strength 0.0017
      // default Young is 1.0
      // typically scaled to real values by multiplying
      // values times 10 GPa
      //------------------------------------------------

      Activate_Lattice(); // construct the lattice

      ChangePoisson(
          0.2); // if it is 0.3333 then bonds break only due to extension

      // change properties of all particles (properties may be reset later,
      // careful)
      WeakenAll(1, 1, 5); // WeakenAll(Young, viscosity (not used), breaking
                          // strength) 1 means nothing changes

      // function from Giulia Fedrizzi to save Csv files for Latte
      getTimeForCsv(experiment_time, 1);

      SetVariationBreakingThreshold(0.6); // Variation on breaking strings

      // change threshold for the relaxation (smaller = more precise but needs
      // more time)
      ChangeRelaxThreshold(0.6);

      MakeGrainBoundaries(1.0, 0.8); // Make grain boundaries weaker

      //=========================================================================
      // Weaken Grain routine uses the Elle grains (Flynns) with integer numbers
      // and applies a given elastic, viscous and breaking strength to all
      // particles in the grain (multiplied by values).
      //
      // WeakenGrain (int nb, float constant, float visc, float break_strength)
      // highest_grain is set in the lattice Activator, its the highest grain
      //==========================================================================

      WeakenGrainAttribute(0, 0.00001, 0,
                           1); // first number is for Flynn Attribute, Youngs
                               // Modulus, viscosity, breaking strength
      // WeakenGrainAttribute (3, 0.00001, 0, 1);
      // WeakenGrain (3,0.00001, 0, 1);

      // Allowing only tensile failure, not shear
      // typically bonds can break in shear and tension following
      // Sachau+Koehn2013

      // Only_Extension(true);

      ChangeRelaxMax(10000); // maximum relaxation loops

      // SetFracturePlot(50,1); 	// plot fractures after 50 bonds broken

      break;
    }
  } else
    cout << "no file open ! " << endl; // no input file in function call, can be
                                       // opened from the interface
}

/******************************************************************
 * A runfunction for the Experiment
 *
 *	In this function each process (again specified by -u in a switch
 *   command) is executed.
 *  	whereas the init function is called only ones, the run function
 * 	is called every time step. This is where the real action takes place
 *
 *  	Latte version 2.0, 2005/6
 *  	Latte version 3.0, 2021
 * 	Latte version 4.0, 2025
 ******************************************************************/

void Experiment::Run() {
  //----------------------------------------------------
  // some local variables
  //----------------------------------------------------

  int i, j, k, kk, ii, iii; // counter
  int time;                 // time
  int process;              // int for process
  UserData udata;           // elle defined user data (-u )

  ElleCheckFiles(); // Check the files

  ElleUserData(udata);     // get the usr data from elle
  process = (int)udata[0]; // process is user data 0

  ElleUpdateDisplay();

  //--------------------------------------------------
  // get the time from the interface
  //--------------------------------------------------

  time = EllemaxStages(); // number of stages

  j = 0;
  kk = 0;
  //--------------------------------------------------
  // loop through the time steps
  //--------------------------------------------------

  for (i = 0; i < time; i++) // cycle through stages
  {
    cout << "time step:" << experiment_time << endl;

    j++; // another time counter

    switch (process) {
    case 0:
      // RunFromFile(experiment_time);
      break;

    case 1: // fracturing

      // function from Giulia Fedrizzi to save Csv files for Latte
      getTimeForCsv(experiment_time, 2);

      //==========================================================
      // For the first part of the simulation defined by
      // experiment time the deformation is uniaxial
      //==========================================================
      if (experiment_time < 84) {
        DeformLattice(0.0001, 0);
      }
      //==========================================================
      // For the rest of the simulation the deformation
      // is for example pure shear, so extending in x
      //==========================================================
      else {
        DeformLatticeExtend(0.0001, 0);
      }
      //==========================================================
      // function that determines the aperture of the cracks
      //==========================================================
      Get_Aperture();

      // and talk to the interface
      UpdateElle();

      break;

    case 2: // ridges

      // function from Giulia Fedrizzi to save Csv files for Latte
      getTimeForCsv(experiment_time, 2);

      //============================================================================================
      // This function deals with the deformation. Here initial geometries are
      // assumed with first the opening relative to the initial box size of 1.0.
      // Then 4 vertical ridge positions are given and the plates will move
      // apart at these positions forcing the ridge opening. Ridges are counting
      // downwards in y from 1.0. numbers are low y1, xpos1, lowy2, xpos2 lowy3,
      // xpos3, xpos4 (last will always be to 0.0). the final number offers a
      // vertical gradient , but that does not work well yet, so if this number
      // is 0.0 there is no gradient.
      //============================================================================================
      DeformLatticedriftmidnew4(-0.001, experiment_time, 0.75, 0.3, 0.5, 0.45,
                                0.25, 0.3, 0.45, 0.0);

      //=============================================================================================
      // this function is looking for empty spots in the repulsion box. Empty
      // spots can be filled with new particles asuming the ridge is intruded
      // with melt (or stretches) the last number is the temperature of new
      // particles. The routine deals with the full insertion of new particles
      // and their connection to the lattice. routine also checks that the empty
      // box is between fractured particles and not just between to stretched
      // particles that are still connected with a long spring.
      //=============================================================================================
      Find_Empty_Box(experiment_time, 1200);

      //=============================================================================================
      // time is relative, you could change the deformation function to have
      // faster spreading however this leads to instabilities in the code and
      // open spaces therefore the time has to be changed in other functions in
      // order to model faster or slower spreading. This is a) the diffusion
      // constant for T, b) the healing probability that is time dependent and
      // c) the viscosity in RelaxSheet and in ViscousStep. ViscousStep is
      // really the viscosity of springs in the elastic lattice and Relax Sheet
      // deals with the second layer that is thought to represent the mantle.
      // space in meters for one box so resolution 200 = 20km for whole box
      // so with 0.001 deformation it opens 20m per step
      // for time = 10 and dif factor = 10 spreading rate is 20cm/year so fast
      // slow dif factor has to be 100
      //==============================================================================================

      //============================================================================
      // solve the heat diffusion
      // this also copies heat back and forth from particles so that the heat
      // moves with the moving particles.
      // the actual action happens in fluid_lattice.cc, which is an object that
      // is reached through the function in lattice.cc
      //============================================================================

      solve_Heat(10000, 10, 0.000004, 1200, 5,
                 10); // space, time (years was 20), dif const; temp, critAge,
                      // diffusion factor

      //==================================================================================
      // Thermomechanical coupling with three prefactors for breaking strength
      // elastic modulus and viscosity. They are multiplied by an exponential of
      // 1/temperature. Prefactors for overall behaviour, exponential for local
      // function of T. Important this function changes everything as a function
      // of a given average background number.
      //==================================================================================
      CoolingDirectExp(0.02, 1.0,
                       1.0); // that is the thermal weakening/strenghtening

      //==================================================================================
      // two vertical boxes on the right and left hand side are created that the
      // sides unbreakable, they are not important for the model. This has to be
      // done because the weakening function weakens everything relativ to an
      // initially determined average. Therefore it makes the sides weak as well
      // so that they might fracture. However they are just supposed to be
      // passive plates.
      //==================================================================================
      WeakenHorizontalBox(0, 1, 0.0, 0.3, 1.0, 1.0, 100, 1);
      WeakenHorizontalBox(0, 1, 0.8, 1.0, 1.0, 1.0, 100, 1);
      //=====================================================================
      // set a variation again, because it was overwritten by the weakening
      //=====================================================================
      SetVariationBreakingThreshold(0.5); // was 0.5

      //=====================================================================
      // relaxation routine that calculates the movements of particles as a
      // function of forces on them and moves them until it finds an equilibrium
      // Once this is achieved bonds can fracture. The bonds that is most
      // probable to break breaks and the model relaxes again. This goes on
      // until no more bonds break, therefore it can take very long. .
      //=====================================================================
      Relaxation(); // relax stresses

      //============================================================================
      // Relaxation of the Maxwell sheet below. This is just an attachement of
      // particles to visco-elastic positions in the third dimension. Its
      // supposed to represent the mantle. Number should be relaxation time.
      //============================================================================
      RelaxSheet(1e2);

      //============================================================================
      // viscous step for the springs. They can also act visco-elastic and their
      // length can increase to allow a stress relaxation.
      //============================================================================
      cout << "viscous" << endl;
      ViscousStep(100000); // time up = less viscous?

      //==============================================================================
      // The second thermomechanical coupling, in this case healing as a
      // function of T. The function test a critical distance for healing (if
      // more than factor times particle size apart dont heal), a time step for
      // a basic probability and potentially a change in elastic constant and
      // breaking strength of new springs the function is basic probability
      // times exponential of -1/T, so that the healing increases with an
      // increase in T.
      //==============================================================================
      cout << "healing" << endl;
      HealingTMol(1.8, 0.01, experiment_time, 1.0, 1.0);

      // giving out some values for particles 10040 to 10160... hopefully across
      // a ridge
      Dump_Young(10040, 10160);

      // and talk to the interface
      UpdateElle();

      break;

    case 3: // ridges random

      // function from Giulia Fedrizzi to save Csv files for Latte
      getTimeForCsv(experiment_time, 2);

      //============================================================================================
      // This function deals with the deformation.
      // here the model is first extended homogenously in the centre and
      // fractures. Works best with a heat seed in the centre. The model then
      // localizes in y-rows once a particle is inserted. It memorizes that
      // particle and opens the plates locally there.
      //============================================================================================
      DeformLatticedrift6(-0.001, experiment_time, 1.0, 0.3, 0.1, 0.7);

      //=============================================================================================
      // this function is looking for empty spots in the repulsion box. Empty
      // spots can be filled with new particles asuming the ridge is intruded
      // with melt (or stretches) the last number is the temperature of new
      // particles. The routine deals with the full insertion of new particles
      // and their connection to the lattice. routine also checks that the empty
      // box is between fractured particles and not just between to stretched
      // particles that are still connected with a long spring.
      //=============================================================================================
      Find_Empty_Box(experiment_time, 1200);

      //=============================================================================================
      // time is relative, you could change the deformation function to have
      // faster spreading however this leads to instabilities in the code and
      // open spaces therefore the time has to be changed in other functions in
      // order to model faster or slower spreading. This is a) the diffusion
      // constant for T, b) the healing probability that is time dependent and
      // c) the viscosity in RelaxSheet and in ViscousStep. ViscousStep is
      // really the viscosity of springs in the elastic lattice and Relax Sheet
      // deals with the second layer that is thought to represent the mantle.
      // space in meters for one box so resolution 200 = 20km for whole box
      // so with 0.001 deformation it opens 20m per step
      // for time = 10 and dif factor = 10 spreading rate is 20cm/year so fast
      // slow dif factor has to be 100
      //==============================================================================================

      //============================================================================
      // solve the heat diffusion
      // this also copies heat back and forth from particles so that the heat
      // moves with the moving particles.
      // the actual action happens in fluid_lattice.cc, which is an object that
      // is reached through the function in lattice.cc
      //============================================================================

      solve_Heat(100, 10, 0.000004, 1200, 5,
                 10); // space, time (years was 20), dif const; 0.000001, temp,
                      // critAge, diffusion factor

      //==================================================================================
      // Thermomechanical coupling with three prefactors for breaking strength
      // elastic modulus and viscosity. They are multiplied by an exponential of
      // 1/temperature. Prefactors for overall behaviour, exponential for local
      // function of T. Important this function changes everything as a function
      // of a given average background number.
      //==================================================================================
      CoolingDirectExp(0.03, 0.50,
                       1.0); // that is the thermal weakening/strenghtening

      //==================================================================================
      // two vertical boxes on the right and left hand side are created that the
      // sides unbreakable, they are not important for the model. This has to be
      // done because the weakening function weakens everything relativ to an
      // initially determined average. Therefore it makes the sides weak as well
      // so that they might fracture. However they are just supposed to be
      // passive plates.
      //==================================================================================
      WeakenHorizontalBox(0, 1, 0.0, 0.3, 1.0, 1.0, 100, 1);
      WeakenHorizontalBox(0, 1, 0.8, 1.0, 1.0, 1.0, 100, 1);
      //=====================================================================
      // set a variation again, because it was overwritten by the weakening
      //=====================================================================
      SetVariationBreakingThreshold(0.5); // was 0.5

      //=====================================================================
      // relaxation routine that calculates the movements of particles as a
      // function of forces on them and moves them until it finds an equilibrium
      // Once this is achieved bonds can fracture. The bonds that is most
      // probable to break breaks and the model relaxes again. This goes on
      // until no more bonds break, therefore it can take very long. .
      //=====================================================================
      Relaxation(); // relax stresses

      //============================================================================
      // Relaxation of the Maxwell sheet below. This is just an attachement of
      // particles to visco-elastic positions in the third dimension. Its
      // supposed to represent the mantle. Number should be relaxation time.
      //============================================================================
      RelaxSheet(1e2);

      //============================================================================
      // viscous step for the springs. They can also act visco-elastic and their
      // length can increase to allow a stress relaxation.
      //============================================================================
      cout << "viscous" << endl;
      ViscousStep(1000000); // time up = less viscous?

      //==============================================================================
      // The second thermomechanical coupling, in this case healing as a
      // function of T. The function test a critical distance for healing (if
      // more than factor times particle size apart dont heal), a time step for
      // a basic probability and potentially a change in elastic constant and
      // breaking strength of new springs the function is basic probability
      // times exponential of -1/T, so that the healing increases with an
      // increase in T.
      //==============================================================================
      cout << "healing" << endl;
      HealingTMol(1.8, 0.01, experiment_time, 1.0, 1.0);

      // giving out some values for particles 10040 to 10160... hopefully across
      // a ridge
      Dump_Young(10040, 10160);

      // and talk to the interface
      UpdateElle();

      break;

    case 4: // fracturing Deformation bands

      // function from Giulia Fedrizzi to save Csv files for Latte
      getTimeForCsv(experiment_time, 2);

      //=============================================
      // here deformation can be added, in this
      // case uniaxial compaction
      //=============================================
      // DeformLattice(0.0001,1);

      // DeformLatticeSimpleShear (0.0001,1);
      DeformLatticePureShear(0.0001, 1);

      // and talk to the interface.
      UpdateElle();

      break;
    }
    experiment_time++;
    Set_TimeFrac(experiment_time);
  }
}
