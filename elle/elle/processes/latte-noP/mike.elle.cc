/*****************************************************************
 * main function for mike
 * calls routines in the lattice class
 * and makes a lattice object
 *
 * book version for Stylolites, Daniel 2005
 *
 * Mike version 2.0 Feb. 2004 Daniel, Jochen and Till
 ****************************************************************/

//-------------------------------------------------------------
// system headers
//-------------------------------------------------------------

#include <cstdio>

#include <cmath>
#include <iostream.h>

//-------------------------------------------------------------
// Elle headers
//-------------------------------------------------------------

#include "attrib.h"
#include "error.h"
#include "file.h"
#include "init.h"
#include "interface.h"
#include "lattice.h" // lattice class of Mike includes particle.h
#include "nodes.h"
#include "phase_lattice.h"
#include "runopts.h"
#include "triattrib.h"
#include "unodes.h" // unodes routines includes the unodeP.h
#include "update.h"

Phase_Lattice *p_mike;
Lattice *mike;

//-----------------------------------------------------------------
// define some local functions
//-----------------------------------------------------------------

int InitSetMike(); // this function is external in the main function
int SetMike();

//----------------------------------------------------------------
// START OF PROGRAM
//----------------------------------------------------------------

/*****************************************************************
 * this function will be run when the application starts,
 * when an elle file is opened or
 * if the user chooses the "Rerun" option
 * InitSetMike is external and declared in the mike.main.cc function
 * where it is connected with Elle
 *
 *****************************************************************/

int InitSetMike() {
  //-------------------------------------------------
  // local variables
  //-------------------------------------------------

  char *infile;   // input
  int err = 0;    // pass errors
  int i;          // counter
  UserData udata; // Elle Structure for data from input
  int process;    // variable for the data (which process)

  //*-----------------------------------------------
  //* clear the data structures
  //*-----------------------------------------------

  ElleReinit();

  //*------------------------------------------------
  //* read the data
  //*-----------------------------------------------

  infile = ElleFile();
  ElleUserData(
      udata); // reads in data from the initial call of the program behind -u
  process = (int)udata[0]; // get this data from Elle

  //----------------------------------------------------
  // processes are defined by simple integers
  //----------------------------------------------------

  if (strlen(infile) > 0) {
    if (err = ElleReadData(infile))
      OnError(infile, err);
  }

  //----------------------------------------------------
  // make the lattice in beginning and put initial
  // boundary conditions on the system.
  // needs an infile at call of mike (-i infile.elle)
  //----------------------------------------------------

  if (process == 0) {
    // this is for the custom-settings from .mike-prefs-file

  } else if (process == 1) // fracturing
  {
    cout << "Fracturing" << endl;
    cout << "Lattice Version 2.0, 2004/5" << endl;

    mike = new Lattice;

    mike->SetPhase(0.0, 0.0, 0.5, 1.2);

    mike->SetGaussianSpringDistribution(1.5, 0.5);

    mike->MakeGrainBoundaries(1.0, 0.3);

    mike->SetFracturePlot(1, 0);
  } else if (process == 2) // fracture boudinage
  {
    cout << "Fracture Boudinage" << endl;
    cout << "Lattice Version 2.0, 2004/5" << endl;

    mike = new Lattice;

    mike->SetPhase(0.0, 0.0, 0.5, 1.2);

    mike->WeakenAll(0.1, 1.0, 1.0);

    mike->WeakenHorizontalParticleLayer(0.20, 0.40, 20.0, 1.0, 1.0);

    mike->WeakenHorizontalParticleLayer(0.85, 0.92, 20.0, 1.0, 1.0);

    mike->SetFracturePlot(1, 0);
  } else if (process == 3) // expanding inclusions
  {
    cout << "Expanding Inclusions" << endl;
    cout << "Lattice Version 2.0, 2004/5" << endl;

    mike = new Lattice;

    mike->SetPhase(0.0, 0.0, 0.2, 1.2);

    mike->SetFracturePlot(1, 0);
  } else if (process == 4) // shrinkage patterns
  {
    cout << "Shrinkage Patterns" << endl;
    cout << "Lattice Version 2.0, 2004/5" << endl;

    mike = new Lattice;

    mike->SetPhase(0.0, 0.0, 0.2, 1.2);

    mike->SetFracturePlot(1, 0);
  } else if (process == 5) // viscoelastic
  {
    mike = new Lattice;

    mike->SetWallBoundaries(1, 1.0);

    mike->SetPhase(0.0, 0.0, 0.5, 1.2);

    mike->WeakenAll(0.1, 1000000000,
                    1.0); // konstante geändert, 0.1 statt 1 als viscosity

    // mike->SetAnisotropyRandom(20,20);

    mike->WeakenHorizontalParticleLayer(
        0.46, 0.53, 100.0, 100.0,
        1.0); //"normal": ...20,20,1 ;;; geringe viskos.: ...100, 100, 1

    mike->SetFracturePlot(1, 0);
  } else if (process == 6) // grooves on free surfaces
  {
    cout << "Dissolution Grooves" << endl;
    cout << "Phase_Lattice Version 2.0, 2004/5" << endl;

    p_mike = new Phase_Lattice;

    p_mike->SetPhase(0.0, 0.0, 500.0, 0.8);

    p_mike->SetGaussianRateDistribution(2.0, 0.5);

    p_mike->WeakenAll(4.0, 1.0, 1.0);

    p_mike->Set_Mineral_Parameters(1);

    p_mike->Set_Absolute_Box_Size(0.0002);

    p_mike->Set_Time(6000.0, 4);

    p_mike->DissolveXRow(0.95, 1.1);

  } else if (process == 7) // Stylolites
  {
    cout << "Stylolite Roughening" << endl;
    cout << "Phase_Lattice Version 2.0, 2004/5" << endl;

    p_mike = new Phase_Lattice;

    //----------------------------------------------------
    // give bonds very high strength (*500) to avoid
    // in this case fracturing during Stylolite growth
    //----------------------------------------------------

    p_mike->SetPhase(0.0, 0.0, 500.0, 0.8);

    //----------------------------------------------------
    // Set a Gaussian distribution on the rate constants
    // of single particles, first mean value, second
    // deviation
    //----------------------------------------------------

    p_mike->SetGaussianRateDistribution(2.0, 0.5);

    //----------------------------------------------------
    // change the Youngs modulus of particles to make
    // them stiffer (*4.0). Second and third number are
    // breaking strength and viscosity, 1.0 means no
    // change of these parameters
    //----------------------------------------------------

    p_mike->WeakenAll(4.0, 1.0, 1.0);

    //----------------------------------------------------
    // Set some mineral paramters, 1 means quartz as
    // mineral, sets the molecular volume and the
    // surface free energy
    //----------------------------------------------------

    p_mike->Set_Mineral_Parameters(1);

    //----------------------------------------------------
    // gives the x dimension of the Elle box in meters
    //----------------------------------------------------

    p_mike->Set_Absolute_Box_Size(0.1);

    //----------------------------------------------------
    // set the time for one deformation step. 6000 years
    // 4 means years
    //----------------------------------------------------

    p_mike->Set_Time(6000.0, 4);

    //----------------------------------------------------
    // dissolve initially one horizontal row of particles
    // in the middle of the Elle box
    // numbers are min and max y value
    //----------------------------------------------------

    p_mike->DissolveYRow(0.49, 0.5);

    p_mike->ChangeRelaxThreshold(100); // factor, 1=normal

  } // end Stylolite

  //-----------------------------------------------------
  // Set the run function
  //-----------------------------------------------------

  ElleSetRunFunction(SetMike);
}

/******************************************************************
 * A runfunction for Mike
 ******************************************************************/

int SetMike() {
  //----------------------------------------------------
  // some local variables
  //----------------------------------------------------

  int i, j; // counter
  int time; // time
  int process;
  UserData udata;

  ElleCheckFiles();

  ElleUserData(udata);
  process = (int)udata[0];

  //--------------------------------------------------
  // get the time from the interface
  //--------------------------------------------------

  time = EllemaxStages(); // number of stages

  //--------------------------------------------------
  // loop through the time steps
  //--------------------------------------------------

  for (i = 0; i < time; i++) // cycle through stages
  {
    cout << "time step" << i << endl;

    if (process == 1) // fracturing
    {
      if (i < 2) {
        mike->DeformLattice(0.001, 1);
      } else {
        mike->DeformLatticePureShear(0.001, 1);
      }

    } else if (process == 2) // fracture boudinage
    {
      mike->DeformLatticePureShear(0.001, 1);
    } else if (process == 3) // expanding inclusions
    {
      mike->ShrinkGrain(9, 1.01, 1);

    } else if (process == 4) // shrinkage patterns
    {
      mike->ShrinkBox(0.001, 1);
    } else if (process == 5) {
      mike->TimeStep = i;

      mike->adjust = false; // for debugging only

      mike->DeformLatticePureShear(0.001, 1);

      mike->ViscousRelax(1, 1000000000.0); // kann auch 100 sein...1000000000.0

    } else if (process == 6) // grooves
    {
      p_mike->DeformLattice(0.002, 1);

      p_mike->Dissolution_Strain(20);

    } else if (process == 7) // Stylolites
    {

      //-------------------------------------------------
      // Deform the lattice from upper and lower
      // boundaries. Upper and lower part of the lattice
      // are now just pressed together assuming
      // there is no resistance. Stresses build up once
      // the two sides meet. Deformation steps have to
      // be very small
      // The side walls are fixed. Movement is vertical
      // in steps of 0.00005 * y size (1.0)
      // second number means make a picture after the
      // movement. 0 means take no picture.
      //-------------------------------------------------

      p_mike->DeformLatticeNewAverage2side(0.00005, 1);

      //---------------------------------------------------
      // Dissolution routine. Particles are just dissolved
      // depending on stress, elastic and surface
      // energies. One particle is dissolved in a time
      // step. 100 means take a picture after 100 particles
      // have dissolved.
      //---------------------------------------------------

      p_mike->Dissolution_Stylos(100, 0, 0);
    }
  }
}
