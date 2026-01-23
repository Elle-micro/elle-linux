/******************************************************
 * Spring Code Lattice 3.0
 *
 * Functions for lattice class in lattice.cc
 *
 * Basic Spring Code for fracturing
 *
 *
 *
 *
 * Daniel Koehn and Jochen Arnold feb. 2002 to feb. 2003
 * Oslo, Mainz
 *
 * Daniel Koehn dec. 2003
 *
 * We thank Anders Malthe-Sörenssen for his enormous help
 * and introduction to these codes
 *
 * Daniel Koehn and Till Sachau 2004/2005
 * koehn_uni-mainz.de
 *
 * Daniel Koehn 2018
 * daniel.koehn@glasgow.ac.uk
 ******************************************************/

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
#include <vector>

// include newmat libs (eigen values and eigen vectors):
// newmat (matrix-operations)
// #include <newmatap.h>
// #include <iomanip>
// #include "newmatio.h"

// ------------------------------------
// elle headers
// ------------------------------------

#include "lattice.h" // include lattice class
#include "unodes.h"  // include unode funct. plus undodesP.h
// for c++
#include "attrib.h" // enth�lt struct coord (hier: xy)
#include "attribarray.h"
#include "attribute.h"
#include "check.h"
#include "convert.h"
#include "display.h"
#include "error.h"
#include "file.h"
#include "fluid_lattice.h"
#include "general.h"
#include "interface.h"
#include "nodes.h"
#include "nodesP.h"
#include "polygon.h"
#include "pressure_lattice.h"
#include "runopts.h"
#include "tripoly.h"
#include "update.h"

// have to define these in the new Elle version

using namespace std;

// newmat (matrix-operations)
// using namespace NEWMAT;
using std::cout;
using std::endl;
using std::vector;
// CONSTRUCTOR
/*******************************************************
 * Constructor for lattice class
 * calls the MakeLattice function and defines some variables
 * at the moment. The Constructor then builds a lattice once
 * the lattice class is called (i.e. by the user in the "elle"
 * function).
 *
 * Daniel spring 2002
 *
 * add security stop for pictdumps
 *
 * Daniel march 2003
 *
 * make lattice is now taken out of the constructor
 * otherwise everything is made in the experiment class
 * has now to be called using Activate Lattice
 ********************************************************/

// ---------------------------------------------------------------
// Constructor of Lattice class
// ---------------------------------------------------------------

Lattice::Lattice()
    :

      // -----------------------------------------------------------------------------
      // default values for some variables. Variables are declared in
      // lattice header
      // -----------------------------------------------------------------------------

      dimx(1.0),      // basic x dimension from Elle
      dimy(1.0),      // basic y dimension from Elle
      particlex(100), // default size for particles in x

      default_Young(1.0),
      // default_Poisson (1.0 / 3.0),
      default_Poisson(0.2),

      nbBreak(0),        // number of broken springs in lattice
      LowSeal_bnd_nb(0), // number of broken springs in lattice under hor_seal
      Seal_bnd_nb(0),    // number of broken springs in seal
      UppSeal_bnd_nb(0), // number of broken springs in lattice above hor_seal

      def_time_p(0), // deformation time zero (pure shear)
      def_time_u(0), // deformation time zero (uniaxial
      // compression)
      grain_counter(0), // initialize grain counter
      visc_rel(0), // flag for viscous relax routines in relaxation (zero is no)
      Pi(3.1415927),         // number pi
      pressure_scale(10000), // default 1.0 is 10 GPa
      pascal_scale(1000000), // to get pascal from megapascal
      set_max_pict(false),   // no security stop for pict dump
      num_pict(0),           // initialize pictcounter
      walls(false),
      // relaxation threshold
      relaxthres(0.0000001), debug_nb(1400000000),
      sheet(false), // no initial viscous sheet

      grav(false), // gravitation for rifting
      grav_Press(false),

      transition(false) {
  // nothing happens
  // except an adaption of the pascal rate for the phase change...
  // pascal_scale *= 100.0;

  alpha = Pi / 3;

  // default value for the wrapping
  wrapping = false;

  P = new Pressure_Lattice(&refParticle);

  cal_time = P->time_a;
  press_cal_time = 0.0;
}

/************************************************************+
 *
 *	Called from the experiment class to construct the
 *	Lattice if it is needed. Connects Elle and Latte
 *
 **************************************************************/

void Lattice::Activate_Lattice() {

  int i, j;

  // --------------------------------------------------------------------
  // call the MakeLattice function and give it a number defining the
  // Lattice. 1 is triangular at the moment
  // --------------------------------------------------------------------

  MakeLattice(1); // Make the Lattice

  Copy_NeigP_List(); // for the viscous relaxation virtual
  // neighbour list

  highest_grain = HighestGrain(); // get highest grain from Elle

  olivine_grain_nb = highest_grain;

  spinel_grain_nb = 0;

  shear = true;

  SetSpringRatio(); // set the spring ratio for viscous part

  visc_flag = false; // debugging

  for (i = 0; i < numParticles;
       i++) // default variables for the viscous routine
  {
    for (j = 0; j < 9; j++)
      runParticle->rep_rad[j] = runParticle->radius;

    runParticle->rot_angle = 0.0;

    runParticle->smallest_rad = runParticle->radius;

    runParticle->merk = 0;

    runParticle->maxis1length = runParticle->radius;
    runParticle->maxis2length = runParticle->radius;

    // this just for security
    runParticle->maxis1angle = 0.0;
    runParticle->maxis2angle = 0.0;

    runParticle = runParticle->nextP;
  }
}

// Copy_NeigP_List
/***************************************************************************************
 * function copies neigP list to NeigP2 list, which isn� subjected to breaking
 *springs
 *
 * function is used for viscous retardation. This virtual neighbour list is
 *important once springs break. The particle should still be able to change its
 *shape. This is a problem since the normal repulsion routine assumes that the
 *particles are round
 *
 * written by Till Sachau 2004
 ****************************************************************************************/
// ------------------------------------------------------------------
// function is called from the constructor of the lattice class
// ------------------------------------------------------------------

void Lattice::Copy_NeigP_List() {
  int i, j;

  for (i = 0; i < numParticles; i++) {
    for (j = 0; j < 8; j++) {
      runParticle->neigP2[j] = runParticle->neigP[j];
    }
    runParticle = runParticle->nextP;
  }
}

// MAKE_LATTICE
/****************************************************************
 * Function that builds the lattice. is given the lattice type.
 * does build its own lattice using particles. Creates particles,
 * defines xy positions, connects all particles by calling connect in
 * the particle class and connects to Unodes (pointer to them)
 * note that some variables are defined in the header lattice.h
 *
 * Note this is a bit simple at the moment. It uses xy coordinates and
 * first finds these even though they are actually allready there in the
 * Unodes. If this function is called and the Unode number and position
 * and the particle number and position do not fit the code goes bananas.
 * Therefore we have to write more advanced functions in the future.
 *
 * The same goes for the connections. It uses also xy position at the
 * moment but would be more clever to draw a circle or sphere around
 * a particle to find its neighbours. This will be especially important
 * once we deal with random lattices.
 *
 * daniel spring and summer 2002
 ******************************************************************/

// --------------------------------------------------------------------
// function MakeLattice of Lattice Class
// input variable int type
//
// function called from the constructor Lattice::Lattice()
// --------------------------------------------------------------------

void Lattice::MakeLattice(int type) {
  // -----------------------------------------------------------------------------
  // local variables
  // ----------------------------------------------------------------------------

  float xx;       // variable for x position in calculation
  float yy;       // variable for y position in calculation
  int irow, icol; // variable for rows in x and y in
  // calculation
  int i, j, jj; // counter
  Coords xy;    // unode position
  double oldy;  // help to check change in y row

  right_lattice_wall_pos = 1.0;
  left_lattice_wall_pos = 0.0;

  // ---------------------------------------------------------------------
  // zero some variables for checks
  // ---------------------------------------------------------------------

  for (i = 0; i < 1000000; i++) {
    nodes_Check[i] = false;
    next_inBox[i] = -1; // no tail node in Box
  }

  // -------------------------------------------------------------------
  // the repulsion box is 2500000 squares wide now. This means that
  // the maximum resolution is a bit more than 1300 particles in
  // the x direction.
  // if larger resolution is wanted this number has to be increased !
  // however if there is not enough memory this will lead
  // to a segmentation fault once the program is started !
  // A segmentation fault even though the program compiles.
  // in addition to the initial definition of the size of the vector
  // repBox[] in the lattice.h header
  // -------------------------------------------------------------------

  for (i = 0; i < 2500000; i++) {
    repBox[i] = 0;
    node_Box[i] = -1;
  }
  local_nbBreak = 0;
  internal_break = 0;

  // talk to the public !

  cout << "H" << endl;
  cout << "E" << endl;
  cout << "L" << endl;
  cout << "L" << endl;
  cout << "O" << endl;
  cout << " " << endl;
  cout << " I am Latte " << endl;
  cout << " Version 3.0 2018 " << endl;
  cout << " " << endl;

  // ------------------------------------------------------------------------------
  // this function can be used in the future in another version of the
  // MakeLattice
  // function to get number of particles.
  // ------------------------------------------------------------------------------

  numParticles = ElleMaxUnodes(); // call elle function to define
  // number of particles

  cout << numParticles << endl; // talk a bit
  cout << "Constructing Lattice\n" << endl;

  // -----------------------------------------------------------------------------------
  // the triangular grid is defined by having rows of particles in x
  // that are shifted
  // by the radius of the particles. The y length is then found as
  // follows:
  // -----------------------------------------------------------------------------------
  // -----------------------------------------------------------------------------------
  // now everything is defined by ppm2elle and read in at the beginning.
  // the following loop finds out the particlex dimension of the
  // lattice. It just
  // runs along x row and checks when the y changes which gives you the
  // length of the
  // row and size of the lattice. Y is then just particlenumber divided
  // by x.
  // ------------------------------------------------------------------------------------

  if (type == 1) // if it is a triangular grid apply ratio
  {
    // --------------------------------------------------------------------
    // ElleGetParticleUnode returns the object Unode with the index
    // that
    // you give the function. Then we can call directly functions
    // within
    // the Unode objects themselves. XY coordinates are not public in
    // the
    // Unode class !
    // ---------------------------------------------------------------------

    ElleGetParticleUnode(0)->getPosition(&xy); // start with
    // first Unode

    oldy = xy.y; // old position

    irow = int(sqrt((double)numParticles)); // x always
    // smaller than y

    for (i = 0; i < irow; i++) // run along first row
    {
      ElleGetParticleUnode(i)->getPosition(&xy); // next particle
      // position

      if (xy.y != oldy) // if it changes we are beginning a new
                        // row
      {
        particlex = i; // i is the x dimension now
        break;
      }
    }
    particley = numParticles / particlex; // and define y
  }

  // cout << " px " << particlex << endl;	// talk a bit again
  // cout << " py " << particley << endl;	// particlex and particley
  //  are dimensions in x,y

  numParticles = particlex * particley; // mike definition of
  // particle number (has to
  // match Unodes)

  cout << " Num Particles " << numParticles << endl;

  /*************************************************************
   * initiate particle list
   * put them in a connected list and give them numbers
   * each particle is connected with the next one in the list
   * and connected backwards so that list is linked two ways
   * list is also circular so that a search can be started
   * anywhere in list by any particle in both directions
   * First and Last list particles are referenced in pointers
   *
   * Note that at the moment a refParticle is used in the lattice class to get
   * the beginning of the list. Pointers did not work so far. And passing
   *pointers from the outside makes things messy. Maybe there is a better trick
   *for this. Note that now once particles may be deleted we have to check not
   *to delete the reference particle.
   *************************************************************/

  runParticle = new Particle; // pointer to first particle, make
  // particle (call particle constructor)
  refParticle = *runParticle; // refParticle is the particle that I
  // point to

  firstParticle = preRunParticle = &refParticle; // set additional
  // pointers

  // -----------------------------------------------------------------------------------
  // This loop does four things: it creates new particles, connects them
  // to the list
  // forwards and backwards and gives particles numbers
  // In the particle objects a pointer points to the next particle
  // (nextP), a pointer
  // points to the previous particle (prevP) and a number is given for
  // particles from
  // 0 to numParticles - 1 (nb)
  // two pointers are used in the loop to define particles (runParticle)
  // and the previous
  // particles (preRunParticle)
  // -----------------------------------------------------------------------------------

  for (i = 1; i < numParticles;) // loop numParticles -1 times
  {
    runParticle = new Particle;          // make new object and point
    preRunParticle->nextP = runParticle; // set pointer in prev
    // object
    runParticle->prevP = preRunParticle; // set pointer to prev
    // object
    preRunParticle->nb = i - 1;   // my number is...start with 0
    preRunParticle = runParticle; // change pointers
    i++;
  }
  lastParticle = runParticle;          // define end of list
  lastParticle->nb = numParticles - 1; // give last particle a
  // number
  lastParticle->nextP = firstParticle; // and connect a circular
  // list forwards
  firstParticle->prevP = lastParticle; // and backwards

  /****************************************************************
   * set x and y positions of particles
   * Run through the list. then define x and y values
   * Note that everything is done asuming the particle radius is 0.5
   * and then rescaled to the elle coordinates (or boxsize)
   *
   ****************************************************************/

  // This routine makes a lattice itself. but the lattice is now
  // actually
  // already made in ppm2elle so that mike just reads the position of
  // Unodes
  // from an elle file. Still left this here, maybe need it some day

  /*
   * runParticle = firstParticle; // go to beginning of list for (i=0;i
   * < numParticles;) // loop through all particles {
   * //---------------------------------------------------------------------
   * // first find the y coordinate by dividing number by xsize and taking
   * the // integer value of the result
   * //---------------------------------------------------------------------
   * irow = (i)/particlex;
   * //---------------------------------------------------------------------
   * // then find x coordinate by subtracting the xrow(100) from the number
   * (i) // . 100 is irow * particlex. This is 0 for the first row and then
   * // minus 100 etc if particlex is 100.
   * //---------------------------------------------------------------------
   * icol = i - (irow * particlex); xx = icol;
   * //----------------------------------------------------------------------
   * // since the lattice is triangular each second row has to be shifted by
   * the // radius of the particles (radius 0.5 in this definition) // ok at
   * the moment since box has still absolute size of particlex
   * //----------------------------------------------------------------------
   * if (irow != (2*(irow/2))) { xx = xx + 0.5; // checks uneven numbers } yy
   * = irow*1.0/xyratio; // scale y-rows for triangular runParticle->xpos =
   * xx/particlex; // scale to elle box runParticle->ypos = yy/particlex;
   * runParticle = runParticle->nextP; // grab next particle i++; }
   */

  // ---------------------------------------------------------------
  // initialization of some Elle attributes that appear in the
  // interface and are written out in the function UpdateElle
  // find them in basecode file.h
  // ----------------------------------------------------------------

  ElleInitUnodeAttribute(U_FRACTURES);
  ElleInitUnodeAttribute(U_PHASE);
  ElleInitUnodeAttribute(U_TEMPERATURE);
  ElleInitUnodeAttribute(U_DIF_STRESS);
  ElleInitUnodeAttribute(U_MEAN_STRESS);
  ElleInitUnodeAttribute(U_DENSITY);
  ElleInitUnodeAttribute(U_YOUNGSMODULUS);
  ElleInitUnodeAttribute(U_ENERGY);
  ElleInitUnodeAttribute(U_DISLOCDEN);
  ElleInitUnodeAttribute(U_VISCOSITY);
  ElleInitUnodeAttribute(U_ATTRIB_A);
  ElleInitUnodeAttribute(U_ATTRIB_B);
  ElleInitUnodeAttribute(U_ATTRIB_C);
  ElleInitUnodeAttribute(U_S_EXPONENT);
  ElleInitUnodeAttribute(U_FINITE_STRAIN);
  ElleInitUnodeAttribute(U_STRAIN);

  /****************************************************************
   * just set the position of the particle equal to its unode
   * position. the positions are made triangular in the right
   * way by ppm2elle now
   ****************************************************************/

  for (i = 0; i < numParticles; i++) {
    ElleGetParticleUnode(runParticle->nb)->getPosition(&xy);
    runParticle->xpos = xy.x;
    runParticle->ypos = xy.y;
    runParticle->oldx = xy.x;
    runParticle->oldy = xy.y;
    runParticle->initialy = xy.y;

    runParticle = runParticle->nextP;
  }

  /*******************************************************************
   * connect the particles
   * This just runs through the list again and calls a function
   * within the particles themselves. (defined in particle.cc)
   * Also set the springs in Particle class by call function SetSpring
   * defined in particle.cc
   ******************************************************************/

  runParticle = firstParticle;       // start
  for (i = 0; i < numParticles; i++) // run
  {
    runParticle->Connect(type, particlex, particley, wrapping); // call
    // Connect
    // and
    // connect
    // (connect
    // in
    // particle.cc)
    runParticle->radius = runParticle->smallest_rad =
        runParticle->radius * 100.0 / particlex; // rescale
    // radius

    runParticle->area =
        runParticle->radius * runParticle->radius * 3.1415927; // set
    // area
    // for
    // reactions

    runParticle = runParticle->nextP; // go to next particle
  }

  // ----------------------------------------------------------------
  // call set Springs in particle.cc after all particles are
  // connected. SetSprings give default spring constants,
  // youngs modulus and viscosity
  // ----------------------------------------------------------------

  for (i = 0; i < numParticles; i++) // run
  {
    runParticle->SetSprings(default_Young, default_Poisson); // call SetSprings

    runParticle = runParticle->nextP; // go to next particle
  }

  // -----------------------------------------------------------------
  // new function that finds the backward pointer of neighbours.
  // springs are always defined in two ways, by the particle
  // itself and by its neighbour. However the particle does not
  // know which spring of its neighbour points backwards.
  // this function finds the backward pointer and writes an
  // integer into neig_spring[] array that is the write counter in the
  // neigP[] array of the neighbour.
  // -----------------------------------------------------------------

  for (i = 0; i < numParticles; i++) // loop
  {
    // a and b for the new radius-function (Alen(int) in particle-class
    // !!!!!!!!!!!!!!!!
    runParticle->a = runParticle->b = runParticle->radius;
    runParticle->fixangle = 0.0;

    for (j = 0; j < 8; j++) // run through neighbours
    {
      if (runParticle->neigP[j]) // if neighbour
      {
        for (jj = 0; jj < 8; jj++) // run through neighbours
                                   // springs
        {
          if (runParticle->neigP[j]->neigP[jj]) // if
                                                // neighbour
          {
            if (runParticle->nb ==
                runParticle->neigP[j]->neigP[jj]->nb) // look
                                                      // backwards
                                                      // and
                                                      // find
                                                      // particle
                                                      // again
            {
              runParticle->neig_spring[j] = jj; // thats
              // the one
            }
          }
        }
      }
    }

    runParticle = runParticle->nextP;
  }

  // ------------------------------------------------------------------------------
  // now finish the MakeLattice call by doing some extra things
  // first call FindUnodeElle defined in Lattice.cc (below)
  // This function runs through the list of particles and calls an Elle
  // function
  // that was written for Mike and returns a pointer to the right Unode.
  // With this we can call functions in the Unode class directly and
  // in addition use all plotting routines. Also do not have to call
  // vector each
  // time we want a specific Unode, each particle has its Unode defined
  // now.
  // I think that makes things faster ?
  // -------------------------------------------------------------------------------

  FindUnodeElle(); // new thing that finds pointer to Unode

  cout << "make RepulsionBox" << endl;

  MakeRepulsionBox(); // fill box for Repulsion

  // ------------------------------------------------------------------------------
  // Now we define whatever boundaries we want for the particles at the
  // edges of
  // the box. Box is not wrapping at the moment !!
  // function SetBoundaries is defined below in lattice.cc
  // ------------------------------------------------------------------------------

  SetBoundaries(); // apply some boundary conditions

  cout << "initiate Grains" << endl;

  InitiateGrainParticle(); // find the grains, grain boundaries and
  // connect nodes

  // -----------------------------------------------------------------------------------------
  // Now we do a basic relaxation at the beginning. Should only run
  // through this
  // once if lattice is in equilibrium (which it should be at the
  // moment)
  // function Relaxation defined below in lattice.cc
  // thats the core of all this ! All the fun is happening there. This
  // we should
  // optimize as best as we can.
  // Relaxation calls FullRelax() in lattice.cc which returns a flag
  // depending on
  // whether or not things are relaxed. If things are not relaxed this
  // routine will
  // loop until they are relaxed. Therefore in this loop the program
  // runs forever if
  // too much happened or something goes wrong.
  // FullRelax calls Relax in the particle (defined in particle.cc)
  // which relaxes each
  // individual particle by looking at its neighbours plus whatever we
  // want (fluidpressure ?)
  // ------------------------------------------------------------------------------------------

  // initialize the neigpos-list of the particles.
  runParticle = &refParticle;
  for (int i = 0; i < numParticles; i++) {
    for (int j = 0; j < 8; j++)
      if (runParticle->neigP[j]) {
        runParticle->neigpos[j][0] =
            runParticle->neigP[j]->xpos - runParticle->xpos;
        runParticle->neigpos[j][1] =
            runParticle->neigP[j]->ypos - runParticle->ypos;
      }
    runParticle = runParticle->nextP;
  }

  Relaxation(); // do a basic relaxation

  // -------------------------------------------------------------------------------
  // and finally we call a routine that updates the interface and plots
  // the
  // Unode stresses (lots still to be done there, we have the full
  // tensor and can
  // do with that whatever we want !
  // UpdateElle also defined below in lattice.cc
  // returns stresses from particles to Unodes for the plotting
  // and calls an Elle routine to update interface
  // --------------------------------------------------------------------------------

  // UpdateElle ();		// talk to Elle
}

// INITIATE_GRAIN_PARTICLE
/********************************************************************************
 * function is called from make mike function in the lattice class.
 *
 * this is a heavy one. Used to take ages of time. Now quite fast.
 * finds the grains from the flynn and gives particles number of grain or flynn
 * that it is in.
 * To do that fast it uses the relaxbox and gets the box position from
 * the Unodes position and then finds the particle. So it only loops
 * through the flynns and their nodes.
 * the second part defines the grain boundary particles by checking which
 * particles have neighbours with different grains. Rather straight forward.
 * The last part finds the nodes in Elle and connects them with the particles.
 * Loops through the particles and defines a region. Then it calls a
 * EllePointInRegion routine with a node coord. In order not to loop through
 * all nodes each time (loops within loops are deadly) we first define a
 * kind of second repulsion box with the same dimensions as the real
 * repulsion box and fill all the nodes in the box. There is a second array
 * for single nodes that contains the tailparticles. Did not want to change
 * the node structure in Elle.
 * Once the box is filled we can find nodes close to a particle (in same or
 * neighbouring boxes) and can easily connect the nodes with particles.
 * just have to make an additional boundary check in case (and that is so)
 * that nodes are not surrounded by particles.
 *
 * daniel summer 2002
 * daniel changed december 2002
 *
 * changed december 2004 -> all Elle nodes found now
 *****************************************************************************/

// -------------------------------------------------------------------------
// function InitiateGrainParticle in lattice class
//
// called from Lattice::MakeLattice()
// -------------------------------------------------------------------------

void Lattice::InitiateGrainParticle() {
  int i, j, ii, jj, iii, jjj; // lots and lots of counter for the
  // thousand loops
  int nb;       // the important number
  int neig1;    // help neighbour
  int maxFlynn; // maxFlynns (max in program not number of
  // flynns !)
  int count; // a counter yeah !
  int max;   // for max of nodes (also not real node
  // number !)
  int x, y;      // conversion for box
  double fx, fy; // more conversion for box

  Coords bbox[3], boundbox[4], xy; // boxes and single coord structures
  // (attrib.h)
  NodeAttrib *p; // Node structure (nodesP.h)

  // ------------------------------------------------------------------------------------
  // the first routine loops through the flynns, gets the Unode list for
  // flynns and
  // then defines a grain for each particle (an int for which flynn the
  // particle is
  // part of.
  // in order to find the particle the position of the Unode is taken
  // and converted to
  // the relaxation box where we look for the right particle
  // a bit around ten corners but very fast and seems to be efficient
  // ------------------------------------------------------------------------------------

  cout << "define grains" << endl;

  vector<int> unodelist; // define unodelist vector for unodes

  maxFlynn = ElleMaxFlynns(); // max number of Flynns (not real max
  // number)

  runParticle = &refParticle; // for particle loops

  // -----------------------------------------------------------------
  // now loop through the flynns and get the Unode lists
  // -----------------------------------------------------------------

  for (j = 0; j < maxFlynn; j++) // loop through flynns
  {
    if (ElleFlynnIsActive(j)) // use only active flynns
    {
      grain_counter = grain_counter + 1;

      ElleGetFlynnUnodeList(j, unodelist); // just gets the
      // Unodelist for
      // flynn j

      count = unodelist.size(); // gets size of Unodelist

      for (i = 0; i < count; i++) // go through Unodelist
      {
        nb = unodelist[i]; // the real number of the Unode

        // -------------------------------------------------------------------------------
        // this routine returns the unode Object if you give it a
        // number so that we then
        // can directly call functions within the unode class like
        // getPosition (the
        // Position of type coord is not public in the Unode
        // class)
        // --------------------------------------------------------------------------------

        ElleGetParticleUnode(nb)->getPosition(&xy);

        // -----------------------------------------------------------------------------------
        // find position in Relax box to find Particle for Unode.
        // Particles point to Unodes
        // but Unodes dont point back. For the box we scale the
        // coord, convert them to ints
        // and make them one-dimensional
        // -----------------------------------------------------------------------------------

        fx = (xy.x * particlex);
        fy = (xy.y * particlex);
        x = int(fx);
        y = int(fy);

        // ---------------------------------------------------------
        // and make one-dimensional
        // ---------------------------------------------------------

        x = (y * particlex * 2) + x;

        // ------------------------------------------------------------
        // now we go to the box and look in it if particle is the
        // right one. if not we loop through the tailparticles
        // until
        // we find the right one
        // -------------------------------------------------------------

        if (nb == repBox[x]->nb) // if the right particle
                                 // is directly in box
        {
          if (repBox[x]->grain == -1) // this strange only works
                                      // like this
          {
            repBox[x]->grain = j; // give particle grain
            // number
          }
        } else // loop through tail
        {
          runParticle = repBox[x]; // redefine pointer for
          // looping

          while (runParticle->next_inBox) // if there is a
                                          // tail (should
                                          // be)
          {
            if (nb == runParticle->next_inBox->nb) // if it
                                                   // is
                                                   // right
                                                   // particle
            {
              if (runParticle->next_inBox->grain == -1) // this
                                                        // strange
                                                        // only
                                                        // works
                                                        // like
                                                        // this
              {
                runParticle->next_inBox->grain = j; // give
                // particle
                // grain
                // number
              }
              break; // jump out of loop, all done
            } else {
              runParticle = runParticle->next_inBox; // go
              // to
              // next
              // in
              // tail
            }
          }
        }
      }
    }
  }

  // -----------------------------------------------------------------------------------------
  // this defines the grain boundary particles by checking whether or
  // not particle has
  // a neighbour that is part of a different grain or flynn.
  // -----------------------------------------------------------------------------------------

  cout << " define grain boundaries " << endl;

  for (i = 0; i < numParticles; i++) // and loop through particles
  {
    for (j = 0; j < 8; j++) // loop through neighbours of particle
    {
      if (runParticle->neigP[j]) // if there is a neighbour
      {
        if (runParticle->neigP[j]->grain != runParticle->grain) // if
                                                                // not
                                                                // of
                                                                // same
                                                                // grain
        {
          runParticle->is_boundary = true; // I am a boundary
          // particle
          break;
        }
      }
    }
    runParticle = runParticle->nextP; // go on looping
  }

  // --------------------------------------------------------------------------------
  // This is the big one. Connects nodes with particles.
  // First put nodes in a Node list similar to the relaxation list so
  // that we can
  // find them easily if we know in what box we want to look
  // --------------------------------------------------------------------------------

  neig1 = -1; // help definition

  // ----------------------------------------------
  // first zero all nodes (= -1)
  // ----------------------------------------------

  max = ElleMaxNodes(); // gets max nodes

  // -----------------------------------------------------------------------------
  // this routine fills the box for nodes. It has the same dimension as
  // the relax box. Therefore it should be easy to find neighbouring
  // nodes
  // and particles using this box. Once it is defined it is rather fast.
  // still have to be careful with negative values.
  // nodes have additional attribute next_inBox (int) that contains tail
  // nodes
  // similar to the tail particles in the relaxbox.
  // the next_inBox array is used as virtual pointer for the nodes.
  // -----------------------------------------------------------------------------

  for (i = 0; i < max; i++) // loop through nodes
  {

    if (ElleNodeIsActive(i)) // if node active
    {
      p = ElleNode(i); // gets structure for the node

      // --------------------------------------------------
      // get position in the box first x and y and make
      // them ints. scale same as particle relax box
      // the box scales with the number of particles
      // --------------------------------------------------

      fx = (p->x * particlex);
      fy = (p->y * particlex);
      x = int(fx);
      y = int(fy);

      // ---------------------------------------------------------
      // and make one-dimensional
      // ---------------------------------------------------------

      x = (y * particlex * 2) + x;

      // ----------------------------------------------------------
      // now fill the node box if it is empty (-1)
      // if it is full loop in the next_inBox array until
      // you find the end of the tail of nodes in that box
      // position
      // -----------------------------------------------------------

      if (node_Box[x] == -1) {
        node_Box[x] = i; // fill the box
      } else             // now we make the tails
      {
        ii = node_Box[x]; // shift to find node with tail

        while (next_inBox[ii] != -1) // loop until end of tail
        {
          ii = next_inBox[ii]; // and shift and shift
        }
        next_inBox[ii] = i; // found tail and fill in node i
      }
    }
  }

  // --------------------------------------------------------------------------------
  // now we come to the connection between particles and nodes
  // loop through particles until you find a boundary particle
  // then check all the neighbours of that particle
  // check the neighbour of the neighbour to find a triangle with the
  // original particle. Then find box position of the particles (plus
  // some
  // extra left and right box otherwise some nodes fall through the
  // net). Then
  // find the nodes with that boxposition in the node box. Then check if
  // that node
  // is really within the triangle that we look for.
  // In the end we take the particles that are boundary of the lattice
  // and
  // define a square around the particle in which we search for nodes on
  // the
  // boundaries.
  // --------------------------------------------------------------------------------

  cout << " and connect nodes to particles " << endl;

  for (i = 0; i < numParticles; i++) // loop through the particles
  {
    neig1 = -1; // flag -1 to find first neighbour

    {
      for (j = 0; j < 8; j++) // loop through its neighbours
      {
        if (runParticle->neigP[j]) // if there is a neighbour
                                   // at that point
        {
          {
            if (neig1 == -1) // if the first neighbour
              // is not yet defined
              neig1 = j; // i am the first

            else // else have to look for the second
                 // neighbour in triangle
            {

              for (ii = 0; ii < 8; ii++) // loop through
                                         // second
                                         // neighbour
                                         // neighbours
              {
                if (runParticle->neigP[neig1]->neigP[ii]) // if
                                                          // neighbour
                                                          // at
                                                          // that
                                                          // point
                {
                  // -------------------------------------------------------
                  // if the number of the first
                  // neighbour is the same as
                  // the number of the neighbour of the
                  // second neighbour
                  // we have defined a closed triangle
                  // --------------------------------------------------------

                  if (runParticle->neigP[neig1]->neigP[ii]->nb ==
                      runParticle->neigP[j]->nb) {
                    // --------------------------------------------------------
                    // the bbox is used for the
                    // definition of a region for the
                    // function EllePointInRegion that
                    // gets a bbox, the dim of
                    // the bbox (here three) and the
                    // structure of type coord
                    // with the coordinates of the
                    // point (in this case the node)
                    // ----------------------------------------------------------

                    // ----------------------------------------------------------
                    // DEFINE bbox[] by using x and y
                    // of the particles of the
                    // found triangle
                    // -----------------------------------------------------------

                    bbox[0].x = runParticle->xpos;
                    bbox[0].y = runParticle->ypos;
                    bbox[1].x = runParticle->neigP[neig1]->xpos;
                    bbox[1].y = runParticle->neigP[neig1]->ypos;
                    bbox[2].x = runParticle->neigP[j]->xpos;
                    bbox[2].y = runParticle->neigP[j]->ypos;

                    // -------------------------------------------------------------
                    // now we loop through 9 different
                    // box positions and all their
                    // tail nodes and check if any of
                    // these nodes fall into the
                    // triangle
                    // the box positions are the
                    // positions from the three
                    // particles
                    // and the left and right box of
                    // each particle. (thats where
                    // we might loose nodes because
                    // box is square lattice)
                    //
                    // x is always the final box
                    // position
                    // --------------------------------------------------------------

                    for (iii = 0; iii < 9; iii++) {
                      if (iii == 0) // box
                                    // position
                                    // for
                                    // first
                                    // particle
                      {
                        // ---------------------------------------------------------
                        // get positions, scale by
                        // lattice and convert to
                        // ints
                        // ---------------------------------------------------------

                        fx = (runParticle->xpos * particlex);
                        fy = (runParticle->ypos * particlex);
                        x = int(fx);
                        y = int(fy);

                        // ---------------------------------------------------------
                        // and make
                        // one-dimensional
                        // ---------------------------------------------------------

                        x = (y * particlex * 2) + x;
                      } else if (iii == 1)
                        x = x + 1; // box to
                      // right
                      else if (iii == 2)
                        x = x - 2; // box to
                      // left
                      else if (iii == 3) {
                        // ---------------------------------------------------------
                        // get positions, scale by
                        // lattice and convert to
                        // ints
                        // ---------------------------------------------------------

                        fx = (runParticle->neigP[neig1]->xpos * particlex);
                        fy = (runParticle->neigP[neig1]->ypos * particlex);
                        x = int(fx);
                        y = int(fy);

                        // ---------------------------------------------------------
                        // and make
                        // one-dimensional
                        // ---------------------------------------------------------

                        x = (y * particlex * 2) + x;
                      } else if (iii == 4)
                        x = x + 1; // box to
                      // right
                      else if (iii == 5)
                        x = x - 2; // box to
                      // left
                      else if (iii == 6) {
                        // ---------------------------------------------------------
                        // get positions, scale by
                        // lattice and convert to
                        // ints
                        // ---------------------------------------------------------

                        fx = (runParticle->neigP[j]->xpos * particlex);
                        fy = (runParticle->neigP[j]->ypos * particlex);
                        x = int(fx);
                        y = int(fy);

                        // ---------------------------------------------------------
                        // and make
                        // one-dimensional
                        // ---------------------------------------------------------

                        x = (y * particlex * 2) + x;
                      } else if (iii == 7)
                        x = x + 1; // box to
                      // right
                      else if (iii == 8)
                        x = x - 2; // box to
                      // left

                      jj = node_Box[x]; // define
                      // node

                      while (jj != -1) // loop
                                       // through
                                       // tails
                      {
                        // ----------------------------------------------------------
                        // if node is active and
                        // has not been checked go
                        // for it
                        // ----------------------------------------------------------

                        if (ElleNodeIsActive(jj) && !nodes_Check[jj]) {
                          p = ElleNode(jj); // get
                          // struct
                          // of
                          // node
                          xy.x = p->x; // convert
                          // x
                          // to
                          // coord
                          xy.y = p->y; // convert
                          // y
                          // to
                          // coord

                          // ----------------------------------------------------
                          // now make the region
                          // check with the
                          // predefined box,
                          // the dimension of
                          // the box and the xy
                          // coord struct
                          // of the node that
                          // you found in the
                          // node box
                          // -----------------------------------------------------

                          if (EllePtInRegion(bbox, 3, &xy)) {
                            // --------------------------------------------------
                            // there are 16
                            // possible places
                            // for nodes in a
                            // particle. check
                            // if one is free
                            // and put node in
                            // do this for all
                            // three bounding
                            // particles
                            // --------------------------------------------------

                            for (jjj = 0; jjj < 16; jjj++) {
                              // main
                              // particle

                              if (runParticle->elle_Node[jjj] == -1) {
                                runParticle->elle_Node[jjj] = jj;
                                break;
                              }
                            }
                            for (jjj = 0; jjj < 16; jjj++) {
                              // first
                              // neighbour

                              if (runParticle->neigP[neig1]->elle_Node[jjj] ==
                                  -1) {
                                runParticle->neigP[neig1]->elle_Node[jjj] = jj;
                                break;
                              }
                            }
                            for (jjj = 0; jjj < 16; jjj++) {
                              // second
                              // neighbour

                              if (runParticle->neigP[j]->elle_Node[jjj] == -1) {
                                runParticle->neigP[j]->elle_Node[jjj] = jj;
                                break;
                              }
                            }
                            // -----------------------------------------------
                            // now define flag
                            // that tells
                            // particles that
                            // this node has
                            // been connected
                            // otherwise we
                            // connect it at
                            // least three
                            // times
                            // to the same
                            // three particles
                            // ------------------------------------------------

                            nodes_Check[jj] = true;
                          }
                        }
                        // ---------------------------------
                        // this is for the while
                        // loop
                        // go to next tail node in
                        // box
                        // ---------------------------------

                        jj = next_inBox[jj];
                      }
                    }
                  }
                }
              }
              // ---------------------------------------------
              // reset neig1 to next neighbour
              // ---------------------------------------------

              neig1 = j;
            }
          }
        }
      }
    }
    // ------------------------------------------------------------------------------------
    // now check the particles at the boundary of the lattice. Some
    // nodes might not be
    // within a triangle of particles. Therefore draw a square around
    // each boundary
    // particle and check if there is an unconnected node in this box.
    // -------------------------------------------------------------------------------------

    if (runParticle->is_lattice_boundary) // if lattice boundary
    {
      // -------------------------------------
      // define a new bbox with 4 dimensions
      // particle is in the centre
      // -------------------------------------

      boundbox[0].x = runParticle->xpos - (runParticle->radius * 4.0); // lower
      // left
      // corner
      boundbox[0].y = runParticle->ypos + (runParticle->radius * 4.0);

      boundbox[1].x = runParticle->xpos + (runParticle->radius * 4.0); // upper
      // left
      // corner
      boundbox[1].y = runParticle->ypos + (runParticle->radius * 4.0);

      boundbox[2].x = runParticle->xpos + (runParticle->radius * 4.0); // upper
      // right
      // corner
      boundbox[2].y = runParticle->ypos - (runParticle->radius * 4.0);

      boundbox[3].x = runParticle->xpos - (runParticle->radius * 4.0); // lower
      // right
      // corner
      boundbox[3].y = runParticle->ypos - (runParticle->radius * 4.0);

      for (jj = 0; jj < max; jj++) // loop through all nodes
                                   // again
      {
        // ----------------------------------------------------------
        // if node is active and has not been connected
        // ----------------------------------------------------------

        if (ElleNodeIsActive(jj) && !nodes_Check[jj]) {
          p = ElleNode(jj); // get node struct and convert
          // coordinates
          xy.x = p->x;
          xy.y = p->y;

          if (EllePtInRegion(boundbox, 4, &xy)) // if in
                                                // the box
          {
            for (jjj = 0; jjj < 16; jjj++) // loop through
                                           // possible nodes
                                           // in particle
            {
              if (runParticle->elle_Node[jjj] == -1) // if
                                                     // free
              {
                runParticle->elle_Node[jjj] = jj; // fill
                nodes_Check[jj] = true;           // flag up
                break;
              } else if (runParticle->elle_Node[jjj] == jj) // in
                                                            // case
                                                            // node
                                                            // is
                                                            // in
                                                            // there
              {
                break; // we can exit
              }
            }
          }
        }
      }
    }

    runParticle = runParticle->nextP; // loop in particles
  }

  cout << " done with grain definitions " << endl;

  // ----------------------------------------------------------------------
  // just a test to see if all Elle nodes are connected now
  // ----------------------------------------------------------------------

  i = 0;

  for (jj = 0; jj < max; jj++) // loop through Elle nodes
  {
    if (ElleNodeIsActive(jj) && !nodes_Check[jj]) // if node is
                                                  // active and has
                                                  // not been
                                                  // checked
    {
      p = ElleNode(jj); // define xy coordinates and pass out to
      // interface

      cout << "node not connected to lattice" << jj << " " << p->x << " "
           << p->y << endl;

      i++;
    }
  }
  if (i != 0)
    cout << i << endl; // how many nodes were not connected
}

/****************************************************************************************
 * new function that is supposed to read in Elle data after things are deformed
 *so that talking back and forth between Elle and mike works better.
 *
 * Function now works. Only problems are fractures. We are now saving the youngs
 *moduli for single grains and also the values of the grainboundaries when they
 *are set in the initialization functions. These are then taken off the
 *particles and springs, the new geometry is read from Elle and the values are
 *again applied to particles. Therefore grain boundaries can now really move
 *including all their properties for the Lattice spring code.
 *
 * You may have to check if all initialization function parameters are saved...
 *
 * Function uses the new Unode attribute Flynn which has directly the grain
 *number
 *
 * Daniel 2005/6
 *****************************************************************************************/

void Lattice::GetNewElleStructure() {
  int i, j, ii, jj, iii, jjj; // lots and lots of counter for the
  // thousand loops
  int nb;       // the important number
  int neig1;    // help neighbour
  int maxFlynn; // maxFlynns (max in program not number of
  // flynns !)
  int count; // a counter yeah !
  int max;   // for max of nodes (also not real node
  // number !)
  int x, y;      // conversion for box
  double fx, fy; // more conversion for box

  Coords bbox[3], boundbox[4], xy; // boxes and single coord structures
  // (attrib.h)
  NodeAttrib *p; // Node structure (nodesP.h)

  // ------------------------------------------------------------------------------------
  // the first routine loops through the flynns, gets the Unode list for
  // flynns and
  // then defines a grain for each particle (an int for which flynn the
  // particle is
  // part of.
  // in order to find the particle the position of the Unode is taken
  // and converted to
  // the relaxation box where we look for the right particle
  // a bit around ten corners but very fast and seems to be efficient
  // ------------------------------------------------------------------------------------

  cout << "define grains" << endl;

  //-------------------------------------------------------
  // This is sooo simple now, Elle does all the work and
  // we just take the flynn number (flynn() function passes
  // that back as an integer)
  //-------------------------------------------------------

  for (j = 0; j < numParticles; j++) {
    runParticle->grain = runParticle->p_Unode->flynn();
    runParticle = runParticle->nextP;
  }

  // -----------------------------------------------------------------------------------------
  // this defines the grain boundary particles by checking whether or
  // not particle has
  // a neighbour that is part of a different grain or flynn.
  // -----------------------------------------------------------------------------------------

  cout << " define grain boundaries " << endl;

  // first we have to kill the old values in order to define new ones
  // The particles have a flag is boundary if they are grain boundaries
  // and also the springs have flags. This is important since we have
  // to track the properties of grains and grain boundaries when these move
  // in Elle due to for example grain growth

  for (i = 0; i < numParticles; i++) {
    if (runParticle->is_boundary == true)
      runParticle->is_boundary = false;

    // grain_young is saving the properties for grains

    runParticle->young = grain_young[runParticle->grain] * 2.0 / sqrt(3.0);
    runParticle->E = grain_young[runParticle->grain] * 2.0 / sqrt(3.0);

    for (j = 0; j < 8; j++) {
      if (runParticle->neigP[j]) // if there is a neighbour
      {
        runParticle->springf[i] =
            sqrt(3.0) * runParticle->E /
            (3.0 * (1.0 - runParticle->poisson_ratio)); // set spring constant
        runParticle->springs[i] = (1.0 - 3.0 * runParticle->poisson_ratio) /
                                  (1.0 + runParticle->poisson_ratio) *
                                  runParticle->springf[i];

        if (runParticle->spring_boundary[j] == true) {
          runParticle->spring_boundary[j] = false;

          runParticle->springf[j] =
              runParticle->springf[j] / boundary_constant; // change
          // spring
          runParticle->break_Str[j] =
              runParticle->break_Str[j] / boundary_strength; // change
        }
      }
    }
    runParticle = runParticle->nextP;
  }

  for (i = 0; i < numParticles; i++) // and loop through particles
  {
    for (j = 0; j < 8; j++) // loop through neighbours of particle
    {
      if (runParticle->neigP[j]) // if there is a neighbour
      {
        if (runParticle->neigP[j]->grain != runParticle->grain) // if
                                                                // not
                                                                // of
                                                                // same
                                                                // grain
        {
          runParticle->is_boundary = true; // I am a boundary
          // particle

          for (j = 0; j < 8; j++) // loop through neighbours
          {
            if (runParticle->neigP[j]) // if neighbour
            {
              if (runParticle->neigP[j]->grain !=
                  runParticle->grain) // if
                                      // neighbour
                                      // is
                                      // in
                                      // different
                                      // grain
              {
                runParticle->spring_boundary[j] = true;

                runParticle->springf[j] =
                    runParticle->springf[j] * boundary_constant; // change
                // spring
                runParticle->break_Str[j] =
                    runParticle->break_Str[j] * boundary_strength; // change
                // breakstrength
              }
            }
          }

          break;
        }
      }
    }
    runParticle = runParticle->nextP; // go on looping
  }

  AdjustConstantGrainBoundaries(); // this is the cleaner to have gradients
                                   // across grains.

  // --------------------------------------------------------------------------------
  // This is the big one. Connects nodes with particles.
  // First put nodes in a Node list similar to the relaxation list so
  // that we can
  // find them easily if we know in what box we want to look
  // --------------------------------------------------------------------------------

  neig1 = -1; // help definition

  // ----------------------------------------------
  // first zero all nodes (= -1)
  // ----------------------------------------------

  for (i = 0; i < 2500000; i++) {

    node_Box[i] = -1;
  }

  max = ElleMaxNodes(); // gets max nodes

  // -----------------------------------------------------------------------------
  // this routine fills the box for nodes. It has the same dimension as
  // the relax box. Therefore it should be easy to find neighbouring
  // nodes
  // and particles using this box. Once it is defined it is rather fast.
  // still have to be careful with negative values.
  // nodes have additional attribute next_inBox (int) that contains tail
  // nodes
  // similar to the tail particles in the relaxbox.
  // -----------------------------------------------------------------------------

  for (i = 0; i < max; i++) // loop through nodes
  {
    next_inBox[i] = -1; // no tail node in Box

    if (ElleNodeIsActive(i)) // if node active
    {
      p = ElleNode(i); // gets structure for the node

      // --------------------------------------------------
      // get position in the box first x and y and make
      // them ints. scale same as particle relax box
      // the box scales with the number of particles
      // --------------------------------------------------

      fx = (p->x * particlex);
      fy = (p->y * particlex);
      x = int(fx);
      y = int(fy);

      // ---------------------------------------------------------
      // and make one-dimensional
      // ---------------------------------------------------------

      x = (y * particlex * 2) + x;

      // ----------------------------------------------------------
      // now fill the node box if it is empty (-1)
      // if it is full loop in the next_inBox array until
      // you find the end of the tail of nodes in that box
      // position
      // -----------------------------------------------------------

      if (node_Box[x] == -1) {
        node_Box[x] = i; // fill the box
      } else             // now we make the tails
      {
        ii = node_Box[x]; // shift to find node with tail

        while (next_inBox[ii] != -1) // loop until end of tail
        {
          ii = next_inBox[ii]; // and shift and shift
        }
        next_inBox[ii] = i; // found tail and fill in node i
      }
    }
  }

  // --------------------------------------------------------------------------------
  // now we come to the connection between particles and nodes
  // loop through particles until you find a boundary particle
  // then check all the neighbours of that particle
  // check the neighbour of the neighbour to find a triangle with the
  // original particle. Then find box position of the particles (plus
  // some
  // extra left and right box otherwise some nodes fall through the
  // net). Then
  // find the nodes with that boxposition in the node box. Then check if
  // that node
  // is really within the triangle that we look for.
  // In the end we take the particles that are boundary of the lattice
  // and
  // define a square around the particle in which we search for nodes on
  // the
  // boundaries.
  // --------------------------------------------------------------------------------

  cout << " and connect nodes to particles " << endl;

  for (i = 0; i < numParticles; i++) // loop through the particles
  {
    neig1 = -1; // flag -1 to find first neighbour

    for (j = 0; j < 32; j++) // maximum of 32 Elle_nodes per
                             // particle at moment
    {
      runParticle->elle_Node[j] = -1;
    }

    {
      for (j = 0; j < 8; j++) // loop through its neighbours
      {
        if (runParticle->neigP[j]) // if there is a neighbour
                                   // at that point
        {
          {
            if (neig1 == -1) // if the first neighbour
              // is not yet defined
              neig1 = j; // i am the first

            else // else have to look for the second
                 // neighbour in triangle
            {

              for (ii = 0; ii < 8; ii++) // loop through
                                         // second
                                         // neighbour
                                         // neighbours
              {
                if (runParticle->neigP[neig1]->neigP[ii]) // if
                                                          // neighbour
                                                          // at
                                                          // that
                                                          // point
                {
                  // -------------------------------------------------------
                  // if the number of the first
                  // neighbour is the same as
                  // the number of the neighbour of the
                  // second neighbour
                  // we have defined a closed triangle
                  // --------------------------------------------------------

                  if (runParticle->neigP[neig1]->neigP[ii]->nb ==
                      runParticle->neigP[j]->nb) {
                    // --------------------------------------------------------
                    // the bbox is used for the
                    // definition of a region for the
                    // function EllePointInRegion that
                    // gets a bbox, the dim of
                    // the bbox (here three) and the
                    // structure of type coord
                    // with the coordinates of the
                    // point (in this case the node)
                    // ----------------------------------------------------------

                    // ----------------------------------------------------------
                    // DEFINE bbox[] by using x and y
                    // of the particles of the
                    // found triangle
                    // -----------------------------------------------------------

                    bbox[0].x = runParticle->xpos;
                    bbox[0].y = runParticle->ypos;
                    bbox[1].x = runParticle->neigP[neig1]->xpos;
                    bbox[1].y = runParticle->neigP[neig1]->ypos;
                    bbox[2].x = runParticle->neigP[j]->xpos;
                    bbox[2].y = runParticle->neigP[j]->ypos;

                    // -------------------------------------------------------------
                    // now we loop through 9 different
                    // box positions and all their
                    // tail nodes and check if any of
                    // these nodes fall into the
                    // triangle
                    // the box positions are the
                    // positions from the three
                    // particles
                    // and the left and right box of
                    // each particle. (thats where
                    // we might loose nodes because
                    // box is square lattice)
                    //
                    // x is always the final box
                    // position
                    // --------------------------------------------------------------

                    for (iii = 0; iii < 9; iii++) {
                      if (iii == 0) // box
                                    // position
                                    // for
                                    // first
                                    // particle
                      {
                        // ---------------------------------------------------------
                        // get positions, scale by
                        // lattice and convert to
                        // ints
                        // ---------------------------------------------------------

                        fx = (runParticle->xpos * particlex);
                        fy = (runParticle->ypos * particlex);
                        x = int(fx);
                        y = int(fy);

                        // ---------------------------------------------------------
                        // and make
                        // one-dimensional
                        // ---------------------------------------------------------

                        x = (y * particlex * 2) + x;
                      } else if (iii == 1)
                        x = x + 1; // box to
                      // right
                      else if (iii == 2)
                        x = x - 2; // box to
                      // left
                      else if (iii == 3) {
                        // ---------------------------------------------------------
                        // get positions, scale by
                        // lattice and convert to
                        // ints
                        // ---------------------------------------------------------

                        fx = (runParticle->neigP[neig1]->xpos * particlex);
                        fy = (runParticle->neigP[neig1]->ypos * particlex);
                        x = int(fx);
                        y = int(fy);

                        // ---------------------------------------------------------
                        // and make
                        // one-dimensional
                        // ---------------------------------------------------------

                        x = (y * particlex * 2) + x;
                      } else if (iii == 4)
                        x = x + 1; // box to
                      // right
                      else if (iii == 5)
                        x = x - 2; // box to
                      // left
                      else if (iii == 6) {
                        // ---------------------------------------------------------
                        // get positions, scale by
                        // lattice and convert to
                        // ints
                        // ---------------------------------------------------------

                        fx = (runParticle->neigP[j]->xpos * particlex);
                        fy = (runParticle->neigP[j]->ypos * particlex);
                        x = int(fx);
                        y = int(fy);

                        // ---------------------------------------------------------
                        // and make
                        // one-dimensional
                        // ---------------------------------------------------------

                        x = (y * particlex * 2) + x;
                      } else if (iii == 7)
                        x = x + 1; // box to
                      // right
                      else if (iii == 8)
                        x = x - 2; // box to
                      // left

                      jj = node_Box[x]; // define
                      // node

                      while (jj != -1) // loop
                                       // through
                                       // tails
                      {
                        // ----------------------------------------------------------
                        // if node is active and
                        // has not been checked go
                        // for it
                        // ----------------------------------------------------------

                        if (ElleNodeIsActive(jj) && !nodes_Check[jj]) {
                          p = ElleNode(jj); // get
                          // struct
                          // of
                          // node
                          xy.x = p->x; // convert
                          // x
                          // to
                          // coord
                          xy.y = p->y; // convert
                          // y
                          // to
                          // coord

                          // ----------------------------------------------------
                          // now make the region
                          // check with the
                          // predefined box,
                          // the dimension of
                          // the box and the xy
                          // coord struct
                          // of the node that
                          // you found in the
                          // node box
                          // -----------------------------------------------------

                          if (EllePtInRegion(bbox, 3, &xy)) {
                            // --------------------------------------------------
                            // there are 16
                            // possible places
                            // for nodes in a
                            // particle. check
                            // if one is free
                            // and put node in
                            // do this for all
                            // three bounding
                            // particles
                            // --------------------------------------------------

                            for (jjj = 0; jjj < 16; jjj++) {
                              // main
                              // particle

                              if (runParticle->elle_Node[jjj] == -1) {
                                runParticle->elle_Node[jjj] = jj;
                                break;
                              }
                            }
                            for (jjj = 0; jjj < 16; jjj++) {
                              // first
                              // neighbour

                              if (runParticle->neigP[neig1]->elle_Node[jjj] ==
                                  -1) {
                                runParticle->neigP[neig1]->elle_Node[jjj] = jj;
                                break;
                              }
                            }
                            for (jjj = 0; jjj < 16; jjj++) {
                              // second
                              // neighbour
                              if (runParticle->neigP[j]->elle_Node[jjj] == -1) {
                                runParticle->neigP[j]->elle_Node[jjj] = jj;
                                break;
                              }
                            }
                            // -----------------------------------------------
                            // now define flag
                            // that tells
                            // particles that
                            // this node has
                            // been connected
                            // otherwise we
                            // connect it at
                            // least three
                            // times
                            // to the same
                            // three particles
                            // ------------------------------------------------

                            nodes_Check[jj] = true;
                          }
                        }
                        // ---------------------------------
                        // this is for the while
                        // loop
                        // go to next tail node in
                        // box
                        // ---------------------------------

                        jj = next_inBox[jj];
                      }
                    }
                  }
                }
              }
              // ---------------------------------------------
              // reset neig1 to next neighbour
              // ---------------------------------------------

              neig1 = j;
            }
          }
        }
      }
    }
    // ------------------------------------------------------------------------------------
    // now check the particles at the boundary of the lattice. Some
    // nodes might not be
    // within a triangle of particles. Therefore draw a square around
    // each boundary
    // particle and check if there is an unconnected node in this box.
    // -------------------------------------------------------------------------------------

    if (runParticle->is_lattice_boundary) // if lattice boundary
    {
      // -------------------------------------
      // define a new bbox with 4 dimensions
      // particle is in the centre
      // -------------------------------------

      boundbox[0].x = runParticle->xpos - (runParticle->radius * 4.0); // lower
      // left
      // corner
      boundbox[0].y = runParticle->ypos + (runParticle->radius * 4.0);

      boundbox[1].x = runParticle->xpos + (runParticle->radius * 4.0); // upper
      // left
      // corner
      boundbox[1].y = runParticle->ypos + (runParticle->radius * 4.0);

      boundbox[2].x = runParticle->xpos + (runParticle->radius * 4.0); // upper
      // right
      // corner
      boundbox[2].y = runParticle->ypos - (runParticle->radius * 4.0);

      boundbox[3].x = runParticle->xpos - (runParticle->radius * 4.0); // lower
      // right
      // corner
      boundbox[3].y = runParticle->ypos - (runParticle->radius * 4.0);

      for (jj = 0; jj < max; jj++) // loop through all nodes
                                   // again
      {
        // ----------------------------------------------------------
        // if node is active and has not been connected
        // ----------------------------------------------------------

        if (ElleNodeIsActive(jj) && !nodes_Check[jj]) {
          p = ElleNode(jj); // get node struct and convert
          // coordinates
          xy.x = p->x;
          xy.y = p->y;

          if (EllePtInRegion(boundbox, 4, &xy)) // if in
                                                // the box
          {
            for (jjj = 0; jjj < 16; jjj++) // loop through
                                           // possible nodes
                                           // in particle
            {
              if (runParticle->elle_Node[jjj] == -1) // if
                                                     // free
              {
                runParticle->elle_Node[jjj] = jj; // fill
                nodes_Check[jj] = true;           // flag up
                break;
              } else if (runParticle->elle_Node[jjj] == jj) // in
                                                            // case
                                                            // node
                                                            // is
                                                            // in
                                                            // there
              {
                break; // we can exit
              }
            }
          }
        }
      }
    }

    runParticle = runParticle->nextP; // loop in particles
  }

  cout << " done with grain definitions " << endl;

  // ----------------------------------------------------------------------
  // just a test to see if all Elle nodes are connected now
  // ----------------------------------------------------------------------

  i = 0;

  for (jj = 0; jj < max; jj++) // loop through Elle nodes
  {
    if (ElleNodeIsActive) {

      /*			void ElleSetNodeAttribute(int node, double val,
      int index)
      {
          int attribindex;
          NodeAttrib *p;

          if ((attribindex=ElleFindNodeAttribIndex(index))==NO_INDX)
              OnError("ElleSetNodeAttribute",RANGE_ERR);
          p = ElleNode(node);
          p->attributes[attribindex] = val;
      }
      double ElleNodeAttribute(int node, int id)
      {
          int attribindex;
          NodeAttrib *p;

          if ((attribindex=ElleFindNodeAttribIndex(id))==NO_INDX)
              OnError("ElleNodeAttribute",RANGE_ERR);
          p = ElleNode(node);
          return(p->attributes[attribindex]);
      }
       */
    }

    if (ElleNodeIsActive(jj) && !nodes_Check[jj]) // if node is
                                                  // active and has
                                                  // not been
                                                  // checked
    {
      p = ElleNode(jj); // define xy coordinates and pass out to
      // interface

      cout << "node not connected to lattice" << jj << " " << p->x << " "
           << p->y << endl;

      i++;
    }
  }
  if (i != 0)
    cout << i << endl; // how many nodes were not connected

  UpdateElle();
}

void Lattice::GetNewElleStructureTwo() {
  int i, j, ii, jj, iii, jjj; // lots and lots of counter for the
  // thousand loops
  int nb;       // the important number
  int neig1;    // help neighbour
  int maxFlynn; // maxFlynns (max in program not number of
  // flynns !)
  int count; // a counter yeah !
  int max;   // for max of nodes (also not real node
  // number !)
  int x, y;      // conversion for box
  double fx, fy; // more conversion for box

  Coords bbox[3], boundbox[4], xy; // boxes and single coord structures
  // (attrib.h)
  NodeAttrib *p; // Node structure (nodesP.h)

  // ------------------------------------------------------------------------------------
  // the first routine loops through the flynns, gets the Unode list for
  // flynns and
  // then defines a grain for each particle (an int for which flynn the
  // particle is
  // part of.
  // in order to find the particle the position of the Unode is taken
  // and converted to
  // the relaxation box where we look for the right particle
  // a bit around ten corners but very fast and seems to be efficient
  // ------------------------------------------------------------------------------------

  cout << "define grains" << endl;

  //-------------------------------------------------------
  // This is sooo simple now, Elle does all the work and
  // we just take the flynn number (flynn() function passes
  // that back as an integer)
  //-------------------------------------------------------

  for (j = 0; j < numParticles; j++) {
    runParticle->grain = runParticle->p_Unode->flynn();
    runParticle = runParticle->nextP;
  }

  // -----------------------------------------------------------------------------------------
  // this defines the grain boundary particles by checking whether or
  // not particle has
  // a neighbour that is part of a different grain or flynn.
  // -----------------------------------------------------------------------------------------

  cout << " define grain boundaries " << endl;

  // first we have to kill the old values in order to define new ones
  // The particles have a flag is boundary if they are grain boundaries
  // and also the springs have flags. This is important since we have
  // to track the properties of grains and grain boundaries when these move
  // in Elle due to for example grain growth

  for (i = 0; i < numParticles; i++) {
    if (runParticle->is_boundary == true)
      runParticle->is_boundary = false;

    // grain_young is saving the properties for grains

    runParticle->young = grain_young[runParticle->grain] * 2.0 / sqrt(3.0);
    runParticle->E = grain_young[runParticle->grain] * 2.0 / sqrt(3.0);

    for (j = 0; j < 8; j++) {
      if (runParticle->neigP[j]) // if there is a neighbour
      {
        runParticle->springf[i] =
            sqrt(3.0) * runParticle->E /
            (3.0 * (1.0 - runParticle->poisson_ratio)); // set spring constant
        runParticle->springs[i] = (1.0 - 3.0 * runParticle->poisson_ratio) /
                                  (1.0 + runParticle->poisson_ratio) *
                                  runParticle->springf[i];

        if (runParticle->spring_boundary[j] == true) {
          runParticle->spring_boundary[j] = false;

          runParticle->springf[j] =
              runParticle->springf[j] / boundary_constant; // change
          // spring
          runParticle->break_Str[j] =
              runParticle->break_Str[j] / boundary_strength; // change
        }
      }
    }
    runParticle = runParticle->nextP;
  }

  for (i = 0; i < numParticles; i++) // and loop through particles
  {
    for (j = 0; j < 8; j++) // loop through neighbours of particle
    {
      if (runParticle->neigP[j]) // if there is a neighbour
      {
        if (runParticle->neigP[j]->grain != runParticle->grain) // if
                                                                // not
                                                                // of
                                                                // same
                                                                // grain
        {
          runParticle->is_boundary = true; // I am a boundary
          // particle

          for (j = 0; j < 8; j++) // loop through neighbours
          {
            if (runParticle->neigP[j]) // if neighbour
            {
              if (runParticle->neigP[j]->grain !=
                  runParticle->grain) // if
                                      // neighbour
                                      // is
                                      // in
                                      // different
                                      // grain
              {
                runParticle->spring_boundary[j] = true;

                runParticle->springf[j] =
                    runParticle->springf[j] * boundary_constant; // change
                // spring
                runParticle->break_Str[j] =
                    runParticle->break_Str[j] * boundary_strength; // change
                // breakstrength
              }
            }
          }

          break;
        }
      }
    }
    runParticle = runParticle->nextP; // go on looping
  }

  AdjustConstantGrainBoundaries(); // this is the cleaner to have gradients
                                   // across grains.

  cout << " done with grain definitions " << endl;

  UpdateElle();
}

// MAKE_REPULSION_BOX
/********************************************************************************
 * function is called from make mike function in the lattice class.
 * It creates a square lattice which lies on top of the triangular one.
 * Uses x and y as number of box. Boxlist is a onedimensional array so that
 * x and y are translated to that array (then x are ones and y are the tens)
 * The box is filled at first with particles. Then a second run checks if there
 * are more than one particle in the box. In the later case the first particle
 * in the box gets a pointer to the next one in the box and the second
 * one a pointer to the third one etc.. (pointer defined in particle.h as
 * next_inBox). Each particle itself has a pointer to the eight boxes around
 * it and its own box (also defined in particle.h as my_Neig_Box).
 *
 * daniel spring and summer 2002
 *****************************************************************************/

// -------------------------------------------------------------------------
// function MakeRepulsionBox in lattice class
//
// called from Lattice::MakeLattice()
// -------------------------------------------------------------------------

void Lattice::MakeRepulsionBox() {
  int i, ii; // counter
  double fx, fy;
  int x, y; // position in box
  int size; // size of box in 1d

  // ---------------------------------------
  // fill main box
  // ---------------------------------------

  size = (particlex * (particlex * 2)) + particlex; // size of box is
  // cubic !

  cout << endl;
  cout << "Nr. of particles along x " << particlex << endl;
  cout << "Total nr. of Repbox " << size << endl;

  // ------------------------------------------------------
  // first we run through all particles and the whole box
  // and just fill the box with particles
  // ------------------------------------------------------

  runParticle = &refParticle;

  for (i = 0; i < numParticles; i++) {
    runParticle->next_inBox = NULL;
    runParticle = runParticle->nextP;
  }

  cout << "fill box" << endl;

  for (ii = 0; ii < numParticles; ii++) // run through
                                        // particlelist
  {

    // -----------------------------------------------------
    // first get integers from position
    // this may be problematic depending on scaling
    // -----------------------------------------------------

    fx = ((runParticle->xpos) * particlex);
    fy = ((runParticle->ypos) * particlex);
    x = int(fx);
    y = int(fy);

    // ----------------------------------------------
    // convert x and y positions to a onedimensional
    // array, x are ones, y are tens
    // -----------------------------------------------

    x = (y * particlex * 2) + x;

    repBox[x] = runParticle;          // pointer to box
    runParticle->box_pos = x;         // position in box for particle
                                      // cout << x << endl;
    runParticle = runParticle->nextP; // run in particle list
  }

  // ------------------------------------------------------------------------
  // now we have to check for secondary and tertiary particles in
  // the boxes. For that we just run through all the particles
  // and check whether they are in the box. If not then another
  // particle is in that box and we add the particle to the local
  // box lists. So in each box there may be one particle and this
  // particle can have endless tails of other particles that are also
  // in the box.
  // ------------------------------------------------------------------------

  cout << "tailparticles" << endl;

  for (i = 0; i < numParticles; i++) // run through list
  {

    // -------------------------------------------------------------
    // again integer x and y
    // -------------------------------------------------------------

    fx = ((runParticle->xpos) * particlex);
    fy = ((runParticle->ypos) * particlex);
    x = int(fx);
    y = int(fy);

    // cout << x << "," << y << endl;

    // ---------------------------------------------------------
    // and make one-dimensional
    // ---------------------------------------------------------

    x = (y * particlex * 2) + x;

    // cout << x << endl;
    // ----------------------------------------------------
    // if the particle is not in the box
    // check the tail of the existing particle
    // ----------------------------------------------------

    if (runParticle->nb != repBox[x]->nb) // if not in box
    {
      // cout << repBox[x]->next_inBox << endl;

      if (repBox[x]->next_inBox != NULL) // does a tailparticle
                                         // exist ?
      {
        // cout << repBox[x]->next_inBox << endl;
        if (repBox[x]->next_inBox->next_inBox != NULL) // does a
                                                       // second
                                                       // tailparticle
                                                       // exist ?
        {

          // ---------------------------------------------------------------------
          // at the moment only four tailparticles can exist
          // a priori. Works so far but might not be scale
          // independent
          // there comes a warning at this point
          // can make this an endless loop then this problem
          // is gone.
          // ---------------------------------------------------------------------

          if (repBox[x]->next_inBox->next_inBox->next_inBox != NULL) // third
          {
            if (repBox[x]->next_inBox->next_inBox->next_inBox->next_inBox !=
                NULL) // four
            {
              cout << "more than four ?" << runParticle->box_pos << endl;
            } else {
              repBox[x]->next_inBox->next_inBox->next_inBox->next_inBox =
                  runParticle;          // tailpointer
              runParticle->box_pos = x; // and position
            }
          } else {
            repBox[x]->next_inBox->next_inBox->next_inBox =
                runParticle;          // tailpointer
            runParticle->box_pos = x; // and position
          }
        } else {

          // ---------------------------------------------------------------
          // in here means we had one tailparticle so that we
          // add our
          // particle to the tail of the first tailparticle
          // ---------------------------------------------------------------

          repBox[x]->next_inBox->next_inBox = runParticle; // tailpointer
          runParticle->box_pos = x;                        // and position
        }
      } else {
        // ------------------------------------------------------------
        // in here we have a particle in the box but no tail so we
        // add our particle to the tail of the one in the box.
        // ------------------------------------------------------------

        repBox[x]->next_inBox = runParticle; // tailpointer
        runParticle->box_pos = x;            // and position
      }
    }
    runParticle = runParticle->nextP; // and next in list
  }
}

// FIND_UNODE_ELLE
/****************************************************************************
 * Function that finds pointer to an Elle Unode.
 * each particle has a pointer to a Unode. Defined in the particle.h class.
 * the function ElleGetParticleUnode takes an int, the Unode number or id and
 * returns the Unode pointer itself to that Unode.
 * function defined in unodes.cc (because need the *Unodes vector)
 *
 * daniel spring 2002
 ***************************************************************************/

// --------------------------------------------------------------------------------
// function FindUnodeElle() of Lattice class
//
// function is called from Lattice::MakeLattice()
// --------------------------------------------------------------------------------

void Lattice::FindUnodeElle() {

  // -----------------------------------------------------------------------------------
  // variables
  // -----------------------------------------------------------------------------------

  int i; // counter variables

  // -----------------------------------------------------------------------------------
  // define starting point
  // -----------------------------------------------------------------------------------

  firstParticle = &refParticle; // get the refParticle, start for
  // particle list

  runParticle = firstParticle; // set pointer for run through
  // list

  // -------------------------------------------------------------------------------------------
  // this loop through all particles calls "Unode
  // *ElleGetParticleUnode(int i)" in unodes.cc of
  // basic Elle functions. The function receives the number of the Unode
  // and returns directly
  // a pointer of type Unode.
  // -------------------------------------------------------------------------------------------

  for (i = 0; i < numParticles; i++) // run through all Unodes and
                                     // particles and connect
  {
    runParticle->p_Unode = ElleGetParticleUnode(i); // get my Unode
    // pointer and
    // define in
    // particle
    runParticle = runParticle->nextP; // next in list
  }

  // ----------------------------------------------------------
  // this is just a test; dont want null pointers
  // ----------------------------------------------------------

  runParticle = firstParticle;

  for (i = 0; i < numParticles; i++) {
    if (!runParticle->p_Unode) // if no Unode in particle
    {
      cout << " no Unode " << runParticle->nb << endl;
    }
    runParticle = runParticle->nextP;
  }
}

/********************************************************************
 * Function for the viscoelastic Deformation, updates the length
 * of springs
 *
 ********************************************************************/

void Lattice::UpdateSpringLength() {
  int j, i, k, rad_nb_buffer;
  float angle, buffer, buffer2, spring_angle_b;

  for (j = 0; j < numParticles; j++) {
    if (!runParticle->is_lattice_boundary) {
      for (i = 0; i < 6; i++) {
        if (runParticle->neigP[i]) {
          buffer = 2 * Pi; // Startwert

          for (k = 0; k < 9; k++) {
            if (k != 4) {
              switch (k) {
              case 0:
                angle = 3.926990817;
                break;
              case 1:
                angle = 4.71238898;
                break;
              case 2:
                angle = 5.497787144;
                break;
              case 3:
                angle = 3.141592654;
                break;
              case 5:
                angle = 0;
                break;
              case 6:
                angle = 2.35619449;
                break;
              case 7:
                angle = 1.570796327;
                break;
              case 8:
                angle = 0.7853981634;
                break;
              }

              spring_angle_b =
                  runParticle->spring_angle[i] + runParticle->rot_angle;

              if (spring_angle_b < 0.0)
                spring_angle_b = 2.0 * 3.14159265 + spring_angle_b;
              else if (spring_angle_b > 2.0 * 3.14159265)
                spring_angle_b = spring_angle_b - 2.0 * 3.14159265;

              // hier suche die spring, die am dichtesten am
              // referenz-winkel liegt, unbeachtet von vz,
              // spruengen,...
              buffer2 = fabs(angle - spring_angle_b);

              if (buffer2 < buffer) {
                buffer = buffer2;
                rad_nb_buffer = k;
              }

              if (k == 5) // dh angle == 0.0
              {
                angle = 2.0 * Pi;

                buffer2 = fabs(angle - spring_angle_b);

                if (buffer2 < buffer) {
                  buffer = buffer2;
                  rad_nb_buffer = k;
                }
              }
            }
          }

          runParticle->rad[i] = runParticle->rep_rad[rad_nb_buffer];
        }
      }
    }

    runParticle = runParticle->nextP;
  }
}

/***********************************************************
 *
 * for the viscoelastic code. Neighbour list is imporant
 * if springs are fractured and the particle is still
 * deforming.
 *
 ************************************************************/
/*
void
Lattice::UpdateNeighList ()
{

  long i, j, k, l;
  int box, pos, neig_nb;
  Particle *neig;

  //runParticle = &refParticle;

  runParticle->done = true;

  // no connected neighbours!


  for (i = 0; i < runParticle->merk; i++)
    {
      runParticle->neigh_list[i]->done = true;
    }

  // update the list of neighbors in the particle itself!
  if (!repBox[runParticle->box_pos])
    {
      cout << " problem replist " << endl;
    }
  else
    {

      // ---------------------------------------------------
      // nine positions in 2d box 3x3
      // ---------------------------------------------------

      for (i = 0; i < 3; i++)
        {
          for (j = 0; j < 3; j++)
            {

              // ------------------------------------------
              // case upper or lower row (+-particlex)
              // ------------------------------------------

              if (j == 0)
                box = 2*particlex * (-1);
              if (j == 1)
                box = 0;
              if (j == 2)
                box = 2*particlex;
              // -----------------------------------------
              // now check this position
              // all particles in this position !!
              // -----------------------------------------

              pos = runParticle->box_pos + (i - 1) + box;

              // -----------------------------------------
              // repBox should always be positive
              // -----------------------------------------

              if (pos > 0)
                {
                  neig = repBox[pos];	// helppointer for while loop

                  while (neig)	// as long as there is somebody
                    {
                      if (!neig->done)
                        {
                          if (neig->xpos)	// just in case
                            {
                              runParticle->neigh_list[runParticle->merk]=neig;
                              runParticle->merk++;
                              neig->done = true;
                            }
                        }

                      if (neig->next_inBox)	// if there is anotherone
                        // in box
                        {
                          neig = neig->next_inBox;	// shift the help
                          // pointer
                        }
                      else
                        {
                          break;	// if not go out to next box
                        }
                    }
                }
            }
        }
    }

  // remove the "done"-flags
  for (i = 0; i < 8; i++)
    {
      if (runParticle->neigP[i])
        runParticle->neigP[i]->done = false;
    }

  for (i = 0; i < runParticle->merk; i++)
    {
      runParticle->neigh_list[i]->done = false;
    }

  runParticle->done = false;

  //if (runParticle->merk > 88)
  //cout << "break" << endl;

}
*/

/**************************************************
 * Clean the neighour list
 *
 ***************************************************/
/*
void
Lattice::EraseNeighList ()
{
  int i, z;

  for (i = 0; i < numParticles; i++)
    {
      for (z = 0; z < 88; z++)
        runParticle->neigh_list[z] = NULL;

      runParticle->merk = 0;

      runParticle = runParticle->nextP;
    }
}
*/

// FULL_RELAX
/********************************************************************
 * This function does a full relaxation. it runs through the particle
 * list and calls the function Relax within each particle (defined in
 * particle.cc. This function will return a true if the particle is moved
 * i.e. is relaxing. Function runs till all particles are relaxed that is
 * they are not moving anymore.
 * Function returns a flag for broken bonds
 *
 * daniel spring 2002
 *********************************************************************/

// -----------------------------------------------------------------------------
// function FullRelax() of Lattice class
// returns flag for broken bonds
//
// called from Lattice::Relaxation()
// -----------------------------------------------------------------------------

bool Lattice::FullRelax() {
  // ---------------------------------------------------
  // local variables
  // ---------------------------------------------------

  bool box_move = false;  // did particle move out of box ?
  bool relaxflag = true;  // flag for relaxation loop
  bool breakBond = false; // second flag for broken bonds
  int i, j = 0;           // loop counters, j counts relaxation
  // steps
  int zaehler = 0;

  int zaehler2 = 0;

  float debug;
  // -----------------------------------------------------------------
  // the relaxation loop, will run until all particles are
  // relaxed. Loops through all particles and calls Relax function
  // in particle objects. Gets a flag from each object indicating
  // whether or not it had to be relaxed.
  // -----------------------------------------------------------------

  // cout << "relax" << endl;

  int preFirst_par_up_row = numParticles - 1 - particlex;
  // cout << right_wall_pos << " " << left_wall_pos << " " << lower_wall_pos <<
  // " " <<
  //  upper_wall_pos << endl;

  runParticle = &refParticle; // go to beginning of list

  debug = 0.0;

  // Tilt ();

  // UpdateSpringLength();
  // EraseNeighList();
  // UpdateNeighList ();
  do {

    relaxflag = false; // new relaxation loop if flag stays false
    // jumps out of loop
    j++; // count relaxation steps

    // cout << j << endl; // use this if you want to check how fast
    // relaxation goes

    // -------------------------------------------------------------------------------
    // loop through all particles and call Relax in the particle class
    // and set flags
    // note that we do not have to redefine the beginning of the list
    // since it is
    // circular, at the end the pointer points again to the first
    // particle
    // -------------------------------------------------------------------------------

    for (i = 0; i < numParticles; i++) // loop through all
                                       // particles
    {
      // ------------------------------------------------------------------------------
      // this does a number of things. It calls the function Relax
      // in runParticle
      // and gives that function the relaxationthreshold defined in
      // the lattice class
      // (see the constructor). It receives a flag that is false if
      // the particle
      // was moved i.e. relaxed. If the flag is false the relaxflag
      // is turned true
      // which means that after all particles are done a new
      // relaxation loop will start
      // -------------------------------------------------------------------------------
      if (!runParticle->isHole) {
        if (!runParticle->Relax(relaxthres, repBox, particlex, walls,
                                right_wall_pos, left_wall_pos, lower_wall_pos,
                                upper_wall_pos, wall_rep, false, debug_nb,
                                visc_rel, sheet, grav, grav_Press, i,
                                preFirst_par_up_row))
          relaxflag = true;

        // ----------------------------------------------------------------------------
        // check if particle moved out of box, if it did call
        // function
        // to go to new position in box
        // ----------------------------------------------------------------------------

        runParticle->ChangeBox(repBox, particlex);
      }

      runParticle = runParticle->nextP; // continue in particle
      // list
    }

    if (j == 500) {
      relaxflag = false;
      cout << "relaxation 500" << endl;
    }

  } while (relaxflag);

  // cout << endl << endl;

  // loop until relaxflag is staying false
  // through a whole relaxation

  // this loop through all particles calculates the stress tensor for
  // each particle and
  // also determines which bonds will break
  // took that out of the relaxroutine so that that is really simple and
  // as fast as it can be.

  for (i = 0; i < numParticles; i++) {
    // call routine in particles that calculates stress and returns a
    // true if
    // a bond is broken, loops only once since the model is relaxed

    if (!runParticle->isHole) {
      runParticle->Relax(repBox, particlex, visc_rel, walls, right_wall_pos,
                         left_wall_pos, lower_wall_pos, upper_wall_pos,
                         wall_rep);
    }

    runParticle = runParticle->nextP;
  }

  // this is new for the new calculation of normal and shear forces and shear
  // failure
  breakBond = Standard_Spring_Failure();

  return breakBond; // and return whether or not a bond was
                    // broken
}

//-----------------------------------------------------------------------
// written by Till 2014/15
//
// includes shear failure, but this still needs a bit of debugging
// daniel 2018
//-----------------------------------------------------------------------

bool Lattice::Standard_Spring_Failure() {
  double break_prob;   // probability of spring to break at node 1
  bool bbreak = false; // does a bond break at all?

  Particle *prtcl = &refParticle;

  for (int j = 0; j < numParticles; j++) {
    prtcl->mbreak = 0.0;

    // connected first neighbours

    for (int i = 0; i < 8; i++) {
      if (prtcl->neigP[i]) {
        break_prob = Standard_Breaking_Probability(prtcl, i);

        if (break_prob > 1.0) {
          if (break_prob > prtcl->mbreak) {
            if (!prtcl->no_break[i]) {
              prtcl->mbreak = break_prob;

              prtcl->neigh = i;
              bbreak = true;
            }
          }
        }
      }
    }
    prtcl = prtcl->nextP;
  }

  return (bbreak);
}

//-------------------------------------------------------------
// this simply turns off the shear component of springs
// makes some of the codes more stable Daniel sept 2018
//
//-------------------------------------------------------------

void Lattice::Only_Extension(bool extension) {
  if (extension)
    shear = false;
  else
    shear = true;
}

//-----------------------------------------------------------------------
// calculation of breaking probabilty after Sachau and Koehn 2013
// check the cohesion, made some changes recently daniel 2018
//-----------------------------------------------------------------------

double Lattice::Standard_Breaking_Probability(Particle *prtcl, int nb) {
  double mc; // shear cohesion

  double break_prob = 0.0; // probability to break

  double E, v, G, ex;
  double usx, usy;
  double unx, uny;
  double dx, dy;
  double ux, uy;
  double nx, ny;
  double alen, dd, ds;
  double f;

  double normal_stress, shear_stress, break_Str;

  break_Str = prtcl->break_Str[nb];
  mc = -(break_Str * (1.0 + sin(prtcl->angle_of_internal_friction)));
  mc = mc / (2.0 * cos(prtcl->angle_of_internal_friction));
  mc = fabs(mc);
  mc = mc * 12.0; // daniel test was 3.0

  //	mc *= model->shear_cohesion_scale;

  dx = prtcl->neigP[nb]->xpos - prtcl->xpos;
  dy = prtcl->neigP[nb]->ypos - prtcl->ypos;

  dd = sqrt((dx * dx) + (dy * dy));

  alen = sqrt(prtcl->neigpos[nb][0] * prtcl->neigpos[nb][0] +
              prtcl->neigpos[nb][1] * prtcl->neigpos[nb][1]);

  // alen = prtcl->radius + prtcl->neigP[nb]->radius;   // daniel back to normal
  // breaking 2018

  if (dd != 0.0) {
    nx = dx / dd;
    ny = dy / dd;
  } else {
    nx = 0;
    ny = 0;
  }

  ux = dx - prtcl->neigpos[nb][0];
  uy = dy - prtcl->neigpos[nb][1];

  f = ux * nx + uy * ny;
  unx = nx * f;
  uny = ny * f;

  usx = ux - unx;
  usy = uy - uny;

  ds = sqrt(usx * usx + usy * usy);

  E = (prtcl->E + prtcl->neigP[nb]->E) / 2.0;
  v = (prtcl->poisson_ratio + prtcl->neigP[nb]->poisson_ratio) / 2.0;

  // Shear modulus
  G = E / (2.0 * (1.0 + v));
  // shear stress
  shear_stress = fabs(G * ds / alen);

  // strain parallel spring
  ex = (alen - dd) / alen;
  // uniaxial stress spring
  normal_stress = ex * E;

  if (normal_stress > 0.0)
    normal_stress = 0.0;

  // tensile failure only in uppermost crust
  //	if (griffith_model)
  //		if (prtcl->diff_stress < 4.0 * fabs(break_Str))
  //			normal_stress = 0.0;

  // if shear failure has been excluded
  //	if (!shear_failure)
  //		shear_stress = 0.0;

  // probability (elliptical fracture criterion):
  if (shear)
    break_prob = (normal_stress / break_Str) * (normal_stress / break_Str) +
                 (shear_stress / mc) * (shear_stress / mc);
  else
    break_prob = (normal_stress / break_Str) * (normal_stress / break_Str);

  return (break_prob);
}

// SET_BOUNDARIES
/*************************************************************************************
 * This function sets whatever boundary conditions we want. At the moment it
 *just fixes the y position of the upper and lower row so that these are not
 * moved during the relaxation routine
 *
 * daniel spring 2002
 *
 * bug fixed for even and uneven lattices Daniel december 2004
 ************************************************************************************/

// -----------------------------------------------------------------------------------
// function SetBoundaries() of Lattice class
//
// called from Lattice::MakeLattice()
// -----------------------------------------------------------------------------------

void Lattice::SetBoundaries() {
  // ------------------------------------
  // local variables
  // ------------------------------------

  int i, j, jj; // loop counter

  runParticle = firstParticle = &refParticle; // set pointers to
  // beginning of list
  runParticle->is_left_lattice_boundary = true;
  runParticle->is_bottom_lattice_boundary = true;
  runParticle->fix_y = true;
  // ---------------------------------------------------------------------------------
  // set flag in particle object telling it not to move its y position
  // set this flag first for the lower row == particles form 0 to
  // particlex-1
  // ---------------------------------------------------------------------------------

  for (i = 0; i < particlex; i++) {
    runParticle->fix_y = true;               // set y
    runParticle->is_lattice_boundary = true; // define as part
    runParticle->is_bottom_lattice_boundary = true;
    // of boundary

    // this is a messy bit that I have to change... at the moment
    // boundary particles
    // cannot have broken bonds... otherwise it crashes.
    // have to define walls some day

    for (j = 0; j < 8; j++) {
      if (runParticle->neigP[j]) {
        runParticle->no_break[j] = true;
        // runParticle->break_Str[j] = 1000000.0;

        for (jj = 0; jj < 8; jj++) {
          if (runParticle->neigP[j]->neigP[jj]) {
            if (runParticle->neigP[j]->neigP[jj]->nb == runParticle->nb) {

              // make it unbreakable

              runParticle->neigP[j]->no_break[jj] = true;
            }
          }
        }
      }
    }
    runParticle = runParticle->nextP;
  }

  // --------------------------------------------------------------------------------
  // then go to the last particle and run backwards particlex times to
  // fix
  // the upper row of particles
  // --------------------------------------------------------------------------------
  firstParticle = &refParticle;
  runParticle = firstParticle->prevP;

  for (i = 0; i < particlex; i++) {
    runParticle->fix_y = true;
    runParticle->is_lattice_boundary = true;
    runParticle->is_top_lattice_boundary = true;

    for (j = 0; j < 8; j++) {
      if (runParticle->neigP[j]) {
        runParticle->no_break[j] = true;

        for (jj = 0; jj < 8; jj++) {
          if (runParticle->neigP[j]->neigP[jj]) {
            if (runParticle->neigP[j]->neigP[jj]->nb == runParticle->nb) {
              runParticle->neigP[j]->no_break[jj] = true;
            }
          }
        }
      }
    }
    runParticle = runParticle->prevP;
  }

  // ------------------------------------------------------------------
  // now do the same for the side boundaries (set x not to move)
  // ------------------------------------------------------------------

  runParticle = firstParticle;             // first do first particle
  runParticle->fix_x = true;               // fix x
  runParticle->is_lattice_boundary = true; // define as part of
  runParticle->is_left_lattice_boundary = true;
  // boundary

  do {
    for (i = 0; i < particlex; i++) // loop two rows
    {
      runParticle = runParticle->nextP;
    }

    // fix left boundary

    runParticle->fix_x = true;
    runParticle->is_lattice_boundary = true;
    runParticle->is_left_lattice_boundary = true;
    for (j = 0; j < 8; j++) {
      if (runParticle->neigP[j]) {
        runParticle->no_break[j] = true;

        for (jj = 0; jj < 8; jj++) {
          if (runParticle->neigP[j]->neigP[jj]) {
            if (runParticle->neigP[j]->neigP[jj]->nb == runParticle->nb) {
              runParticle->neigP[j]->no_break[jj] = true;
            }
          }
        }
      }
    }

    // fix right boundary

    runParticle->prevP->fix_x = true;
    runParticle->prevP->is_lattice_boundary = true;
    runParticle->prevP->is_right_lattice_boundary = true;
    for (j = 0; j < 8; j++) {
      if (runParticle->prevP->neigP[j]) {
        runParticle->prevP->no_break[j] = true;

        for (jj = 0; jj < 8; jj++) {
          if (runParticle->prevP->neigP[j]->neigP[jj]) {
            if (runParticle->prevP->neigP[j]->neigP[jj]->nb ==
                runParticle->prevP->nb) {
              runParticle->prevP->neigP[j]->no_break[jj] = true;
            }
          }
        }
      }
    }

  } while (runParticle != firstParticle); // loop once through all
}

/************************************************************
 * Function that sets "wall boundaries" around the particle
 * box. Walls are straight elastic lines at the moment
 *
 * Walls are only repulsive at the moment, therefore
 * they cannot build up tensile strain !
 *
 * Daniel 2005/6
 ************************************************************/

void Lattice::SetWallBoundaries(int both, float constant) {
  // ------------------------------------
  // local variables
  // ------------------------------------

  int i, j; // loop counter

  if (both == 0) {
    for (i = 0; i < numParticles; i++) {
      runParticle->fix_y = false; // set y
      runParticle->fix_x = false;
      // runParticle->is_lattice_boundary = false;	// define as part
      //  of boundary

      for (j = 0; j < 8; j++) {
        if (runParticle->neigP[j])
          runParticle->no_break[j] = false;
      }

      runParticle = runParticle->nextP;
    }
  }

  walls = true;

  // positions of walls

  cout << "nr. & ypos of 1st particle in lower 1st row" << " "
       << firstParticle->nb << " " << firstParticle->ypos << endl;
  // lower_wall_pos = firstParticle->nextP->ypos - firstParticle->ypos;
  lower_wall_pos = firstParticle->nextP->ypos - firstParticle->radius;

  cout << "nr. & xpos of 1st particle in lower 1st row" << " "
       << firstParticle->nb << " " << firstParticle->xpos << endl;
  for (i = 0; i < particlex; i++) {
    firstParticle = firstParticle->nextP;
  }
  cout << "nr. & xpos of 1st particle in lower 2nd row" << " "
       << firstParticle->nb << " " << firstParticle->xpos << endl;

  left_wall_pos = firstParticle->xpos - firstParticle->radius;

  runParticle = firstParticle->prevP;

  cout << "nr. & ypos of last particle in upper last row" << " "
       << runParticle->nb << " " << runParticle->ypos << endl;
  upper_wall_pos = runParticle->ypos + runParticle->radius;

  cout << "nr. & xpos of last particle in upper last row" << " "
       << runParticle->nb << " " << runParticle->xpos << endl;
  right_wall_pos = runParticle->xpos + runParticle->radius;

  for (i = 0; i < particlex; i++) {
    runParticle = runParticle->prevP;
  }
  if ((runParticle->xpos + runParticle->radius) > right_wall_pos)
    right_wall_pos = runParticle->xpos + runParticle->radius;
  cout << "nr. & xpos of last particle in upper 2nd-last row" << " "
       << runParticle->nb << " " << runParticle->xpos << endl;

  cout << "wallpos   " << lower_wall_pos << " lower " << upper_wall_pos
       << "upper " << endl;
  cout << left_wall_pos << " left " << right_wall_pos << " right " << endl;

  wall_rep = constant; // elastic constant of wall
}

void Lattice::SetWallBoundariesX(int release, int fracture, int walls,
                                 float constant) {
  // ------------------------------------
  // local variables
  // ------------------------------------

  int i, j;

  for (i = 0; i < numParticles; i++) {
    if (runParticle->fix_y == true) {
      if (release)
        runParticle->fix_y = false;
      if (runParticle->fix_x == false) {
        if (fracture)

        {
          runParticle->is_lattice_boundary = false; // define as part
                                                    // of boundary

          for (j = 0; j < 8; j++) {
            if (runParticle->neigP[j])
              runParticle->no_break[j] = false;
          }
        }
      }
    }

    runParticle = runParticle->nextP;
  }

  if (walls) {
    walls = true;

    // positions of walls

    cout << "nr. & ypos of 1st particle in lower 1st row" << " "
         << firstParticle->nb << " " << firstParticle->ypos << endl;
    // lower_wall_pos = firstParticle->nextP->ypos - firstParticle->ypos;
    lower_wall_pos = firstParticle->nextP->ypos - firstParticle->radius;

    cout << "nr. & xpos of 1st particle in lower 1st row" << " "
         << firstParticle->nb << " " << firstParticle->xpos << endl;
    for (i = 0; i < particlex; i++) {
      firstParticle = firstParticle->nextP;
    }
    cout << "nr. & xpos of 1st particle in lower 2nd row" << " "
         << firstParticle->nb << " " << firstParticle->xpos << endl;

    left_wall_pos = firstParticle->xpos - firstParticle->radius;

    runParticle = firstParticle->prevP;

    cout << "nr. & ypos of last particle in upper last row" << " "
         << runParticle->nb << " " << runParticle->ypos << endl;
    upper_wall_pos = runParticle->ypos + runParticle->radius;

    cout << "nr. & xpos of last particle in upper last row" << " "
         << runParticle->nb << " " << runParticle->xpos << endl;
    right_wall_pos = runParticle->xpos + runParticle->radius;

    for (i = 0; i < particlex; i++) {
      runParticle = runParticle->prevP;
    }
    if ((runParticle->xpos + runParticle->radius) > right_wall_pos)
      right_wall_pos = runParticle->xpos + runParticle->radius;
    cout << "nr. & xpos of last particle in upper 2nd-last row" << " "
         << runParticle->nb << " " << runParticle->xpos << endl;

    cout << "wallpos   " << lower_wall_pos << " lower " << upper_wall_pos
         << "upper " << endl;
    cout << left_wall_pos << " left " << right_wall_pos << " right " << endl;

    wall_rep = constant; // elastic constant of wall
  }
}

// RELAXATION
/*********************************************************************************
 * function that does a full relaxation
 * calls FullRelax() in Lattice class. That function does a full relaxatiion
 * until the lattice is relaxed and returns a flag wether or not a new bond is
 * broken. If that is the case the function searches for the broken bond which
 * is most probable to break, breaks it and starts the relaxation again until
 *the lattice is fully relaxed and all possible bonds are broken
 *
 * daniel spring 2002
 **********************************************************************************/

// ----------------------------------------------------------------------------------
// function Relaxation() of Lattice class
//
// called in Lattice:MakeLattice() and Lattice:DeformLattice()
// ----------------------------------------------------------------------------------

void Lattice::Relaxation() {
  // -------------------------------------------------------
  // some local variables
  // -------------------------------------------------------

  float maxbreak; // probability of bond to break
  int i, z;       // loop counter,z for debugging

  BoxRad(); // function writes splined radius-values

  if (internal_break == 1)
    local_nbBreak = 0;

  // -------------------------------------------------------------------------------------
  // loop through the full relaxation routine until no more bonds break
  // -------------------------------------------------------------------------------------

  nbBreak_last = nbBreak;

  while (FullRelax()) {
    // -----------------------------------------------------------------
    // if in here a bond wants to break
    // -----------------------------------------------------------------
    nbBreak++;       // one more broken bond
    local_nbBreak++; // counter for fracture plot

    // --------------------------------------------------------------------
    // This plots a picture in between timesteps if the user defines
    // a number larger than zero. Each time the code fractures
    // internally
    // that is within a time step n times it plots a picture.
    // --------------------------------------------------------------------

    if (fract_Plot_Numb != 0) // in this case (default) no
                              // internal picts taken
    {
      if (local_nbBreak == fract_Plot_Numb) // otherwise if
                                            // more fractures
                                            // than n
      {
        UpdateElle();
        local_nbBreak = 0; // and reset pointer
      }
    }

    cout << "Broken Bonds " << nbBreak << " , "; // tell the world

    // -----------------------------------------------------------------------
    // now run through particles and find particle that contains the
    // bond
    // with the highest probability to break. Only break that bond.
    // What the other do that would break also will be seen after
    // the next relaxation routine. Want to relax lattice after we
    // break a bond
    // ------------------------------------------------------------------------

    maxbreak = 0.0; // set to zero

    runParticle = &refParticle; // beginning of list

    for (i = 0; i < numParticles;) // run through list
    {
      if (runParticle->mbreak > maxbreak) // if the particle max
                                          // breaking probablity is
                                          // higher
      {
        maxbreak = runParticle->mbreak; // thats the new number
        preRunParticle = runParticle;   // and catch the particle
      }
      runParticle = runParticle->nextP;
      i++;
    }

    preRunParticle->BreakBond(); // call BreakBond() in the object
    preRunParticle->fracture_time = def_time;

    if (visc_flag == true)
      preRunParticle->visc_flag = true; // debugging
  }

  nbBreak_current = nbBreak - nbBreak_last;
}

void Lattice::Set_TimeFrac(int time) {
  def_time = time;
  // cout<<time<<endl;
}

void Lattice::Zonally_brkBond_Cal() {
  int limit;

  if (preRunParticle->ypos <
      P->LowSeal_y) // broken bond's particle y_pos w.r.t
                    // lower y_limit of hor_seal in pressure_lattice class
    limit = 1;
  else if (preRunParticle->ypos >= P->LowSeal_y &&
           preRunParticle->ypos <= P->UppSeal_y) // w.r.t lower and
                                                 // upper y_limits of hor_seal
                                                 // in pressure_lattice class
    limit = 2;
  else
    limit = 3;

  switch (limit) {
  case 1: // cout << "its case 1" << endl;
    LowSeal_bnd_nb++;
    cout << " in layer 1" << endl;
    break;

  case 2: // cout << "its case 2" << endl;
    Seal_bnd_nb++;
    cout << " in layer 2" << endl;
    break;

  case 3: // cout << "its case 3" << endl;
    UppSeal_bnd_nb++;
    cout << " in layer 3" << endl;
    break;

  default:
    cout << "no case present" << endl;
  }
}

void Lattice::HalfMod_brkBond_Cal() {
  if (preRunParticle->ypos <
      P->ModHalf_y) // broken bond's particle y_pos w.r.t
                    // lower part of model in pressure_lattice class
  {
    ModLowr_bnd_nb++;
    cout << " in Lower layer" << endl;
  } else if (preRunParticle->ypos >= P->ModHalf_y) // w.r.t upper part of model
                                                   // in pressure_lattice class
  {
    ModUppr_bnd_nb++;
    cout << " in Upper layer" << endl;
  } else
    cout << "no case present" << endl;
}

void Lattice::SetRightLatticeWallPos() {

  int i; // counter for loop

  for (i = 0; i < numParticles; i++) {

    runParticle->right_lattice_wall_pos = right_lattice_wall_pos;
    // runParticle->left_lattice_wall_pos = left_lattice_wall_pos;

    runParticle = runParticle->nextP;
  }
}

// UPDATE_ELLE
/**************************************************************************************
 * this function writes the new values to the Unode objects for elle plotting
 * and also adjusts positions of Unodes
 * Calls ElleUpdate to talk to the Elle interface
 *
 * daniel summer and spring 2002
 * Rianne, Daniel and Jochen December 2002, Jan 2003
 *
 * add security stop
 *
 * Daniel March 2003
 **************************************************************************************/

// --------------------------------------------------------------------------------
// function UpdateElle() of Lattice class
//
// called from Lattice::MakeLattice() and Lattice::DeformLattice()
// --------------------------------------------------------------------------------

void Lattice::UpdateElle() {
  // ----------------------------------------------------------
  // local variables
  // ----------------------------------------------------------

  double pressure; // variable for the sum of the trace of
  // the strain tensor sxx+syy
  Coords xy; // coordinate structure defined in
  // attrib.h
  int i, ii, max; // counter
  NodeAttrib *p;  // for nodes (x and y)
  float smax, smin, sxx, syy, sxy,
      energy; // some local variables for stresses etc.
  float uxy, uxx, uyy, lame1, lame2, pois;
  double frac;

  bool writeCSV = false;
  // flag to set if csv will be written or not
  int timestep_from_melt;
  // start counting timesteps from the moment I add melt

  pois = 0.33333;

  runParticle = &refParticle; // jump to beginning of list

  timestep_from_melt = experimentTime;
  // start the time from when I first produce melt

  // check if I need to write CSV output files

  if (timestep_from_melt % saveCSVInterval == 0 and timestep_from_melt >= 0)
  // start the time from when I first produce melt
  {
    // if (experimentTime>1) // the first two files are useless - don't need
    // this condition anymore
    //{
    writeCSV = true;
    // flag to say: save csv file!
    string fileroot = "my_experiment";
    // set the name of the files
    char csvFileName[1024];
    // initialise with some size

    // Name of the csv file:
    sprintf(csvFileName, "%s%05d%s", fileroot.c_str(), timestep_from_melt,
            ".csv");
    // take the root name, add the number of zeros you need to fill 5 spaces,
    // put the timestep, add the extension .csv
    myExpCsv.open(csvFileName, ios::out);
    cout << "Saving file: " << csvFileName << endl;
    //}
  }

  max = ElleMaxNodes();

  for (ii = 0; ii < max; ii++) {
    nodes_Check[ii] = false;
  }

  // --------------------------------------------------------------------------------------
  // loop through the list of particles and change attributes and
  // positions in the Undoes
  // has to be called everytime something is changed and we talk to elle
  // again
  // this should include all kinds of plotting things, stresses,
  // energies, fractures etc.
  // --------------------------------------------------------------------------------------

  for (i = 0; i < numParticles; i++) // and loop again
  {
    if (runParticle->p_Unode) {
      // otherwise stuff gets too long

      sxx = runParticle->sxx;
      sxy = runParticle->sxy;
      syy = runParticle->syy;

      //----------------------------------------------------------------------
      // setting elastic constants
      //----------------------------------------------------------------------

      lame1 = runParticle->young * pois / (1.0 - (2.0 * pois) * (1.0 + pois));
      lame2 = runParticle->young / (2.0 + 2.0 * pois);

      //-----------------------------------------------------------------------
      // getting strain tensor
      //-----------------------------------------------------------------------

      uxy = ((1 + pois) / runParticle->young) * runParticle->sxy;
      uxx = (1 / runParticle->young) *
            (runParticle->sxx - (runParticle->syy * pois));
      uyy = (1 / runParticle->young) *
            (runParticle->syy - (runParticle->sxx * pois));

      //-----------------------------------------------------------------------
      // calculate elastic energy for particle
      //-----------------------------------------------------------------------

      runParticle->eEl = 0.5 * lame1 * (uxx + uyy) * (uxx + uyy);
      runParticle->eEl =
          runParticle->eEl +
          (lame2 * ((uxx * uxx) + (uyy * uyy) + 2 * (uxy * uxy)));

      //---------------------------------------------------------
      // scaling energy
      //---------------------------------------------------------

      runParticle->eEl =
          runParticle->eEl * runParticle->mV * pressure_scale * pascal_scale;

      // runParticle->p_Unode->setAttribute(U_ENERGY, runParticle->eEl);

      // runParticle->p_Unode->setAttribute(U_ENERGY,
      // 1000000*sqrt((runParticle->v_darcyy*runParticle->v_darcyy)+(runParticle->v_darcyx*runParticle->v_darcyx)));
      runParticle->p_Unode->setAttribute(U_ENERGY,
                                         runParticle->v_darcyy * 1000000000);

      // runParticle->p_Unode->setAttribute(U_ENERGY, runParticle->rate_factor);

      // plot healing time

      // runParticle->p_Unode->setAttribute(U_ATTRIB_A, runParticle->healing);

      runParticle->p_Unode->setAttribute(U_ATTRIB_A, runParticle->Mg);

      // plot the fracture time
      runParticle->p_Unode->setAttribute(U_ATTRIB_B, runParticle->SheetT);
      // runParticle->p_Unode->setAttribute(U_ATTRIB_B,
      // runParticle->fracture_time);
      // runParticle->p_Unode->setAttribute(U_ATTRIB_B,
      // runParticle->fluid_particles);

      runParticle->p_Unode->setAttribute(U_ATTRIB_C, runParticle->Mg_fluid);

      runParticle->p_Unode->setAttribute(U_DISLOCDEN, runParticle->Fe_fluid);
      // -----------------------------------------------------------------------
      // stress tensor is set in the particle object themselves at
      // relaxation
      // therefore call particles, get stresses and calculate whatever
      // you want, then pass to Unodes.
      // here now pressure = trace of tensor sxx + syy
      // -----------------------------------------------------------------------

      pressure = (sxx + syy) / 2.0; // calc

      // ----------------------------------------------------------------------
      // this function setAttribute is defined in unodeP.h which
      // includes the
      // definition for the class plus all the functions. There is no
      // unodesP.cc in Elle
      // this function wants a reference int for which value we want to
      // change that is CONC_A at the moment (normally used for
      // concentration)
      // and also wants the value here as a double.
      //
      // this call just calls the function setAttribute in the Unode
      // connected to the particle
      // -----------------------------------------------------------------------
      runParticle->p_Unode->setAttribute(U_MEAN_STRESS, pressure);

      // runParticle->p_Unode->setAttribute(U_DENSITY,
      // runParticle->real_density);
      runParticle->p_Unode->setAttribute(U_DENSITY, runParticle->particle_age);

      runParticle->p_Unode->setAttribute(U_YOUNGSMODULUS, runParticle->young);

      // runParticle->p_Unode->setAttribute
      // (U_VISCOSITY,runParticle->crack_seal_counter);

      runParticle->p_Unode->setAttribute(U_VISCOSITY, runParticle->healing);

      // -------------------------------------------------------------------------
      // determine smax - smin (Eigenvalues)
      // -------------------------------------------------------------------------

      smax = ((sxx + syy) / 2.0) +
             sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) + sxy * sxy);

      smin = ((sxx + syy) / 2.0) -
             sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) + sxy * sxy);

      pressure = smax - smin;

      //~ runParticle->p_Unode->setAttribute(U_MEAN_STRESS, pressure);
      runParticle->p_Unode->setAttribute(U_DIF_STRESS, pressure);

      frac = 0.0;
      if (runParticle->draw_break) {
        if (runParticle->no_spring_flag)
          frac = -0.3;
        else if (runParticle->visc_flag)
          frac = -0.6;
        else
          frac = -1.0;

        runParticle->p_Unode->setAttribute(U_FRACTURES, frac); // something
                                                               // to
                                                               // visualize
                                                               // fractures
      } else
        runParticle->p_Unode->setAttribute(U_FRACTURES, frac);

      pressure = 0.0;

      // if (!transition) {

      //  if (runParticle->isHoleBoundary)
      //   pressure = 1.0;

      //  if (runParticle->isHole)
      // pressure = runParticle->conc;
      //   pressure = -1;

      // runParticle->p_Unode->setAttribute(U_PHASE, pressure);
      //   }
      // else {
      //  pressure = (runParticle->mineral - 1);
      // runParticle->p_Unode->setAttribute(U_PHASE, pressure);
      //   }

      if (runParticle->isHoleBoundary)
        pressure = 1.0;

      if (runParticle->isHole)
        // pressure = runParticle->conc;
        pressure = -1;

      // runParticle->p_Unode->setAttribute (U_DENSITY,
      // runParticle->crack_seal_counter);
      runParticle->p_Unode->setAttribute(U_ATTRIB_C, pressure);

      runParticle->p_Unode->setAttribute(U_TEMPERATURE, runParticle->conc);

      // runParticle->p_Unode->setAttribute(U_TEMPERATURE,
      // runParticle->rate_factor);

      runParticle->p_Unode->setAttribute(U_ATTRIB_A, runParticle->abreak);

      // runParticle->p_Unode->setAttribute(U_ATTRIB_B, runParticle->conc);

      if (runParticle->draw_break) {
        smin = smin * 10.0; // something to visualize fractures
      }

      pressure = smax;

      // runParticle->p_Unode->setAttribute(U_CONC_A,
      // runParticle->crack_seal_counter);

      // runParticle->p_Unode->setAttribute(ENERGY,pressure);

      // pressure = arctan ((smax-sxx)/sxx);

      // runParticle->p_Unode->setAttribute(DISLOCATION_DENSITY,pressure);

      // -------------------------------------------------------------------------
      // and now also change the Unode position because the particles
      // have moved
      // I am not exactly shure what the consequences of that for Elle
      // are, maybe
      // something crashes at some point.
      // since the function setPosition in the Unode class wants a Coord
      // structure defined in attrib.h in Elle we create one and give it
      // the
      // x and y positions in our particle and then pass it to the Unode
      // function
      // -------------------------------------------------------------------------

      xy.x = runParticle->xpos; // pass x
      xy.y = runParticle->ypos; // pass y

      runParticle->p_Unode->setPosition(&xy); // and pass that to the
      // Unode

      // ------------------------------------------------------------------------------
      // This routine moves Elle nodes according to movements in the
      // particle lattice
      // At the moment a bit rough, node just moved on top of particle
      // depends a bit on resolution, make triangulation of unodes crash
      // !
      // ------------------------------------------------------------------------------

      for (ii = 0; ii < 32; ii++) {
        if (runParticle->elle_Node[ii] != -1) {
          if (!nodes_Check[runParticle->elle_Node[ii]]) {
            ElleNodePosition(runParticle->elle_Node[ii], &xy);

            nodes_Check[runParticle->elle_Node[ii]] = true;

            xy.y = xy.y + (runParticle->ypos - runParticle->oldy);
            xy.x = xy.x + (runParticle->xpos - runParticle->oldx);
            ElleSetPosition(runParticle->elle_Node[ii], &xy);
          }
        }
      }
      runParticle->oldx = runParticle->xpos;
      runParticle->oldy = runParticle->ypos;

      if (writeCSV) {
        // save all of this in a Map that will be written in a csv file
        setMapForCSV(i);
      }
    }
    runParticle = runParticle->nextP; // and the next particle
  }

  // Now that it finished with all of the particle:
  if (writeCSV) {
    myExpCsv.close();
  }
  doneWithHeaders = 0;

  // ---------------------------------------------------------------------------
  // call ElleUpdate() an Elle function that updates the interface
  // then the new values will be plotted
  // if security stop is set and max number of picts is reached dont
  // plot
  // pict anymore.
  // ---------------------------------------------------------------------------

  if (!set_max_pict) {
    // cout << "interface" << endl;
    ElleUpdate();
  } else if (num_pict < max_pict) {
    // cout << "interface" << endl;
    ElleUpdate();
    max_pict = max_pict + 1;
  }
}

void Lattice::getTimeForCsv(int expTime, int saveInt) {
  experimentTime = expTime;
  saveCSVInterval = saveInt;
}

//----------------------------------------------
// Write output in a csv file
// Giulia Apr 2022

// save unode values in a Map
// key = name of attribute
// value = value of attribute
//------------------------------------------------

void Lattice::setMapForCSV(int i) {
  std::map<std::string, double> attribMap;
  // map of name (string) and attribute value (double)
  float smax, smin, sxx, syy, sxy;

  // rename some stuff coming from particle.cc so that it's shorter to use
  sxx = runParticle->sxx;
  syy = runParticle->syy;
  sxy = runParticle->sxy;

  // get stress components

  // assign a key (name, will be the header of the column) and its value to the
  // map
  attribMap["id"] = i;
  // save the particle number (id)
  // attribMap["x velocity"] = runParticle->v_darcyx*1000000000;
  // attribMap["y velocity"] = runParticle->v_darcyy*1000000000;
  // attribMap["Fracture normal stress"] =
  // runParticle->frac_stress_n*pressure_scale * pascal_scale;
  // attribMap["Fracture shear stress"] = runParticle->frac_shear*pressure_scale
  // * pascal_scale; attribMap["Gravity"] = runParticle->fg* pressure_scale *
  // pascal_scale; attribMap["Porosity"] = runParticle->average_porosity;

  double meanStress = (sxx + syy) / 2.0;
  // mean stress
  attribMap["Mean Stress"] = meanStress * pressure_scale * pascal_scale;

  attribMap["Youngs Modulus"] =
      runParticle->young * pressure_scale * pascal_scale;
  // attribMap["Healing"] = runParticle->healing;

  // Eigenvalues:
  smax = ((sxx + syy) / 2.0) -
         sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) + sxy * sxy);

  smin = ((sxx + syy) / 2.0) +
         sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) + sxy * sxy);

  double diffStress = sqrt((smax - smin) * (smax - smin));
  attribMap["Differential Stress"] = diffStress * pressure_scale * pascal_scale;
  attribMap["Sigma_1"] = smax * pressure_scale * pascal_scale;
  // eignevalue
  attribMap["Sigma_2"] = smin * pressure_scale * pascal_scale;
  // eignevalue

  double frac = 0.0;
  // set initially to zero, change if there is a fracture
  if (runParticle->draw_break) {
    if (runParticle->no_spring_flag)
      frac = -0.3;
    else if (runParticle->visc_flag)
      frac = -0.6;
    else
      frac = -1.0;

    runParticle->p_Unode->setAttribute(U_FRACTURES, frac);
    // something
  } else {
    runParticle->p_Unode->setAttribute(U_FRACTURES, frac);
  }
  attribMap["Fractures"] = (double)frac;
  // 0 = no fractures, -1 = there's at least a fracure
  // attribMap["Broken Bonds"] = runParticle->crack_seal_counter;
  // attribMap["Permeability"] = (double)runParticle->average_permeability;
  attribMap["Aperture"] = (double)runParticle->aperture;
  // fluid pressure (already in Pa)
  // attribMap["Tensile Strength"] = runParticle->average_break_str; // average
  // of break_Str of neighbours attribMap["Fracture Time"] =
  // runParticle->fracture_time;
  // test this one, I haven't checked if it works

  attribMap["x coord"] = runParticle->xpos;
  // x and y coordinates of the particle
  attribMap["y coord"] = runParticle->ypos;
  attribMap["z coord"] = 0.0;
  // paraview wants the 3rd dimension, so I set z to be always zero

  // This is the end of the map.

  // Now go to the function that actually writes the file.
  writeCSVfile(attribMap);
}

void Lattice::writeCSVfile(std::map<std::string, double> attribMap) {
  // writes the variables that were set in setMapForCSV to a csv file.

  vector<int> nameAttr;
  // initialise some variables
  vector<int> valueAttr;

  // create an iterator to go through the map elements: names and values
  map<string, double>::iterator itMap = attribMap.begin();

  if (doneWithHeaders == 0) // = haven't written headers yet
  {
    // Write HEADERS :
    while (itMap != attribMap.end()) {
      string aName = itMap->first;
      // get the attribute name
      myExpCsv << aName;
      // write the headers name by name
      if (itMap != --attribMap.end()) {
        myExpCsv << ",";
        // if it's not the last one, write a comma
      }
      itMap++;
    }
    doneWithHeaders = 1;
    // wrote the headers = true

    myExpCsv << endl;
  }

  while (itMap != attribMap.end())
  // write values
  {
    double aVal = itMap->second;
    // get the attribute value
    myExpCsv << aVal;
    // write the value.
    if (itMap != --attribMap.end())
    // Then, if it's not the last one,
    {
      myExpCsv << ",";
      // write a comma
    }
    itMap++;
    // increase the counter
  }
  myExpCsv << endl;
  // we're done with this particle, end the current line
}

void Lattice::UpdateElle(bool movenode) {
  // ----------------------------------------------------------
  // local variables
  // ----------------------------------------------------------

  double pressure; // variable for the sum of the trace of
  // the strain tensor sxx+syy
  Coords xy; // coordinate structure defined in
  // attrib.h
  int i, ii, max; // counter
  NodeAttrib *p;  // for nodes (x and y)
  float smax, smin, sxx, syy, sxy,
      energy; // some local variables for stresses etc.
  float uxy, uxx, uyy, lame1, lame2, pois;
  double frac;

  pois = 0.33333;

  runParticle = &refParticle; // jump to beginning of list

  max = ElleMaxNodes();

  for (ii = 0; ii < max; ii++) {
    nodes_Check[ii] = false;
  }

  // --------------------------------------------------------------------------------------
  // loop through the list of particles and change attributes and
  // positions in the Undoes
  // has to be called everytime something is changed and we talk to elle
  // again
  // this should include all kinds of plotting things, stresses,
  // energies, fractures etc.
  // --------------------------------------------------------------------------------------

  for (i = 0; i < numParticles; i++) // and loop again
  {

    // otherwise stuff gets too long

    sxx = runParticle->sxx;
    sxy = runParticle->sxy;
    syy = runParticle->syy;

    //----------------------------------------------------------------------
    // setting elastic constants
    //----------------------------------------------------------------------

    lame1 = runParticle->young * pois / (1.0 - (2.0 * pois) * (1.0 + pois));
    lame2 = runParticle->young / (2.0 + 2.0 * pois);

    //-----------------------------------------------------------------------
    // getting strain tensor
    //-----------------------------------------------------------------------

    uxy = ((1 + pois) / runParticle->young) * runParticle->sxy;
    uxx = (1 / runParticle->young) *
          (runParticle->sxx - (runParticle->syy * pois));
    uyy = (1 / runParticle->young) *
          (runParticle->syy - (runParticle->sxx * pois));

    //-----------------------------------------------------------------------
    // calculate elastic energy for particle
    //-----------------------------------------------------------------------

    runParticle->eEl = 0.5 * lame1 * (uxx + uyy) * (uxx + uyy);
    runParticle->eEl = runParticle->eEl +
                       (lame2 * ((uxx * uxx) + (uyy * uyy) + 2 * (uxy * uxy)));

    //---------------------------------------------------------
    // scaling energy
    //---------------------------------------------------------

    runParticle->eEl =
        runParticle->eEl * runParticle->mV * pressure_scale * pascal_scale;

    runParticle->p_Unode->setAttribute(U_ENERGY, runParticle->eEl);

    // -----------------------------------------------------------------------
    // stress tensor is set in the particle object themselves at
    // relaxation
    // therefore call particles, get stresses and calculate whatever
    // you want, then pass to Unodes.
    // here now pressure = trace of tensor sxx + syy
    // -----------------------------------------------------------------------

    pressure = (sxx + syy) / 2.0; // calc

    // ----------------------------------------------------------------------
    // this function setAttribute is defined in unodeP.h which
    // includes the
    // definition for the class plus all the functions. There is no
    // unodesP.cc in Elle
    // this function wants a reference int for which value we want to
    // change that is CONC_A at the moment (normally used for
    // concentration)
    // and also wants the value here as a double.
    //
    // this call just calls the function setAttribute in the Unode
    // connected to the particle
    // -----------------------------------------------------------------------
    runParticle->p_Unode->setAttribute(U_MEAN_STRESS, pressure);

    // runParticle->p_Unode->setAttribute(U_DENSITY,runParticle->fluid_particles);

    // runParticle->p_Unode->setAttribute(U_YOUNGSMODULUS,runParticle->rad_par_fluid);
    // runParticle->p_Unode->setAttribute(U_ATTRIB_B,runParticle->grain);

    runParticle->p_Unode->setAttribute(U_YOUNGSMODULUS, runParticle->young);

    // ------------------------------------------------------------------------
    // And pass some other values like the full stress tensor
    // in N_Attrib_A N_Attrib_B and N_Attrib_C
    // ------------------------------------------------------------------------

    // runParticle->p_Unode->setAttribute(CAXIS,runParticle->sxx);

    // runParticle->p_Unode->setAttribute(EULER_3,runParticle->syy);

    // runParticle->p_Unode->setAttribute(VISCOSITY,runParticle->sxy);

    // -------------------------------------------------------------------------
    // determine smax - smin (Eigenvalues)
    // -------------------------------------------------------------------------

    smax = ((sxx + syy) / 2.0) +
           sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) + sxy * sxy);

    smin = ((sxx + syy) / 2.0) -
           sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) + sxy * sxy);

    pressure = smax - smin;

    runParticle->p_Unode->setAttribute(U_DIF_STRESS, pressure);

    frac = 0.0;
    if (runParticle->draw_break) {
      if (runParticle->no_spring_flag)
        frac = -0.3;
      else if (runParticle->visc_flag)
        frac = -0.6;
      else
        frac = -1.0;

      runParticle->p_Unode->setAttribute(U_FRACTURES, frac); // something
      // to
      // visualize
      // fractures
    } else
      runParticle->p_Unode->setAttribute(U_FRACTURES, frac);

    pressure = 0.0;

    if (runParticle->isHoleBoundary)
      pressure = 1.0;

    if (runParticle->isHole)
      // pressure = runParticle->conc;
      pressure = -1;

    runParticle->p_Unode->setAttribute(U_PHASE, pressure);

    runParticle->p_Unode->setAttribute(U_TEMPERATURE, runParticle->temperature);

    if (runParticle->draw_break) {
      smin = smin * 10.0; // something to visualize fractures
    }

    pressure = smax;

    // runParticle->p_Unode->setAttribute(ENERGY,pressure);

    // pressure = arctan ((smax-sxx)/sxx);

    // runParticle->p_Unode->setAttribute(DISLOCATION_DENSITY,pressure);

    // -------------------------------------------------------------------------
    // and now also change the Unode position because the particles
    // have moved
    // I am not exactly shure what the consequences of that for Elle
    // are, maybe
    // something crashes at some point.
    // since the function setPosition in the Unode class wants a Coord
    // structure defined in attrib.h in Elle we create one and give it
    // the
    // x and y positions in our particle and then pass it to the Unode
    // function
    // -------------------------------------------------------------------------

    xy.x = runParticle->xpos; // pass x
    xy.y = runParticle->ypos; // pass y

    runParticle->p_Unode->setPosition(&xy); // and pass that to the
    // Unode

    // ------------------------------------------------------------------------------
    // This routine moves Elle nodes according to movements in the
    // particle lattice
    // At the moment a bit rough, node just moved on top of particle
    // depends a bit on resolution, make triangulation of unodes crash
    // !
    // ------------------------------------------------------------------------------

    if (movenode || !runParticle->isHole) {
      for (ii = 0; ii < 32; ii++) {
        if (runParticle->elle_Node[ii] != -1) {
          if (!nodes_Check[runParticle->elle_Node[ii]]) {
            ElleNodePosition(runParticle->elle_Node[ii], &xy);

            nodes_Check[runParticle->elle_Node[ii]] = true;

            xy.y = xy.y + (runParticle->ypos - runParticle->oldy);
            xy.x = xy.x + (runParticle->xpos - runParticle->oldx);
            ElleSetPosition(runParticle->elle_Node[ii], &xy);
          }
        }
      }
      runParticle->oldx = runParticle->xpos;
      runParticle->oldy = runParticle->ypos;
    }

    runParticle = runParticle->nextP; // and the next particle
  }

  // ---------------------------------------------------------------------------
  // call ElleUpdate() an Elle function that updates the interface
  // then the new values will be plotted
  // if security stop is set and max number of picts is reached dont
  // plot
  // pict anymore.
  // ---------------------------------------------------------------------------

  if (!set_max_pict) {
    // cout << "interface" << endl;
    ElleUpdate();
  } else if (num_pict < max_pict) {
    // cout << "interface" << endl;
    ElleUpdate();
    max_pict = max_pict + 1;
  }
}

//---------------------------------------------------------------------
// deformation routines, first simple shear
//---------------------------------------------------------------------

float Lattice::DeformLatticeSimpleShear(float move, int plot) {
  // ------------------------------------------------------
  // and some local variables
  // ------------------------------------------------------

  int i;          // counter for loop
  float wall_pos; // wall position

  cout << " Simple Shear " << endl;

  def_time_u = def_time_u + 1; // set deformation time

  runParticle = &refParticle; // start of list

  runParticle = runParticle->prevP; // go back one

  wall_pos = runParticle->ypos + runParticle->radius; // define wall position

  // ------------------------------------------------------------------------
  // run through particles and move their y position starting from upper
  // row
  // Upper row is pushed inwards. Remember that the yposition of upper
  // row
  // dont change in the relaxation
  // routine, they are fixed.
  // Other particles are moved as if all deformation happens in y
  // direction
  // This is a first approximation and makes the relaxation a lot faster
  // and
  // more accurate.
  // ------------------------------------------------------------------------

  for (i = 0; i < numParticles; i++) // go through all particles
                                     // starting from back
  {
    runParticle->xpos =
        runParticle->xpos + (move * (runParticle->ypos / wall_pos));

    runParticle->ChangeBox(repBox, particlex);

    runParticle = runParticle->prevP;
  }

  // upper_wall_pos = upper_wall_pos - move;

  // for wrapping: the slope of the boundary on the sides
  // lattice_wall_slope is preset to 1, thus the "if" is necessary
  // lattice_wall_slope += move;

  Relaxation(); // and do a full relaxation

  if (plot == 1)
    UpdateElle(); // and update the interface of Elle

  SetRightLatticeWallPos();

  // DumpStatisticStressBox ( 0.1, 0.9, 0.1, 0.9, 0.001 );
}

//-------------------------------------------------------------------------------
// Apply gravity for large scale and hydrofracture dynamics
//-------------------------------------------------------------------------------

void Lattice::ApplyGravity(float height) {
  int i;
  float average_density;
  float average_young;
  float move;

  average_density = average_young = 0.0;

  for (i = 0; i < numParticles; i++) {
    average_density = average_density + runParticle->real_density;
    average_young = average_young + runParticle->real_young;

    runParticle = runParticle->prevP;
  }
  average_density = average_density / i;
  average_young = average_young / i;

  move = 1.5 * 10 * height * average_density /
         average_young; // dl = (stress*L)/E = (F/A)*(L/E) = (m*g/A)*(L/E) =
                        // (rho*V*g/A)*(L/E)

  cout << average_density << endl;
  cout << average_young << endl;
  cout << move << endl;

  for (i = 0; i < numParticles; i++) {
    runParticle->ypos = runParticle->ypos - (runParticle->ypos * move);

    runParticle->ChangeBox(repBox, particlex);

    runParticle = runParticle->prevP;
  }
  runParticle = &refParticle;
  surf_pos = runParticle->prevP->ypos;
  volume = surf_pos * (1.0 - (1 / (2 * particlex)));
  cout << volume << endl;
}

// simple extension in x

float Lattice::DeformLatticeExtend(float move, int plot) {
  int i;
  float wall_pos, wall_Xpos;
  Coords newbbox[4];

  // def_time _u = def_time_u + 1;

  runParticle = &refParticle;
  runParticle = runParticle->prevP;

  wall_pos = wall_Xpos = runParticle->xpos;

  for (i = 0; i < particlex; i++) {
    runParticle = runParticle->prevP;
  }

  if (runParticle->xpos > wall_Xpos)
    wall_Xpos = runParticle->xpos;

  wall_Xpos = wall_Xpos + (0.5 / particlex);

  for (i = 0; i < numParticles; i++) {
    runParticle->xpos =
        runParticle->xpos + (runParticle->xpos * move / wall_Xpos);
    runParticle->ChangeBox(repBox, particlex);

    runParticle = runParticle->prevP;
  }
  if (walls) {
    right_wall_pos += move;
  }
  Relaxation();

  newbbox[0].x = 0.0;
  newbbox[0].y = 0.0;

  newbbox[1].x = wall_Xpos + move;
  newbbox[1].y = 0.0;

  newbbox[2].x = wall_Xpos + move;
  newbbox[2].y = 1.0;

  newbbox[3].x = 0.0;
  newbbox[3].y = 1.0;

  // and set the new boxsize

  ElleSetCellBBox(&newbbox[0], &newbbox[1], &newbbox[2], &newbbox[3]);
  if (plot == 1)
    UpdateElle();
}

// more complicated extension that can also shear

void Lattice::DeformLatticeExtension(float move, int plot, float dx,
                                     float shear) {
  // ------------------------------------------------------
  // and some local variables
  // ------------------------------------------------------

  int i;                                  // counter for loop
  float wall_pos, dy, shearfactor;        // wall position
  float movex, position_Xwall, wall_Xpos; // stuff for x wall

  cout << " Pure Extension in x " << endl;

  dy = 1.0 - dx; // dx is component of dx direction, if 1.0 all extension
                 // horizontal, if 0.0 all extension vertical

  Coords newbbox[4]; // have to change Elle box size in file

  def_time_p = def_time_p + 1; // set deformation time

  runParticle = &refParticle; // start of list

  runParticle = runParticle->prevP; // go back one

  wall_pos = runParticle->ypos; // get wall position y

  wall_Xpos = runParticle->xpos; // get wall position x

  // have to check if particle on edge is furthest in x or particle
  // below it
  // therefore we run one x row below the highest particle and check
  // that p as well

  for (i = 0; i < particlex; i++) {
    runParticle = runParticle->prevP;
  }

  if (runParticle->xpos > wall_Xpos) {
    wall_Xpos = runParticle->xpos;
  }
  // now we can calculate the initial area of the box at time 1

  /* if (def_time_p == 1)
     {

       area_box = wall_pos * wall_Xpos;
       cout << area_box << endl;

     }*/
  // ------------------------------------------------------------------------
  // run through particles and move their y position starting from upper
  // row
  // Upper row is pushed inwards. Remember that the yposition of upper
  // row
  // dont change in the relaxation
  // routine, they are fixed.
  // Other particles are moved as if all deformation happens in y
  // direction
  // This is a first approximation and makes the relaxation a lot faster
  // and
  // more accurate.
  // ------------------------------------------------------------------------

  position_Xwall = wall_Xpos + (move * dx);

  for (i = 0; i < numParticles; i++) // go through all particles
                                     // starting from back
  {
    // average the deformation in x..
    shearfactor = shear * (runParticle->ypos - 1);

    runParticle->xpos =
        runParticle->xpos + (runParticle->xpos * (move * dx) / wall_Xpos) +
        (runParticle->xpos * (move * dx) / wall_Xpos) * shearfactor;
    // runParticle->xpos = runParticle->xpos +
    // runParticle->xpos*(-(runParticle->ypos-1)*shear);
    runParticle->ypos =
        runParticle->ypos + (runParticle->ypos * (move * dy) / wall_pos);
    if (sheet) {
      runParticle->sheet_x =
          runParticle->sheet_x +
          (runParticle->sheet_x * (move * dx) / wall_Xpos) +
          (runParticle->sheet_x * (move * dx) / wall_Xpos) * shearfactor;
      // runParticle->sheet_x = runParticle->sheet_x +
      // runParticle->sheet_x*(-(runParticle->sheet_y-1)*shear);
      runParticle->sheet_y = runParticle->sheet_y +
                             (runParticle->sheet_y * (move * dy) / wall_pos);
    }

    runParticle->ChangeBox(repBox, particlex);

    runParticle = runParticle->prevP;
  }

  if (walls) {
    right_wall_pos = right_wall_pos + (move * dx);
    upper_wall_pos = upper_wall_pos + (move * dy);
  }

  // now define the new boxsize for elle

  newbbox[0].x = 0.0;
  newbbox[0].y = 0.0;

  newbbox[1].x = position_Xwall;
  newbbox[1].y = 0.0;

  newbbox[2].x = position_Xwall;
  newbbox[2].y = wall_pos + (move * dy);

  newbbox[3].x = 0.0;
  newbbox[3].y = wall_pos + (move * dy);

  // and set the new boxsize

  ElleSetCellBBox(&newbbox[0], &newbbox[1], &newbbox[2], &newbbox[3]);

  /*for (i = 0; i < numParticles; i++)
    {
      runParticle->ChangeBox (repBox, particlex);

      runParticle = runParticle->nextP;
    }*/

  // Relaxation ();		// and do a full relaxation

  if (plot == 1)
    UpdateElle(); // and update the interface of Elle if it
  // is wanted
}
// DEFORM_LATTICE
/*******************************************************************************
 * This function does whatever deformation we want to perform with the lattice
 * It is called from the outside by the User and does all that is needed for
 * a timestep including all the relaxation etc.
 * at the moment it receives a float for how much the walls are moved
 * as a deformation it moves upper and lower row of particles inwards if
 * positive value
 *
 * Here we can include all kinds of stuff in the future
 *
 * daniel spring and summer 2002
 ********************************************************************************/

// ---------------------------------------------------------------------------------
// function DeformLattice(float move) in Lattice class
// receives a number for deformation
// dont know why I return something ?
//
// called from the outside ! i.e. mike_elle whatever routine
// ---------------------------------------------------------------------------------

float Lattice::DeformLattice(float move, int plot) {
  // ------------------------------------------------------
  // and some local variables
  // ------------------------------------------------------

  int i;          // counter for loop
  float wall_pos; // wall position

  cout << " Uniaxial Compression" << endl;

  def_time_u = def_time_u + 1; // set deformation time

  runParticle = &refParticle; // start of list

  runParticle = runParticle->prevP; // go back one

  wall_pos = runParticle->ypos + runParticle->radius; // define wall position

  // ------------------------------------------------------------------------
  // run through particles and move their y position starting from upper
  // row
  // Upper row is pushed inwards. Remember that the yposition of upper
  // row
  // dont change in the relaxation
  // routine, they are fixed.
  // Other particles are moved as if all deformation happens in y
  // direction
  // This is a first approximation and makes the relaxation a lot faster
  // and
  // more accurate.
  // ------------------------------------------------------------------------

  for (i = 0; i < numParticles; i++) // go through all particles
                                     // starting from back
  {
    runParticle->ypos =
        runParticle->ypos - (runParticle->ypos * move /
                             wall_pos); // strain = dl/l_u = (l_d - l_u)/l_u

    runParticle->ChangeBox(repBox, particlex);

    runParticle = runParticle->prevP;
  }

  upper_wall_pos = upper_wall_pos - move;

  Relaxation(); // and do a full relaxation

  if (plot == 1)
    UpdateElle(); // and update the interface of Elle

  SetRightLatticeWallPos();
}

float Lattice::DeformLatticeX(float move, int plot) {
  // ------------------------------------------------------
  // and some local variables
  // ------------------------------------------------------

  int i;                         // counter for loop
  float wall_pos, wall_pos_left; // wall position
  Coords newbbox[4];

  cout << " Uniaxial Compression in X" << endl;

  def_time_u = def_time_u + 1; // set deformation time

  //  runParticle = &refParticle;	// start of list

  //	wall_pos_left = runParticle->xpos;
  //
  //	for (i = 0; i < particlex; i++)
  //    {
  //      runParticle = runParticle->nextP;
  //    }

  //  if (runParticle->xpos < wall_pos_left)
  //    {
  //      wall_pos_left = runParticle->xpos;
  //    }

  runParticle = &refParticle; // start of list

  runParticle = runParticle->prevP; // go back one

  wall_pos = runParticle->xpos; // define wall position

  for (i = 0; i < particlex; i++) {
    runParticle = runParticle->prevP;
  }

  if (runParticle->xpos > wall_pos) {
    wall_pos = runParticle->xpos;
  }

  // ------------------------------------------------------------------------
  // run through particles and move their y position starting from upper
  // row
  // Upper row is pushed inwards. Remember that the yposition of upper
  // row
  // dont change in the relaxation
  // routine, they are fixed.
  // Other particles are moved as if all deformation happens in y
  // direction
  // This is a first approximation and makes the relaxation a lot faster
  // and
  // more accurate.
  // ------------------------------------------------------------------------

  for (i = 0; i < numParticles; i++) // go through all particles
                                     // starting from back
  {
    runParticle->xpos =
        runParticle->xpos - (runParticle->xpos * move / wall_pos);

    runParticle->ypos =
        runParticle->ypos -
        (runParticle->ypos * (surf_pos - (volume / (wall_pos + (move * -1)))) /
         surf_pos);

    runParticle->ChangeBox(repBox, particlex);

    runParticle = runParticle->prevP;
  }

  right_wall_pos = right_wall_pos - move;
  right_lattice_wall_pos -= move;
  surf_pos = surf_pos - (surf_pos - (volume / (wall_pos + (move * -1))));
  cout << surf_pos << endl;
  //  right_lattice_wall_pos = wall_pos + refParticle.radius;
  //  left_lattice_wall_pos = wall_pos_left - refParticle.radius;

  Relaxation(); // and do a full relaxation

  if (plot == 1)
    UpdateElle(); // and update the interface of Elle

  SetRightLatticeWallPos();

  newbbox[0].x = 0.0;
  newbbox[0].y = 0.0;

  newbbox[1].x = wall_pos;
  newbbox[1].y = 0.0;

  newbbox[2].x = wall_pos;
  newbbox[2].y = 1.0;

  newbbox[3].x = 0.0;
  newbbox[3].y = 1.0;

  // and set the new boxsize

  ElleSetCellBBox(&newbbox[0], &newbbox[1], &newbbox[2], &newbbox[3]);
}

//~ float
//~ Lattice::DeformLatticeXMidTwo(float move, int plot)
//~ {
//~ // ------------------------------------------------------
//~ // and some local variables
//~ // ------------------------------------------------------
//~
//~ int i;			// counter for loop
//~ float wall_pos, wall_pos_left;		// wall position
//~ Coords newbbox[4];
//~
//~ cout << " Uniaxial Compression in X" << endl;
//~
//~ def_time_u = def_time_u + 1;	// set deformation time
//~
//~ //  runParticle = &refParticle;	// start of list
//~
//~ //	wall_pos_left = runParticle->xpos;
//~ //
//~ //	for (i = 0; i < particlex; i++)
//~ //    {
//~ //      runParticle = runParticle->nextP;
//~ //    }
//~
//~ //  if (runParticle->xpos < wall_pos_left)
//~ //    {
//~ //      wall_pos_left = runParticle->xpos;
//~ //    }
//~
//~ runParticle = &refParticle;	// start of list
//~
//~ runParticle = runParticle->prevP;	// go back one
//~
//~ wall_pos = runParticle->xpos;	// define wall position
//~
//~ for (i = 0; i < particlex; i++)
//~ {
//~ runParticle = runParticle->prevP;
//~ }
//~
//~ if (runParticle->xpos > wall_pos)
//~ {
//~ wall_pos = runParticle->xpos;
//~ }
//~
//~
//~ // ------------------------------------------------------------------------
//~ // run through particles and move their y position starting from upper
//~ // row
//~ // Upper row is pushed inwards. Remember that the yposition of upper
//~ // row
//~ // dont change in the relaxation
//~ // routine, they are fixed.
//~ // Other particles are moved as if all deformation happens in y
//~ // direction
//~ // This is a first approximation and makes the relaxation a lot faster
//~ // and
//~ // more accurate.
//~ // ------------------------------------------------------------------------
//~
//~ for (i = 0; i < numParticles; i++)	// go through all particles
//~ // starting from back
//~ {
//~ if (runParticle->xpos > 0.25 && runParticle->xpos < 0.4 && runParticle->ypos
//> 0.7) ~ {
//~ runParticle->xpos =
//~ runParticle->xpos -
//~ ((runParticle->xpos-0.25) * move / (wall_pos-1.2));
//~ }
//~ if (runParticle->xpos >= 0.4 && runParticle->xpos < 0.65 &&
//runParticle->ypos > 0.7) ~ { ~ runParticle->xpos = runParticle->xpos -
//move/2.0;
//~ }
//~ if (runParticle->xpos > 0.63 && runParticle->xpos < 0.65 &&
//runParticle->ypos > 0.7) ~ {
//~ runParticle->xpos =
//~ runParticle->xpos -
//~ ((runParticle->xpos-0.63) * move / (wall_pos-1.2));
//~ }
//~ if (runParticle->xpos > 0.65 && runParticle->ypos > 0.7)
//~ runParticle->xpos = runParticle->xpos - move;
//~ if (runParticle->ypos <= 0.7)
//~ runParticle->xpos =
//~ runParticle->xpos -
//~ (runParticle->xpos * move / wall_pos);
//~
//~ runParticle->ypos = runParticle->ypos - (runParticle->ypos *
//(surf_pos-(volume/(wall_pos+(move*-1))))/surf_pos);
//~
//~ runParticle->ChangeBox (repBox, particlex);
//~
//~ runParticle = runParticle->prevP;
//~ }
//~
//~ right_wall_pos = right_wall_pos - move;
//~ right_lattice_wall_pos -= move;
//~ surf_pos = surf_pos - (surf_pos-(volume/(wall_pos+(move*-1))));
//~ cout << surf_pos << endl;
//~ //  right_lattice_wall_pos = wall_pos + refParticle.radius;
//~ //  left_lattice_wall_pos = wall_pos_left - refParticle.radius;
//~
//~ Relaxation ();		// and do a full relaxation
//~
//~ if (plot == 1)
//~ UpdateElle ();	// and update the interface of Elle
//~
//~ SetRightLatticeWallPos();
//~
//~ newbbox[0].x = 0.0;
//~ newbbox[0].y = 0.0;
//~
//~ newbbox[1].x = wall_pos;
//~ newbbox[1].y = 0.0;
//~
//~ newbbox[2].x = wall_pos;
//~ newbbox[2].y = 1.0;
//~
//~ newbbox[3].x = 0.0;
//~ newbbox[3].y = 1.0;
//~
//~ // and set the new boxsize
//~
//~ ElleSetCellBBox (&newbbox[0], &newbbox[1], &newbbox[2], &newbbox[3]);
//~ }

float Lattice::DeformLatticeXMidTwo_b(float move, float ymov, int plot) {
  // ------------------------------------------------------
  // and some local variables
  // ------------------------------------------------------

  int i;                                   // counter for loop
  float wall_pos, wall_pos_left, ymovsurf; // wall position
  Coords newbbox[4];

  cout << " Uniaxial Compression in X" << endl;

  ymovsurf = ymov * surf_pos;

  def_time_u = def_time_u + 1; // set deformation time

  runParticle = &refParticle; // start of list

  runParticle = runParticle->prevP; // go back one

  wall_pos = runParticle->xpos; // define wall position

  for (i = 0; i < particlex; i++) {
    runParticle = runParticle->prevP;
  }

  if (runParticle->xpos > wall_pos) {
    wall_pos = runParticle->xpos;
  }

  // ------------------------------------------------------------------------
  // run through particles and move their y position starting from upper
  // row
  // Upper row is pushed inwards. Remember that the yposition of upper
  // row
  // dont change in the relaxation
  // routine, they are fixed.
  // Other particles are moved as if all deformation happens in y
  // direction
  // This is a first approximation and makes the relaxation a lot faster
  // and
  // more accurate.
  // ------------------------------------------------------------------------

  for (i = 0; i < numParticles; i++) // go through all particles
                                     // starting from back
  {
    if (runParticle->xpos > 0.25 && runParticle->xpos < 0.4 &&
        runParticle->ypos > ymovsurf) {
      runParticle->xpos =
          runParticle->xpos -
          ((runParticle->xpos - 0.25) * (move / wall_pos) * (0.5 / 0.15));
    }
    if (runParticle->xpos >= 0.4 && runParticle->xpos < 0.75 &&
        runParticle->ypos > ymovsurf) {
      runParticle->xpos = runParticle->xpos - move / 2.0;
    }
    if (runParticle->xpos > 0.6 && runParticle->xpos < 0.75 &&
        runParticle->ypos > ymovsurf) {
      runParticle->xpos =
          runParticle->xpos -
          ((runParticle->xpos - 0.6) * (move / wall_pos) * (0.5 / 0.15));
    }
    if (runParticle->xpos > 0.75 && runParticle->ypos > ymovsurf)
      runParticle->xpos = runParticle->xpos - move;
    if (runParticle->ypos <= ymovsurf)
      runParticle->xpos =
          runParticle->xpos - (runParticle->xpos * move / wall_pos);

    // now average for y movement due to stretching

    if (runParticle->ypos <= ymovsurf)

      runParticle->ypos =
          runParticle->ypos -
          (runParticle->ypos *
           (surf_pos - (volume / (wall_pos + (move * -1)))) / surf_pos);
    else {
      if (runParticle->xpos > 0.25 && runParticle->xpos < 0.4)
        runParticle->ypos =
            runParticle->ypos -
            (runParticle->ypos *
             (surf_pos - (volume / (wall_pos + (move * -1)))) / surf_pos);
      else if (runParticle->xpos > 0.6 && runParticle->xpos < 0.75)
        runParticle->ypos =
            runParticle->ypos -
            (runParticle->ypos *
             (surf_pos - (volume / (wall_pos + (move * -1)))) / surf_pos);
      else
        runParticle->ypos =
            runParticle->ypos -
            (ymovsurf * (surf_pos - (volume / (wall_pos + (move * -1)))) /
             surf_pos);
    }

    runParticle->ChangeBox(repBox, particlex);

    runParticle = runParticle->prevP;
  }

  right_wall_pos = right_wall_pos - move;
  right_lattice_wall_pos -= move;
  surf_pos = surf_pos - (surf_pos - (volume / (wall_pos + (move * -1))));
  cout << surf_pos << endl;
  //  right_lattice_wall_pos = wall_pos + refParticle.radius;
  //  left_lattice_wall_pos = wall_pos_left - refParticle.radius;

  Relaxation(); // and do a full relaxation

  if (plot == 1)
    UpdateElle(); // and update the interface of Elle

  SetRightLatticeWallPos();

  newbbox[0].x = 0.0;
  newbbox[0].y = 0.0;

  newbbox[1].x = wall_pos;
  newbbox[1].y = 0.0;

  newbbox[2].x = wall_pos;
  newbbox[2].y = 1.0;

  newbbox[3].x = 0.0;
  newbbox[3].y = 1.0;

  // and set the new boxsize

  ElleSetCellBBox(&newbbox[0], &newbbox[1], &newbbox[2], &newbbox[3]);
}

float Lattice::DeformLatticeXMidTwo(float move, float ymov, int plot,
                                    float left_margin_width,
                                    float right_margin_width,
                                    float left_rift_width,
                                    float right_rift_width) {

  int i;
  float wall_pos, wall_pos_left, ymovsurf;
  Coords newbbox[4];

  ymovsurf = ymov * surf_pos;

  def_time_u = def_time_u + 1;

  runParticle = &refParticle;
  runParticle = runParticle->prevP;

  wall_pos = runParticle->xpos;

  for (i = 0; i < particlex; i++) {
    runParticle = runParticle->prevP;
  }

  if (runParticle->xpos > wall_pos) {
    wall_pos = runParticle->xpos;
  }

  left_margin_width *= wall_pos;
  right_margin_width *= wall_pos;
  left_rift_width *= wall_pos;
  right_rift_width *= wall_pos;

  for (i = 0; i < numParticles; i++) {
    if (runParticle->xpos > left_margin_width &&
        runParticle->xpos < left_margin_width + left_rift_width &&
        runParticle->ypos > ymovsurf) {
      runParticle->xpos =
          runParticle->xpos -
          ((runParticle->xpos - left_margin_width) * move /
           (wall_pos - (1 - right_margin_width + left_rift_width)));
    }

    if (runParticle->xpos >= left_margin_width + left_rift_width &&
        runParticle->xpos < wall_pos - right_margin_width &&
        runParticle->ypos > ymovsurf) {
      runParticle->xpos = runParticle->xpos - move / 2.0;
    }

    if (runParticle->xpos > wall_pos - right_margin_width - right_rift_width &&
        runParticle->xpos < wall_pos - right_margin_width &&
        runParticle->ypos > ymovsurf) {
      runParticle->xpos =
          runParticle->xpos -
          ((runParticle->xpos - right_margin_width - right_rift_width) * move /
           (wall_pos - (1 - right_margin_width + left_rift_width)));
    }

    if (runParticle->xpos > wall_pos - right_margin_width &&
        runParticle->ypos > ymovsurf)
      runParticle->xpos = runParticle->xpos - move;

    if (runParticle->ypos <= ymovsurf)
      runParticle->xpos =
          runParticle->xpos - (runParticle->xpos * move / wall_pos);

    if (runParticle->ypos <= ymovsurf)
      runParticle->ypos =
          runParticle->ypos -
          (runParticle->ypos *
           (surf_pos - (volume / (wall_pos + (move * -1)))) / surf_pos);
    else {
      if (runParticle->xpos > left_margin_width &&
          runParticle->xpos < left_margin_width + left_rift_width)
        runParticle->ypos =
            runParticle->ypos -
            (runParticle->ypos *
             (surf_pos - (volume / (wall_pos + (move * -1)))) / surf_pos);
      else if (runParticle->xpos >
                   wall_pos - right_margin_width - right_rift_width &&
               runParticle->xpos < wall_pos - right_margin_width)
        runParticle->ypos =
            runParticle->ypos -
            (runParticle->ypos *
             (surf_pos - (volume / (wall_pos + (move * -1)))) / surf_pos);
      else
        runParticle->ypos =
            runParticle->ypos -
            (ymovsurf * (surf_pos - (volume / (wall_pos + (move * -1)))) /
             surf_pos);
    }

    runParticle->ChangeBox(repBox, particlex);

    runParticle = runParticle->prevP;
  }

  right_wall_pos = right_wall_pos - move;
  right_lattice_wall_pos -= move;
  surf_pos = surf_pos - (surf_pos - (volume / (wall_pos + (move * -1))));
  cout << surf_pos << endl;

  Relaxation();

  if (plot == 1)
    UpdateElle();

  SetRightLatticeWallPos();

  newbbox[0].x = 0.0;
  newbbox[0].y = 0.0;

  newbbox[1].x = wall_pos;
  newbbox[1].y = 0.0;

  newbbox[2].x = wall_pos;
  newbbox[2].y = 1.0;

  newbbox[3].x = 0.0;
  newbbox[3].y = 1.0;

  // and set the new boxsize

  ElleSetCellBBox(&newbbox[0], &newbbox[1], &newbbox[2], &newbbox[3]);
}

float Lattice::DeformLatticeXMid(float move, int plot, float neck_height,
                                 float left, float right) {
  // ------------------------------------------------------
  // and some local variables
  // ------------------------------------------------------

  int i;                         // counter for loop
  float wall_pos, wall_pos_left; // wall position
  Coords newbbox[4];

  cout << " Uniaxial Compression in X" << endl;

  def_time_u = def_time_u + 1; // set deformation time

  //  runParticle = &refParticle;	// start of list

  //	wall_pos_left = runParticle->xpos;
  //
  //	for (i = 0; i < particlex; i++)
  //    {
  //      runParticle = runParticle->nextP;
  //    }

  //  if (runParticle->xpos < wall_pos_left)
  //    {
  //      wall_pos_left = runParticle->xpos;
  //    }

  runParticle = &refParticle; // start of list

  runParticle = runParticle->prevP; // go back one

  wall_pos = runParticle->xpos; // define wall position

  for (i = 0; i < particlex; i++) {
    runParticle = runParticle->prevP;
  }

  if (runParticle->xpos > wall_pos) {
    wall_pos = runParticle->xpos;
  }

  // ------------------------------------------------------------------------
  // run through particles and move their y position starting from upper
  // row
  // Upper row is pushed inwards. Remember that the yposition of upper
  // row
  // dont change in the relaxation
  // routine, they are fixed.
  // Other particles are moved as if all deformation happens in y
  // direction
  // This is a first approximation and makes the relaxation a lot faster
  // and
  // more accurate.
  // ------------------------------------------------------------------------

  for (i = 0; i < numParticles; i++) // go through all particles
                                     // starting from back
  {
    if (runParticle->xpos > left && runParticle->xpos < right &&
        runParticle->ypos > neck_height) {
      runParticle->xpos =
          runParticle->xpos -
          ((runParticle->xpos - left) * move / (wall_pos - (1 - right + left)));
    }

    if (runParticle->xpos > right && runParticle->ypos > neck_height)
      runParticle->xpos = runParticle->xpos - move;
    if (runParticle->ypos <= neck_height)
      runParticle->xpos =
          runParticle->xpos - (runParticle->xpos * move / wall_pos);

    // runParticle->ypos = runParticle->ypos - (runParticle->ypos *
    // (surf_pos-(volume/(wall_pos+(move*-1))))/surf_pos);

    runParticle->ChangeBox(repBox, particlex);

    runParticle->oldy = runParticle->ypos;

    runParticle = runParticle->prevP;
  }

  right_wall_pos = right_wall_pos - move;
  right_lattice_wall_pos -= move;
  surf_pos = surf_pos - (surf_pos - (volume / (wall_pos + (move * -1))));
  cout << surf_pos << endl;
  //  right_lattice_wall_pos = wall_pos + refParticle.radius;
  //  left_lattice_wall_pos = wall_pos_left - refParticle.radius;

  // Relaxation ();		// and do a full relaxation

  // if (plot == 1)
  //  UpdateElle ();	// and update the interface of Elle

  // SetRightLatticeWallPos();

  newbbox[0].x = 0.0;
  newbbox[0].y = 0.0;

  newbbox[1].x = wall_pos;
  newbbox[1].y = 0.0;

  newbbox[2].x = wall_pos;
  newbbox[2].y = 1.0;

  newbbox[3].x = 0.0;
  newbbox[3].y = 1.0;

  // and set the new boxsize

  ElleSetCellBBox(&newbbox[0], &newbbox[1], &newbbox[2], &newbbox[3]);
}

float Lattice::DeformLatticedriftmid(float move, int plot, float change,
                                     float leftup, float rightup,
                                     float leftdown, float rightdown) {
  // ------------------------------------------------------
  // and some local variables
  // ------------------------------------------------------

  int i;                         // counter for loop
  float wall_pos, wall_pos_left; // wall position
  Coords newbbox[4];

  cout << "extension in the middle" << endl;

  def_time_u = def_time_u + 1; // set deformation time

  runParticle = &refParticle; // start of list

  runParticle = runParticle->prevP; // go back one

  wall_pos = runParticle->xpos; // define wall position

  // ------------------------------------------------------------------------
  // run through particles and move their y position starting from upper
  // row
  // Upper row is pushed inwards. Remember that the yposition of upper
  // row
  // dont change in the relaxation
  // routine, they are fixed.
  // Other particles are moved as if all deformation happens in y
  // direction
  // This is a first approximation and makes the relaxation a lot faster
  // and
  // more accurate.
  // ------------------------------------------------------------------------

  for (i = 0; i < numParticles; i++) // go through all particles
  // starting from back
  {
    // if (runParticle->xpos > leftup && runParticle->xpos < (wall_pos -
    // rightup) && runParticle->ypos > change)
    //{
    //	runParticle->xpos = runParticle->xpos - ((runParticle->xpos-leftup) *
    //move / (wall_pos-(rightup+leftup)));
    // }

    // if (runParticle->xpos > (wall_pos - rightup) && runParticle->ypos >
    // change) 	runParticle->xpos = runParticle->xpos - move;

    // if (runParticle->xpos > leftdown && runParticle->xpos <= (wall_pos -
    // rightdown) && runParticle->ypos < change)
    //{
    //	runParticle->xpos = runParticle->xpos - ((runParticle->xpos-leftdown) *
    //move / (wall_pos-(rightdown+leftdown)));
    // }

    // if (runParticle->xpos > (wall_pos - rightdown) && runParticle->ypos <=
    // change) 	runParticle->xpos = runParticle->xpos - move;

    if (runParticle->xpos > leftup && runParticle->xpos < rightup &&
        runParticle->ypos > change) {
      runParticle->xpos =
          runParticle->xpos - ((runParticle->xpos - leftup) * move /
                               (wall_pos - (1 - rightup + leftup)));
    }

    if (runParticle->xpos > rightup && runParticle->ypos > change)
      runParticle->xpos = runParticle->xpos - move;

    if (runParticle->xpos > leftdown && runParticle->xpos <= rightdown &&
        runParticle->ypos < change) {
      runParticle->xpos =
          runParticle->xpos - ((runParticle->xpos - leftdown) * move /
                               (wall_pos - (1 - rightdown + leftdown)));
    }

    if (runParticle->xpos > rightdown && runParticle->ypos <= change)
      runParticle->xpos = runParticle->xpos - move;

    // if (runParticle->ypos <= neck_height)
    //	runParticle->xpos =
    // runParticle->xpos -
    //(runParticle->xpos * move / wall_pos);

    // runParticle->ypos = runParticle->ypos - (runParticle->ypos *
    // (surf_pos-(volume/(wall_pos+(move*-1))))/surf_pos);

    runParticle->ChangeBox(repBox, particlex);

    runParticle->oldy = runParticle->ypos;

    runParticle = runParticle->prevP;
  }

  right_wall_pos = right_wall_pos - move;
  right_lattice_wall_pos -= move;
  // surf_pos = surf_pos - (surf_pos-(volume/(wall_pos+(move*-1))));
  // cout << surf_pos << endl;
  //  right_lattice_wall_pos = wall_pos + refParticle.radius;
  //  left_lattice_wall_pos = wall_pos_left - refParticle.radius;

  // Relaxation ();		// and do a full relaxation

  // if (plot == 1)
  //  UpdateElle ();	// and update the interface of Elle

  // SetRightLatticeWallPos();

  newbbox[0].x = 0.0;
  newbbox[0].y = 0.0;

  newbbox[1].x = wall_pos;
  newbbox[1].y = 0.0;

  newbbox[2].x = wall_pos;
  newbbox[2].y = 1.0;

  newbbox[3].x = 0.0;
  newbbox[3].y = 1.0;

  // and set the new boxsize

  ElleSetCellBBox(&newbbox[0], &newbbox[1], &newbbox[2], &newbbox[3]);
}

float Lattice::DeformLatticedriftmidnew(float move, int plot, float change,
                                        float leftup, float rightup,
                                        float leftdown, float rightdown) {
  // ------------------------------------------------------
  // and some local variables
  // ------------------------------------------------------

  int i;                         // counter for loop
  float wall_pos, wall_pos_left; // wall position
  Coords newbbox[4];

  cout << "extension in the middle" << endl;

  def_time_u = def_time_u + 1; // set deformation time

  runParticle = &refParticle; // start of list

  runParticle = runParticle->prevP; // go back one

  wall_pos = runParticle->xpos; // define wall position

  // ------------------------------------------------------------------------
  // run through particles and move their y position starting from upper
  // row
  // Upper row is pushed inwards. Remember that the yposition of upper
  // row
  // dont change in the relaxation
  // routine, they are fixed.
  // Other particles are moved as if all deformation happens in y
  // direction
  // This is a first approximation and makes the relaxation a lot faster
  // and
  // more accurate.
  // ------------------------------------------------------------------------

  for (i = 0; i < numParticles; i++) // go through all particles
  // starting from back
  {
    if (runParticle->xpos > leftup &&
        runParticle->xpos < (wall_pos - rightup) &&
        runParticle->ypos > change) {
      runParticle->xpos =
          runParticle->xpos - ((runParticle->xpos - leftup) * move * wall_pos /
                               (wall_pos - (rightup + leftup)));
    }

    if (runParticle->xpos > (wall_pos - rightup) && runParticle->ypos > change)
      runParticle->xpos = runParticle->xpos - move;

    if (runParticle->xpos > leftdown &&
        runParticle->xpos <= (wall_pos - rightdown) &&
        runParticle->ypos < change) {
      runParticle->xpos =
          runParticle->xpos - ((runParticle->xpos - leftdown) * move *
                               wall_pos / (wall_pos - (rightdown + leftdown)));
    }

    if (runParticle->xpos > (wall_pos - rightdown) &&
        runParticle->ypos <= change)
      runParticle->xpos = runParticle->xpos - move;

    // if (runParticle->ypos <= neck_height)
    //	runParticle->xpos =
    // runParticle->xpos -
    //(runParticle->xpos * move / wall_pos);

    // runParticle->ypos = runParticle->ypos - (runParticle->ypos *
    // (surf_pos-(volume/(wall_pos+(move*-1))))/surf_pos);

    runParticle->ChangeBox(repBox, particlex);

    runParticle->oldy = runParticle->ypos;

    runParticle = runParticle->prevP;
  }

  right_wall_pos = right_wall_pos - move;
  right_lattice_wall_pos -= move;
  // surf_pos = surf_pos - (surf_pos-(volume/(wall_pos+(move*-1))));
  // cout << surf_pos << endl;
  //  right_lattice_wall_pos = wall_pos + refParticle.radius;
  //  left_lattice_wall_pos = wall_pos_left - refParticle.radius;

  // Relaxation ();		// and do a full relaxation

  // if (plot == 1)
  //  UpdateElle ();	// and update the interface of Elle

  // SetRightLatticeWallPos();

  newbbox[0].x = 0.0;
  newbbox[0].y = 0.0;

  newbbox[1].x = wall_pos;
  newbbox[1].y = 0.0;

  newbbox[2].x = wall_pos;
  newbbox[2].y = 1.0;

  newbbox[3].x = 0.0;
  newbbox[3].y = 1.0;

  // and set the new boxsize

  ElleSetCellBBox(&newbbox[0], &newbbox[1], &newbbox[2], &newbbox[3]);
}

float Lattice::DeformLatticedriftmidnew2(float move, int plot, float change,
                                         float leftup, float rightup,
                                         float leftdown, float rightdown) {
  // ------------------------------------------------------
  // and some local variables
  // ------------------------------------------------------

  int i;                              // counter for loop
  float wall_pos, wall_pos_left, cor; // wall position
  Coords newbbox[4];

  cout << "extension in the middle" << endl;

  def_time_u = def_time_u + 1; // set deformation time

  runParticle = &refParticle; // start of list

  runParticle = runParticle->prevP; // go back one

  wall_pos = runParticle->xpos; // define wall position

  cor = (wall_pos - 0.99) / 2.0;

  // ------------------------------------------------------------------------
  // run through particles and move their y position starting from upper
  // row
  // Upper row is pushed inwards. Remember that the yposition of upper
  // row
  // dont change in the relaxation
  // routine, they are fixed.
  // Other particles are moved as if all deformation happens in y
  // direction
  // This is a first approximation and makes the relaxation a lot faster
  // and
  // more accurate.
  // ------------------------------------------------------------------------

  for (i = 0; i < numParticles; i++) // go through all particles
  // starting from back
  {
    if (runParticle->xpos > leftup + cor &&
        runParticle->xpos < (wall_pos - (rightup + cor)) &&
        runParticle->ypos > change) {
      runParticle->xpos =
          runParticle->xpos -
          ((runParticle->xpos - (leftup + cor)) * move * wall_pos /
           (wall_pos - (rightup + leftup + cor + cor)));
    }

    if (runParticle->xpos > (wall_pos - (rightup + cor)) &&
        runParticle->ypos > change)
      runParticle->xpos = runParticle->xpos - move;

    if (runParticle->xpos > leftdown + cor &&
        runParticle->xpos <= (wall_pos - (rightdown + cor)) &&
        runParticle->ypos < change) {
      runParticle->xpos =
          runParticle->xpos -
          ((runParticle->xpos - leftdown + cor) * move * wall_pos /
           (wall_pos - (rightdown + leftdown + cor + cor)));
    }

    if (runParticle->xpos > (wall_pos - rightdown + cor) &&
        runParticle->ypos <= change)
      runParticle->xpos = runParticle->xpos - move;

    // if (runParticle->ypos <= neck_height)
    //	runParticle->xpos =
    // runParticle->xpos -
    //(runParticle->xpos * move / wall_pos);

    // runParticle->ypos = runParticle->ypos - (runParticle->ypos *
    // (surf_pos-(volume/(wall_pos+(move*-1))))/surf_pos);

    runParticle->ChangeBox(repBox, particlex);

    runParticle->oldy = runParticle->ypos;

    runParticle = runParticle->prevP;
  }

  right_wall_pos = right_wall_pos - move;
  right_lattice_wall_pos -= move;
  // surf_pos = surf_pos - (surf_pos-(volume/(wall_pos+(move*-1))));
  // cout << surf_pos << endl;
  //  right_lattice_wall_pos = wall_pos + refParticle.radius;
  //  left_lattice_wall_pos = wall_pos_left - refParticle.radius;

  // Relaxation ();		// and do a full relaxation

  // if (plot == 1)
  //  UpdateElle ();	// and update the interface of Elle

  // SetRightLatticeWallPos();

  newbbox[0].x = 0.0;
  newbbox[0].y = 0.0;

  newbbox[1].x = wall_pos;
  newbbox[1].y = 0.0;

  newbbox[2].x = wall_pos;
  newbbox[2].y = 1.0;

  newbbox[3].x = 0.0;
  newbbox[3].y = 1.0;

  // and set the new boxsize

  ElleSetCellBBox(&newbbox[0], &newbbox[1], &newbbox[2], &newbbox[3]);
}

float Lattice::DeformLatticedriftmidnew3(float move, int plot, float change,
                                         float leftup, float rightup,
                                         float leftdown, float rightdown) {
  // ------------------------------------------------------
  // and some local variables
  // ------------------------------------------------------

  int i;                                            // counter for loop
  float wall_pos, wall_pos_left, cor, upper, lower; // wall position
  Coords newbbox[4];

  cout << "extension in the middle" << endl;

  def_time_u = def_time_u + 1; // set deformation time

  runParticle = &refParticle; // start of list

  runParticle = runParticle->prevP; // go back one

  wall_pos = runParticle->xpos; // define wall position

  cor = (wall_pos - 0.99) / 2.0;

  upper = leftup + (move * plot * -0.5);

  lower = leftdown + (move * plot * -0.5);

  // ------------------------------------------------------------------------
  // run through particles and move their y position starting from upper
  // row
  // Upper row is pushed inwards. Remember that the yposition of upper
  // row
  // dont change in the relaxation
  // routine, they are fixed.
  // Other particles are moved as if all deformation happens in y
  // direction
  // This is a first approximation and makes the relaxation a lot faster
  // and
  // more accurate.
  // ------------------------------------------------------------------------

  for (i = 0; i < numParticles; i++) // go through all particles
  // starting from back
  {
    if (runParticle->xpos > upper && runParticle->ypos > change) {
      runParticle->xpos = runParticle->xpos - move;
    }

    if (runParticle->xpos > lower && runParticle->ypos <= change)
      runParticle->xpos = runParticle->xpos - move;

    // if (runParticle->ypos <= neck_height)
    //	runParticle->xpos =
    // runParticle->xpos -
    //(runParticle->xpos * move / wall_pos);

    // runParticle->ypos = runParticle->ypos - (runParticle->ypos *
    // (surf_pos-(volume/(wall_pos+(move*-1))))/surf_pos);

    runParticle->ChangeBox(repBox, particlex);

    runParticle->oldy = runParticle->ypos;

    runParticle = runParticle->prevP;
  }

  right_wall_pos = right_wall_pos - move;
  right_lattice_wall_pos -= move;
  // surf_pos = surf_pos - (surf_pos-(volume/(wall_pos+(move*-1))));
  // cout << surf_pos << endl;
  //  right_lattice_wall_pos = wall_pos + refParticle.radius;
  //  left_lattice_wall_pos = wall_pos_left - refParticle.radius;

  // Relaxation ();		// and do a full relaxation

  // if (plot == 1)
  //  UpdateElle ();	// and update the interface of Elle

  // SetRightLatticeWallPos();

  newbbox[0].x = 0.0;
  newbbox[0].y = 0.0;

  newbbox[1].x = wall_pos;
  newbbox[1].y = 0.0;

  newbbox[2].x = wall_pos;
  newbbox[2].y = 1.0;

  newbbox[3].x = 0.0;
  newbbox[3].y = 1.0;

  // and set the new boxsize

  ElleSetCellBBox(&newbbox[0], &newbbox[1], &newbbox[2], &newbbox[3]);
}

float Lattice::DeformLatticedriftmidnew4(float move, int time, float changeone,
                                         float leftup, float changetwo,
                                         float lefttwo, float changethree,
                                         float leftthree, float leftdown,
                                         float grad) {
  // ------------------------------------------------------
  // and some local variables
  // ------------------------------------------------------

  int i; // counter for loop
  float wall_pos, wall_pos_left, cor, upper, uppertwo, upperthree,
      lower; // wall position
  Coords newbbox[4];

  cout << "extension in the middle" << endl;

  def_time_u = def_time_u + 1; // set deformation time

  runParticle = &refParticle; // start of list

  runParticle = runParticle->prevP; // go back one

  wall_pos = runParticle->xpos; // define wall position

  cor = (wall_pos - 0.99) / 2.0;

  upper = leftup + (move * time * -0.5);

  uppertwo = lefttwo + (move * time * -0.5);

  upperthree = leftthree + (move * time * -0.5);

  lower = leftdown + (move * time * -0.5);

  // ------------------------------------------------------------------------
  // run through particles and move their y position starting from upper
  // row
  // Upper row is pushed inwards. Remember that the yposition of upper
  // row
  // dont change in the relaxation
  // routine, they are fixed.
  // Other particles are moved as if all deformation happens in y
  // direction
  // This is a first approximation and makes the relaxation a lot faster
  // and
  // more accurate.
  // ------------------------------------------------------------------------

  for (i = 0; i < numParticles; i++) // go through all particles
  // starting from back
  {
    if (runParticle->xpos > upper && runParticle->ypos > changeone) {
      if (grad) {
        runParticle->xpos =
            runParticle->xpos - move * (grad * runParticle->ypos);
        runParticle->sheet_x =
            runParticle->sheet_x - move * (grad * runParticle->ypos);
      } else {
        runParticle->xpos = runParticle->xpos - move;
        runParticle->sheet_x = runParticle->sheet_x - move;
      }
    }

    if (runParticle->xpos > uppertwo && runParticle->ypos <= changeone &&
        runParticle->ypos > changetwo) {
      if (grad) {
        runParticle->xpos =
            runParticle->xpos - move * (grad * runParticle->ypos);
        runParticle->sheet_x =
            runParticle->sheet_x - move * (grad * runParticle->ypos);
      } else {
        runParticle->xpos = runParticle->xpos - move;
        runParticle->sheet_x = runParticle->sheet_x - move;
      }
    }

    if (runParticle->xpos > upperthree && runParticle->ypos <= changetwo &&
        runParticle->ypos > changethree) {
      if (grad) {
        runParticle->xpos =
            runParticle->xpos - move * (grad * runParticle->ypos);
        runParticle->sheet_x =
            runParticle->sheet_x - move * (grad * runParticle->ypos);
      } else {
        runParticle->xpos = runParticle->xpos - move;
        runParticle->sheet_x = runParticle->sheet_x - move;
      }
    }

    if (runParticle->xpos > lower && runParticle->ypos <= changethree) {
      if (grad) {
        runParticle->xpos =
            runParticle->xpos - move * (grad * runParticle->ypos);
        runParticle->sheet_x =
            runParticle->sheet_x - move * (grad * runParticle->ypos);
      } else {
        runParticle->xpos = runParticle->xpos - move;
        runParticle->sheet_x = runParticle->sheet_x - move;
      }
    }

    runParticle->ChangeBox(repBox, particlex);

    runParticle->oldy = runParticle->ypos;

    runParticle = runParticle->prevP;
  }

  right_wall_pos = right_wall_pos - move;
  right_lattice_wall_pos -= move;

  newbbox[0].x = 0.0;
  newbbox[0].y = 0.0;

  newbbox[1].x = wall_pos;
  newbbox[1].y = 0.0;

  newbbox[2].x = wall_pos;
  newbbox[2].y = 1.0;

  newbbox[3].x = 0.0;
  newbbox[3].y = 1.0;

  // and set the new boxsize

  ElleSetCellBBox(&newbbox[0], &newbbox[1], &newbbox[2], &newbbox[3]);
}

float Lattice::DeformLatticedriftmidnew5(float move, int time, float changeone,
                                         float leftup, float changetwo,
                                         float lefttwo, float changethree,
                                         float leftthree, float leftdown,
                                         float grad) {
  // ------------------------------------------------------
  // and some local variables
  // ------------------------------------------------------

  int i; // counter for loop
  float wall_pos, wall_pos_left, cor, upper, uppertwo, upperthree,
      lower; // wall position
  Coords newbbox[4];

  cout << "extension in the middle" << endl;

  def_time_u = def_time_u + 1; // set deformation time

  runParticle = &refParticle; // start of list

  runParticle = runParticle->prevP; // go back one

  wall_pos = runParticle->xpos; // define wall position

  cor = (wall_pos - 0.99) / 2.0;

  upper = leftup + (move * time * -0.5);

  uppertwo = lefttwo + (move * time * -0.5);

  upperthree = leftthree + (move * time * -0.5);

  lower = leftdown + (move * time * -0.5);

  // ------------------------------------------------------------------------
  // run through particles and move their y position starting from upper
  // row
  // Upper row is pushed inwards. Remember that the yposition of upper
  // row
  // dont change in the relaxation
  // routine, they are fixed.
  // Other particles are moved as if all deformation happens in y
  // direction
  // This is a first approximation and makes the relaxation a lot faster
  // and
  // more accurate.
  // ------------------------------------------------------------------------

  for (i = 0; i < numParticles; i++) // go through all particles
  // starting from back
  {
    if (runParticle->xpos > upper && runParticle->ypos > changeone) {
      if (grad) {
        runParticle->xpos =
            runParticle->xpos - move * (grad * runParticle->ypos);
        runParticle->sheet_x =
            runParticle->sheet_x - move * (grad * runParticle->ypos);
      } else {
        runParticle->xpos = runParticle->xpos - move;
        runParticle->sheet_x = runParticle->sheet_x - move;
      }
    }

    if (runParticle->xpos > uppertwo + cor && runParticle->ypos <= changeone &&
        runParticle->ypos > changetwo) {
      if (runParticle->xpos < upperthree) {

        runParticle->xpos =
            runParticle->xpos -
            ((runParticle->xpos - (uppertwo + cor)) * move * wall_pos /
             (wall_pos - (upperthree + uppertwo + cor + cor)));
        runParticle->sheet_x =
            runParticle->sheet_x -
            ((runParticle->sheet_x - (uppertwo + cor)) * move * wall_pos /
             (wall_pos - (upperthree + uppertwo + cor + cor)));
      } else {
        runParticle->xpos = runParticle->xpos - move;
        runParticle->sheet_x = runParticle->sheet_x - move;
      }
    }

    if (runParticle->xpos > lower && runParticle->ypos <= changethree) {
      if (grad) {
        runParticle->xpos =
            runParticle->xpos - move * (grad * runParticle->ypos);
        runParticle->sheet_x =
            runParticle->sheet_x - move * (grad * runParticle->ypos);
      } else {
        runParticle->xpos = runParticle->xpos - move;
        runParticle->sheet_x = runParticle->sheet_x - move;
      }
    }

    runParticle->ChangeBox(repBox, particlex);

    runParticle->oldy = runParticle->ypos;

    runParticle = runParticle->prevP;
  }

  right_wall_pos = right_wall_pos - move;
  right_lattice_wall_pos -= move;

  newbbox[0].x = 0.0;
  newbbox[0].y = 0.0;

  newbbox[1].x = wall_pos;
  newbbox[1].y = 0.0;

  newbbox[2].x = wall_pos;
  newbbox[2].y = 1.0;

  newbbox[3].x = 0.0;
  newbbox[3].y = 1.0;

  // and set the new boxsize

  ElleSetCellBBox(&newbbox[0], &newbbox[1], &newbbox[2], &newbbox[3]);
}

float Lattice::InitDriftList(int pos) {
  int i;
  int box;

  for (i = 0; i < 200; i++) {
    box = (i * particlex * 2.0) + pos;

    driftlist[i] = repBox[box];
  }
}

float Lattice::DeformLatticedrift6(float move, int time, float changeone,
                                   float leftup, float changetwo,
                                   float lefttwo) {
  // ------------------------------------------------------
  // and some local variables
  // ------------------------------------------------------

  int i; // counter for loop
  float wall_pos, wall_pos_left, cor, upper, uppertwo, upperthree, lower,
      ran; // wall position
  Coords newbbox[4];

  cout << "extension in the middle" << endl;

  def_time_u = def_time_u + 1; // set deformation time

  runParticle = &refParticle; // start of list

  runParticle = runParticle->prevP; // go back one

  wall_pos = runParticle->xpos; // define wall position

  cor = (wall_pos - 0.99) / 2.0;

  upper = leftup + (move * time * -0.5);

  uppertwo = lefttwo + (move * time * -0.5);

  // ------------------------------------------------------------------------
  // run through particles and move their y position starting from upper
  // row
  // Upper row is pushed inwards. Remember that the yposition of upper
  // row
  // dont change in the relaxation
  // routine, they are fixed.
  // Other particles are moved as if all deformation happens in y
  // direction
  // This is a first approximation and makes the relaxation a lot faster
  // and
  // more accurate.
  // ------------------------------------------------------------------------

  for (i = 0; i < numParticles; i++) // go through all particles
  // starting from back
  {
    ran = rand() / (float)RAND_MAX;
    // if (ran < 0.99) //0.9, 0.7
    // an = -1;
    // else
    // ran = 1;

    // if (runParticle->xpos > upper  && runParticle->ypos > changeone) //
    // everything right of the plate, fixed move to the right
    //{

    //{
    // runParticle->xpos = runParticle->xpos - move;
    // runParticle->sheet_x = runParticle->sheet_x- move;
    //}
    //}

    if (runParticle->xpos >
        0.25) // everything right of the always fixed position
    {
      if (driftlist[int(
              (runParticle->ypos) *
              particlex)]) // localize this is the list with particles in y
      {
        if (runParticle->xpos >
            (driftlist[int((runParticle->ypos) * particlex)]
                 ->xpos)) // all particles right of the ridge move right
        {
          runParticle->xpos = runParticle->xpos - move;
          runParticle->sheet_x = runParticle->sheet_x - move;
        } else if (runParticle->xpos ==
                   (driftlist[int((runParticle->ypos) * particlex)]->xpos)) {
          if (ran < 0.5) {
            runParticle->xpos = runParticle->xpos - move;
            runParticle->sheet_x = runParticle->sheet_x - move;
          }
        }
      } else // no drift list for this
      {
        if (runParticle->xpos < 0.75) // average
        {
          cout << "average" << endl;
          // runParticle->xpos = runParticle->xpos - ((runParticle->xpos-0.25)
          // * move *2.0); runParticle->sheet_x = runParticle->sheet_x-
          // ((runParticle->sheet_x-0.25) * move*2.0 );
          runParticle->xpos = runParticle->xpos - move;
          runParticle->sheet_x = runParticle->sheet_x - move;
        } else {
          runParticle->xpos = runParticle->xpos - move;
          runParticle->sheet_x = runParticle->sheet_x - move;
        }
      }
    }

    runParticle->ChangeBox(repBox, particlex);

    runParticle->oldy = runParticle->ypos;

    runParticle = runParticle->prevP;
  }

  right_wall_pos = right_wall_pos - move;
  right_lattice_wall_pos -= move;

  newbbox[0].x = 0.0;
  newbbox[0].y = 0.0;

  newbbox[1].x = wall_pos;
  newbbox[1].y = 0.0;

  newbbox[2].x = wall_pos;
  newbbox[2].y = 1.0;

  newbbox[3].x = 0.0;
  newbbox[3].y = 1.0;

  // and set the new boxsize

  ElleSetCellBBox(&newbbox[0], &newbbox[1], &newbbox[2], &newbbox[3]);
}

float Lattice::DeformLatticedrift7(float move, int time, float changeone,
                                   float leftup, float changetwo,
                                   float lefttwo) {
  // ------------------------------------------------------
  // and some local variables
  // ------------------------------------------------------

  int i; // counter for loop
  float wall_pos, wall_pos_left, cor, upper, uppertwo, upperthree, lower,
      ran; // wall position
  Coords newbbox[4];

  cout << "extension in the middle" << endl;

  def_time_u = def_time_u + 1; // set deformation time

  runParticle = &refParticle; // start of list

  runParticle = runParticle->prevP; // go back one

  wall_pos = runParticle->xpos; // define wall position

  cor = (wall_pos - 0.99) / 2.0;

  upper = leftup + (move * time * -0.5);

  uppertwo = lefttwo + (move * time * -0.5);

  // ------------------------------------------------------------------------
  // run through particles and move their y position starting from upper
  // row
  // Upper row is pushed inwards. Remember that the yposition of upper
  // row
  // dont change in the relaxation
  // routine, they are fixed.
  // Other particles are moved as if all deformation happens in y
  // direction
  // This is a first approximation and makes the relaxation a lot faster
  // and
  // more accurate.
  // ------------------------------------------------------------------------

  for (i = 0; i < numParticles; i++) // go through all particles
  // starting from back
  {
    ran = rand() / (float)RAND_MAX;
    // if (ran < 0.99) //0.9, 0.7
    // an = -1;
    // else
    // ran = 1;

    // if (runParticle->xpos > upper  && runParticle->ypos > changeone) //
    // everything right of the plate, fixed move to the right
    //{

    //{
    // runParticle->xpos = runParticle->xpos - move;
    // runParticle->sheet_x = runParticle->sheet_x- move;
    //}
    //}

    if (runParticle->xpos >
        0.25) // everything right of the always fixed position
    {
      if (driftlist[int(
              (runParticle->ypos) *
              particlex)]) // localize this is the list with particles in y
      {
        if (runParticle->xpos >
            (driftlist[int((runParticle->ypos) * particlex)]
                 ->xpos)) // all particles right of the ridge move right
        {
          runParticle->xpos = runParticle->xpos - move;
          runParticle->sheet_x = runParticle->sheet_x - move;
        } else if (runParticle->xpos ==
                   (driftlist[int((runParticle->ypos) * particlex)]->xpos)) {
          if (ran < 0.5) {
            runParticle->xpos = runParticle->xpos - move;
            runParticle->sheet_x = runParticle->sheet_x - move;
          }
        }
      } else // no drift list for this
      {
        // if (runParticle->xpos < 0.75)  // average
        {
          // cout << "average" << endl;
          cout << "average" << runParticle->ypos << "two" << endl;
          // runParticle->xpos = runParticle->xpos - ((runParticle->xpos-0.25)
          // * move *2.0); runParticle->sheet_x = runParticle->sheet_x-
          // ((runParticle->sheet_x-0.25) * move*2.0 );
        }
        // else
        {
          // runParticle->xpos = runParticle->xpos - move;
          // runParticle->sheet_x = runParticle->sheet_x- move;
        }
      }
    }

    runParticle->ChangeBox(repBox, particlex);

    runParticle->oldy = runParticle->ypos;

    runParticle = runParticle->prevP;
  }

  right_wall_pos = right_wall_pos - move;
  right_lattice_wall_pos -= move;

  newbbox[0].x = 0.0;
  newbbox[0].y = 0.0;

  newbbox[1].x = wall_pos;
  newbbox[1].y = 0.0;

  newbbox[2].x = wall_pos;
  newbbox[2].y = 1.0;

  newbbox[3].x = 0.0;
  newbbox[3].y = 1.0;

  // and set the new boxsize

  ElleSetCellBBox(&newbbox[0], &newbbox[1], &newbbox[2], &newbbox[3]);
}

//-----------------------------------------------------------------------
// this simply measures movement in y, written for faults initially
// daniel 2018
//-----------------------------------------------------------------------

float Lattice::Save_Particle_Velocity() {
  for (i = 0; i < numParticles; i++) {

    runParticle->fracture_reaction = runParticle->initialy - runParticle->ypos;
    // runParticle->oldy = runParticle->ypos;

    runParticle = runParticle->prevP;
  }
}

// DEFORM_LATTICE_PURESHEAR
/*******************************************************************************
 * This function does whatever deformation we want to perform with the lattice
 * It is called from the outside by the User and does all that is needed for
 * a timestep including all the relaxation etc.
 * at the moment it receives a float for how much the walls are moved
 * as a deformation it moves upper and lower row of particles inwards if
 * positive value. then it uses area conservation to calculate movement in
 * x direction of the right wall. Note that at the moment elle still scales
 * everything to the window so that it looks as if the x dimension does not
 * change and the volume shrinks...
 *
 * Here we can include all kinds of stuff in the future
 *
 * Daniel and Rianne December 2002,
 * area bug fixed Daniel and Jochen Feb. 2003
 ********************************************************************************/

// ---------------------------------------------------------------------------------
// function DeformLatticePUREShear in Lattice class
// receives a number for deformation
// dont know why I return something ?
//
// called from the outside ! i.e. mike_elle whatever routine
// ---------------------------------------------------------------------------------

float Lattice::DeformLatticePureShear(float move, int plot) {
  // ------------------------------------------------------
  // and some local variables
  // ------------------------------------------------------

  int i;                                  // counter for loop
  float wall_pos;                         // wall position
  float movex, position_Xwall, wall_Xpos; // stuff for x wall

  cout << " Pure Shear Deformation " << endl;

  Coords newbbox[4]; // have to change Elle box size in file

  def_time_p = def_time_p + 1; // set deformation time

  runParticle = &refParticle; // start of list

  runParticle = runParticle->prevP; // go back one

  wall_pos = runParticle->ypos; // get wall position y

  wall_Xpos = runParticle->xpos; // get wall position x

  // have to check if particle on edge is furthest in x or particle
  // below it
  // therefore we run one x row below the highest particle and check
  // that p as well

  for (i = 0; i < particlex; i++) {
    runParticle = runParticle->prevP;
  }

  if (runParticle->xpos > wall_Xpos) {
    wall_Xpos = runParticle->xpos;
  }
  // now we can calculate the initial area of the box at time 1

  // if (def_time_p == 1)
  {

    area_box = wall_pos * wall_Xpos;
    cout << "area of box: " << area_box << endl;
  }
  // ------------------------------------------------------------------------
  // run through particles and move their y position starting from upper
  // row
  // Upper row is pushed inwards. Remember that the yposition of upper
  // row
  // dont change in the relaxation
  // routine, they are fixed.
  // Other particles are moved as if all deformation happens in y
  // direction
  // This is a first approximation and makes the relaxation a lot faster
  // and
  // more accurate.
  // ------------------------------------------------------------------------

  for (i = 0; i < numParticles; i++) // go through all particles
                                     // starting from back
  {
    runParticle->ypos =
        runParticle->ypos - (runParticle->ypos * move / wall_pos);

    runParticle->ChangeBox(repBox, particlex); // always call
    // Repbox if
    // positions
    // changes

    runParticle = runParticle->prevP;
  }

  if (walls)
    upper_wall_pos = upper_wall_pos - (upper_wall_pos * move / wall_pos);

  position_Xwall = area_box / (wall_pos - move);

  movex = wall_Xpos - position_Xwall; // now calculate x position
  // asuming area remains constant

  movex = -move * wall_Xpos / (wall_pos - move);

  for (i = 0; i < numParticles; i++) // go through all particles
                                     // starting from back
  {
    // average the deformation in x..

    runParticle->xpos =
        runParticle->xpos - (runParticle->xpos * movex / wall_Xpos);

    runParticle->ChangeBox(repBox, particlex);

    runParticle = runParticle->prevP;
  }

  right_lattice_wall_pos =
      right_lattice_wall_pos -
      (right_lattice_wall_pos * movex / right_lattice_wall_pos);

  if (walls)
    right_wall_pos = right_wall_pos - (right_wall_pos * movex / wall_Xpos);

  // now define the new boxsize for elle

  newbbox[0].x = 0.0;
  newbbox[0].y = 0.0;

  newbbox[1].x = position_Xwall;
  newbbox[1].y = 0.0;

  newbbox[2].x = position_Xwall;
  newbbox[2].y = wall_pos - move;

  newbbox[3].x = 0.0;
  newbbox[3].y = wall_pos - move;

  // and set the new boxsize

  ElleSetCellBBox(&newbbox[0], &newbbox[1], &newbbox[2], &newbbox[3]);

  for (i = 0; i < numParticles; i++) {
    runParticle->ChangeBox(repBox, particlex);

    runParticle = runParticle->nextP;
  }

  Relaxation(); // and do a full relaxation

  if (plot == 1)
    UpdateElle(); // and update the interface of Elle if it
                  // is wanted

  SetRightLatticeWallPos();
}

// DEFORM_LATTICE_PURESHEAR_AREACHANGE
/*******************************************************************************
 * This function does a coaxial deformation including area change
 * It is called from the outside by the User and does all that is needed for
 * a timestep including all the relaxation etc.
 * at the moment it receives a float for how much the walls are moved
 * as a deformation it moves upper and lower row of particles inwards if
 * positive value. then it uses area change given to calculate movement in
 * x direction of the right wall. Note that at the moment elle still scales
 * everything to the window so that it looks as if the x dimension does not
 * change and the volume shrinks...
 *
 * Here we can include all kinds of stuff in the future
 *
 * Daniel July 2003
 ********************************************************************************/

// ---------------------------------------------------------------------------------
// function DeformLatticePUREShearAreaChange in Lattice class
// receives a number for deformation
// dont know why I return something ?
//
// called from the outside ! i.e. mike_elle whatever routine
// ---------------------------------------------------------------------------------

float Lattice::DeformLatticePureShearAreaChange(float move, float area_change,
                                                int plot) {
  // ------------------------------------------------------
  // and some local variables
  // ------------------------------------------------------

  int i;                                  // counter for loop
  float wall_pos;                         // wall position
  float movex, position_Xwall, wall_Xpos; // stuff for x wall

  Coords newbbox[4]; // have to change Elle box size in file

  def_time_p = def_time_p + 1; // set deformation time

  runParticle = &refParticle; // start of list

  runParticle = runParticle->prevP; // go back one

  wall_pos = runParticle->ypos; // get wall position y

  wall_Xpos = runParticle->xpos; // get wall position x

  // have to check if particle on edge is furthest in x or particle
  // below it
  // therefore we run one x row below the highest particle and check
  // that p as well

  for (i = 0; i < particlex; i++) {
    runParticle = runParticle->prevP;
  }

  if (runParticle->xpos > wall_Xpos) {
    wall_Xpos = runParticle->xpos;
  }
  // now we can calculate the initial area of the box at time 1

  if (def_time_p == 1) {
    area_box = wall_pos * wall_Xpos;
    cout << area_box << endl;
  }

  area_box = area_box + (area_box * area_change);

  // ------------------------------------------------------------------------
  // run through particles and move their y position starting from upper
  // row
  // Upper row is pushed inwards. Remember that the yposition of upper
  // row
  // dont change in the relaxation
  // routine, they are fixed.
  // Other particles are moved as if all deformation happens in y
  // direction
  // This is a first approximation and makes the relaxation a lot faster
  // and
  // more accurate.
  // ------------------------------------------------------------------------

  for (i = 0; i < numParticles; i++) // go through all particles
                                     // starting from back
  {
    runParticle->ypos =
        runParticle->ypos - (runParticle->ypos * move / wall_pos);

    runParticle->ChangeBox(repBox, particlex); // always call
    // Repbox if
    // positions
    // changes

    runParticle = runParticle->prevP;
  }

  // area is not exactly 1.0, but I dont know how this changes when
  // amount of particles changes
  // is now for 400 times something particles.

  position_Xwall = area_box / (wall_pos - move);

  movex = wall_Xpos - position_Xwall; // now calculate x position
  // asuming area remains constant

  for (i = 0; i < numParticles; i++) // go through all particles
                                     // starting from back
  {
    // average the deformation in x..

    runParticle->xpos =
        runParticle->xpos - (runParticle->xpos * movex / wall_Xpos);

    runParticle->ChangeBox(repBox, particlex);

    runParticle = runParticle->prevP;
  }

  // now define the new boxsize for elle

  newbbox[0].x = 0.0;
  newbbox[0].y = 0.0;

  newbbox[1].x = position_Xwall;
  newbbox[1].y = 0.0;

  newbbox[2].x = position_Xwall;
  newbbox[2].y = wall_pos - move;

  newbbox[3].x = 0.0;
  newbbox[3].y = wall_pos - move;

  // and set the new boxsize

  ElleSetCellBBox(&newbbox[0], &newbbox[1], &newbbox[2], &newbbox[3]);

  Relaxation(); // and do a full relaxation

  if (plot == 1)
    UpdateElle(); // and update the interface of Elle if it
  // is wanted

  SetRightLatticeWallPos();
}

/***********************************************************************
 * Function that deform uniaxially by pushing down upper wall
 * no average is taken so that only the upper wall particles
 * are moved. This will result in serious gradients in deformation
 * if the relaxation threshold is not seriously reduced.
 * Very slow in relaxation
 ***********************************************************************/

float Lattice::DeformLatticeNoAverage(float move, int plot) {
  // ------------------------------------------------------
  // and some local variables
  // ------------------------------------------------------

  int i;          // counter for loop
  float wall_pos; // wall position

  def_time_u = def_time_u + 1; // set deformation time

  runParticle = &refParticle; // start of list

  runParticle = runParticle->prevP; // go back one

  wall_pos = runParticle->ypos; // define wall position

  // ------------------------------------------------------------------------
  // run through particles and move their y position starting from upper
  // row
  // Upper row is pushed inwards. Remember that the yposition of upper
  // row
  // dont change in the relaxation
  // routine, they are fixed.
  // Other particles are moved as if all deformation happens in y
  // direction
  // This is a first approximation and makes the relaxation a lot faster
  // and
  // more accurate.
  // ------------------------------------------------------------------------

  for (i = 0; i < particlex; i++) // go through all particles
                                  // starting from back
  {
    runParticle->ypos =
        runParticle->ypos - (runParticle->ypos * move / wall_pos);

    runParticle->ChangeBox(repBox, particlex);

    runParticle = runParticle->prevP;
  }

  Relaxation(); // and do a full relaxation

  if (plot == 1)
    UpdateElle(); // and update the interface of Elle

  SetRightLatticeWallPos();
}

/***********************************************************************
 *	function that deforms uniaxially but now moves the upper and
 * lower wall particles inwards (or outwards if move is negative)
 * no average is taken of the deformation within the box. This
 * will result in gradients in deformation if the relaxation threshold
 * is not lowered.
 *
 * variable move indicates vertical strain per deformation step
 * variable plot indicates whether or not a picture is taken
 * after the deformation step (0 = no picture, 1 = picture)
 ***********************************************************************/

float Lattice::DeformLatticeNoAverage2side(float move, int plot) {
  // ------------------------------------------------------
  // and some local variables
  // ------------------------------------------------------

  int i;          // counter for loop
  float wall_pos; // wall position

  def_time_u = def_time_u + 1; // set deformation time

  runParticle = &refParticle; // start of list

  runParticle = runParticle->prevP; // go back one

  // ------------------------------------------------------------------------
  // run through particles and move their y position starting from upper
  // row
  // Upper row is pushed inwards. Remember that the yposition of upper
  // row
  // dont change in the relaxation
  // routine, they are fixed.
  // Other particles are moved as if all deformation happens in y
  // direction
  // This is a first approximation and makes the relaxation a lot faster
  // and
  // more accurate.
  // ------------------------------------------------------------------------

  for (i = 0; i < particlex; i++) // go through all particles
                                  // starting from back, only upper
                                  // row
  {
    runParticle->ypos = runParticle->ypos - move;

    runParticle->ChangeBox(repBox, particlex);

    runParticle = runParticle->prevP;
  }

  runParticle = &refParticle; // start of list

  // ------------------------------------------------------------------------
  // run through particles and move their y position starting from upper
  // row
  // Upper row is pushed inwards. Remember that the yposition of upper
  // row
  // dont change in the relaxation
  // routine, they are fixed.
  // Other particles are moved as if all deformation happens in y
  // direction
  // This is a first approximation and makes the relaxation a lot faster
  // and
  // more accurate.
  // ------------------------------------------------------------------------

  for (i = 0; i < particlex; i++) // go through all particles
                                  // starting from back lower row
  {
    runParticle->ypos = runParticle->ypos + move;

    runParticle->ChangeBox(repBox, particlex);

    runParticle = runParticle->nextP;
  }

  Relaxation(); // and do a full relaxation

  if (plot == 1)
    UpdateElle(); // and update the interface of Elle

  SetRightLatticeWallPos();
}

/*************************************************************************************
 * function that deforms uniaxially. Upper and lower row are pushed a factor
 *"move" inwards. Average in this case means that all particles are just moved
 *by the factor move. Only works if there is an open surface in the middle of
 *the model and two materials are pressed together. Function dissolve x row
 *gives the input for which part of the crystal is upper and which lower.
 *
 * if plot is 0 no picture is taken, if it is 1 a picture is taken
 *************************************************************************************/

float Lattice::DeformLatticeNewAverage2side(float move, int plot) {
  // ------------------------------------------------------
  // and some local variables
  // ------------------------------------------------------

  int i, j, count; // counter for loop
  float wall_pos;  // wall position

  def_time_u = def_time_u + 1; // set deformation time

  runParticle = &refParticle; // start of list

  // ------------------------------------------------------------------------
  // run through particles and move their y position starting from upper
  // row
  // Upper row is pushed inwards. Remember that the yposition of upper
  // row
  // dont change in the relaxation
  // routine, they are fixed.
  // Other particles are moved as if all deformation happens in y
  // direction
  // This is a first approximation and makes the relaxation a lot faster
  // and
  // more accurate.
  // ------------------------------------------------------------------------

  for (i = 0; i < numParticles; i++) // go through all particles
                                     // starting from back
  {
    if (!runParticle->isHole) {
      if (!runParticle->isHoleBoundary) {
        if (runParticle->isUpperBox) // upper part move downwards

          runParticle->ypos = runParticle->ypos - move;

        else // lower part move upwards

          runParticle->ypos = runParticle->ypos + move;

        runParticle->ChangeBox(repBox, particlex);
      } else {
        count = 0;
        for (j = 0; j < 8; j++) {
          if (runParticle->neigP[j]) {
            if (!runParticle->neigP[j]->isHole)
              count++;
          }
        }
        if (count > 1) {
          if (runParticle->isUpperBox) // upper part move downwards

            runParticle->ypos = runParticle->ypos - move;

          else // lower part move upwards

            runParticle->ypos = runParticle->ypos + move;

          runParticle->ChangeBox(repBox, particlex);
        } else
          runParticle->rate_factor = runParticle->rate_factor * 20.0;
      }
    }

    runParticle = runParticle->nextP;
  }

  runParticle = &refParticle; // start of list

  Relaxation(); // and do a full relaxation

  if (plot == 1)
    UpdateElle(); // and update the interface of Elle

  SetRightLatticeWallPos();
}

/*****************************************************************************
 * just a dummy function that looks for the highest grain number
 *****************************************************************************/

int Lattice::HighestGrain() {
  int i;
  int highest;

  highest = 0;

  runParticle = &refParticle;

  for (i = 0; i < numParticles; i++) {
    if (runParticle->grain > highest) {
      highest = runParticle->grain;
    }
    runParticle = runParticle->nextP;
  }
  return (highest);
}

// CHANGEPARTICLE (int number, float shrink)
/*********************************************************************************
 * simple function to shrink particle of number x by a certain given amount
 *
 * daniel spring 2002
 ********************************************************************************/

// ----------------------------------------------------------------------------
// function ChangeParticle in lattice class
//
// receives number of particle and amount to change radius of particle
// shrink in percent 0.9 is 90%
//
// called from outside, i.e. mike.elle by user
// ----------------------------------------------------------------------------

void Lattice::ChangeParticle(int number, float shrink, int plot) {
  int i, j; // just a counter

  runParticle = &refParticle;

  for (i = 0; i < numParticles; i++) {
    // ------------------------------------------------------------------
    // run through all particles, if I find my particle I shrink it
    // ------------------------------------------------------------------

    if (runParticle->nb == number) {
      runParticle->radius = runParticle->radius * shrink;

      runParticle->a = runParticle->a - runParticle->a * shrink;
      runParticle->b = runParticle->b - runParticle->b * shrink;

      for (j = 0; j < 8; j++) {
        if (runParticle->neigP[j]) {
          runParticle->rad[j] = runParticle->radius;
        }
      }
    }

    runParticle = runParticle->nextP;
  }

  Relaxation(); // and a full relax
  if (plot == 1)
    UpdateElle(); // and some elle update for plotting and
  // stress
}

// CHANGEPARTICLESTRESS(float shrink)
/*********************************************************************
 * function that shrinks the particle that has the highest pressure
 * where pressure is sigma one plus sigma two
 *
 * daniel spring 2002
 ********************************************************************/

// ----------------------------------------------------------------------
// function ChangeParticleStress in Lattice class
//
// receives amount to shrink particle
// shrink in percent 0.9 is 90%
//
// called from the outside by the user in a deformation loop for example
// in mike.elle
// -----------------------------------------------------------------------

void Lattice::ChangeParticleStress(float shrink, int plot) {
  int i, j;                     // just a counter
  float pressure, max_pressure; // for pressure

  runParticle = preRunParticle = &refParticle;

  // -------------------------------------------------------------------
  // first define pressure (sxx plus syy) memorize first particle
  // pressure and then loop through all particles and find particle
  // with highest pressure
  // --------------------------------------------------------------------

  pressure = max_pressure = runParticle->sxx + runParticle->syy;

  for (i = 0; i < numParticles; i++) // loop loop
  {
    pressure = runParticle->sxx + runParticle->syy; // find pressure

    if (pressure < max_pressure) // if higher
    {
      max_pressure = pressure;      // define new
      preRunParticle = runParticle; // and set pointer
    }
    runParticle = runParticle->nextP;
  }
  // -------------------------------------------------------------
  // shrink the particle with highest pressure
  // -------------------------------------------------------------

  preRunParticle->radius = preRunParticle->radius * shrink;

  preRunParticle->a = preRunParticle->a - preRunParticle->a * shrink;
  preRunParticle->b = preRunParticle->b - preRunParticle->b * shrink;

  for (j = 0; j < 8; j++) {
    if (preRunParticle->neigP[j]) {
      preRunParticle->rad[j] = preRunParticle->radius;
    }
  }

  Relaxation(); // and a full relax
  if (plot == 1)
    UpdateElle(); // and some elle update
}

/**********************************************************
 * routine shrinks the grain nb by an amount shrink.
 * shrink is given in percent (0.0 is no change,
 * smaller 1.0 is growth and larger shrinking
 * routine called from the outside, normally experiment.cc
 *
 * Daniel summer 2002
 **********************************************************/

void Lattice::ShrinkGrain(int nb, float shrink, int plot) {
  int i, j; // run variable
  double ratio;

  runParticle = &refParticle; // start

  for (i = 0; i < numParticles; i++) // loop through particles
  {
    if (runParticle->grain == nb) // if particle is in the right
                                  // grain
    {

      ;
      runParticle->radius =
          runParticle->radius + runParticle->radius * shrink; // shrink

      ratio = 1 - shrink;
      for (j = 0; j < 8; j++) {
        if (runParticle->neigP[j]) {
          runParticle->neigpos[j][0] *= ratio; // change
          runParticle->neigpos[j][1] *= ratio;
        }
      }

      // runParticle->radius = runParticle->radius-runParticle->radius * shrink;
      // // shrink
      // runParticle->a = runParticle->a - runParticle->a * shrink;
      // runParticle->b = runParticle->b - runParticle->b * shrink;

      // for (j = 0; j < 8; j++)
      //{

      // if (runParticle->neigP[j])
      //  {
      // runParticle->neigpos[j][0] *= shrink;	// change
      //  runParticle->neigpos[j][1] *= shrink;
      //  individual
      //  radius
      //}

      //}
    }

    runParticle = runParticle->nextP; // go to next
  }

  // Relaxation ();		// and relax
  // if (plot == 1)
  // UpdateElle ();	// and update elle (plot picture on
  // interface)
}

/*************************************************************
 * function used to shrink the whole box
 * all particles shrink or expand by the same amount
 * good for mudcrack structures !
 *
 * Daniel 2004/5
 *************************************************************/

void Lattice::ShrinkBox(float shrink, int plot, int yes) {
  int i, j;
  float percent;

  percent = 1 - shrink;

  for (i = 0; i < numParticles; i++) {
    if (!runParticle->is_lattice_boundary) // dont want to shrink
                                           // boundaries
    {
      if (!yes) {
        runParticle->radius =
            runParticle->radius - runParticle->radius * shrink; // shrink
        // runParticle->a = runParticle->a - runParticle->a * shrink;
        // runParticle->b = runParticle->b - runParticle->b * shrink;
        for (j = 0; j < 8; j++) {
          if (runParticle->neigP[j]) {
            runParticle->neigpos[j][0] *= percent; // change
            runParticle->neigpos[j][1] *= percent;
            // individual
            // radius
          }
        }
      } else if (yes == 1) {
        runParticle->radius =
            runParticle->radius -
            runParticle->radius * shrink / runParticle->young; // shrink
        runParticle->a =
            runParticle->a - runParticle->a * shrink / runParticle->young;
        runParticle->b =
            runParticle->b - runParticle->b * shrink / runParticle->young;
      } else if (yes == 2) {
        runParticle->radius =
            runParticle->radius -
            runParticle->radius * shrink * runParticle->conc; // shrink
        runParticle->a =
            runParticle->a - runParticle->a * shrink / runParticle->young;
        runParticle->b =
            runParticle->b - runParticle->b * shrink / runParticle->young;
      } else if (yes == 3) {
        if (runParticle->conc != 0.0)
          runParticle->radius =
              runParticle->radius - runParticle->radius * shrink -
              runParticle->radius * shrink * runParticle->conc; // shrink
        // runParticle->a = runParticle->a - runParticle->a *
        // shrink/runParticle->young; runParticle->b = runParticle->b -
        // runParticle->b * shrink/runParticle->young;
      } else if (yes == 4) {
        // if (runParticle->conc != 0.0)
        if (runParticle->xpos > 0.45 & runParticle->xpos < 0.55)
          runParticle->radius = runParticle->radius -
                                runParticle->radius * shrink * 2.0; // shrink
        // runParticle->a = runParticle->a - runParticle->a *
        // shrink/runParticle->young; runParticle->b = runParticle->b -
        // runParticle->b * shrink/runParticle->young;
      }
      // !!

      for (j = 0; j < 8; j++) {
        if (runParticle->neigP[j]) {
          runParticle->rad[j] = runParticle->radius; // change
          // individual
          // radius
        }
      }
    }
    runParticle = runParticle->nextP; // and go to next !
  }
  Relaxation(); // and a full relax
  if (plot == 1)
    UpdateElle(0); // and some elle update for plotting and
  // stress
}

/*********************************************************************
 * routine to define the properties of the grain boundaries
 * first property is changing the spring constant ( in percent, 1.0 is
 * no change, larger than 1.0 is stronger, smaller is weaker).
 * second number changes the tensile strength of the grainboundary
 * springs (in percent, larger 1.0 is stronger, smaller is break
 * easier.
 * function is called from mike.elle. normally after lattice is build
 *
 * daniel autumn 2002
 *********************************************************************/

void Lattice::MakeGrainBoundaries(float constant, float break_strength) {
  int i, j; // counters

  boundary_strength = break_strength;
  boundary_constant = constant;

  runParticle = &refParticle; // start

  for (i = 0; i < numParticles; i++) // run along
  {
    if (runParticle->is_boundary) // if boundaryparticle
    {
      for (j = 0; j < 8; j++) // loop through neighbours
      {
        if (runParticle->neigP[j]) // if neighbour
        {
          if (runParticle->neigP[j]->grain != runParticle->grain) // if
                                                                  // neighbour
                                                                  // is
                                                                  // in
                                                                  // different
                                                                  // grain
          {
            runParticle->spring_boundary[j] = true;

            runParticle->springf[j] =
                runParticle->springf[j] * constant; // change
                                                    // spring
            runParticle->springs[j] = runParticle->springs[j] * constant;
            runParticle->break_Str[j] =
                runParticle->break_Str[j] * break_strength; // change
            // breakstrength
          }
        }
      }
    }
    runParticle = runParticle->nextP; // and loop
  }

  // Relaxation ();		// and relax
}

/************************************************************************
 * just a clean function that takes an average of springs that have
 * different constants on both sides. Uses the youngs modulus for that
 *
 * Daniel 2004/5
 ************************************************************************/

void Lattice::AdjustParticleConstants() {
  int i, j;

  double E;
  double poisson;

  for (i = 0; i < numParticles; i++) // run along
  {
    for (j = 0; j < 8; j++) {
      if (runParticle->neigP[j]) {
        E = (runParticle->E + runParticle->neigP[j]->E) / 2.0;
        poisson = (runParticle->poisson_ratio +
                   runParticle->neigP[j]->poisson_ratio) /
                  2.0;

        runParticle->springf[j] = sqrt(3.0) * E / (3.0 * (1.0 - poisson));
        runParticle->springs[j] =
            (1.0 - 3.0 * poisson) / (1.0 + poisson) * runParticle->springf[j];
        runParticle->springv[j] =
            ((runParticle->viscosity + runParticle->neigP[j]->viscosity) *
             sqrt(3.0) / 2.0) /
            2.0;
      }
    }
    runParticle = runParticle->nextP;
  }
}

/******************************************************************************
 * this function adjusts the spring constant of grain boundary springs if
 * the two particles have different spring constants. weakening grains does
 * not affect the spring constant of boundary springs !
 * Function should be called after constants of all grains are defined
 *
 * if you want to have different strength for boundary springs than the grains
 * you should call this function first and then afterwards change the spring
 * strength with the MakeGrainBoundaries function
 *
 * daniel december 2002
 ******************************************************************************/

void Lattice::AdjustConstantGrainBoundaries() {
  int i, j, jj, jjj; // counters

  float grain_one_strength, grain_two_strength; // strength of the two grains

  runParticle = &refParticle; // start

  for (i = 0; i < numParticles; i++) // run along
  {
    if (runParticle->is_boundary) // if boundaryparticle
    {
      // first find your own spring constant
      // so find a spring within the grain

      for (jjj = 0; jjj < 8; jjj++) {
        if (runParticle->neigP[jjj]) {
          if (runParticle->neigP[jjj]->grain == runParticle->grain) {
            grain_one_strength = runParticle->springf[jjj];
          }
        }
      }

      // then find the neighbour spring and strength of neighbours
      // may be more than one....

      for (j = 0; j < 8; j++) // loop through neighbours
      {
        if (runParticle->neigP[j]) // if neighbour
        {
          // have to find springs that connect boundaries and
          // then calculate
          // mean of spring constant of both grains.

          if (runParticle->neigP[j]->grain != runParticle->grain) // if
                                                                  // neighbour
                                                                  // is
                                                                  // in
                                                                  // different
                                                                  // grain
          {
            for (jj = 0; jj < 8; jj++) {
              // have to find constant of particle(grain) so
              // have to find
              // spring that is not attached to particle of
              // other grain

              if (runParticle->neigP[j]->neigP[jj]) // find
                                                    // spring
                                                    // with
                                                    // right
                                                    // constant
              {
                // if the spring is connected to right
                // particle

                if (runParticle->neigP[j]->neigP[jj]->grain ==
                    runParticle->neigP[j]->grain) {
                  // should be strength of neighbour
                  // grain

                  grain_two_strength = runParticle->neigP[j]->springf[jj];
                  // apply average
                  // if the strength is zero (hole) the
                  // boundary springs
                  // should also be zero

                  if (grain_two_strength == 0.0) {
                    runParticle->springf[j] = 0.0;
                  } else if (grain_one_strength == 0.0) {
                    runParticle->springf[j] = 0.0;
                  } else {
                    runParticle->springf[j] =
                        (grain_one_strength + grain_two_strength) / 2.0;
                    runParticle->springs[j] =
                        (grain_one_strength + grain_two_strength) / 2.0;
                  }
                }
              }
            }
          }
        }
      }
    }
    runParticle = runParticle->nextP; // and loop
  }

  Relaxation(); // and relax
}

/***********************************************************************************
 * this routine weakens or strenghtens a grain. gets the number of the grain
 * we are talking about and the change of the spring constant in percent (1.0 is
 *zero change, smaller is weaker and more is harder) and then a number for
 *changing the breaking strength in percent (again 1.0 is nothing and more is
 *harder or breaks less easy and less is break more easy). normally called after
 *lattice is constructed by mike.elle but can also be called from within a
 *routine to have some time dependent weakening maybe ?
 *
 * daniel summer 2002
 *
 * add youngs modulus,
 *
 * Daniel, march 2003
 ***********************************************************************************/

void Lattice::WeakenGrain(int nb, float constant, float visc,
                          float break_strength) {
  int i, j; // counters

  runParticle = &refParticle; // start

  for (i = 0; i < numParticles; i++) // loop through particles
  {
    if (runParticle->grain == nb) // if the right grain
    {
      // --------------------------------------------------------
      // adjust youngs modulus of particle, default is for
      // spring constant of 1.0
      // --------------------------------------------------------

      runParticle->young = runParticle->young * constant;
      runParticle->E = runParticle->E * constant;
      runParticle->viscosity = runParticle->viscosity * visc;

      for (j = 0; j < 8; j++) // loop through neighbours
      {

        if (runParticle->neigP[j]) // if neighbour
        {
          runParticle->break_Str[j] =
              runParticle->break_Str[j] * break_strength; // change
        }
      }
    }
    runParticle = runParticle->nextP; // go on looping
  }

  AdjustParticleConstants();
}

/*******************************************************************************************
 * Routine to make a horizontal layer weak or hard. gets as input y minimum and
 *maximum position of layer (taken for whole grains) and a constant for the
 *elastic modulus and a constant for the spring. Changes in percent, 1.0 is no
 *change, > 1.0 stronger and < 1.0 softer.
 *
 * routine can be called several times.
 *
 *
 * daniel, December 2002
 *
 * adjust youngs modulus,
 *
 * daniel, March 2003
 **********************************************************************************************/

void Lattice::WeakenHorizontalLayer(double y_min, double y_max, float constant,
                                    float break_strength) {
  int i, j, jj; // counters
  bool isgrain; // is grain in list ?

  // first set the vector layer to -2 (have zero grain number)

  for (j = 0; j < 10000; j++) {
    layer[j] = -2;
  }

  runParticle = &refParticle; // start

  for (i = 0; i < numParticles; i++) // run through particles
  {
    if (runParticle->ypos > y_min) // larger than y min
    {
      if (runParticle->ypos < y_max) // smaller than y max ->
                                     // in layer
      {
        isgrain = false; // counter is false

        for (j = 0; j < 10000; j++) // run through the list
        {
          if (layer[j] == runParticle->grain) // if already
                                              // defined
          {
            isgrain = true; // grain is already in list
            break;
          } else if (layer[j] == -2) // end of list break,
                                     // grain not in list
          {
            break;
          }
        }
        if (!isgrain) // grain in layer but not yet in the list
        {
          layer[j] = runParticle->grain; // put it in the
          // list
        }
      }
    }
    runParticle = runParticle->nextP;
  }

  for (i = 0; i < numParticles; i++) // loop through particles
  {
    for (jj = 0; jj < 10000; jj++) // loop through list
    {
      if (layer[jj] == -2)
        break; // end of list exit

      else if (runParticle->grain == layer[jj]) // if the right
                                                // grain
      {
        // -----------------------------------------------------
        // adjust youngs modulus of particle
        // -----------------------------------------------------

        runParticle->young = runParticle->young * constant;

        for (j = 0; j < 8; j++) // loop through neighbours
        {
          // -------------------------------------------------------------------------
          // we do not want to change the grain boundaries as
          // well, only springs
          // within the grain, therefore check which springs are
          // between particles
          // of the same grain.
          //
          // call function adjustgrainboundaries later to get
          // mean constants for boundaries
          // --------------------------------------------------------------------------

          if (runParticle->neigP[j]) // if neighbour
          {
            if (runParticle->neigP[j]->grain == runParticle->grain) // if
                                                                    // the
                                                                    // same
                                                                    // grain
            {
              runParticle->springf[j] =
                  runParticle->springf[j] * constant; // change
              // constant
              runParticle->break_Str[j] =
                  runParticle->break_Str[j] * break_strength; // change
              // tensile
            }
          }
        }
      }
    }
    runParticle = runParticle->nextP; // go on looping
  }
}

//---------------------------------------------------------------------------------
// function that sets a distribution on the porosity for the fluid codes
// does not affect elasticity. Changes the rad par fluid radius of particles
// that is only used to calculate solid fraction for fluid codes
//
// daniel 2018
//--------------------------------------------------------------------------------

void Lattice::SetDistributionPorosity(float dis) {
  int i;
  float ran_nb, montecarlo;

  srand(time(0));

  for (i = 0; i < numParticles; i++) {
    montecarlo = rand() / (float)RAND_MAX;

    montecarlo = dis * montecarlo / 100.0;
    // cout << montecarlo << endl;

    runParticle->rad_par_fluid =
        runParticle->rad_par_fluid + runParticle->rad_par_fluid * montecarlo;

    runParticle->area_par_fluid =
        Pi * pow(runParticle->rad_par_fluid, 2.0); // area of particle radius

    runParticle = runParticle->nextP;
  }
}

/*******************************************************************************************
 * Routine to make a horizontal layer weak or hard (particles). gets as input y
 *minimum and maximum position of layer  and a constant for the elastic modulus
 *and a constant for the spring. Changes in percent, 1.0 is no change, > 1.0
 *stronger and < 1.0 softer.
 *
 * routine can be called several times.
 *
 *
 * daniel, January 2004
 *
 **********************************************************************************************/

void Lattice::WeakenHorizontalBox(double y_min, double y_max, double x_min,
                                  double x_max, float constant, float vis,
                                  float break_strength, float rad_frac) {
  int i, j; // counters

  runParticle = &refParticle;

  for (i = 0; i < numParticles; i++) // loop through particles
  {
    if (runParticle->ypos >= y_min) {
      if (runParticle->ypos <= y_max) {
        if (runParticle->xpos >= x_min) {
          if (runParticle->xpos <= x_max) {
            // -----------------------------------------------------
            // adjust youngs modulus of particle
            // -----------------------------------------------------

            runParticle->mineral = 0;

            runParticle->young = runParticle->young * constant;
            runParticle->E = runParticle->E * constant;
            runParticle->springconst = runParticle->springconst * constant;
            runParticle->viscosity *= vis;
            runParticle->rad_par_fluid *= rad_frac;
            runParticle->area_par_fluid *= pow(rad_frac, 2.0);

            for (j = 0; j < 8; j++) // loop through neighbours
            {
              // -------------------------------------------------------------------------
              // we do not want to change the grain boundaries as
              // well, only springs
              // within the grain, therefore check which springs are
              // between particles
              // of the same grain.
              //
              // call function adjustgrainboundaries later to get
              // mean constants for boundaries
              // --------------------------------------------------------------------------

              if (runParticle->neigP[j]) // if neighbour
              {
                if (runParticle->neigP[j]->ypos >= y_min) {
                  if (runParticle->neigP[j]->ypos <= y_max) {

                    runParticle->springf[j] =
                        runParticle->springf[j] * constant; // change
                    // constant
                    runParticle->springs[j] =
                        runParticle->springs[j] * constant;
                    runParticle->springv[j] = runParticle->springv[j] * vis;
                    runParticle->break_Str[j] =
                        runParticle->break_Str[j] * break_strength; // change
                    // tensile
                  }
                }
              }
            }
          }
        }
      }
    }

    // else
    //{
    // runParticle->young = runParticle->young * 0.5;

    // for (j = 0; j < 8; j++)	// loop through neighbours
    //{
    // if (runParticle->neigP[j])	// if neighbour
    //{
    // if (runParticle->neigP[j]->ypos < y_min)
    //{
    // runParticle->springf[j] = runParticle->springf[j] * 0.5;	// change
    //// constant
    //}
    //}
    //}
    //}
    runParticle = runParticle->nextP; // go on looping
  }

  AdjustParticleConstants();
}

void Lattice::WeakenHorizontalParticleLayer(double x_min, double x_max,
                                            double y_min, double y_max,
                                            float constant, float vis,
                                            float break_strength,
                                            float rad_frac) {
  int i, j; // counters

  runParticle = &refParticle;

  for (i = 0; i < numParticles; i++) // loop through particles
  {
    if (runParticle->ypos >= y_min) {
      if (runParticle->ypos <= y_max) {
        if (runParticle->xpos >= x_min) {
          if (runParticle->xpos <= x_max) {
            // -----------------------------------------------------
            // adjust youngs modulus of particle
            // -----------------------------------------------------

            runParticle->mineral = 0;

            runParticle->young = runParticle->young * constant;
            runParticle->E = runParticle->E * constant;
            runParticle->springconst = runParticle->springconst * constant;
            runParticle->viscosity *= vis;
            runParticle->rad_par_fluid *= rad_frac;
            runParticle->area_par_fluid *= pow(rad_frac, 2.0);

            for (j = 0; j < 8; j++) // loop through neighbours
            {
              // -------------------------------------------------------------------------
              // we do not want to change the grain boundaries as
              // well, only springs
              // within the grain, therefore check which springs are
              // between particles
              // of the same grain.
              //
              // call function adjustgrainboundaries later to get
              // mean constants for boundaries
              // --------------------------------------------------------------------------

              if (runParticle->neigP[j]) // if neighbour
              {
                if (runParticle->neigP[j]->ypos >= y_min) {
                  if (runParticle->neigP[j]->ypos <= y_max) {

                    runParticle->springf[j] =
                        runParticle->springf[j] * constant; // change
                    // constant
                    runParticle->springs[j] =
                        runParticle->springs[j] * constant;
                    runParticle->springv[j] = runParticle->springv[j] * vis;
                    runParticle->break_Str[j] =
                        runParticle->break_Str[j] * break_strength; // change
                    // tensile
                  }
                }
              }
            }
          }
        }
      }
    }

    // else
    //{
    // runParticle->young = runParticle->young * 0.5;

    // for (j = 0; j < 8; j++)	// loop through neighbours
    //{
    // if (runParticle->neigP[j])	// if neighbour
    //{
    // if (runParticle->neigP[j]->ypos < y_min)
    //{
    // runParticle->springf[j] = runParticle->springf[j] * 0.5;	// change
    //// constant
    //}
    //}
    //}
    //}
    runParticle = runParticle->nextP; // go on looping
  }

  AdjustParticleConstants();
}

// this was from Irfan

void Lattice::PoreDepnttBrkStr_Seal(float y_min, float y_max,
                                    double parRad_modl) {
  int i, j; // counters
  double rad_m, del_rad, Brk_prct;

  runParticle = &refParticle;

  for (i = 0; i < numParticles; i++) // loop through particles
  {
    if (runParticle->ypos >= y_min) {
      if (runParticle->ypos <= y_max) {
        for (j = 0; j < 8; j++) // loop through neighbours
        {
          if (runParticle->neigP[j]) // if neighbour
          {
            if (runParticle->neigP[j]->ypos >= y_min) {
              if (runParticle->neigP[j]->ypos <= y_max) {
                rad_m = 0;
                del_rad = 0;
                Brk_prct = 0;

                rad_m = (runParticle->rad_par_fluid +
                         runParticle->neigP[j]->rad_par_fluid) /
                        2.0;
                del_rad = (rad_m - parRad_modl) / parRad_modl;
                Brk_prct = runParticle->break_Str[j] * del_rad;
                runParticle->break_Str[j] =
                    runParticle->break_Str[j] - Brk_prct;

                // cout << runParticle->rad_par_fluid << ","
                // <<runParticle->neigP[j]->rad_par_fluid << ","
                //<< rad_m << "," << del_rad << "," << Brk_prct << "," <<
                //runParticle->break_Str[j] << " "; cout <<
                // runParticle->break_Str[j] << " "; cout << "i : " << i <<
                // endl; cout << " neigP[i]->nb : " << neigP[i]->nb << endl;
                // cout << " neigP[i]->rad_par_fluid : " <<
                // runParticle->neigP[j]->rad_par_fluid << " ";
              }
            }
          }
          // cout << endl;
        }
        // cout << endl<<endl;
      }
    }

    runParticle = runParticle->nextP; // go on looping
  }

  AdjustParticleConstants();

  // cout << endl << endl << endl;
  // runParticle = &refParticle;
  // for (i = 0; i < numParticles; i++)	// loop through particles
  //{
  // for (j = 0; j < 8; j++)	// loop through neighbours
  //{
  // if (runParticle->neigP[j])	// if neighbour
  //{
  // cout << runParticle->break_Str[j] << " ";
  //}
  //}
  // cout << endl<<endl;

  // runParticle = runParticle->nextP;
  //}
}

// this was fro irfan

void Lattice::WeakenHorizontal_Slot(double y_min, double y_max, double x_min,
                                    double x_max, float constant, float vis,
                                    float break_strength) {
  int i, j; // counters

  runParticle = &refParticle;

  for (i = 0; i < numParticles; i++) // loop through particles
  {
    if ((runParticle->ypos >= y_min) && (runParticle->ypos <= y_max)) {
      if (runParticle->xpos <= x_min || runParticle->xpos >= x_max)
      // if ((runParticle->xpos < x_min-0.25) || (runParticle->xpos > x_min-0.20
      // && runParticle->xpos < x_min) || (runParticle->xpos > x_max &&
      //runParticle->xpos < x_max+0.20) || (runParticle->xpos > x_max+0.25))
      {
        // -----------------------------------------------------
        // adjust youngs modulus of particle
        // -----------------------------------------------------

        runParticle->mineral = 0;

        runParticle->young = runParticle->young * constant;
        runParticle->viscosity *= vis;

        for (j = 0; j < 8; j++) // loop through neighbours
        {
          if (runParticle->neigP[j]) // if neighbour
          {
            if (runParticle->neigP[j]->ypos > y_min) {
              if (runParticle->neigP[j]->ypos < y_max) {
                runParticle->springf[j] =
                    runParticle->springf[j] * constant; // change
                // constant
                runParticle->springv[j] = runParticle->springv[j] * vis;
                runParticle->break_Str[j] =
                    runParticle->break_Str[j] * break_strength; // change
                // tensile
              }
            }
          }
        }
      }
    }
    runParticle = runParticle->nextP; // go on looping
  }

  AdjustParticleConstants();
}

//-----------------------------------------------------------------------
// horizontal layer that goes from xmin to xmax
//-----------------------------------------------------------------------

void Lattice::WeakenHorizontalParticleLayerX(double y_min, double y_max,
                                             double x_min, double x_max,
                                             float constant, float vis,
                                             float break_strength) {
  int i, j; // counters

  for (i = 0; i < numParticles; i++) // loop through particles
  {
    if (runParticle->ypos > y_min) {
      if (runParticle->ypos < y_max) {
        if (runParticle->xpos > x_min) {
          if (runParticle->xpos < x_max) {
            // -----------------------------------------------------
            // adjust youngs modulus of particle
            // -----------------------------------------------------

            runParticle->mineral = 0;

            runParticle->young = runParticle->young * constant;
            runParticle->viscosity *= vis;

            for (j = 0; j < 8; j++) // loop through neighbours
            {
              // -------------------------------------------------------------------------
              // we do not want to change the grain boundaries as
              // well, only springs
              // within the grain, therefore check which springs are
              // between particles
              // of the same grain.
              //
              // call function adjustgrainboundaries later to get
              // mean constants for boundaries
              // --------------------------------------------------------------------------

              if (runParticle->neigP[j]) // if neighbour
              {
                if (runParticle->neigP[j]->ypos > y_min) {
                  if (runParticle->neigP[j]->ypos < y_max) {
                    runParticle->springf[j] =
                        runParticle->springf[j] * constant; // change
                    // constant
                    runParticle->springv[j] = runParticle->springv[j] * vis;
                    runParticle->break_Str[j] =
                        runParticle->break_Str[j] * break_strength; // change
                    // tensile
                  }
                }
              }
            }
          }
        }
      }
    }
    runParticle = runParticle->nextP; // go on looping
  }

  AdjustParticleConstants();
}

/*******************************************************************************************
 * Routine to make a horizontal layer weak or hard (particles). gets as input y
 *minimum and maximum position of layer  and a constant for the elastic modulus
 *and a constant for the spring. Changes in percent, 1.0 is no change, > 1.0
 *stronger and < 1.0 softer.
 *
 * routine can be called several times.
 *
 *
 * daniel, October 2005
 *
 **********************************************************************************************/

void Lattice::WeakenTiltedParticleLayer(double y_min, double y_max,
                                        double shift, float constant, float vis,
                                        float break_strength, float rad_frac) {
  int i, j; // counters

  for (i = 0; i < numParticles; i++) // loop through particles
  {
    if (runParticle->ypos > (y_min + (shift * runParticle->xpos))) {
      if (runParticle->ypos < (y_max + (shift * runParticle->xpos))) {
        // -----------------------------------------------------
        // adjust youngs modulus of particle
        // -----------------------------------------------------

        runParticle->young = runParticle->young * constant;
        runParticle->viscosity *= vis;

        runParticle->rad_par_fluid *= rad_frac;
        runParticle->area_par_fluid *= pow(rad_frac, 2.0);

        for (j = 0; j < 8; j++) // loop through neighbours
        {
          // -------------------------------------------------------------------------
          // we do not want to change the grain boundaries as
          // well, only springs
          // within the grain, therefore check which springs are
          // between particles
          // of the same grain.
          //
          // call function adjustgrainboundaries later to get
          // mean constants for boundaries
          // --------------------------------------------------------------------------

          if (runParticle->neigP[j]) // if neighbour
          {
            if (runParticle->neigP[j]->ypos > y_min) {
              if (runParticle->neigP[j]->ypos < y_max) {
                runParticle->springf[j] =
                    runParticle->springf[j] * constant; // change
                // constant
                runParticle->springv[j] = runParticle->springv[j] * vis;
                runParticle->break_Str[j] =
                    runParticle->break_Str[j] * break_strength; // change
                // tensile
              }
            }
          }
        }
      }
    }
    runParticle = runParticle->nextP; // go on looping
  }

  AdjustParticleConstants();
}

void Lattice::WeakenTiltedParticleLayerX(double x_min, double x_max,
                                         double shift, float constant,
                                         float vis) {
  int i, j; // counters

  for (i = 0; i < numParticles; i++) // loop through particles
  {
    if (runParticle->xpos > (x_min + (shift * runParticle->ypos))) {
      if (runParticle->xpos < (x_max + (shift * runParticle->ypos))) {
        // -----------------------------------------------------
        // adjust youngs modulus of particle
        // -----------------------------------------------------

        runParticle->young = runParticle->young * constant;
        runParticle->viscosity *= vis;
      }
    }
    runParticle = runParticle->nextP; // go on looping
  }

  AdjustParticleConstants();
}

//---------------------------------------------------------------------------------
// routine that overprints existing porosity, rad/frac then has to have the
// right value for the lattice resolution. daniel 2018
//----------------------------------------------------------------------------------

void Lattice::WeakenTiltedParticleLayerabsolute(double y_min, double y_max,
                                                double shift, float constant,
                                                float vis, float break_strength,
                                                float rad_frac) {
  int i, j; // counters

  for (i = 0; i < numParticles; i++) // loop through particles
  {
    if (runParticle->xpos > (y_min + (shift * runParticle->ypos))) {
      if (runParticle->xpos < (y_max + (shift * runParticle->ypos))) {
        // -----------------------------------------------------
        // adjust youngs modulus of particle
        // -----------------------------------------------------

        runParticle->young = constant;
        runParticle->viscosity = vis;

        runParticle->rad_par_fluid = rad_frac;
        runParticle->area_par_fluid = 3.14159 * pow(rad_frac, 2.0);

        for (j = 0; j < 8; j++) // loop through neighbours
        {
          // -------------------------------------------------------------------------
          // we do not want to change the grain boundaries as
          // well, only springs
          // within the grain, therefore check which springs are
          // between particles
          // of the same grain.
          //
          // call function adjustgrainboundaries later to get
          // mean constants for boundaries
          // --------------------------------------------------------------------------

          if (runParticle->neigP[j]) // if neighbour
          {

            runParticle->break_Str[j] = break_strength; // change
            // tensile
          }
        }
      }
    }
    runParticle = runParticle->nextP; // go on looping
  }

  AdjustParticleConstants();
}

/*********************************************************************
 * function in lattice class used for multi layer boudins
 * function creates layers of two sorts. Input is the spring constant
 * of layer a and b, the viscosity of layer a and b, the breaking strength of
 *layer a and b and the thickness of layer a and b. layers are arranged parallel
 *to the x direction (horizontal) starting from the bottom. Elle grains are
 *ignored, only particles are taken
 *
 * arranged like : 	b
 *                       a
 *                       b
 *                       a
 *                       ....
 * last input (integer av) determines whether or not an average
 * of the two layer constants is taken (0 = no average, 1 = average)
 * Becomes important when layers are only as thick as a particle
 *
 * Daniel and Arzu 2004
 **********************************************************************/

void Lattice::MakeHorizontalLayers(float young_a, float viscous_a,
                                   float break_a, float thick_a, float young_b,
                                   float viscous_b, float break_b,
                                   float thick_b, int av) {
  int layernumber;
  int i, j, jj, k;
  float average;

  layernumber = 1.0 / (thick_a + thick_b); // determine how many
  // layers fit in the model

  cout << "layernumber " << layernumber << endl; // give out layer
  // number

  for (i = 0; i < numParticles; i++) // loop through all particles
  {
    // algorithm runs through the layers in sets of two, always a
    // sequence,
    // then checks first the sequence for the particle
    // and afterwards determines if its layer a or b in that sequence.

    for (j = 1; j < (layernumber + 1); j++) // loop through the layers
                                            // to find layer of
                                            // particle
    {
      if (runParticle->ypos < ((thick_a + thick_b) * j)) // in this
                                                         // layer
                                                         // sequence
                                                         // (defined
                                                         // by j)
      {
        if (runParticle->ypos < ((thick_a + thick_b) * (j - 1)) + thick_a) // is
                                                                           // it
                                                                           // a
                                                                           // or
                                                                           // b
                                                                           // ?
        {
          for (k = 0; k < 8; k++) // its layer a, give it a
                                  // constants !
          {

            runParticle->springv[k] = viscous_a; // viscosity
            // of
            // spring
            // (absolute)
            runParticle->young = young_a * 2.0 / sqrt(3.0); // young
            // modulus
            // (absolute)
            runParticle->E = young_a * 2.0 / sqrt(3.0); // young

            runParticle->springf[k] =
                sqrt(3.0) * runParticle->E /
                (3.0 *
                 (1.0 - runParticle->poisson_ratio)); // set spring constant
            runParticle->springs[k] = (1.0 - 3.0 * runParticle->poisson_ratio) /
                                      (1.0 + runParticle->poisson_ratio) *
                                      runParticle->springf[k];

            runParticle->break_Str[k] =
                runParticle->break_Str[k] * break_a; // breaking
            // strength
            // (factor)
          }
          break; // and exist, we have found the layer
        } else   // oh, this is layer b -> give it the b
               // constants !
        {
          for (k = 0; k < 8; k++) // loop through springs
          {

            runParticle->springv[k] = viscous_b; // viscosity
            // (absolute)
            runParticle->young = young_b * 2.0 / sqrt(3.0); // young
                                                            // modulus
                                                            // (absolute)
            runParticle->E = young_b * 2.0 / sqrt(3.0);     // young
            // modulus
            // (absolute)
            runParticle->springf[k] =
                sqrt(3.0) * runParticle->E /
                (3.0 *
                 (1.0 - runParticle->poisson_ratio)); // set spring constant
            runParticle->springs[k] = (1.0 - 3.0 * runParticle->poisson_ratio) /
                                      (1.0 + runParticle->poisson_ratio) *
                                      runParticle->springf[k];

            runParticle->break_Str[k] =
                runParticle->break_Str[k] * break_b; // breaking
            // strength
            // (factor)
          }
          break; // get out !
        }
      }
      if (runParticle->ypos > ((thick_a + thick_b) * layernumber)) // last
                                                                   // layer...
      {
        for (k = 0; k < 8; k++) {

          runParticle->springv[k] = viscous_a;
          runParticle->young = young_a;
          runParticle->E = young_a;

          runParticle->springf[k] =
              sqrt(3.0) * runParticle->E /
              (3.0 * (1.0 - runParticle->poisson_ratio)); // set spring constant
          runParticle->springs[k] = (1.0 - 3.0 * runParticle->poisson_ratio) /
                                    (1.0 + runParticle->poisson_ratio) *
                                    runParticle->springf[k];

          runParticle->break_Str[k] = runParticle->break_Str[k] * break_a;
        }
      }
    }
    runParticle = runParticle->nextP; // next particle please...
  }

  // now we have to adjust the spring constants, spring constants on
  // layer boundaries are now different
  // event though springs are the same. Now a spring has a constant b if
  // looking from layer b and one of
  // a if we look from layer a... this is unstable and physically
  // wrong...

  for (i = 0; i < numParticles; i++) // run along
  {

    // first find your own spring constant
    // so find a spring within the grain

    for (j = 0; j < 8; j++) // loop through springs
    {
      if (runParticle->neigP[j]) // if there is a spring
      {
        for (jj = 0; jj < 8; jj++) // loop second time
        {
          if (runParticle->neigP[j]->neigP[jj]) // find my
                                                // neighbour
                                                // (spring
                                                // from
                                                // both
                                                // sides..)
          {
            if (runParticle->neigP[j]->neigP[jj]->nb == runParticle->nb) {
              if (av == 1) // take the average of the
                           // two values
              {
                average = (runParticle->springf[j] +
                           runParticle->neigP[j]->springf[jj]) /
                          2.0;
                runParticle->springf[j] = average;
                runParticle->neigP[j]->springf[jj] = average;
                average = (runParticle->springs[j] +
                           runParticle->neigP[j]->springs[jj]) /
                          2.0;
                runParticle->springs[j] = average;
                runParticle->neigP[j]->springs[jj] = average;
                average = (runParticle->springv[j] +
                           runParticle->neigP[j]->springv[jj]) /
                          2.0;
                runParticle->springv[j] = average;
                runParticle->neigP[j]->springv[jj] = average;
              } else // no average, take always the
                     // lower constant
              {
                if (runParticle->springf[j] >
                    runParticle->neigP[j]->springf[jj]) {
                  runParticle->springf[j] = runParticle->neigP[j]->springf[jj];
                  runParticle->springs[j] = runParticle->neigP[j]->springs[jj];
                } else {
                  runParticle->neigP[j]->springf[jj] = runParticle->springf[j];
                  runParticle->neigP[j]->springs[jj] = runParticle->springs[j];
                }
                if (runParticle->springv[j] >
                    runParticle->neigP[j]->springv[jj]) {
                  runParticle->springv[j] = runParticle->neigP[j]->springv[jj];
                } else {
                  runParticle->neigP[j]->springv[jj] = runParticle->springv[j];
                }
              }
            }
          }
        }
      }
    }

    runParticle = runParticle->nextP; // and loop
  }
}

void Lattice::ChangeYoung(float change) {
  int i;

  for (i = 0; i < numParticles; i++) {
    runParticle->young = runParticle->young * change;
    runParticle->E = runParticle->E * change;
    runParticle->springconst = runParticle->springconst * change;

    runParticle = runParticle->nextP;
  }
}

void Lattice::DamageModel(float change, float stress, int damage) {
  int i;
  float yield_stress;

  for (i = 0; i < numParticles; i++) {
    yield_stress = runParticle->sxx - runParticle->syy;
    yield_stress = sqrt(yield_stress * yield_stress);

    if (yield_stress > stress && runParticle->damage < damage &&
        runParticle->ypos > 0.5) {
      cout << yield_stress << endl;
      runParticle->viscosity = runParticle->viscosity * change;
      runParticle->damage = runParticle->damage + 1;
      cout << "change viscosity" << endl;
    }

    runParticle = runParticle->nextP;
  }
}

// was used for large scale simulations

void Lattice::OutputTopParticles(double time) {
  ofstream file;
  file.open("Top_Particles_XY", ios::app);

  int i;
  double size;

  runParticle = firstParticle->prevP;

  for (i = 0; i < numParticles; i++) {
    if (runParticle->is_top_lattice_boundary && runParticle->ypos > 0.5) {
      size = runParticle->real_radius / runParticle->radius;

      file << runParticle->xpos * size / 1000.0 << " " << time << " "
           << runParticle->ypos * size / 1000.0 << endl;
    }
    runParticle = runParticle->nextP;
  }
  file << "\n";

  file.close();
}

/***********************************************************************
 * function changes elastic constant, viscosity and breaking strength of all
 * Particles.
 * all variables are factors for the change -> so that initial set
 * distributions are not changed, just shifted
 *
 * Daniel January 2004
 *************************************************************************/

void Lattice::WeakenAll(float constant, float viscosity, float break_strength) {
  int i, j; // loop counter

  for (i = 0; i < numParticles; i++) {
    runParticle->young = runParticle->young * constant;
    runParticle->real_young *= constant;
    runParticle->E = runParticle->E * constant;

    runParticle->springconst = runParticle->springconst * constant;

    runParticle->viscosity *= viscosity;

    for (j = 0; j < 8; j++) {
      if (runParticle->neigP[j]) {
        runParticle->springf[j] = runParticle->springf[j] * constant;
        runParticle->springs[j] = runParticle->springs[j] * constant;
        runParticle->springv[j] = runParticle->springv[j] * viscosity;
        runParticle->break_Str[j] = runParticle->break_Str[j] * break_strength;
      }
    }
    runParticle = runParticle->nextP;
  }
}

// debug

void Lattice::SetSinAnisotropy(float nb, float ampl) {
  int i;

  for (i = 0; i < numParticles; i++) {
    runParticle->young =
        runParticle->young * ((sin(runParticle->xpos * nb) + 1.1) * ampl);

    runParticle = runParticle->nextP;
  }

  AdjustParticleConstants();
}

// Cosinus function on the horizontal distribution of Sheet Youngs Modulus
// factor changes the width of the function
// factor b the height

void Lattice::SetCosViscSheetMod(float factor, float factorb) {
  int i;
  float x, y;

  for (i = 0; i < numParticles; i++) {
    if (!runParticle->is_lattice_boundary) {
      y = (runParticle->xpos - 0.5) * factorb + 0.5;
      x = cos(2 * Pi * y);
      x = x + 1.000001;

      x = x * factor;
      runParticle->sheet_young = runParticle->sheet_young * x;
    }
    runParticle = runParticle->nextP;
  }
}

/******************************************************************************************
 *  function that sets an initial anisotropy by setting a pseudorandom
 *distribution of grains that consist of one to three particles aligned parallel
 *to x. breaking strength, youngs modulus and viscosity of these particles can
 *be changed. length determines how long they are in x , at the moment only 1-3
 *particles ratio determines how many of these particles there are (percent from
 *0.0 to 1.0) these new particles will have no direct neighbours, there is
 *always supposed to be at least one other particle in between.
 *******************************************************************************************/

void Lattice::SetAnisotropy(float break_strength, float youngs_mod,
                            float viscosity, float ratio, float length) {
  int i, j, jj, nb;    // counter
  float ran_nb;        // pseudorandom number
  float str_nb, br_nb; // number for spring constant, number for
  // breaking strength
  float str_m, br_m; // helpers for the mean
  bool no_neighbour; // check if there is a matrix
  // particle in between

  runParticle = &refParticle; // and go

  ran_nb = rand() / (float)RAND_MAX; // pic pseudorandom from 0 to 1.0
  // for grains

  for (i = 0; i < numParticles; i++) // loop through particles
  {
    if (runParticle->xpos > 0.02) {
      if (runParticle->xpos < 0.98) {
        if (runParticle->ypos > 0.02) {
          if (runParticle->ypos < 0.98) // leave a matrix
                                        // boundary at the
                                        // rims of the
                                        // box.
          {
            ran_nb = rand() / (float)RAND_MAX; // pic
            // pseudorandom
            // from 0
            // to 1.0
            // for
            // grains
            no_neighbour = true; // assume I have only
            // matrix neighbours

            if (ran_nb < ratio / length) // rannumber used
                                         // to determine
                                         // propability for
                                         // particles to be
                                         // there
            {
              if (!runParticle->done) // if I have not
                                      // been converted
                                      // before
              {
                if (!runParticle->is_lattice_boundary) // if
                                                       // I
                                                       // am
                                                       // not
                                                       // at
                                                       // the
                                                       // lattice
                                                       // boundary
                                                       // (just
                                                       // in
                                                       // case)
                {
                  for (j = 0; j < 8; j++) // loop
                                          // through
                                          // springs
                  {
                    if (runParticle->neigP[j]) // if
                                               // I
                                               // have
                                               // a
                                               // neighbour
                    {
                      if (runParticle->neigP[j]->done)
                        no_neighbour = false; // if
                      // this
                      // neighbour
                      // is
                      // already
                      // transformed
                      // get
                      // out
                      // !
                    }
                  }
                  if (no_neighbour) // if I have no
                                    // neighbour that
                                    // is transformed
                  {
                    no_neighbour = true;                                  // ???
                    runParticle->young = runParticle->young * youngs_mod; // now
                    // apply
                    // the
                    // factors
                    runParticle->viscosity *= viscosity;

                    for (j = 0; j < 8; j++) // loop
                                            // through
                                            // neighbours
                    {
                      if (runParticle->neigP[j]) // if
                                                 // neighbour
                      {
                        runParticle->springf[j] =
                            runParticle->springf[j] * youngs_mod;
                        runParticle->break_Str[j] =
                            runParticle->break_Str[j] * break_strength;
                        runParticle->springv[j] =
                            runParticle->springv[j] * viscosity;
                      }
                    }
                    runParticle->done = true; // ok
                    // this
                    // one
                    // is
                    // done
                    if (length > 1) // if length is
                                    // larger than one
                                    // than try to
                                    // transform
                                    // neighbours in x
                                    // direction as
                                    // well
                    {
                      for (j = 0; j < 8; j++) // loop
                                              // through
                                              // neighbours
                      {
                        if (runParticle->neigP[j]) // if
                                                   // neighbour
                        {
                          if (runParticle->neigP[j]->ypos ==
                              runParticle->ypos) // if
                                                 // neighbour
                                                 // in
                                                 // same
                                                 // row
                          {
                            if (!runParticle->neigP[j]->done) // if
                                                              // not
                                                              // done
                            {
                              if (!runParticle->neigP[j]
                                       ->is_lattice_boundary) // if
                                                              // not
                                                              // lattice
                              {
                                for (jj = 0; jj < 8; jj++) // loop
                                                           // through
                                                           // neighbours
                                {
                                  if (runParticle->neigP[j]
                                          ->neigP[jj]) // if
                                                       // neighbour
                                  {
                                    if (runParticle->neigP[j]
                                            ->neigP[jj]
                                            ->done) // if
                                                    // that
                                                    // one
                                                    // is
                                                    // transformed
                                    {
                                      if (runParticle->nb !=
                                          runParticle->neigP[j]
                                              ->neigP[jj]
                                              ->nb) // if
                                        // its
                                        // not
                                        // myself
                                        no_neighbour = false; // dont
                                      // transform
                                    }
                                  }
                                }
                                if (no_neighbour) // if
                                                  // I
                                                  // can
                                                  // transform
                                {
                                  runParticle->neigP[j]->young =
                                      runParticle->neigP[j]->young * youngs_mod;
                                  runParticle->neigP[j]->done = true;
                                  for (jj = 0; jj < 8; jj++) // loop
                                                             // through
                                                             // neighbours
                                  {
                                    if (runParticle->neigP[j]
                                            ->neigP[jj]) // if
                                                         // neighbour
                                    {
                                      runParticle->neigP[j]->springf[jj] =
                                          runParticle->neigP[j]->springf[jj] *
                                          youngs_mod;
                                      runParticle->neigP[j]->break_Str[jj] =
                                          runParticle->neigP[j]->break_Str[jj] *
                                          break_strength;
                                      runParticle->neigP[j]->springv[jj] =
                                          runParticle->neigP[j]->springv[jj] *
                                          viscosity;
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                        if (length == 2) // if
                          // length
                          // is
                          // two
                          // stop
                          // otherwise
                          // try
                          // to
                          // do
                          // the
                          // second
                          // neibhour
                          // in
                          // x
                          // as
                          // well
                          break;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    runParticle = runParticle->nextP;
  }
  for (i = 0; i < numParticles; i++) // loop through particles to do
                                     // some averaging.
  {
    runParticle->done = false;
    for (j = 0; j < 8; j++) // loop through neighbours
    {
      if (runParticle->neigP[j]) // if neighbour
      {
        runParticle->springf[j] =
            (runParticle->neigP[j]->young + runParticle->young) / 2.0;
      }
    }

    runParticle = runParticle->nextP;
  }
}

/*******************************************************************+
 * Function that set a horizontal anisotropy. Litte horizontally
 * aligned mica flakes are generated with varying length (3-7) and
 * a randomly chosen Youngs modulus between 1 and max.
 *
 * Daniel 2005/6
 *********************************************************************/

void Lattice::SetAnisotropyRandom(float max, float boundary) {
  int size, i, j, ii;
  float ran_nb;

  runParticle = firstParticle;

  srand(std::time(0)); // pick random seed

  for (ii = 0; ii < numParticles;) {

    ran_nb = rand() / (float)RAND_MAX; // draw a number

    // distribute length of micas randomly bewteen 3 and 7

    if (ran_nb < 0.2)
      size = 3;
    else if (ran_nb < 0.4)
      size = 4;
    else if (ran_nb < 0.6)
      size = 5;
    else if (ran_nb < 0.8)
      size = 6;
    else
      size = 7;

    ran_nb = rand() / (float)RAND_MAX; // draw second number

    ran_nb = (ran_nb * (max - 1)) + 1; // base for youngs modulus variation

    for (i = 0; i < size; i++) {
      if (runParticle->is_lattice_boundary) {
        runParticle->young = runParticle->young * boundary;
        runParticle->viscosity *= boundary;

        for (j = 0; j < 8; j++) {
          if (runParticle->neigP[j]) {
            runParticle->springf[j] = runParticle->springf[j] * boundary;
            runParticle->springv[j] = runParticle->springv[j] * boundary;
          }
        }
      } else {
        runParticle->young = runParticle->young * ran_nb;
        runParticle->viscosity *= ran_nb;

        for (j = 0; j < 8; j++) {
          if (runParticle->neigP[j]) {
            runParticle->springf[j] = runParticle->springf[j] * ran_nb;
            runParticle->springv[j] = runParticle->springv[j] * ran_nb;
          }
        }
      }

      ii++;
      if (ii >= numParticles)
        break;

      runParticle = runParticle->nextP;
    }
  }

  AdjustParticleConstants();
}

/***********************************************************
 * function to set a max for pict dumps so that the
 * harddisk is not filled up completely
 *
 * Daniel March 2003
 ************************************************************/

void Lattice::Set_Max_Pict(int max) {
  set_max_pict = true;
  max_pict = max;
}

/***************************************************************************
 * This little thing can be called at the beginning after the lattice
 * was made (in mike.elle). You give it a number. It makes a plot of the
 * stuff each time that n fractures developed within a time step. Is
 * useful sometimes when a lot of fractures develop after a deformation
 * step and you want to see whats happening. Or you want to watch the
 * fracture develop. Default is 0 so if you dont call this function the
 * program will only make plots after time steps.
 *
 * daniel summer 2002
 **************************************************************************/

void Lattice::SetFracturePlot(int numbfrac, int local) {
  internal_break = local;     // default 0
  fract_Plot_Numb = numbfrac; // define the number (default is 0)
}

/********************************************************************************
 * Function that is used to dump a statistic file. This file gets the strain as
 * input from the main function (the strain that you put in the deformation
 *file) and an int for the grain that you want the statistic for. Dumps mean
 *values for stress tensor and mean and differential stress for this grain.
 *Called normally by the user from the mike.elle.cc file within the Run section
 *after deformation.
 *
 * file is called statisticgrain.txt
 *
 * dumps strain, mean stress smax + smin/2, sxx, syy, sxy and differential
 *stress smax - smin. Note that all compressive stresses are by definition
 *negative !
 *
 * can be directly importet by all kinds of programs including Matlab and Excel
 *
 * Daniel and Jochen Feb. 2003
 * thanks for some help from Jens
 **********************************************************************************/

void Lattice::DumpStatisticStressGrain(double strain, int nb) {
  FILE *stat;               // file pointer
  int i, p_counter;         // counters
  float mean, differential; // mean stress, differential
  // stress, pressure
  float m_sxx, m_syy, m_sxy;       // mean values of tensor
  float finite_strain;             // finite strain on the box
  float sxx, syy, sxy, smax, smin; // stress tensor and eigenvalues

  // zero out some values

  p_counter = 0;

  mean = 0.0;
  differential = 0.0;
  m_sxx = 0.0;
  m_syy = 0.0;
  m_sxy = 0.0;

  finite_strain = strain * (def_time_p + def_time_u); // calculate
  // finite strain

  runParticle = &refParticle; // start to loop through particles

  for (i = 0; i < numParticles; i++) // loop
  {
    if (runParticle->grain == nb) // if its the grain
    {
      p_counter = p_counter + 1; // counts number of particles

      // get the tensor from particles, just used to save space

      sxx = runParticle->sxx;
      syy = runParticle->syy;
      sxy = runParticle->sxy;

      // calculate eigenvalues

      smax = ((sxx + syy) / 2.0) +
             sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) + sxy * sxy);

      smin = ((sxx + syy) / 2.0) -
             sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) + sxy * sxy);

      // calculate mean stress

      mean = mean + (smax + smin) / 2.0;

      // calculate tensor means

      m_sxx = m_sxx + sxx;

      m_syy = m_syy + syy;

      m_sxy = m_sxy + sxy;

      // calculate differential stress

      differential = differential + (smax - smin);
    }
    runParticle = runParticle->nextP;
  }

  // divide by number of particles

  mean = mean / p_counter;

  m_sxx = m_sxx / p_counter;

  m_syy = m_syy / p_counter;

  m_sxy = m_sxy / p_counter;

  differential = differential / p_counter;

  // and dump the data for this step

  // append "a" means he is just adding values directly to the file
  // if the program crashes the file is still there with data in it
  // but if you start the program twice he will add data to the file
  // that is originally there so you have several steps of loading !!!!
  // just delete or move or rename old file and he makes a new one.

  stat = fopen("statisticgrain.txt", "a"); // open statistic output
  // append file

  // and write the data in one line with a semicolon ; between the
  // numbers
  // and jump to next line for the next deformation step

  fprintf(stat, "%d", nb);
  fprintf(stat, ";%f", finite_strain);
  fprintf(stat, ";%f", mean);
  fprintf(stat, ";%f", m_sxx);
  fprintf(stat, ";%f", m_syy);
  fprintf(stat, ";%f", m_sxy);
  fprintf(stat, ";%f\n", differential);

  fclose(stat); // close file
}

/*************************************************************************
 * ok, this is almost the same as the above code (to dump statistic for just
 * one grain). This is now for two grains. So the function gets the strain
 * from the deformation you set and two grain numbers. This is useful to
 * compare stresses in a hard layer and a soft layer. The function is an
 * append file again, so the program writes the data directly to the file
 * each time and the file is there even if the code crashes (thanks, Jens)
 * dumps number of grain,
 * strain on box, mean stress, differential stress and sxx, syy and sxy
 *
 * Jochen and Daniel Feb 2003
 *******************************************************************************/

void Lattice::DumpStatisticStressTwoGrains(double strain, int nb1, int nb2) {
  FILE *stat;               // file pointer
  int i, j, p_counter;      // counters
  float mean, differential; // mean stress, differential
  // stress, pressure
  float m_sxx, m_syy, m_sxy;       // mean stresses for tensor
  float finite_strain;             // finite strain on the box
  float sxx, syy, sxy, smax, smin; // stress tensor and eigenvalues

  // and zero values

  p_counter = 0;

  mean = 0.0;
  differential = 0.0;
  m_sxx = 0.0;
  m_syy = 0.0;
  m_sxy = 0.0;

  finite_strain = strain * (def_time_p + def_time_u); // calculate
  // finite strain

  runParticle = &refParticle;

  for (i = 0; i < numParticles; i++) // and loop
  {
    if (runParticle->grain == nb1) // if grain one
    {
      p_counter = p_counter + 1; // counter for particles

      // get tensor

      sxx = runParticle->sxx;
      syy = runParticle->syy;
      sxy = runParticle->sxy;

      // eigenvalues

      smax = ((sxx + syy) / 2.0) +
             sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) + sxy * sxy);

      smin = ((sxx + syy) / 2.0) -
             sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) + sxy * sxy);

      // mean stress

      mean = mean + (smax + smin) / 2.0;

      // mean of tensor

      m_sxx = m_sxx + sxx;

      m_syy = m_syy + syy;

      m_sxy = m_sxy + sxy;

      // differential stress

      differential = differential + (smax - smin);
    }
    runParticle = runParticle->nextP;
  }

  // and divide by number of particles

  mean = mean / p_counter;

  m_sxx = m_sxx / p_counter;

  m_syy = m_syy / p_counter;

  m_sxy = m_sxy / p_counter;

  differential = differential / p_counter;

  stat = fopen("statistictwograins.txt", "a"); // open statistic
  // output append
  // file

  // and dump data for first grain

  fprintf(stat, "%d", nb1);
  fprintf(stat, ";%f", finite_strain);
  fprintf(stat, ";%f", mean);
  fprintf(stat, ";%f", m_sxx);
  fprintf(stat, ";%f", m_syy);
  fprintf(stat, ";%f", m_sxy);
  fprintf(stat, ";%f", differential);

  // and do second grain
  // reset counters

  p_counter = 0;

  mean = 0.0;
  differential = 0.0;
  m_sxx = 0.0;
  m_syy = 0.0;
  m_sxy = 0.0;

  for (j = 0; j < numParticles; j++) // and loop again
  {
    if (runParticle->grain == nb2) // if second grain
    {
      p_counter = p_counter + 1; // count particles

      // stress tensor

      sxx = runParticle->sxx;
      syy = runParticle->syy;
      sxy = runParticle->sxy;

      // eigenvalues

      smax = ((sxx + syy) / 2.0) +
             sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) + sxy * sxy);

      smin = ((sxx + syy) / 2.0) -
             sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) + sxy * sxy);

      // mean stress

      mean = mean + (smax + smin) / 2.0;

      // mean stress tensor

      m_sxx = m_sxx + sxx;

      m_syy = m_syy + syy;

      m_sxy = m_sxy + sxy;

      // differential stress

      differential = differential + (smax - smin);
    }
    runParticle = runParticle->nextP;
  }

  // and divide by number of particles

  mean = mean / p_counter;

  m_sxx = m_sxx / p_counter;

  m_syy = m_syy / p_counter;

  m_sxy = m_sxy / p_counter;

  differential = differential / p_counter;

  // and put in the file

  fprintf(stat, "%d", nb2);
  fprintf(stat, ";%f", mean);
  fprintf(stat, ";%f", m_sxx);
  fprintf(stat, ";%f", m_syy);
  fprintf(stat, ";%f", m_sxy);
  fprintf(stat, ";%f\n", differential);

  fclose(stat); // and close
}

/*************************************************************************
 * function that plots the x and y coordinates of surface particles
 * and the strain. First row is strain than x and y coordinates
 * of one particle. Note that this is a scan along the particles so that
 * an additional function is needed to sort the data so that a
 * continuous function of the surface is attained.
 ************************************************************************/

void Lattice::DumpStatisticSurface(double strain) {
  FILE *stat;               // file pointer
  int i, j, jj, p_counter;  // counters
  float mean, differential; // mean stress, differential
  // stress, pressure
  float m_sxx, m_syy, m_sxy; // stress tensor
  float finite_strain;       // finite strain within the box of
  // interest
  float sxx, syy, sxy, smax, smin; // stress tensor and eigenvalues

  // set counters to zero

  finite_strain = strain * (def_time_u + def_time_p); // calculate
  // finite strain

  stat = fopen("statisticsurface.txt", "a"); // open statistic output
  // append file

  runParticle = &refParticle; // start

  for (i = 0; i < numParticles; i++) // look which particles are in the
                                     // box
  {
    if (runParticle->isHoleBoundary) // larger ymin
    {

      // and dump the data
      fprintf(stat, "%f", finite_strain);
      fprintf(stat, " %f", runParticle->xpos);
      fprintf(stat, " %f\n", runParticle->ypos);
    }
    runParticle = runParticle->nextP;
  }

  fclose(stat); // close file
}
void Lattice::Dumplayer(double strain) {
  FILE *stat;               // file pointer
  int i, j, jj, p_counter;  // counters
  float mean, differential; // mean stress, differential
  // stress, pressure
  float m_sxx, m_syy, m_sxy; // stress tensor
  float finite_strain;       // finite strain within the box of
  // interest
  float sxx, syy, sxy, smax, smin; // stress tensor and eigenvalues

  // set counters to zero

  finite_strain = strain * (def_time_u + def_time_p); // calculate
  // finite strain

  stat = fopen("layer.txt", "a"); // open statistic output
  // append file
  p_counter = 0;
  mean = 0;

  runParticle = &refParticle; // start

  for (i = 0; i < numParticles; i++) // look which particles are in the
                                     // box
  {
    if (runParticle->young > 0.0) // larger ymin
    {
      p_counter++;
      mean = runParticle->young + mean;
      // and dump the data
    }
    runParticle = runParticle->nextP;
  }

  mean = mean / p_counter;

  // fprintf (stat, "%f", finite_strain);
  // fprintf (stat, " %f", runParticle->xpos);
  fprintf(stat, " %f\n", mean);

  fclose(stat); // close file
}

void Lattice::Dump_concentration_growth(double strain) {
  FILE *stat;               // file pointer
  int i, j, jj, p_counter;  // counters
  float mean, differential; // mean stress, differential
  // stress, pressure
  float m_sxx, m_syy, m_sxy; // stress tensor
  float finite_strain;       // finite strain within the box of
  // interest
  float sxx, syy, sxy, smax, smin; // stress tensor and eigenvalues

  // set counters to zero

  finite_strain = strain * (def_time_u + def_time_p); // calculate
  // finite strain

  stat = fopen("growth.txt", "a"); // open statistic output
  // append file
  p_counter = 0;
  mean = 0;

  runParticle = &refParticle; // start

  for (i = 0; i < numParticles; i++) // look which particles are in the
                                     // box
  {
    if (runParticle->nb >= 20000 && runParticle->nb < 20200) // larger ymin
    {

      fprintf(stat, " %f", runParticle->conc);
      fprintf(stat, " %f\n", runParticle->calcite);
    }
    runParticle = runParticle->nextP;
  }

  // fprintf (stat, "%f", finite_strain);
  // fprintf (stat, " %f", runParticle->xpos);

  fclose(stat); // close file
}

/**************************************************************************************
 * function that dumps the y coordinate of surface particles within a region
 * xmin < x < xmax. end of line is the strain
 **************************************************************************************/

void Lattice::DumpTimeStatisticSurface(double strain, double xmin,
                                       double xmax) {
  FILE *stat;               // file pointer
  int i, j, jj, p_counter;  // counters
  float mean, differential; // mean stress, differential
  // stress, pressure
  float m_sxx, m_syy, m_sxy; // stress tensor
  float finite_strain;       // finite strain within the box of
  // interest
  float sxx, syy, sxy, smax, smin; // stress tensor and eigenvalues

  // set counters to zero

  finite_strain = strain * (def_time_u + def_time_p); // calculate
  // finite strain

  stat = fopen("timestatsurf.txt", "a"); // open statistic output
  // append file

  runParticle = &refParticle; // start

  for (i = 0; i < numParticles; i++) // look which particles are in the
                                     // box
  {
    if (runParticle->isHoleBoundary) // larger ymin
    {
      if (runParticle->xpos < xmax) {
        if (runParticle->xpos > xmin) {
          // and dump the data

          fprintf(stat, " %f", runParticle->ypos); // one row
          // of y
          // coordinates
        }
      }
    }
    runParticle = runParticle->nextP;
  }

  fprintf(stat, " %f\n", finite_strain); // give out strain

  fclose(stat); // close file
}

// SetGaussianStrengthDistribution
/*******************************************************************
 * puts a quenched noise on the breaking strength of the springs
 * in the box. This is a gaussian distribution of breaking strength
 * with a mean g_mean and a deviation of g_sigma
 *
 * The distribution is determined as follows: Pic a pseudorandom number (from
 *mean plus/minus 2 times sigma)  to give a particle a dissolution constant Then
 *determine the probability for that constant to appear from the Gauss function.
 *Then pic a second pseudo random number (from 0 to sigma) to determine whether
 *or not the grain will actually get this spring constant. If the random number
 *is below the probability the spring contant is accepted, if it is higher it is
 *rejected and a new one is drawn.
 *
 * Function writes also a text file with the set spring strengths distribution
 * calld gauss.txt
 *******************************************************************/

void Lattice::SetGaussianStrengthDistribution(double g_mean, double g_sigma) {
  FILE *stat; // for the outfile

  float prob;     // probability from gauss
  float k_spring; // pict spring
  int i, j;       // counters
  float ran_nb;   // rand number

  stat = fopen("gauss.txt", "a"); // open statistic output append
  // file

  // grain_counter counter how many grain there were intially
  // times two is for security, not all grains may be there
  // this can certainly be made nicer

  runParticle = &refParticle;

  srand(time(0));

  for (i = 0; i < numParticles; i++) // loop through the particles now
  {
    for (j = 0; j < 8; j++) {
      if (runParticle->neigP[j]) {
        do {
          k_spring = rand() / (float)RAND_MAX; // pic a
          // pseudorandom
          // float
          // between
          // 0 and 1

          // ----------------------------------------------------------------
          // now convert this distribution to a distribution
          // from
          // zero minus 8 times sigma to plus 8 times sigma
          // -----------------------------------------------------------------

          k_spring = (k_spring - 0.5) * 8.0 * (8.0 * g_sigma);

          // --------------------------------------------------------------
          // and shift it to mean plus minus 2 times sigma
          // --------------------------------------------------------------

          k_spring = k_spring + g_mean;

          // ---------------------------------------------------------------
          // now apply gauss function to determine a probability
          // for this
          // reaction constant to occur
          // ---------------------------------------------------------------

          prob = (1 / (g_sigma * sqrt(2.0 * Pi))); // part
          // one
          // (right
          // part)
          // of
          // function

          prob = prob * (exp(-0.5 * (((k_spring)-g_mean) / g_sigma) *
                             (((k_spring)-g_mean) / g_sigma))); // rest

          // --------------------------------------------------------------
          // now adjust probablity to run from 0 to 1.0
          // --------------------------------------------------------------

          prob = prob * sqrt(2.0 * Pi) * g_sigma;

          // ---------------------------------------------------------------
          // pic the second pseudorandom number
          // ---------------------------------------------------------------

          ran_nb = rand() / (float)RAND_MAX;

          // --------------------------------------------------------------
          // if the number picted is smaller or the same as the
          // probability
          // from the gauss function accept the rate constant,
          // if not go on
          // looping.
          // -------------------------------------------------------------------

          if (ran_nb <= prob) // if smaller
          {
            if (k_spring <= 0.0) // if not 0.0 or smaller
              // (dont want that)
              ;
            else {

              // scale to spring default strength

              k_spring = k_spring * 0.0125;

              // set the constant in the list !
              runParticle->break_Str[j] = k_spring;

              // and dump some data

              fprintf(stat, "%f", k_spring);
              fprintf(stat, "\n");

              // and break the while loop !

              break;
            }
          }
        } while (1); // loop until break
      }
    }

    runParticle = runParticle->nextP;
  }
}

// SetGaussianYoungDistribution
/*******************************************************************
 * puts a quenched noise on the breaking strength of the springs
 * in the box. This is a gaussian distribution of breaking strength
 * with a mean g_mean and a deviation of g_sigma
 *
 * The distribution is determined as follows: Pic a pseudorandom number (from
 *mean plus/minus 2 times sigma)  to give a particle a dissolution constant Then
 *determine the probability for that constant to appear from the Gauss function.
 *Then pic a second pseudo random number (from 0 to sigma) to determine whether
 *or not the grain will actually get this spring constant. If the random number
 *is below the probability the spring contant is accepted, if it is higher it is
 *rejected and a new one is drawn.
 *
 * Function writes also a text file with the set spring strengths distribution
 * calld gauss.txt
 *******************************************************************/

void Lattice::SetGaussianYoungDistribution(double g_mean, double g_sigma) {
  FILE *stat; // for the outfile

  float prob;     // probability from gauss
  float k_spring; // pict spring
  int i, j;       // counters
  float ran_nb;   // rand number

  stat = fopen("gauss.txt", "a"); // open statistic output append file

  // grain_counter counter how many grain there were intially
  // times two is for security, not all grains may be there
  // this can certainly be made nicer

  runParticle = &refParticle;

  srand(time(0));

  for (i = 0; i < numParticles; i++) // loop through the particles now
  {

    do {
      k_spring = rand() / (float)RAND_MAX; // pic a
      // pseudorandom float between 0 and 1

      // ----------------------------------------------------------------
      // now convert this distribution to a distribution
      // from zero minus 8 times sigma to plus 8 times sigma
      // -----------------------------------------------------------------

      k_spring = (k_spring - 0.5) * 8.0 * (8.0 * g_sigma);

      // --------------------------------------------------------------
      // and shift it to mean plus minus 2 times sigma
      // --------------------------------------------------------------

      k_spring = k_spring + g_mean;

      // ---------------------------------------------------------------
      // now apply gauss function to determine a probability
      // for this reaction constant to occur
      // ---------------------------------------------------------------

      prob = (1 / (g_sigma * sqrt(2.0 * 3.1415927))); // part
      // one (right part) of function

      prob = prob * (exp(-0.5 * (((k_spring)-g_mean) / g_sigma) *
                         (((k_spring)-g_mean) / g_sigma))); // rest

      // --------------------------------------------------------------
      // now adjust probablity to run from 0 to 1.0
      // --------------------------------------------------------------

      prob = prob * sqrt(2.0 * 3.1415927) * g_sigma;

      // ---------------------------------------------------------------
      // pic the second pseudorandom number
      // ---------------------------------------------------------------

      ran_nb = rand() / (float)RAND_MAX;

      // --------------------------------------------------------------
      // if the number picted is smaller or the same as the
      // probability from the gauss function accept the rate constant,
      // if not go on looping.
      // -------------------------------------------------------------------

      if (ran_nb <= prob) // if smaller
      {
        if (k_spring <= 0.0) // if not 0.0 or smaller
          // (dont want that)
          ;
        else {
          // scale to spring default strength
          k_spring = k_spring * 1;

          // set the constant in the list !
          runParticle->young = k_spring;

          // and dump some data
          fprintf(stat, "%f", k_spring);
          fprintf(stat, "\n");

          // and break the while loop !
          break;
        }
      }
    } while (1); // loop until break

    runParticle = runParticle->nextP;
  }
  AdjustParticleConstants();
}

/***************************************************************************
 * This function generates a gaussian distribution of springconstants k_gauss
 *for a a given grain population. This means that the number of grains with a
 *certain function sprinconstant f(k) within a desired intervall of values
 *follows a density function of the form f(k)=(sqrt(2*pi)*g_sigma)**(-1)
 *exp((k_gauss-g_mean)/g_sigma)**2 where g_sigma is the standard variation,
 *g_mean is the mean value. Remember that g_mean determines the centre of the
 *distribution and sigma the sharpness and height of the distribution.
 *
 * the function receives the mean and the sigma for the distribution.
 * AS mentioned above the mean will be the mean of the spring constants of the
 * grains so that a mean statistic of the whole box should give this constant as
 * long as the system behaves linear elastic. The distibution is set for grains
 * so that all particles in a grain will have the same spring constants (i.e.
 *their springs will). mean will give absolute spring constants so that this
 *function should be called first before anything else is set, otherwise it will
 *be wiped out. No problem however in setting a gauss distribution and
 *afterwards making a layer hard or soft. Default spring constant in the model
 *is 1.0 so that a mean of 1.0 will produce values clustering around 1.0 and 2.0
 *values with a mean spring constant of 2.0 etc.. We determine distribution
 *within the interval mean plus minus two times sigma !
 *
 * The distribution is determined as follows: Pic a pseudorandom number (from
 *mean plus/minus 2 times sigma)  to give a grain a spring constant. Then
 *determine the probability for that constant to appear from the Gauss function.
 *Then pic a second pseudo random number (from 0 to sigma) to determine whether
 *or not the grain will actually get this spring constant. If the random number
 *is below the probability the spring contant is accepted, if it is higher it is
 *rejected and a new one is drawn.
 *
 * Function writes also a text file with the set spring constant distribution
 * calld gauss.txt
 *
 * Daniel and Jochen Feb. 2003
 *
 * add adjustment to youngs modulus of particles,
 *
 * Daniel, March 2003
 **************************************************************************/

void Lattice::SetGaussianSpringDistribution(double g_mean, double g_sigma) {
  FILE *stat;               // for the outfile
  float spring_list[10000]; // distribution list for grains
  float prob;               // probability from gauss
  float k_spring;           // pict spring
  int i, j;                 // counters
  float ran_nb;             // rand number

  stat = fopen("gauss.txt", "a"); // open statistic output append
  // file

  // grain_counter counter how many grain there were intially
  // times two is for security, not all grains may be there
  // this can certainly be made nicer

  for (i = 0; i < (grain_counter * 2); i++) // loop through grains
  {
    do {
      k_spring = rand() / (float)RAND_MAX; // pic a
      // pseudorandom
      // float between 0
      // and 1

      // ----------------------------------------------------------------
      // now convert this distribution to a distribution from
      // zero minus 2 times sigma to plus 2 times sigma
      // -----------------------------------------------------------------

      k_spring = (k_spring - 0.5) * 2.0 * (2.0 * g_sigma);

      // --------------------------------------------------------------
      // and shift it to mean plus minus 2 times sigma
      // --------------------------------------------------------------

      k_spring = k_spring + g_mean;

      // ---------------------------------------------------------------
      // now apply gauss function to determine a probability for
      // this
      // spring constant to occur
      // ---------------------------------------------------------------

      prob = (1 / (g_sigma * sqrt(2.0 * 3.1415927))); // part
      // one
      // (right
      // part)
      // of
      // function

      prob = prob * (exp(-0.5 * (((k_spring)-g_mean) / g_sigma) *
                         (((k_spring)-g_mean) / g_sigma))); // rest

      // --------------------------------------------------------------
      // now adjust probablity to run from 0 to 1.0
      // --------------------------------------------------------------

      prob = prob * sqrt(2.0 * 3.1415927) * g_sigma;

      // ---------------------------------------------------------------
      // pic the second pseudorandom number
      // ---------------------------------------------------------------

      ran_nb = rand() / (float)RAND_MAX;

      // --------------------------------------------------------------
      // if the number picted is smaller or the same as the
      // probability
      // from the gauss function accept the spring constant, if not
      // go on
      // looping.
      // -------------------------------------------------------------------

      if (ran_nb <= prob) // if smaller
      {
        if (k_spring <= 0.0) // if not 0.0 or smaller (dont
          // want that)
          ;
        else {
          // set the constant in the list !

          spring_list[i] = k_spring;
          grain_young[i] = k_spring;

          // and dump some data

          fprintf(stat, "%f", k_spring);
          fprintf(stat, "\n");

          // and break the while loop !

          break;
        }
      }
    } while (1); // loop until break
  }

  runParticle = &refParticle;

  for (i = 0; i < numParticles; i++) // loop through the particles now
  {
    for (j = 0; j < 8; j++) // loop through neighbours
    {

      if (runParticle->neigP[j]) // if neighbour
      {

        // ------------------------------------------------------------------------
        // we use the place in the spring_list as number of the
        // grain
        // then we set all particles in the grain to that spring
        // constant
        // including all springs of particles.
        // ------------------------------------------------------------------------

        runParticle->springf[j] = spring_list[runParticle->grain];

        // ---------------------------------------------------
        // and set the new youngs modulus
        // ---------------------------------------------------

        runParticle->young = runParticle->springf[j] * 2.0 / sqrt(3.0);
      }
    }
    runParticle = runParticle->nextP;
  }

  fclose(stat); // and close file

  AdjustConstantGrainBoundaries(); // clean up (make mean of grain
  // boundary springs)
}

/************************************************************************************************
 * Another file to dump statistics in a text file. This is for a whole box. It
 *gets the input values box position with is box min in y box max in y, box min
 *in x and box max in x and also the strain applied on the whole box (the elle
 *box now) from your deformation file. The function takes a mean of all
 *particles in the box. It calculates means for the mean stress, the
 *differential stress and the stress tensor sxx, syy, sxy.
 *
 * output is finite strain on box, mean stress, sxx,syy,sxy and differential
 *stress. Function can be called from whereever you want it.
 *
 * daniel and jochen Feb. 2003
 *****************************************************************************************************/

void Lattice::DumpStatisticStressBox(double y_box_min, double y_box_max,
                                     double x_box_min, double x_box_max,
                                     double strain) {
  FILE *stat;       // file pointer
  FILE *scaledstat; // another file pointer where the stresses are scaled by the
                    // pressure_scale variable
  int i, j, jj, p_counter;  // counters
  float mean, differential; // mean stress, differential
  // stress, pressure
  float m_sxx, m_syy, m_sxy; // stress tensor
  float finite_strain;       // finite strain within the box of
  // interest
  float sxx, syy, sxy, smax, smin; // stress tensor and eigenvalues

  float wall_pos, wall_Xpos;

  runParticle = &refParticle;       // start of list
  runParticle = runParticle->prevP; // go back one
  wall_pos = runParticle->ypos;     // get wall position y
  wall_Xpos = runParticle->xpos;    // get wall position x

  for (i = 0; i < particlex; i++) {
    runParticle = runParticle->prevP;
  }

  if (runParticle->xpos > wall_Xpos) {
    wall_Xpos = runParticle->xpos;
  }

  p_counter = 0;

  mean = 0.0;
  differential = 0.0;
  m_sxx = 0.0;
  m_syy = 0.0;
  m_sxy = 0.0;

  finite_strain = strain * (def_time_u + def_time_p); // calculate
  // finite_strain = def_time_p;

  runParticle = &refParticle; // start

  for (i = 0; i < numParticles; i++) // look which particles are in the
                                     // box
  {
    if (runParticle->ypos > y_box_min * wall_pos) // larger ymin
    {
      if (runParticle->ypos < y_box_max * wall_pos) // smaller ymax
      {
        if (runParticle->xpos > x_box_min * wall_Xpos) // larger xmin
        {
          if (runParticle->xpos < x_box_max * wall_Xpos) // smaller xmax
          {
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

            // values added for particles

            m_sxx = m_sxx + sxx;

            m_syy = m_syy + syy;

            m_sxy = m_sxy + sxy;

            mean = mean + (smax + smin) / 2.0;

            differential = differential + (smax - smin);
          }
        }
      }
    }
    runParticle = runParticle->nextP;
  }

  // and divide by number of particles

  mean = mean / p_counter;

  differential = differential / p_counter;

  m_sxx = m_sxx / p_counter;

  m_syy = m_syy / p_counter;

  m_sxy = m_sxy / p_counter;

  stat = fopen("statisticbox.txt", "a"); // open statistic output
  // append file

  // and dump the data

  fprintf(stat, "%f", finite_strain);
  fprintf(stat, ";%f", mean);
  fprintf(stat, ";%f", m_sxx);
  fprintf(stat, ";%f", m_syy);
  fprintf(stat, ";%f", m_sxy);
  fprintf(stat, ";%f\n", differential);

  fclose(stat); // close file

  scaledstat = fopen("statisticbox-scaled.txt", "a"); // open statistic output
  // append file

  // and dump the data

  fprintf(scaledstat, "%f", finite_strain);
  fprintf(scaledstat, ";%f", mean * pressure_scale);
  fprintf(scaledstat, ";%f", m_sxx * pressure_scale);
  fprintf(scaledstat, ";%f", m_syy * pressure_scale);
  fprintf(scaledstat, ";%f", m_sxy * pressure_scale);
  fprintf(scaledstat, ";%f\n", differential * pressure_scale);

  fclose(scaledstat); // close file
}

void Lattice::DumpStatisticGravityStress(int particle_row) {
  FILE *stat; // file pointer

  int i, j; // counters

  runParticle = &refParticle; // start of list

  for (i = 0; i < particle_row; i++) {
    runParticle = runParticle->nextP;
  }

  stat = fopen("Gravitystatistic.txt", "a"); // open statistic output

  for (j = 0; j < particlex / 2; j++) {
    fprintf(stat, "%i", runParticle->nb);
    fprintf(stat, ";%f", runParticle->sxx);
    fprintf(stat, ";%f", runParticle->syy);
    fprintf(stat, ";%f\n", runParticle->sxy);

    for (i = 0; i < particlex * 2; i++) {
      runParticle = runParticle->nextP;
    }
  }
  fclose(stat); // close file
}

/*************************************************************************************************
 * Function that sets a pseudorandom distribution on spring constants and/or on
 *breaking- strength of springs. Spring constants are set for grains, breaking
 *strength for the whole model. First number is for the mean spring constant.
 *This number is multiplied by the set values (default 1.0) and the result
 *substracted or added so that 0.0 does not change the mean. Therefore
 *anisotropies or distributions in the model will not be distroyed but
 *redistributed depending on what they are. The second number gives the
 *distribution size. The distribution is original mean times str_mean plus/minus
 *half of str_size. The same applies for the breaking strength of springs. if
 *second or fourth numbers are 0.0 no distributions are set for these values
 *
 * Daniel and Jochen, Feb. 2003
 *
 * adjust youngs modulus of particle, bug in mean, add helpers
 *
 * Daniel, March 2003
 *************************************************************************************************/

void Lattice::SetPhase(float str_mean, float str_size, float br_mean,
                       float br_size) {
  int i, j, nb;        // counter
  float ran_nb;        // pseudorandom number
  float str_nb, br_nb; // number for spring constant, number for
  // breaking strength
  float str_m, br_m; // helpers for the mean

  srand(time(0));

  runParticle = &refParticle; // and go

  for (nb = 0; nb < (grain_counter * 2); nb++) // loop through
                                               // grains
  {
    if (str_size != 0.0) // in case a distribution is wanted
    {
      ran_nb = rand() / (float)RAND_MAX; // pic pseudorandom from 0
      // to 1.0 for grains
      str_nb = ran_nb * str_size; // now from 0 to distribution size
    }

    for (i = 0; i < numParticles; i++) // loop through particles
    {
      if (runParticle->grain == nb) // if particle of grain
      {
        if (str_size != 0.0) // if distribution
        {
          if (!runParticle->is_lattice_boundary) {
            for (j = 0; j < 8; j++) // loop through neighbours
            {

              if (runParticle->neigP[j]) // if neighbour
              {
                // ------------------------------------------------
                // multiply mean by springconstant
                // ------------------------------------------------

                str_m = runParticle->springf[j] * str_mean;

                // --------------------------------------------------
                // apply random distribution to spring
                // --------------------------------------------------

                runParticle->springf[j] =
                    runParticle->springf[j] * (1.0 + str_nb - (str_size / 2.0));

                // ---------------------------------------------------
                // and shift by the str_mean
                // ---------------------------------------------------

                runParticle->springf[j] = runParticle->springf[j] + str_m;
              }
            }
            // ------------------------------------------
            // and set youngs modulus of particle
            // ------------------------------------------

            str_m = runParticle->young * str_mean;
            runParticle->young =
                runParticle->young * (1.0 + str_nb - (str_size / 2.0));
            runParticle->young = runParticle->young + str_m;
          }
        }
        // --------------------------------------------------------
        // now do the same for all springs
        // --------------------------------------------------------

        if (br_size != 0.0) // if distribution is wanted
        {
          // cout << runParticle->rad_par_fluid << " ";

          for (j = 0; j < 8; j++) // loop through neighbours
          {

            if (runParticle->neigP[j]) // if neighbour
            {
              // cout << runParticle->neigP[j]->rad_par_fluid << " ";

              // ------------------------------------------------------
              // we pic a different pseudorandom number for
              // each spring
              // --------------------------------------------------------

              ran_nb = rand() / (float)RAND_MAX; // pic
              // number between 0.0 and 1.0

              br_nb = ran_nb * br_size; // now from 0.0 to br_size

              // br_m = runParticle->break_Str[j] * br_mean;	// set br_mean

              // --------------------------------------------------------
              // and apply distribution to the spring
              // breaking strength
              // -------------------------------------------------------

              runParticle->break_Str[j] =
                  runParticle->break_Str[j] +
                  (runParticle->break_Str[j] * (br_nb - (br_size / 2.0)));

              // ----------------------------------------------------------
              // and shift the distribution
              // ---------------------------------------------------------

              runParticle->break_Str[j] = runParticle->break_Str[j] * br_mean;

              //~ cout << runParticle->break_Str[j] << " " ;
            }
          }
        }
      }

      runParticle = runParticle->nextP;
    }
  }
  //~ cout << endl;
  AdjustConstantGrainBoundaries();
}

void Lattice::SetVariationBreakingThreshold(float var) {
  int i, j;
  float ran_nb;

  srand(time(0));

  for (i = 0; i < numParticles; i++) {
    if (!runParticle->is_lattice_boundary) {
      for (j = 0; j < 8; j++) // loop through neighbours
      {
        if (runParticle->neigP[j]) {
          ran_nb = rand() / (float)RAND_MAX;
          ran_nb = ((ran_nb * 2.0) - 1.0) * var;

          runParticle->break_Str[j] =
              runParticle->break_Str[j] + (runParticle->break_Str[j] * ran_nb);
          // cout << runParticle->break_Str[j] << endl;
        }
      }
    }

    runParticle = runParticle->nextP;
  }

  // AdjustConstantGrainBoundaries ();
}

void Lattice::SetViscous(float max) {
  int i;

  for (i = 0; i < numParticles; i++) {
    if (runParticle->ypos < max)
      runParticle->novisc = true;
    else
      runParticle->novisc = false;
    runParticle = runParticle->nextP;
  }
}

/**********************************************************************************
 * viscous retardation after the elastic deformation step.
 * the equilibrium lengths of springs is adjusted following an exponential
 * stress retardation depending on the time step and the youngs modulus to
 * viscosity ratio.
 * input is an integer that determines whether or not a picture is taken
 * (0 is no picture, 1 is take a picture).
 * the second variable is the time step for the viscous relaxation.
 * It does determine how fast the viscous step relaxes stresses before the next
 * elastic deformation step.
 ***********************************************************************************/

/*
 * kleine fehlerliste/problem-list/komments: in
 *
 * bool Particle::Relax(Particle **list,int size)
 *
 * wrd noch immer der radius ( statt rad[j] fuer die sxx/syy-berrechnung
 * benutzt, in 2. Relax-funktion sicher auch.
 *
 * * runParticle->springv ==> viscosity of spring reales verhaeltnis der
 * springs statt 2.0? formel f. ldl: haeh?
 *
 */

void Lattice::ViscousRelax(int dump, float timestep) {
  int i, j, ii, jj, k, nbts;
  double sigma;
  double angle, n_force; // repulsion: relevant angle with regard to str.tensor
                         // of neig, n_force normal force in that dir.
  double phi, sigma_angle, tau, mean_stress, merker1, merker2, r, sigma_new, G;

  double x1, y1, x2, y2;

  double lowest_ay, lowest_bx;

  double buffer, sxx, syy, sxy, sigma1, sigma2, x, y;

  double eq_springs[6][2];

  // new variables for new function:
  double c11, c12, c21, c22, K, k1, k2, dex, dey, dgamma, a11, a12, a21, a22,
      exx, eyy, exy, R_t, zeta, delta, winkel, xdiff, ydiff, v11, v12, v21, v22;

  double ang, viscosity = 1e22; // for angular calculations, buffer

  double debugf1 = 0.0, debugf2 = 1.0;
  int debugi1, debugi2, debugi3, debugi4;

  visc_rel = 1;

  cout << "in viscous relax" << endl;

  for (i = 0; i < numParticles; i++) // loop through particles
  {
    viscosity = runParticle->viscosity;

    if (!runParticle->novisc && !runParticle->is_bottom_lattice_boundary &&
        !runParticle->is_left_lattice_boundary &&
        !runParticle->is_right_lattice_boundary) // only olivine
    {

      /*
         stresses
         */

      mean_stress = (runParticle->sxx + runParticle->syy) / 2.0;

      sxx = (runParticle->sxx);
      syy = (runParticle->syy);
      sxy = (runParticle->sxy);

      // normal stresses
      sigma_angle = 0.5 * atan((2.0 * sxy) / (sxx - sxy));

      sigma1 = 0.5 * (sxx + syy) + 0.5 * (sxx - syy) * cos(2.0 * sigma_angle) +
               sxy * sin(2.0 * sigma_angle);
      sigma2 = 0.5 * (sxx + syy) +
               0.5 * (sxx - syy) * cos(2.0 * (sigma_angle + Pi / 2.0)) +
               sxy * sin(2.0 * (sigma_angle + Pi / 2.0));

      if (sigma2 > sigma1) {
        sigma_angle += Pi / 2.0;
        buffer = sigma1;
        sigma1 = sigma2;
        sigma2 = buffer;
      }

      sxx -= mean_stress;
      syy -= mean_stress;
      sigma1 -= mean_stress;
      sigma2 -= mean_stress;

      sxx *= pressure_scale * pascal_scale;
      syy *= pressure_scale * pascal_scale;
      sigma1 *= pressure_scale * pascal_scale;
      sigma2 *= pressure_scale * pascal_scale;

      /*
         tensor calculations and settings.
         compare the ramberg paper
         */

      zeta =
          (sigma_angle + (Pi / 4.0)); // zeta == angle of shear-force (korrekt?)

      tau = 0.5 * (sigma1 - sigma2);    // maximum shear stress
      dgamma = tau / (2.0 * viscosity); // now deformation-rate

      dex = sxx / (2.0 * viscosity); // negativ?????
      dey = syy / (2.0 * viscosity);

      // rate-of-displacement tensor (==deformation-tensor)
      a11 = dex - dgamma * sin(zeta) * cos(zeta);
      a12 = dgamma * cos(zeta) * cos(zeta);
      a21 = -dgamma * sin(zeta) * sin(zeta);
      a22 = dey + dgamma * sin(zeta) * cos(zeta);

      // von weijermars. für den deviatorischen stresstensor gilt: xx=-yy
      a11 = dex;
      a12 = sxy / (2.0 * viscosity);
      a21 = sxy / (2.0 * viscosity);
      a22 = -dex;
      // a22 = dey;

      // eigenvalues of r-o-d-tensor
      k1 = 0.5 * sqrt(dex * dex - 2.0 * dex * dgamma * sin(zeta) * cos(zeta));
      k2 = -0.5 * sqrt(dex * dex - 2.0 * dex * dgamma * sin(zeta) * cos(zeta));
      k1 = 0.5 * sqrt((a11 - a22) * (a11 - a22) + 4.0 * a12 * a21);
      k2 = -0.5 * sqrt((a11 - a22) * (a11 - a22) + 4.0 * a12 * a21);

      // move the points stored in xy[6][2]
      // these are the EQUILIBRIUM end-points of the springs
      for (j = 0; j < 6; j++) {

        // get eq_angle of spring to pos. x-axis
        // find springs

        // for every spring you will find the value of the spring-endpoint (x,y)
        // in xy[6][2]
        x = runParticle->neigpos[j][0] * runParticle->real_radius;
        y = runParticle->neigpos[j][1] * runParticle->real_radius;

        // coefficients of particle-path tensor
        c11 = ((k2 - a11) * x - a12 * y) / (k2 - k1);
        c12 = ((a11 - k1) * x + a12 * y) / (k2 - k1);
        c21 = ((k1 - a11) * (k2 - a11) * (1.0 / a12) * x - (k1 - a11) * y) /
              (k2 - k1);
        c22 = ((a11 - k1) * (k2 - a11) * (1.0 / a12) * x + (k2 - a11) * y) /
              (k2 - k1);

        runParticle->neigpos[j][0] =
            (c11 * exp(k1 * timestep) + c12 * exp(k2 * timestep)) /
            runParticle->real_radius;
        runParticle->neigpos[j][1] =
            (c21 * exp(k1 * timestep) + c22 * exp(k2 * timestep)) /
            runParticle->real_radius;
      }

      // now spring-length, as stored in rad[6]
      for (j = 0; j < 6; j++) {
        runParticle->rad[j] =
            sqrt((runParticle->xy[j][0] * runParticle->xy[j][0]) +
                 (runParticle->xy[j][1] * runParticle->xy[j][1]));

        // scale
        runParticle->rad[j] *= runParticle->radius;
        if (runParticle->nb == 5050 || runParticle->nb == 0 ||
            runParticle->nb == 99 || runParticle->nb == 10000 ||
            runParticle->nb == 10099)
          cout << runParticle->nb << " rad[" << j
               << "]: " << runParticle->rad[j] << endl;
      }

      // and finally the conjugate radii for the hull of the particles
      // (repulsion) this is the same idea as behind the spring-xy calculaitons
      x = runParticle->ax * runParticle->real_radius;
      y = runParticle->ay * runParticle->real_radius;

      // coefficients of particle-path tensor
      c11 = ((k2 - a11) * x - a12 * y) / (k2 - k1);
      c12 = ((a11 - k1) * x + a12 * y) / (k2 - k1);
      c21 = ((k1 - a11) * (k2 - a11) * (1.0 / a12) * x - (k1 - a11) * y) /
            (k2 - k1);
      c22 = ((a11 - k1) * (k2 - a11) * (1.0 / a12) * x + (k2 - a11) * y) /
            (k2 - k1);

      runParticle->ax = c11 * exp(k1 * timestep) + c12 * exp(k2 * timestep);
      runParticle->ay = c21 * exp(k1 * timestep) + c22 * exp(k2 * timestep);

      x = runParticle->bx * runParticle->real_radius;
      y = runParticle->by * runParticle->real_radius;

      // coefficients of particle-path tensor
      c11 = ((k2 - a11) * x - a12 * y) / (k2 - k1);
      c12 = ((a11 - k1) * x + a12 * y) / (k2 - k1);
      c21 = ((k1 - a11) * (k2 - a11) * (1.0 / a12) * x - (k1 - a11) * y) /
            (k2 - k1);
      c22 = ((a11 - k1) * (k2 - a11) * (1.0 / a12) * x + (k2 - a11) * y) /
            (k2 - k1);

      runParticle->bx = c11 * exp(k1 * timestep) +
                        c12 * exp(k2 * timestep) / runParticle->real_radius;
      runParticle->by = c21 * exp(k1 * timestep) +
                        c22 * exp(k2 * timestep) / runParticle->real_radius;
    }
    runParticle = runParticle->nextP;
  }

  cout << "out viscous relax" << endl;

  Relaxation();
}

// if no neighbour exists in the Adjust_Radii-function (or in the
// SetSpringRatio-function), then go and get x-/y-coordinates HERE
// returns x-value, if xy==0;;;y-value, if xy==1

double Lattice::VirtualNeighbour(int spring, int xy) {
  double coord, alpha, beta, total_angle, numerator, denominator, ratio,
      upper_angle, lower_angle, upper_xdiff, lower_xdiff, upper_ydiff,
      lower_ydiff, buffer, buffer2, buffer3, buffer4;
  int i, j, k, l, n, upper_empty_springs, lower_empty_springs,
      total_empty_springs, lower_real_spring, upper_real_spring, sign;
  bool upper_atan_exists, lower_atan_exists;

  // first find out, if there exist more virtual neighbours, and if so,
  // count them
  // while, da 5/0 grenze ueberschritten werden kann
  upper_empty_springs = lower_empty_springs = 0;
  // nach oben
  i = (spring < 5) ? spring + 1 : 0;
  while (!runParticle->neigP[i]) {
    upper_empty_springs++;
    i++;
    if (i == 6) {
      i = 0;
    }
  }
  // nach unten
  i = (spring > 0) ? spring - 1 : 5;
  while (!runParticle->neigP[i]) {
    lower_empty_springs++;
    i--;
    if (i == -1) {
      i = 5;
    }
  }

  total_empty_springs = upper_empty_springs + lower_empty_springs + 1;

  if (total_empty_springs != 6) // achtung: s.u.
  {
    lower_real_spring = spring - lower_empty_springs - 1;
    lower_real_spring =
        (lower_real_spring < 0) ? 6 + lower_real_spring : lower_real_spring;
    upper_real_spring = spring + upper_empty_springs + 1;
    upper_real_spring =
        (upper_real_spring > 5) ? upper_real_spring - 6 : upper_real_spring;
  } else // dh wenn nur noch eine spring exitiert
         // (kreis->keine 5 total_empty_springs!!!)
  { // dieser workaround sollte bei
    // bestimmungsgemaessem ablauf NICHT zur
    // anwendung kommen!
    cout << "fucked up in VirtualNeighbour()" << endl;
    i = 0;
    while (!runParticle->neigP[i]) {
      i++;
    }
    lower_real_spring = i;
    upper_real_spring = i;
  }

  // differenzen festlegen, um richtung zu indizieren
  upper_xdiff = runParticle->xpos - runParticle->neigP[upper_real_spring]->xpos;
  upper_ydiff = runParticle->ypos - runParticle->neigP[upper_real_spring]->ypos;
  lower_xdiff = runParticle->xpos - runParticle->neigP[lower_real_spring]->xpos;
  lower_ydiff = runParticle->ypos - runParticle->neigP[lower_real_spring]->ypos;

  // falls atan nicht loesbar...
  if (lower_xdiff == 0) {
    lower_atan_exists = false;
  } else {
    lower_angle =
        atan((runParticle->neigP[lower_real_spring]->ypos - runParticle->ypos) /
             (runParticle->neigP[lower_real_spring]->xpos - runParticle->xpos));
    lower_atan_exists = true;
  }

  if (upper_xdiff == 0) {
    upper_atan_exists = false;
  } else {
    upper_angle =
        atan((runParticle->neigP[upper_real_spring]->ypos - runParticle->ypos) /
             (runParticle->neigP[upper_real_spring]->xpos - runParticle->xpos));
    upper_atan_exists = true;
  }

  if (lower_xdiff > 0) {
    lower_angle = lower_angle + Pi;
  } else if (lower_xdiff < 0 && lower_ydiff > 0) {
    lower_angle = lower_angle + 2 * Pi;
  }

  if (upper_xdiff > 0) {
    upper_angle = upper_angle + Pi;
  } else if (upper_xdiff < 0 && upper_ydiff > 0) {
    upper_angle = upper_angle + 2 * Pi;
  }

  if (lower_real_spring < upper_real_spring) {
    if (lower_angle < upper_angle)
      total_angle = upper_angle - lower_angle;
    else
      total_angle = 2 * Pi - lower_angle + upper_angle;
  } else {
    if (lower_angle > upper_angle)
      total_angle = 2 * Pi - lower_angle + upper_angle;
    else
      total_angle = upper_angle - lower_angle;
  }

  if (total_empty_springs == 5)
    total_angle = 2 * Pi;

  // hier erstmal ratio ueber gesamtheit der kaputten (nachbar-)springs
  // kalkulieren
  // variablen: ratio, buffer, i, j u.v.m.
  // ursprngliche version

  // dritte version

  ratio = 0;

  i = lower_real_spring + 1;
  if (i > 5)
    i = 0;

  n = spring + 1;
  if (n > 5)
    n = 0;

  while (i != n) {
    j = i;
    buffer = 1.0;

    while (j != n) {
      buffer /= runParticle->ratio[j];

      j++;
      if (j > 5)
        j = 0;
    }

    buffer2 = buffer;

    j = spring + 1;
    if (j > 5)
      j = 0;

    while (j != upper_real_spring) {
      buffer2 += buffer2 / runParticle->ratio[j];

      j++;
      if (j > 5)
        j = 0;
    }

    ratio += 1.0 / buffer2;

    i++;
    if (i > 5)
      i = 0;
  }

  // winkel berechnen
  alpha = total_angle / ((1.0 / ratio) + 1); // angle OVER missing
  // spring (in terms of
  // indices)
  beta = total_angle / (ratio + 1); // abgle UNDER missing spring (in
  // terms of indices)

  if (xy == 0) // x-wert berechnen und zurueckgeben
  {
    // hier vermutlich noch jede MENGE SPEZIALFAELLE (vielleicht auch
    // nicht...)
    coord =
        cos(alpha + lower_angle) * runParticle->rad[spring] + runParticle->xpos;
  } else // y-wert berechnen und zurueckgeben
  {
    // hier vermutlich noch jede MENGE SPEZIALFAELLE (hier wohl auch
    // nicht...)
    coord =
        sin(alpha + lower_angle) * runParticle->rad[spring] + runParticle->ypos;
  }

  return (coord);
}

void Lattice::SetSpringAngle() {
  int i, j;
  float spring_angle[6];
  float ref_xpos, ref_ypos, xdiff, ydiff;
  bool around_flag;

  for (i = 0; i < numParticles; i++) {
    around_flag = false;
    // einmal durch alle springs des partikels loopen und jeder spring
    // einen absolut-winkelwert zuweisen
    for (j = 0; j < 6; j++) {
      if (runParticle->neigP[j]) {
        ref_xpos = runParticle->neigP[j]->xpos;
        ref_ypos = runParticle->neigP[j]->ypos;

        spring_angle[j] = atan((runParticle->ypos - ref_ypos) /
                               (runParticle->xpos - ref_xpos));

        xdiff = runParticle->xpos - ref_xpos;
        ydiff = runParticle->ypos - ref_ypos;

        if (xdiff > 0) {
          spring_angle[j] = spring_angle[j] + Pi;
        } else if (xdiff < 0 && ydiff > 0) {
          spring_angle[j] = spring_angle[j] + 2 * Pi;
        }

        if (xdiff == 0.0) {
          if (ydiff < 0.0) {
            spring_angle[j] = 3.14159265 / 2.0;
          } else {
            spring_angle[j] = 3.14159265 * 1.5;
          }
        }

        runParticle->spring_angle[j] = spring_angle[j];
      }
    }
    runParticle = runParticle->nextP;
  }
}

// setze nach jeder deformation neue ratio fuer noch AKTIVE springs
// aufruf am ende von AdjustRadii + konstruktor
void Lattice::SetSpringRatio() {
  int i, j, k, l, m, o;
  double ratio, angle[3], xpos, ypos, ref_xpos, ref_ypos, ref_angle, xdiff,
      ydiff;
  bool springs_present, atan_exists, ref_atan_exists, zero_jump_up,
      zero_jump_down;

  runParticle = &refParticle;

  for (j = 0; j < numParticles; j++) {
    // find out, if there exist ANY springs (otherwise seg-fault in
    // VirtualNeighbour)
    o = 0;
    springs_present = true;
    for (m = 0; m < 6; m++) {
      if (!runParticle->neigP[m]) {
        o++;
      }
    }
    if (o > 5)
      springs_present = false;

    // if there are springs, continue
    if (!runParticle->is_lattice_boundary) {
      if (springs_present) {

        runParticle->no_spring_flag = false;

        for (i = 0; i < 6; i++) {

          if (runParticle->neigP[i]) {

            // for neigP i, angle = angle from positive x-axis
            ref_xpos = runParticle->neigP[i]->xpos;
            ref_ypos = runParticle->neigP[i]->ypos;
            ref_angle = atan((ref_ypos - runParticle->ypos) /
                             (ref_xpos - runParticle->xpos));
            ref_atan_exists = true;

            xdiff = runParticle->xpos - ref_xpos;
            ydiff = runParticle->ypos - ref_ypos;

            if (xdiff > 0) {
              ref_angle += Pi;
            } else if (xdiff < 0 && ydiff > 0) {
              ref_angle += 2 * Pi;
            }

            for (k = i - 1; k <= i + 1; k += 2) {

              l = k;

              zero_jump_up = false;
              zero_jump_down = false;

              if (l == -1) {
                l = 5;
                zero_jump_down = true;
              }
              if (l == 6) {
                l = 0;
                zero_jump_up = true;
              }
              if (l > 6)
                cout << "autsch" << endl;

              if (runParticle->neigP[l]) {
                xpos = runParticle->neigP[l]->xpos;
                ypos = runParticle->neigP[l]->ypos;
              } else {
                xpos = VirtualNeighbour(l, 0);
                ypos = VirtualNeighbour(l, 1);
              }

              xdiff = runParticle->xpos - xpos;
              ydiff = runParticle->ypos - ypos;

              angle[k - i + 1] = atan(ydiff / xdiff);

              if (xdiff > 0) {
                angle[k - i + 1] += Pi;
              } else if (xdiff < 0 && ydiff > 0) {
                angle[k - i + 1] += 2 * Pi;
              }

              if (zero_jump_down) {
                if (ref_angle < angle[k - i + 1]) {
                  angle[k - i + 1] = 2 * Pi - angle[k - i + 1] + ref_angle;
                } else {
                  angle[k - i + 1] = ref_angle - angle[k - i + 1];
                }
              } else if (zero_jump_up) {
                if (ref_angle < angle[k - i + 1]) {
                  angle[k - i + 1] = angle[k - i + 1] - ref_angle;
                } else {
                  angle[k - i + 1] = 2 * Pi + angle[k - i + 1] - ref_angle;
                }
              } else {
                if (k < i) {
                  if (ref_angle > angle[k - i + 1]) {
                    angle[k - i + 1] = ref_angle - angle[k - i + 1];
                  } else {
                    angle[k - i + 1] = 2 * Pi + ref_angle - angle[k - i + 1];
                  }
                } else {
                  if (angle[k - i + 1] > ref_angle) {
                    angle[k - i + 1] = angle[k - i + 1] - ref_angle;
                  } else {
                    angle[k - i + 1] = 2 * Pi - ref_angle + angle[k - i + 1];
                  }
                }
              }
            }
            // achtung: so rum oder anders rum???
            runParticle->ratio[i] = fabs(angle[2] / angle[0]);
          } else {
            // dummy else for debugging
          }
        }
      } else // no springs present
      {
        runParticle->no_spring_flag = true;
      }
    }
    runParticle = runParticle->nextP;
  }
}

// function, that writes the averaged radius into a field, depending on
// the direction
// directions associated to values 1..9 depend on box in particle-class
// (for the repulsion)

// result will be written in particle->rep_rad[9]

// new stuff: now calculates length of the rep_rad using quadratic
// interpolation, based on the slopes on the edges and the x/y-positions
// of the edges themselves
// variable/name definition: y = a + bx + cx^2
void Lattice::BoxRad() {}

/*************************************************************
 * function that changes the relaxation threshold for
 * the relaxation routine. number < 1.0 is better relaxation
 * but relaxation takes longer.
 ************************************************************/

void Lattice::ChangeRelaxThreshold(float change) {
  relaxthres = relaxthres * change;
}

/***************************************************
 * function tilts elliptical particles
 *****************************************************/
/*
void
Lattice::Tilt ()
{
  int i;
  float buffer;

  buffer = 0.0;

  EraseNeighList ();

  //SetSpringAngle ();

  runParticle = &refParticle;

  for (i = 0; i < numParticles; i++)
    {
      // runParticle->Tilt();
      UpdateNeighList ();

      runParticle = runParticle->nextP;
    }

  //SetSpringAngle ();

  //UpdateSpringLength ();
}

*/

void Lattice::WeakenVerticalParticleLayer(double x_min, double x_max,
                                          double y_min, double y_max,
                                          float constant, float vis,
                                          float break_strength) {
  int i, j; // counters

  for (i = 0; i < numParticles; i++) // loop through particles
  {
    if (runParticle->xpos > x_min && runParticle->ypos > y_min) {
      if (runParticle->xpos < x_max && runParticle->ypos < y_max) {
        // -----------------------------------------------------
        // adjust youngs modulus of particle
        // -----------------------------------------------------

        runParticle->mineral = 2;

        runParticle->young = runParticle->young * constant;
        runParticle->viscosity *= vis;

        for (j = 0; j < 8; j++) // loop through neighbours
        {
          // -------------------------------------------------------------------------
          // we do not want to change the grain boundaries as
          // well, only springs
          // within the grain, therefore check which springs are
          // between particles
          // of the same grain.
          //
          // call function adjustgrainboundaries later to get
          // mean constants for boundaries
          // --------------------------------------------------------------------------

          if (runParticle->neigP[j]) // if neighbour
          {
            if (runParticle->neigP[j]->xpos > x_min &&
                runParticle->neigP[j]->ypos > y_min) {
              if (runParticle->neigP[j]->xpos < x_max &&
                  runParticle->neigP[j]->ypos < y_max) {
                runParticle->springf[j] =
                    runParticle->springf[j] * constant; // change
                // constant
                runParticle->springv[j] = runParticle->springv[j] * vis;
                runParticle->break_Str[j] =
                    runParticle->break_Str[j] * break_strength; // change
                                                                // tensile
              }
            }
          }
        }
      }
    }
    runParticle = runParticle->nextP; // go on looping
  }

  AdjustParticleConstants();
}

void Lattice::SetRealDensityYoung(double x_max, double x_min, double y_max,
                                  double y_min, double den, double you) {
  int i;

  runParticle = &refParticle;

  for (i = 0; i < numParticles; i++) {

    if (runParticle->xpos < x_max && runParticle->xpos > x_min) {

      if (runParticle->ypos < y_max && runParticle->ypos > y_min) {

        runParticle->real_density = den;
        runParticle->real_young = you;
      }
    }

    runParticle = runParticle->nextP;
  }
}

void Lattice::Set_Absolute_Box(double size, int gr) {
  int i; // counter

  runParticle = &refParticle; // start at beginning

  for (i = 0; i < numParticles; i++) // and loop
  {
    runParticle->real_radius =
        runParticle->radius * size; // set a real radius for particles

    runParticle = runParticle->nextP; // loop to next particle
  }
  if (gr)
    grav = true; // turn gravity on for rift models
}

//-------------------------------------------------------------
// routine sets a fully scaled youngs modulus in mega pascal
// and a density for the particles as default, used for rifts
//-------------------------------------------------------------

void Lattice::InitRealDensityYoung(double den, double you) {
  int i;

  runParticle = &refParticle;

  for (i = 0; i < numParticles; i++) {

    runParticle->real_density = den;
    runParticle->real_young = you;

    runParticle = runParticle->nextP;
  }
}

//------------------------------------------------------
// Important routine that converts the set real scale
// young back to model scale constants
//------------------------------------------------------

void Lattice::InitInternalYoung(double young_standard) {
  int i, j;

  runParticle = &refParticle;

  for (i = 0; i < numParticles; i++) {

    for (j = 0; j < 6; j++) {
      if (runParticle->neigP[j]) {
        runParticle->springf[j] = sqrt(3.0) * runParticle->E /
                                  (3.0 * (1.0 - runParticle->poisson_ratio));
        runParticle->springs[j] = (1.0 - 3.0 * runParticle->poisson_ratio) /
                                  (1.0 + runParticle->poisson_ratio) *
                                  runParticle->springf[j];
      }
    }

    runParticle->young *= runParticle->real_young / young_standard;
    runParticle->E *= runParticle->E / young_standard;

    runParticle = runParticle->nextP;
  }
}

double Lattice::SpreadingStrain(double m, double box_size, double years) {

  double right_pos, absolute_size;

  runParticle = &refParticle;
  runParticle = runParticle->prevP;
  right_pos = runParticle->xpos;

  right_pos *= box_size;

  return ((m * years) / right_pos);
}

void Lattice::InitViscosity(double visc) {

  int i;

  runParticle = &refParticle;

  for (i = 0; i < numParticles; i++) {

    runParticle->viscosity = visc;

    runParticle = runParticle->nextP;
  }
}

void Lattice::SetViscosity(double x_max, double x_min, double y_max,
                           double y_min, double visc) {
  int i;

  runParticle = &refParticle;

  for (i = 0; i < numParticles; i++) {

    if (runParticle->xpos < x_max && runParticle->xpos > x_min) {

      if (runParticle->ypos < y_max && runParticle->ypos > y_min) {

        runParticle->viscosity = visc;
      }
    }

    runParticle = runParticle->nextP;
  }
}

void Lattice::SetBreakingStrength(double x_max, double x_min, double y_max,
                                  double y_min, double factor) {

  int i, j;

  runParticle = &refParticle;

  for (i = 0; i < numParticles; i++) {

    if (runParticle->xpos <= x_max && runParticle->xpos >= x_min) {

      if (runParticle->ypos <= y_max && runParticle->ypos >= y_min) {

        for (j = 0; j < 6; j++) {
          runParticle->break_Str[j] *= factor;
        }
      }
    }

    runParticle = runParticle->nextP;
  }
}

void Lattice::UnsetNobreakUpperRow() {

  int i, j;

  runParticle = &refParticle;

  for (i = 0; i < numParticles; i++) {
    if (runParticle->is_top_lattice_boundary) {
      for (j = 0; j < 6; j++) {
        if (runParticle->neigP[j]) {
          runParticle->no_break[j] = false;
        }
      }
      runParticle->fix_y = false;
      cout << runParticle->ypos << endl;
    }

    runParticle = runParticle->nextP;
  }
}

// viscous sheet for pseudo 3D simulations

void Lattice::ActivateSheet(float visc, float young) {
  int i;

  sheet = true;
  cout << visc << endl;
  cout << young << endl;
  for (i = 0; i < numParticles; i++) // run
  {
    runParticle->sheet_x = runParticle->xpos; // call SetSprings
    runParticle->sheet_y = runParticle->ypos;
    runParticle->sheet_young = young;
    runParticle->sheet_visc = visc;

    runParticle = runParticle->nextP; // go to next particle
  }
}
//--------------------------------------------------------------------
// This is for boundary particles, for the upper and lower boundary
// these are typically fixed in y so they dont move out of the box
// here they can be released
// In addtion normally they are not allowed to fracture,
// This is also allowed here
//---------------------------------------------------------------------

void Lattice::ReleaseBoundaryParticlesY(int fracture, int release) {
  int i, j, jj;

  for (i = 0; i < numParticles; i++) {
    if (runParticle->fix_y == true) // if fixed in y
    {
      if (release) // turn of fix y
        runParticle->fix_y = false;

      if (runParticle->fix_x ==
          false) // in case its not a left or right boundary particle as well
      {
        if (fracture) // allow it to fracture

        {
          runParticle->is_lattice_boundary = false; // not a boundary any more

          for (j = 0; j < 8; j++) {
            if (runParticle->neigP[j]) // if there is a neighbour
            {
              runParticle->no_break[j] = false; // make it breakable

              for (jj = 0; jj < 8; jj++) {
                if (runParticle->neigP[j]->neigP[jj]) {
                  if (runParticle->neigP[j]->neigP[jj]->nb == runParticle->nb) {

                    // make it breakable, this is the link back from the
                    // neighbour

                    runParticle->neigP[j]->no_break[jj] = false;
                  }
                }
              }
            }
          }
        }
      }
    }

    runParticle = runParticle->nextP;
  }
}

// Relaxation routine for the viscous sheet
// each particle has a sheet position sheet_x and sheet_y
//

void Lattice::RelaxSheet(float t_step) {
  int i;
  double difx, dify, dx, dy, ddx, ddy;

  for (i = 0; i < numParticles; i++) // loop
  {
    difx = runParticle->xpos - runParticle->sheet_x;
    dify = runParticle->ypos - runParticle->sheet_y;

    dx = sqrt(difx * difx);
    dy = sqrt(dify * dify);

    // exponential relaxation

    ddx = dx - dx * exp(-runParticle->sheet_young * t_step /
                        runParticle->sheet_visc);
    ddy = dy - dy * exp(-runParticle->sheet_young * t_step /
                        runParticle->sheet_visc);

    if (dx != 0.0)
      difx = ddx * difx / dx;
    if (dy != 0.0)
      dify = ddy * dify / dy;

    // move the sheet points

    runParticle->sheet_x = runParticle->sheet_x + difx;
    runParticle->sheet_y = runParticle->sheet_y + dify;

    // cout << runParticle->sheet_x << difx << endl;

    runParticle = runParticle->nextP;
  }
}

void Lattice::SetBreakingStrengthGradient(double x_max, double x_min,
                                          double y_max, double y_min,
                                          double c) {
  int i, j;

  runParticle = &refParticle;

  for (i = 0; i < numParticles; i++) {

    for (j = 0; j < 6; j++) {

      runParticle->break_Str[j] *=
          (1.0 + (1.0 - runParticle->ypos) * (c - 1.0));
    }

    runParticle = runParticle->nextP;
  }
}

void Lattice::Initialize_Pressure_Lattice(double Press_Background, int scale) {
  int i;

  runParticle = &refParticle;
  for (i = 0; i < numParticles; i++) {
    runParticle->velx = 0;
    runParticle->vely = 0;

    runParticle->oldx = runParticle->xpos;
    runParticle->oldy = runParticle->ypos;

    runParticle->rad_par_fluid = runParticle->radius;

    runParticle->area_par_fluid =
        Pi * pow(runParticle->rad_par_fluid, 2.0); // area of particle radius

    runParticle = runParticle->nextP;
  }

  P->d_pressure_node_pos_n_area(scale);

  P->pressure_background(Press_Background, 1);
}

void Lattice::Make_Seal(float seal_y_min, float seal_y_max, float seal_x_min,
                        float seal_x_max, float rad_mean, float rad_sigma) {
  float ran_nb;   // random number
  float rad_frac; // radius fraction constant
  float prob;

  int i;

  // runParticle = &refParticle;
  // for (i = 0; i < numParticles; i++)
  //{
  // if ( runParticle->ypos >= seal_y_min )
  //{
  // if (runParticle->ypos <= seal_y_max)
  //{
  // if (runParticle->xpos >= seal_x_min)
  //{
  // if (runParticle->xpos <= seal_x_max)
  //{
  // runParticle->rad_par_fluid *= rad_mean;
  // runParticle->area_par_fluid *= pow(rad_mean,2.0);
  //}
  //}
  //}
  //}
  // runParticle = runParticle->nextP;
  //}

  runParticle = &refParticle;
  for (i = 0; i < numParticles; i++) {
    if (runParticle->ypos >= seal_y_min) {
      if (runParticle->ypos <= seal_y_max) {
        // if (runParticle->xpos >=seal_x_min)
        {
          // if (runParticle->xpos <= seal_x_max)
          {
            ran_nb = rand() / (float)RAND_MAX;

            rad_frac = ran_nb * rad_sigma; // now from 0.0 to rad_sigma

            runParticle->rad_par_fluid =
                runParticle->rad_par_fluid +
                (runParticle->rad_par_fluid * (rad_frac - (rad_sigma / 2.0)));

            runParticle->rad_par_fluid = runParticle->rad_par_fluid * rad_mean;

            runParticle->area_par_fluid =
                Pi * pow(runParticle->rad_par_fluid, 2.0);
          }
        }
      }
    }
    runParticle = runParticle->nextP;
  }

  //~ runParticle = &refParticle;
  //~
  //~ for (i = 0; i < numParticles; i++)
  //~ {
  //~ if ( runParticle->ypos >= seal_y_min )
  //~ {
  //~ if (runParticle->ypos <= seal_y_max)
  //~ {
  //~ //if (runParticle->xpos >=seal_x_min)
  //~ {
  //~ //if (runParticle->xpos <= seal_x_max)
  //~ {
  //~ do
  //~ {
  //~ // pic a pseudorandom float between 0 and 1
  //~ rad_frac = rand () / (float) RAND_MAX;
  //~
  //~ // ----------------------------------------------------------------
  //~ // now convert this distribution to a distribution
  //~ // from zero minus 2 times sigma to plus 2 times sigma
  //~ // -----------------------------------------------------------------
  //~ rad_frac = (rad_frac - 0.5) * 2.0 * (2.0 * rad_sigma);
  //~
  //~ // --------------------------------------------------------------
  //~ // and shift it to mean plus minus 2 times sigma
  //~ // --------------------------------------------------------------
  //~ rad_frac = rad_frac + rad_mean;
  //~
  //~ // ---------------------------------------------------------------
  //~ // now apply gauss function to determine a probability
  //~ // for this reaction constant to occur
  //~ // ---------------------------------------------------------------
  //~ prob = (1 / (rad_sigma * sqrt (2.0 * 3.1415927))); // part one (right
  //part) of function
  //~
  //~ prob = prob * (exp (-0.5 * (((rad_frac) - rad_mean) / rad_sigma)
  //~ * (((rad_frac) - rad_mean) / rad_sigma)));	//rest of Gauss function
  //~
  //~ // --------------------------------------------------------------
  //~ // now adjust probablity to run from 0 to 1.0
  //~ // --------------------------------------------------------------
  //~ prob = prob * sqrt (2.0 * 3.1415927) * rad_sigma;
  //~
  //~ // ---------------------------------------------------------------
  //~ // pic the second pseudorandom number
  //~ // ---------------------------------------------------------------
  //~ ran_nb = rand () / (float) RAND_MAX;
  //~
  //~ // --------------------------------------------------------------
  //~ // if the number picted is smaller or the same as the
  //~ // probability from the gauss function accept the rate constant,
  //~ // if not go on looping.
  //~ // -------------------------------------------------------------------
  //~
  //~ if (ran_nb <= prob)	// if smaller
  //~ {
  //~ if (rad_frac <= 0.0) // if 0.0 or smaller (dont want that)
  //~ ;
  //~ else
  //~ {
  //~ // scale to particle default radius
  //~ rad_frac = rad_frac * 1.0;
  //~
  //~ // set the constant to change particle radius and area !
  //~ runParticle->rad_par_fluid *= rad_frac;
  //~ runParticle->area_par_fluid *= pow(rad_frac,2.0);
  //~
  //~ //cout << runParticle->rad_par_fluid << " ";
  //~
  //~ // and break the while loop !
  //~ break;
  //~ }
  //~ }
  //~ } while (1);	// loop until break
  //~
  //~ //if (runParticle->xpos >= seal_x_max-0.02 && runParticle->xpos <=
  //seal_x_max) ~ //cout << endl << endl ;
  //~ }
  //~ }
  //~ }
  //~ }
  //~ runParticle = runParticle->nextP;
  //~ }
}

void Lattice::Make_Seal_Slot(float seal_y_min, float seal_y_max,
                             float seal_x_min, float seal_x_max,
                             float rad_frac) {
  int i;

  runParticle = &refParticle;

  for (i = 0; i < numParticles; i++) {
    if ((runParticle->ypos >= seal_y_min) &&
        (runParticle->ypos <= seal_y_max)) {
      // if ( (runParticle->xpos <= seal_x_min) || (runParticle->xpos >=
      // seal_x_max))
      {
        runParticle->rad_par_fluid *= rad_frac;
        runParticle->area_par_fluid *= pow(rad_frac, 2.0);
      }
    }
    runParticle = runParticle->nextP;
  }
}

void Lattice::Pressure_init_Boundary_cond(int LboundrySourc, double boundry_val,
                                          int vertical_z_depth) {
  P->boundary_fixed(LboundrySourc, boundry_val, vertical_z_depth);

  // P->background_press();		//need in case of using
  // Calculate_Pressure() in Run.
  P->source_pressure(); // part of pressure source term: approximated pressure
                        // value at the end of last time step
}

void Lattice::Pressure_initialize_Hydrostatic_gradient(int box_size,
                                                       int depth) {
  // P->boundary_fixed(LboundrySourc, boundry_val, vertical_z_depth);
  P->boundary_set_hydrostatic_gradient(box_size, depth);

  // P->background_press();		//need in case of using
  // Calculate_Pressure() in Run.
  P->source_pressure(); // part of pressure source term: approximated pressure
                        // value at the end of last time step
}

void Lattice::Pressure_initialize_Hydrostatic(int box_size, int depth,
                                              double add_pressure) {
  // P->boundary_fixed(LboundrySourc, boundry_val, vertical_z_depth);
  P->boundary_set_hydrostatic(box_size, depth, add_pressure);

  // P->background_press();		//need in case of using
  // Calculate_Pressure() in Run.
  P->source_pressure(); // part of pressure source term: approximated pressure
                        // value at the end of last time step
}

void Lattice::pressure_for_run(float seal_y_min) {
  P->grad_pressure(seal_y_min);
  // P->background_press();	//need in case of using Calculate_Pressure() in
  // Run
  P->source_pressure(); // part of pressure source term: approximated pressure
                        // value at the end of last time step
}

void Lattice::Pressure_init_source(int press_node_y, int press_node_x,
                                   double press_increment) {
  P->pressure_incre_source(press_node_y, press_node_x, press_increment);
}

void Lattice::Pressure_init_RandSource(int rand_min_y, int rand_max_y,
                                       int rand_min_x, int rand_max_x,
                                       double pressure_min,
                                       double pressure_max) {
  srand((unsigned)time(0));
  int rand_integer_y, rand_integer_x;
  double rand_press_increment;

  rand_integer_y = rand_min_y + rand() % (rand_max_y + 1 - rand_min_y);
  rand_integer_x = rand_min_x + rand() % (rand_max_x + 1 - rand_min_x);

  // rand_press_increment = pressure_min + rand() % (pressure_max+1 -
  // pressure_min);
  rand_press_increment = (rand() / (static_cast<double>(RAND_MAX) + 1.0)) *
                             (pressure_max - pressure_min) +
                         pressure_min;
  /********************************************************
   * double temp;
   * temp = (rand() / (static_cast<double>(RAND_MAX) + 1.0))* (5 - 1) + 1;
   * *****************************************************/

  P->rand_source_PressIncrement(rand_integer_y, rand_integer_x,
                                rand_press_increment);
}

void Lattice::Pressure_init_RandSource_hor_lay(int hor_lay, int rand_min_x,
                                               int rand_max_x,
                                               double pressure_min,
                                               double pressure_max) {
  srand((unsigned)time(0));
  int rand_integer_x;
  double rand_press_increment;

  rand_integer_x = rand_min_x + rand() % (rand_max_x + 1 - rand_min_x);

  // rand_press_increment = pressure_min + rand() % (pressure_max+1 -
  // pressure_min);
  rand_press_increment = (rand() / (static_cast<double>(RAND_MAX) + 1.0)) *
                             (pressure_max - pressure_min) +
                         pressure_min;
  /********************************************************
   * double temp;
   * temp = (rand() / (static_cast<double>(RAND_MAX) + 1.0))* (5 - 1) + 1;
   * *****************************************************/

  P->rand_source_PressIncrement(hor_lay, rand_integer_x, rand_press_increment);
}

// void Lattice::Pressure_init_twoXsources(int y_press_node, int
// x_press_node_diff, double increment_press) //(y_position, difference b/w
// x_positions from central node)
//{
// P->Pressure_incre_twoXsources(y_press_node, x_press_node_diff,
// increment_press); P->source_pressure(); 	// pressure source term:
// approximated pressure value at the end of last time step
//}

// void Lattice::Pressure_init_twoYsources(int x_press_node, int
// y_press_node_diff, double increment_press) //(x_position, difference b/w
// y_positions from central node)
//{
// P->Pressure_incre_twoYsources(x_press_node, y_press_node_diff,
// increment_press); P->source_pressure(); 	// pressure source term:
// approximated pressure value at the end of last time step
//}

// void Lattice::Pressure_init_source_n_sink(int press_node_y_a, int
// press_node_x_a, int press_node_y_b, int press_node_x_b, double
// press_decrement, double press_increment)
//{
// P->Pressure_source_n_sink(press_node_y_a, press_node_x_a, press_node_y_b,
// press_node_x_b, press_decrement, press_increment); P->source_pressure();
// // pressure source term: approximated pressure value at the end of last time
// step
//}

// void Lattice::Pressure_init_sink(int press_node_y, int press_node_x, double
// press_decrement)
//{
// P->Pressure_point_sink(press_node_y, press_node_x, press_decrement);
// P->source_pressure();	// pressure source term: approximated pressure
// value at the end of last time step
//}

// void Lattice::Pressure_constt_SinkArea(float y_box_min, float y_box_max,
// float x_box_min, float x_box_max, double press_decrement)
//{
// P->Pressure_sink(y_box_min, y_box_max, x_box_min, x_box_max,
// press_decrement); P->source_pressure();	// pressure source term:
// approximated pressure value at the end of last time step
//}

// void Lattice::Pressure_SinkArea_cont(float y_box_min, float y_box_max, float
// x_box_min, float x_box_max, double press_decrement)
//{
// P->Pressure_sink_cont(y_box_min, y_box_max, x_box_min, x_box_max,
// press_decrement); P->source_pressure();	// pressure source term:
// approximated pressure value at the end of last time step
//}

void Lattice::Pressure_force_range_text(double y_box_min, double y_box_max,
                                        double x_box_min, double x_box_max) {
  FILE *pressure;
  int i;

  pressure = fopen("pressure_force.txt", "a");

  for (i = 0; i < numParticles; i++) {
    if (runParticle->xpos > x_box_min) {
      if (runParticle->xpos < x_box_max) {
        if (runParticle->ypos > y_box_min) {
          if (runParticle->ypos < y_box_max) {
            fprintf(pressure, "%d:", i);
            fprintf(pressure, " %20lf", runParticle->F_P_x);
            fprintf(pressure, " %24lf\n", runParticle->F_P_y);
          }
        }
      }
    }

    runParticle = runParticle->nextP;
  }

  fprintf(pressure, "\n");
  fprintf(pressure, "\t");

  fclose(pressure);
}

void Lattice::Pressure_RandSource_fix(int rand_min_y, int rand_max_y,
                                      int rand_min_x, int rand_max_x,
                                      double press_increment) {
  int rand_integer_y, rand_integer_x;

  rand_integer_y = rand_min_y + rand() % (rand_max_y + 1 - rand_min_y);
  rand_integer_x = rand_min_x + rand() % (rand_max_x + 1 - rand_min_x);

  //~ cout << P->pressure_dist[rand_integer_y][rand_integer_x]<< " ";
  //~ if (P->pressure_dist[rand_integer_y][rand_integer_x] < 1.3e+07)
  {
    P->rand_source_PressIncrement(rand_integer_y, rand_integer_x,
                                  press_increment);
    cout << P->pressure_dist[rand_integer_y][rand_integer_x] << endl;
  }
  //~ else
  //~ ;
}

void Lattice::Pressure_RandSource(int rand_min_y, int rand_max_y,
                                  int rand_min_x, int rand_max_x,
                                  double pressure_min, double pressure_max) {
  int rand_integer_y, rand_integer_x;
  double rand_press_increment;

  rand_integer_y = rand_min_y + rand() % (rand_max_y + 1 - rand_min_y);
  rand_integer_x = rand_min_x + rand() % (rand_max_x + 1 - rand_min_x);

  // rand_press_increment = pressure_min + rand() % (pressure_max+1 -
  // pressure_min);
  rand_press_increment = (rand() / (static_cast<double>(RAND_MAX) + 1.0)) *
                             (pressure_max - pressure_min) +
                         pressure_min;
  /********************************************************
   * double temp;
   * temp = (rand() / (static_cast<double>(RAND_MAX) + 1.0))* (5 - 1) + 1;
   * *****************************************************/

  P->rand_source_PressIncrement(rand_integer_y, rand_integer_x,
                                rand_press_increment);
}

void Lattice::Pressure_RandFixSource_hor_lay(int lay_y, int rand_min_x,
                                             int rand_max_x, double press_fix) {
  int rand_integer_x;

  rand_integer_x = rand_min_x + rand() % (rand_max_x + 1 - rand_min_x);

  P->rand_source_PressIncrement(lay_y, rand_integer_x, press_fix);
}

void Lattice::Pressure_RandSource_hor_lay(int lay_y, int rand_min_x,
                                          int rand_max_x, double pressure_min,
                                          double pressure_max) {
  int rand_integer_x;
  double rand_press_increment;

  rand_integer_x = rand_min_x + rand() % (rand_max_x + 1 - rand_min_x);

  // rand_press_increment = pressure_min + rand() % (pressure_max+1 -
  // pressure_min);
  rand_press_increment = (rand() / (static_cast<double>(RAND_MAX) + 1.0)) *
                             (pressure_max - pressure_min) +
                         pressure_min;
  /********************************************************
   * double temp;
   * temp = (rand() / (static_cast<double>(RAND_MAX) + 1.0))* (5 - 1) + 1;
   * *****************************************************/

  P->rand_source_PressIncrement(lay_y, rand_integer_x, rand_press_increment);
}

void Lattice::Pressure_RandSource_ver_lay(int lay_x, int rand_min_y,
                                          int rand_max_y, double pressure_min,
                                          double pressure_max) {
  int rand_integer_y;
  double rand_press_increment;

  rand_integer_y = rand_min_y + rand() % (rand_max_y + 1 - rand_min_y);

  // rand_press_increment = pressure_min + rand() % (pressure_max+1 -
  // pressure_min);
  rand_press_increment = (rand() / (static_cast<double>(RAND_MAX) + 1.0)) *
                             (pressure_max - pressure_min) +
                         pressure_min;
  /********************************************************
   * double temp;
   * temp = (rand() / (static_cast<double>(RAND_MAX) + 1.0))* (5 - 1) + 1;
   * *****************************************************/

  P->rand_source_PressIncrement(rand_integer_y, lay_x, rand_press_increment);
}

void Lattice::Pressure_xGauss_RandSource_area(int y_rand_min, int y_rand_max,
                                              int x_rand_min, int x_rand_max,
                                              double press_increment,
                                              float P_sigma) {
  int rand_integer_y;
  int pos_x;
  float ran_nb;
  float prob;

  P_x_mean = 0.5 * P->max_size;

  sourc_Ymin = y_rand_min * P->max_size;
  sourc_Ymax = y_rand_max * P->max_size;
  sourc_Xmin = x_rand_min * P->max_size;
  sourc_Xmax = x_rand_max * P->max_size;

  // rand_integer_y = sourc_Ymin + rand() % (sourc_Ymax+1 - sourc_Ymin);
  rand_integer_y = (rand() / (static_cast<float>(RAND_MAX) + 1.0)) *
                       (sourc_Ymax - sourc_Ymin) +
                   sourc_Ymin;

  do {
    ran_nb =
        rand() / (float)RAND_MAX; // get a float pseudorandom number b/w 0 and 1
    ran_nb =
        (ran_nb - 0.5) * 3.0 *
        (3.0 * P_sigma); // get a distribution from plus/minus 3 times sigma
    ran_nb = ran_nb + P_x_mean; // shift it around x_mean

    if (ran_nb > sourc_Xmin && ran_nb < sourc_Xmax) {
      ran_nb = ran_nb + 0.5;
      pos_x = ran_nb;

      break;
    }
  } while (1);

  prob =
      (1.0 / (P_sigma *
              sqrt(2.0 * Pi))); // first part of Gaussian distribution function
  prob =
      prob * (exp(-0.5 * ((pos_x - P_x_mean) / P_sigma) *
                  ((pos_x - P_x_mean) / P_sigma))); // rest of Gaussian function
  prob = prob * sqrt(2.0 * Pi) *
         P_sigma; // adjust probability to run from 1 to 1.0

  // prob = prob * 100.0;	// get percentage of probability
  // press_increment = (press_increment * prob)/100.0;  // get fraction of given
  // pressure
  press_increment = press_increment * prob;

  if (P->pressure_dist[rand_integer_y][pos_x] < press_increment) {
    P->rand_source_PressIncrement(rand_integer_y, pos_x, press_increment);
  } else
    ;
}

void Lattice::Pressure_yGauss_RandSource_area(int y_rand_min, int y_rand_max,
                                              int x_rand_min, int x_rand_max,
                                              double press_increment,
                                              float P_y_mean, float P_sigma) {
  int rand_integer_x;
  int pos_y;
  float ran_nb;
  float prob;

  sourc_Ymin = y_rand_min * P->max_size;
  sourc_Ymax = y_rand_max * P->max_size;
  sourc_Xmin = x_rand_min * P->max_size;
  sourc_Xmax = x_rand_max * P->max_size;

  // rand_integer_x = sourc_Xmin + rand() % (sourc_Xmax+1 - sourc_Xmin);
  rand_integer_x = (rand() / (static_cast<float>(RAND_MAX) + 1.0)) *
                       (sourc_Xmax - sourc_Xmin) +
                   sourc_Xmin;

  do {
    ran_nb =
        rand() / (float)RAND_MAX; // get a float pseudorandom number b/w 0 and 1
    ran_nb =
        (ran_nb - 0.5) * 3.0 *
        (3.0 * P_sigma); // get a distribution from plus/minus 3 times sigma
    ran_nb = ran_nb + P_y_mean; // shift it around y_mean

    if (ran_nb > sourc_Ymin && ran_nb < sourc_Ymax) {
      ran_nb = ran_nb + 0.5;
      pos_y = ran_nb;

      break;
    }
  } while (1);

  prob =
      (1.0 / (P_sigma *
              sqrt(2.0 * Pi))); // first part of Gaussian distribution function
  prob =
      prob * (exp(-0.5 * ((pos_y - P_y_mean) / P_sigma) *
                  ((pos_y - P_y_mean) / P_sigma))); // rest of Gaussian function
  prob = prob * sqrt(2.0 * Pi) *
         P_sigma; // adjust probability to run from 0 to 1.0

  // prob = prob * 100.0;	// get percentage of probability
  // press_increment = (press_increment * prob)/100.0;  // get fraction of given
  // pressure
  press_increment = press_increment * prob;

  if (P->pressure_dist[pos_y][rand_integer_x] < press_increment) {
    P->rand_source_PressIncrement(pos_y, rand_integer_x, press_increment);
  } else
    ;
}

void Lattice::Melt_Parameters(double viscosity, double compressibility) {
  P->Make_Melt(viscosity, compressibility);
}

void Lattice::Fracture_Reaction(float fluid_threshold, float prob,
                                float max_young) {
  int i;
  float ran_nb;
  // float prob;

  srand(time(0));
  for (i = 0; i < numParticles; i++) {
    if (runParticle->draw_break) // is bond fractured? runParticle->temperature
    {
      if (runParticle->temperature > fluid_threshold) {
        ran_nb = rand() / (float)RAND_MAX;

        runParticle->fracture_reaction =
            runParticle->fracture_reaction +
            (1 / (1 + runParticle->fracture_reaction)) * prob;
        if (runParticle->young < max_young) {
          runParticle->young =
              runParticle->young + runParticle->fracture_reaction;
          runParticle->E = runParticle->E + runParticle->fracture_reaction;
        } else {
          runParticle->young = max_young;
          runParticle->E = max_young;
        }

        for (int j = 0; j < 8; j++) {
          runParticle->springf[j] =
              runParticle->springf[j] * runParticle->fracture_reaction;
          runParticle->springs[j] =
              runParticle->springs[j] * runParticle->fracture_reaction;
        }
      }
    }
    runParticle = runParticle->nextP;
  }
  AdjustParticleConstants();
}

void Lattice::Pressure_2dCircularGauss_RandSource_area(
    float y_rand_min, float y_rand_max, float x_rand_min, float x_rand_max,
    double press_increment, float P_y_mean, float P_sigma) {
  float ran_nb;
  float prob;
  int pos_y, pos_x;

  sourc_Ymin = y_rand_min * P->max_size;
  sourc_Ymax = y_rand_max * P->max_size;
  sourc_Xmin = x_rand_min * P->max_size;
  sourc_Xmax = x_rand_max * P->max_size;

  P_x_mean = 0.5 * P->max_size;

  do {
    ran_nb =
        rand() / (float)RAND_MAX; // get a float pseudorandom number b/w 0 and 1
    ran_nb =
        (ran_nb - 0.5) * 3.0 *
        (3.0 * P_sigma); // get a distribution from plus/minus 3 times sigma
    ran_nb = ran_nb + P_y_mean; // shift it around y_mean

    if (ran_nb > sourc_Ymin && ran_nb < sourc_Ymax) {
      ran_nb = ran_nb + 0.5;
      pos_y = ran_nb;

      break;
    }
  } while (1);

  do {
    ran_nb =
        rand() / (float)RAND_MAX; // get a float pseudorandom number b/w 0 and 1
    ran_nb =
        (ran_nb - 0.5) * 3.0 *
        (3.0 * P_sigma); // get a distribution from plus/minus 3 times sigma
    ran_nb = ran_nb + P_x_mean; // shift it around x_mean

    if (ran_nb > sourc_Xmin && ran_nb < sourc_Xmax) {
      ran_nb = ran_nb + 0.5;
      pos_x = ran_nb;

      break;
    }
  } while (1);

  prob = exp(-0.5 * ((pos_x - P_x_mean) / P_sigma) *
             ((pos_x - P_x_mean) / P_sigma)); // exponential part with x
  prob =
      prob * (exp(-0.5 * ((pos_y - P_y_mean) / P_sigma) *
                  ((pos_y - P_y_mean) / P_sigma))); // exponential part with y
  prob = prob *
         (1.0 / (2.0 * Pi *
                 pow(P_sigma,
                     2.0))); // first part of 2d_Gaussian distribution function
  prob = prob * 2.0 * Pi *
         pow(P_sigma, 2.0); // adjust probability to run from 0 to 1.0

  press_increment = press_increment * prob;

  cout << "P->pressure_dist[pos_y][pos_x] : " << P->pressure_dist[pos_y][pos_x]
       << endl;

  if (P->pressure_dist[pos_y][pos_x] < press_increment) {
    P->rand_source_PressIncrement(pos_y, pos_x, press_increment);

    cout << "pos_x : " << pos_x << ", pos_y : " << pos_y << endl;
    cout << "press_increment : " << press_increment << endl;
  } else
    ;
}

void Lattice::Pressure_2dEllipticGauss_RandSource_area(
    float y_rand_min, float y_rand_max, float x_rand_min, float x_rand_max,
    double press_increment, float P_y_mean, float P_y_sigma, float P_x_sigma) {
  float ran_nb = 0;
  float prob = 0;
  int pos_y = 0, pos_x = 0;

  // srand (time (0));

  sourc_Ymin = y_rand_min * P->max_size;
  sourc_Ymax = y_rand_max * P->max_size;
  sourc_Xmin = x_rand_min * P->max_size;
  sourc_Xmax = x_rand_max * P->max_size;

  P_x_mean = 0.5 * P->max_size;

  // cout << "whats going on here in starting of pressure assignment!" << endl;
  do {
    ran_nb =
        rand() / (float)RAND_MAX; // get a float pseudorandom number b/w 0 and 1
    ran_nb =
        (ran_nb - 0.5) * 2.0 *
        (2.0 * P_y_sigma); // get a distribution from plus/minus 2 times sigma
    ran_nb = ran_nb + P_y_mean; // shift it around y_mean

    if (ran_nb > sourc_Ymin && ran_nb < sourc_Ymax) {
      ran_nb = ran_nb + 0.5;
      pos_y = ran_nb;

      break;
    }
  } while (1);

  ran_nb = 0;
  do {
    ran_nb =
        rand() / (float)RAND_MAX; // get a float pseudorandom number b/w 0 and 1
    ran_nb =
        (ran_nb - 0.5) * 2.0 *
        (2.0 * P_x_sigma); // get a distribution from plus/minus 2 times sigma
    ran_nb = ran_nb + P_x_mean; // shift it around x_mean

    if (ran_nb > sourc_Xmin && ran_nb < sourc_Xmax) {
      ran_nb = ran_nb + 0.5;
      pos_x = ran_nb;

      break;
    }
  } while (1);

  prob = exp(-0.5 * ((pos_x - P_x_mean) / P_x_sigma) *
             ((pos_x - P_x_mean) / P_x_sigma)); // exponential part with x
  prob =
      prob * (exp(-0.5 * ((pos_y - P_y_mean) / P_y_sigma) *
                  ((pos_y - P_y_mean) / P_y_sigma))); // exponential part with y
  prob = prob *
         (1.0 / (2.0 * Pi * P_x_sigma *
                 P_y_sigma)); // first part of 2d_Gaussian distribution function
  prob = prob * 2.0 * Pi * P_x_sigma *
         P_y_sigma; // adjust probability to run from 0 to 1.0

  // prob = prob * 100.0;	// get percentage of probability
  // press_increment = (press_increment * prob)/100.0;  // get fraction of given
  // pressure
  press_increment = press_increment * prob;
  // cout << "porb : " << prob << endl;
  // cout << "pos_x : " << pos_x << ", pos_y : " << pos_y << endl;
  // cout << "press_increment : " << press_increment << endl;
  cout << "P->pressure_dist[pos_y][pos_x] : " << P->pressure_dist[pos_y][pos_x]
       << endl;

  if (pos_y >= 5) {
    if (pos_y <= 20) {
      if (pos_x >= 5) {
        if (pos_x <= 45) {
          if (P->pressure_dist[pos_y][pos_x] < press_increment) {
            P->rand_source_PressIncrement(pos_y, pos_x, press_increment);

            cout << "pos_x : " << pos_x << ", pos_y : " << pos_y << endl;
            cout << "press_increment : " << press_increment << endl;
          }
        }
      }
    }
  } else
    ;
}

void Lattice::Pressure_HorSource(int lay_y, int min_x, int max_x,
                                 double press_increment) {
  P->source_hor_PressIncrement(lay_y, min_x, max_x, press_increment);
}

void Lattice::Pressure_VerSource(int lay_x, int min_y, int max_y,
                                 double press_increment) {
  P->source_ver_PressIncrement(lay_x, min_y, max_y, press_increment);
}

void Lattice::Pressure_AreaSource(int min_y, int max_y, int min_x, int max_x,
                                  double press_increment) {
  P->source_area_PressIncrement(min_y, max_y, min_x, max_x, press_increment);
}

void Lattice::Pressure_BuildUp(int press_node_y, int press_node_x,
                               double press_increment) {
  cout << P->pressure_dist[press_node_y][press_node_x] << endl;

  P->pressure_incre_source(press_node_y, press_node_x, press_increment);

  cout << P->pressure_dist[press_node_y][press_node_x] << endl;
}
void Lattice::Pressure_Insert_Random_Nodey(double ymin, double ymax, int res,
                                           double press_increment) {
  float ran_pos_x, ran_pos_y, ran_pressure; // rand number
  double press_increment_final;
  int node_x, node_y;

  srand(time(0));

  ran_pos_x = rand() / (float)RAND_MAX;
  ran_pos_y = rand() / (float)RAND_MAX;
  ran_pressure = rand() / (float)RAND_MAX;

  node_x = int((((ran_pos_x * (ymax - ymin)) + ymin) * (res - 3)) + 2);
  node_y = int((ran_pos_y * (res - 3)) + 2);

  press_increment_final = press_increment * ran_pressure;

  P->pressure_incre_source(node_x, node_y, press_increment_final);
}
void Lattice::Pressure_Insert_Random_Node(int res, double press_increment) {
  float ran_pos_x, ran_pos_y, ran_pressure; // rand number
  double press_increment_final;
  int node_x, node_y;

  srand(time(0));

  ran_pos_x = rand() / (float)RAND_MAX;
  ran_pos_y = rand() / (float)RAND_MAX;
  ran_pressure = rand() / (float)RAND_MAX;

  node_x = int((ran_pos_x * (res - 3)) + 2);
  node_y = int((ran_pos_y * (res - 3)) + 2);
  press_increment_final = press_increment * ran_pressure;

  P->pressure_incre_source(node_x, node_y, press_increment_final);
}

void Lattice::Pressure_BuildUp_twoXsources(int y_press_node,
                                           int x_press_node_diff,
                                           double increment_press) {
  P->Pressure_incre_twoXsources(y_press_node, x_press_node_diff,
                                increment_press);
}

void Lattice::Pressure_BuildUp_twoYsources(int x_press_node,
                                           int y_press_node_diff,
                                           double increment_press) {
  P->Pressure_incre_twoYsources(x_press_node, y_press_node_diff,
                                increment_press);
}

void Lattice::Pressure_build_Up_n_Down(int press_node_y_a, int press_node_x_a,
                                       int press_node_y_b, int press_node_x_b,
                                       double press_increment,
                                       double press_decrement) {
  P->Pressure_source_n_sink(press_node_y_a, press_node_x_a, press_node_y_b,
                            press_node_x_b, press_increment, press_decrement);
}

void Lattice::Pressure_BuildDown_point(int press_node_y, int press_node_x,
                                       double press_decrement) {
  P->pressure_point_sink(press_node_y, press_node_x, press_decrement);
}

void Lattice::Pressure_BuildDown(float y_box_min, float y_box_max,
                                 float x_box_min, float x_box_max,
                                 double press_decrement) {
  P->pressure_sink(y_box_min, y_box_max, x_box_min, x_box_max, press_decrement);
}

void Lattice::Pressure_cont_BuildDown(float y_box_min, float y_box_max,
                                      float x_box_min, float x_box_max,
                                      double press_decrement) {
  P->pressure_sink_cont(y_box_min, y_box_max, x_box_min, x_box_max,
                        press_decrement);
}

void Lattice::Pressure_Boundary_cond(int LboundrySourc, double boundry_val,
                                     int vertical_z_depth) {
  P->boundary_fixed(LboundrySourc, boundry_val, vertical_z_depth);

  P->source_pressure(); // pressure source term: approximated pressure value at
                        // the end of last time step
}

void Lattice::Particle_Fluid_parameters() {
  int i;

  runParticle = &refParticle;
  for (i = 0; i < numParticles; i++) {
    runParticle->velx = runParticle->xpos - runParticle->oldx;
    runParticle->vely = runParticle->ypos - runParticle->oldy;

    runParticle->oldx = runParticle->xpos;
    runParticle->oldy = runParticle->ypos;

    runParticle = runParticle->nextP;
  }

  P->d_solid_frac_n_vel(repBox); // function extnsion _d, dimension (53,53), low
                                 // resolution res100.elle fil

  P->numberDensity_porosity_kappa();
}

void Lattice::Calculate_Pressure_1st_time(int stable, double step_time) {
  double sub_time = step_time / 10.0;
  double l;

  {

    P->alpha_beta_set(l);

    double **P_last = NULL;

    P_last = new double *[P->max_size];

    if (P_last == NULL)
      printf(
          "Error: Dynamic Memory Allocation for %d rows of P_last matrix Fails",
          P->max_size);

    for (i = P->min_size; i < P->max_size; i++) {
      P_last[i] = new double[P->max_size];

      if (P_last[i] == NULL)
        printf("Error: Dynamic Memory Allocation for %d columns of P_last "
               "matrix Fails",
               P->max_size);
    }

    for (i = P->min_size; i < P->max_size; i++) {
      for (j = P->min_size; j < P->max_size; j++)

        P_last[i][j] = P->pressure_dist[i][j];
    }

    diff = 2.0;

    // while(diff > 1.0)
    {
      diff = 0;

      /************************
       * Scheme - test
       *
       * this setupt is new with corrected functions of source term
       * and pressure calculation an dassignment to particles
       *
       * *********************/
      {
        P->dan_matrices_set_1();
        P->dan_multiplication_1st();
        P->test_pressure_with_source(l); // adding the source term (test please)
        P->transpose_trans_P();
        P->dan_invert();
        P->dan_multiplication_inv();

        P->dan_matrices_set_2();
        P->dan_multiplication_1st();
        P->transpose_trans_P();
        P->test_pressure_with_source(l); // adding the source term (test please)
        P->dan_invert();
        P->dan_multiplication_inv();
      }

      for (i = P->min_size + 1; i < P->max_size - 1; i++) {
        for (j = P->min_size + 1; j < P->max_size - 1; j++)

          diff = diff + fabs((P_last[i][j] - P->pressure_dist[i][j]));
      }

      for (i = P->min_size; i < P->max_size; i++) {
        for (j = P->min_size; j < P->max_size; j++)

          P_last[i][j] = P->pressure_dist[i][j];
      }
    }

    for (i = P->min_size; i < P->max_size; i++)
      delete[] P_last[i];

    delete[] P_last;

    if (stable)

      P->d_pressure_particles(repBox);

    else
      ;
  }
}

// void Lattice::Calculate_Pressure()
//{
////int i,j;

////runParticle = &refParticle;
////for (i = 0; i < numParticles; i++)
////{
////runParticle->velx = runParticle->xpos - runParticle->oldx;
////runParticle->vely = runParticle->ypos - runParticle->oldy;

////runParticle->oldx = runParticle->xpos;
////runParticle->oldy = runParticle->ypos;

////runParticle = runParticle->nextP;
////}

////P->solid_frac_n_vel(repBox);

////P->solid_frac_circled(repBox);
////P->solid_frac_n_vel_b(repBox);

////P->numberDensity_porosity_kappa(); 	// Note 1: In case of using
///solid_frac_circled(repBox) along with soild_frac_n_vel(repBox), /
///only for the calculation of rho[i][j], rho_s[i][j] should be turned into
///rho_s_circle[i][j]. / Note 2: while using solid_frac_n_vel_b(repBox), the
///rho_s[i][j] in calculation of all three parameters /         should be turned
///into rho_s_circle[i][j].
// P->alpha_beta_set();

// P->ohne_background_press();

////double **P_last = NULL;

////P_last = new double *[P->max_size];

////if (P_last == NULL)
////printf("Error: Dynamic Memory Allocation for %d rows of P_last matrix
///Fails", P->max_size);

////for(i = P->min_size; i < P->max_size; i++)
////{
////P_last[i] = new double[P->max_size];

////if (P_last[i] == NULL)
////printf("Error: Dynamic Memory Allocation for %d columns of P_last matrix
///Fails", P->max_size);
////}

////for(i = P->min_size; i < P->max_size; i++)
////{
////for (j = P->min_size; j < P->max_size; j++)

////P_last[i][j] = P->pressure_dist[i][j];
////}

////diff = 2.0;

//////while(diff > 1.0)
////{
////diff = 0.0;

// P->dan_matrices_set_1();
// P->dan_multiplication_1st();
// P->test_pressure_with_source();
// P->transpose_trans_P();
// P->dan_invert();
// P->dan_multiplication_inv();

// P->dan_matrices_set_2();
// P->dan_multiplication_1st();
// P->transpose_trans_P();
// P->test_pressure_with_source();
// P->dan_invert();
// P->dan_multiplication_inv();

////for(i = P->min_size+1; i < P->max_size-1; i++)
////{
////for(j = P->min_size+1; j < P->max_size-1; j++)

////diff = diff + fabs((P_last[i][j] - P->pressure_dist[i][j]));
////}

//////cout << diff << endl;

////for(i = P->min_size; i < P->max_size; i++)
////{
////for(j = P->min_size; j < P->max_size; j++)

////P_last[i][j] = P->pressure_dist[i][j];
////}
////}

////for(i = P->min_size; i < P->max_size; i++)
////delete [] P_last[i];

////delete [] P_last;

// P->mit_background_press();

// P->pressure_with_source();

////P->a_pressure_particles(repBox);
////P->circled_pressure_particles(repBox);
////P->b_pressure_particles(repBox);
////P->c_pressure_particles(repBox);
////P->c_dasmal_pressure_particles(repBox);
// P->d_pressure_particles(repBox);
////P->cir_d_pressure_particles(repBox);

////P->pressure_particles(repBox);

// UpdateElle();
//}

void Lattice::DumpStat_Boxes_Press_brkBond(int press_node_y, int press_node_x) {
  FILE *stat;       // file pointer
  int i, p_counter; // counters
  double pressure;  // mean continuum pressure, elapsed time;

  float mid_node_x, mid_node_y;
  float x_box_min, x_box_max, y_box_min, y_box_max;
  float j, x, y, z;

  press_cal_time += cal_time;
  // cout << press_cal_time << endl;

  P->press_nodes(press_node_y, press_node_x);

  mid_node_x = P->press_pos_x - (P->delta_x / 2.0);
  mid_node_y = P->press_pos_y - (P->delta_y / 2.0);

  stat = fopen("Boxes_Press_brkBond.txt", "a");

  j = 1.0;
  x = P->delta_x, y = P->delta_y;

  while (x <= mid_node_x && y <= mid_node_y) {
    x_box_min = mid_node_x - x;
    x_box_max = mid_node_x + x;

    y_box_min = mid_node_y - y;
    y_box_max = mid_node_y + y;

    // cout <<y_box_min << " :::" << endl;

    p_counter = 0;
    pressure = 0;

    runParticle = &refParticle; // start

    for (i = 0; i < numParticles; i++) // look which particles are in the box
    {
      if (runParticle->ypos >= y_box_min) {
        if (runParticle->ypos <= y_box_max) {
          if (runParticle->xpos >= x_box_min) {
            if (runParticle->xpos <= x_box_max) {
              p_counter = p_counter + 1; // count particles

              // set pressure tensor

              pressure =
                  pressure + sqrt((runParticle->F_P_x * runParticle->F_P_x) +
                                  (runParticle->F_P_y * runParticle->F_P_y));
            }
          }
        }
      }
      runParticle = runParticle->nextP;
    }
    pressure =
        pressure /
        p_counter; // and divide by number of particles to get average pressure

    fprintf(stat, "%lf", pressure);
    fprintf(stat, "\t");

    // if(j < 10.0)
    {
      z = 5.0 * j;
      j += 1.0;
    }
    // else
    // break;

    x = P->delta_x * z;
    y = P->delta_y * z;
  }

  fprintf(stat, "%d\t", nbBreak);
  fprintf(stat, "%d\t", nbBreak_last);
  fprintf(stat, "%d\t", nbBreak_current);
  fprintf(stat, "%le", press_cal_time);

  fprintf(stat, "\n");

  fclose(stat); // close file
}

void Lattice::DumpStat_Boxes_Press_brkBond_LHRes_c_d(int press_node_y,
                                                     int press_node_x) {
  FILE *stat;       // file pointer
  int i, p_counter; // counters
  double pressure;  // mean continuum pressure, elapsed time;

  float x_box_min, x_box_max, y_box_min, y_box_max;
  float j, x, y, z;

  press_cal_time += cal_time;

  P->press_nodes(press_node_y, press_node_x);

  stat = fopen("Boxes_Press_brkBond.txt", "a");

  j = 1.0;
  x = P->delta_x, y = P->delta_y;

  while (x <= P->press_pos_x && y <= P->press_pos_y) {
    x_box_min = P->press_pos_x - x;
    x_box_max = P->press_pos_x + x;

    y_box_min = P->press_pos_y - y;
    y_box_max = P->press_pos_y + y;

    p_counter = 0;
    pressure = 0;

    runParticle = &refParticle; // start

    for (i = 0; i < numParticles; i++) // look which particles are in the box
    {
      if (runParticle->ypos >= y_box_min) {
        if (runParticle->ypos <= y_box_max) {
          if (runParticle->xpos >= x_box_min) {
            if (runParticle->xpos <= x_box_max) {
              p_counter = p_counter + 1; // count particles

              pressure =
                  pressure + sqrt((runParticle->F_P_x * runParticle->F_P_x) +
                                  (runParticle->F_P_y * runParticle->F_P_y));
            }
          }
        }
      }
      runParticle = runParticle->nextP;
    }
    pressure =
        pressure /
        p_counter; // and divide by number of particles to get average pressure

    fprintf(stat, "%lf", pressure);
    fprintf(stat, "\t");

    // if(j > 5.0)
    // break;

    z = 5.0 * j;
    j += 1.0;

    x = P->delta_x * z;
    y = P->delta_y * z;
  }

  fprintf(stat, "%d\t", nbBreak);
  fprintf(stat, "%d\t", nbBreak_last);
  fprintf(stat, "%d\t", nbBreak_current);
  fprintf(stat, "%le", press_cal_time);

  fprintf(stat, "\n");

  fclose(stat); // close file
}

void Lattice::DumpStat_UnderSealUpper_Press_brkBond_LHRes_c_d(float sealy_min,
                                                              float sealy_max) {
  FILE *stat;       // file pointer
  int i, p_counter; // counters
  double pressure;  // mean continuum pressure, elapsed time;

  press_cal_time += cal_time;

  stat = fopen("UnderSealUpper_Press_brkBond.txt", "a");

  p_counter = 0;
  pressure = 0;
  runParticle = &refParticle;        // start
  for (i = 0; i < numParticles; i++) // look which particles are in the box
  {
    if (runParticle->ypos > 0.1) {
      if (runParticle->ypos < sealy_min) {
        if (runParticle->xpos > 0.1) {
          if (runParticle->xpos < 0.9) {
            p_counter = p_counter + 1; // count particles

            pressure =
                pressure + sqrt((runParticle->F_P_x * runParticle->F_P_x) +
                                (runParticle->F_P_y * runParticle->F_P_y));
          }
        }
      }
    }
    runParticle = runParticle->nextP;
  }
  pressure =
      pressure /
      p_counter; // and divide by number of particles to get average pressure
  fprintf(stat, "%lf", pressure);
  fprintf(stat, "\t");
  fprintf(stat, "%d\t", LowSeal_bnd_nb);
  fprintf(stat, "%d\t", nbLowSeal_last);
  fprintf(stat, "%d\t", nbLowSeal_current);

  p_counter = 0;
  pressure = 0;
  runParticle = &refParticle;        // start
  for (i = 0; i < numParticles; i++) // look which particles are in the box
  {
    if (runParticle->ypos >= sealy_min) {
      if (runParticle->ypos <= sealy_max) {
        if (runParticle->xpos > 0.1) {
          if (runParticle->xpos < 0.9) {
            p_counter = p_counter + 1; // count particles

            pressure =
                pressure + sqrt((runParticle->F_P_x * runParticle->F_P_x) +
                                (runParticle->F_P_y * runParticle->F_P_y));
          }
        }
      }
    }
    runParticle = runParticle->nextP;
  }
  pressure =
      pressure /
      p_counter; // and divide by number of particles to get average pressure
  fprintf(stat, "%lf", pressure);
  fprintf(stat, "\t");
  fprintf(stat, "%d\t", Seal_bnd_nb);
  fprintf(stat, "%d\t", nbSeal_last);
  fprintf(stat, "%d\t", nbSeal_current);

  p_counter = 0;
  pressure = 0;
  runParticle = &refParticle;        // start
  for (i = 0; i < numParticles; i++) // look which particles are in the box
  {
    if (runParticle->ypos > sealy_max) {
      if (runParticle->ypos < 0.9) {
        if (runParticle->xpos > 0.1) {
          if (runParticle->xpos < 0.9) {
            p_counter = p_counter + 1; // count particles

            pressure =
                pressure + sqrt((runParticle->F_P_x * runParticle->F_P_x) +
                                (runParticle->F_P_y * runParticle->F_P_y));
          }
        }
      }
    }
    runParticle = runParticle->nextP;
  }
  pressure =
      pressure /
      p_counter; // and divide by number of particles to get average pressure
  fprintf(stat, "%lf", pressure);
  fprintf(stat, "\t");
  fprintf(stat, "%d\t", UppSeal_bnd_nb);
  fprintf(stat, "%d\t", nbUppSeal_last);
  fprintf(stat, "%d\t", nbUppSeal_current);

  fprintf(stat, "%le", press_cal_time);

  fprintf(stat, "\n");

  fclose(stat); // close file
}
void Lattice::DumpStat_EdgeSeal_Press_brkBond_LHRes_c_d(float sealy_min,
                                                        float sealy_max) {
  FILE *stat;       // file pointer
  int i, p_counter; // counters
  double pressure;  // mean continuum pressure, elapsed time;

  press_cal_time += cal_time;

  stat = fopen("EdgeSeal_Press_brkBond.txt", "a");

  p_counter = 0;
  pressure = 0;
  runParticle = &refParticle;        // start
  for (i = 0; i < numParticles; i++) // look which particles are in the box
  {
    if (runParticle->ypos >= sealy_min - 0.05) {
      if (runParticle->ypos <= sealy_min + 0.05) {
        if (runParticle->xpos > 0.1) {
          if (runParticle->xpos < 0.9) {
            p_counter = p_counter + 1; // count particles

            pressure =
                pressure + sqrt((runParticle->F_P_x * runParticle->F_P_x) +
                                (runParticle->F_P_y * runParticle->F_P_y));
          }
        }
      }
    }
    runParticle = runParticle->nextP;
  }
  pressure =
      pressure /
      p_counter; // and divide by number of particles to get average pressure
  fprintf(stat, "%lf", pressure);
  fprintf(stat, "\t");
  fprintf(stat, "%d\t", LowSeal_bnd_nb);
  fprintf(stat, "%d\t", nbLowSeal_last);
  fprintf(stat, "%d\t", nbLowSeal_current);

  p_counter = 0;
  pressure = 0;
  runParticle = &refParticle;        // start
  for (i = 0; i < numParticles; i++) // look which particles are in the box
  {
    if (runParticle->ypos >= sealy_max - 0.05) {
      if (runParticle->ypos <= sealy_max + 0.05) {
        if (runParticle->xpos > 0.1) {
          if (runParticle->xpos < 0.9) {
            p_counter = p_counter + 1; // count particles

            pressure =
                pressure + sqrt((runParticle->F_P_x * runParticle->F_P_x) +
                                (runParticle->F_P_y * runParticle->F_P_y));
          }
        }
      }
    }
    runParticle = runParticle->nextP;
  }
  pressure =
      pressure /
      p_counter; // and divide by number of particles to get average pressure
  fprintf(stat, "%lf", pressure);
  fprintf(stat, "\t");
  fprintf(stat, "%d\t", Seal_bnd_nb);
  fprintf(stat, "%d\t", nbSeal_last);
  fprintf(stat, "%d\t", nbSeal_current);

  fprintf(stat, "%le", press_cal_time);

  fprintf(stat, "\n");

  fclose(stat); // close file
}

void Lattice::DumpStat_Boxes_Press_brkBond_HRes_c_d(int press_node_y,
                                                    int press_node_x) {
  FILE *stat;       // file pointer
  int i, p_counter; // counters
  double pressure;  // mean continuum pressure, elapsed time;

  float x_box_min, x_box_max, y_box_min, y_box_max;
  float j, x, y, z;

  press_cal_time += cal_time;

  P->press_nodes(press_node_y, press_node_x);

  stat = fopen("Boxes_Press_brkBond.txt", "a");

  j = 1.0;
  x = P->delta_x, y = P->delta_y;

  while (x <= P->press_pos_x && y <= P->press_pos_y) {
    x_box_min = P->press_pos_x - x;
    x_box_max = P->press_pos_x + x;

    y_box_min = P->press_pos_y - y;
    y_box_max = P->press_pos_y + y;

    // cout << x_box_min << " :::" << x_box_max << endl;

    p_counter = 0;
    pressure = 0;

    runParticle = &refParticle; // start

    for (i = 0; i < numParticles; i++) // look which particles are in the box
    {
      if (runParticle->ypos >= y_box_min) {
        if (runParticle->ypos <= y_box_max) {
          if (runParticle->xpos >= x_box_min) {
            if (runParticle->xpos <= x_box_max) {
              p_counter = p_counter + 1; // count particles

              // set pressure tensor

              pressure =
                  pressure + sqrt((runParticle->F_P_x * runParticle->F_P_x) +
                                  (runParticle->F_P_y * runParticle->F_P_y));
            }
          }
        }
      }
      runParticle = runParticle->nextP;
    }
    pressure =
        pressure /
        p_counter; // and divide by number of particles to get average pressure

    fprintf(stat, "%lf", pressure);
    fprintf(stat, "\t");

    // if(j > 10.0)
    // break;

    z = 10.0 * j;
    j += 1.0;

    x = P->delta_x * z;
    y = P->delta_y * z;
  }

  fprintf(stat, "%d\t", nbBreak);
  fprintf(stat, "%d\t", nbBreak_last);
  fprintf(stat, "%d\t", nbBreak_current);
  fprintf(stat, "%le", press_cal_time);

  fprintf(stat, "\n");

  fclose(stat); // close file
}

void Lattice::DumpStat_Boxes_Stress(int press_node_y, int press_node_x) {
  FILE *stat;                       // file pointer
  int i, p_counter;                 // counters
  double sxx, syy, sxy, smax, smin; // stress tensor and eigenvalues
  double sigma_1, sigma_3, mean,
      differential; // mean stress, differential stress, pressure

  // float mid_node_x, mid_node_y;
  float x_box_min, x_box_max, y_box_min, y_box_max;
  float j, x, y, z;

  P->press_nodes(press_node_y, press_node_x);

  // mid_node_x = P->press_pos_x - (P->delta_x / 2.0);
  // mid_node_y = P->press_pos_y - (P->delta_y / 2.0);

  // cout << P->press_pos_x  << " ::: "  << P->press_pos_y  << endl;

  stat = fopen("Boxes_Stresses.txt", "a"); // open statistic output
                                           // append file
                                           // and dump the data

  j = 1.0;
  x = P->delta_x, y = P->delta_y;

  while (x <= P->press_pos_x && y <= P->press_pos_y) {
    x_box_min = P->press_pos_x - x;
    x_box_max = P->press_pos_x + x;

    y_box_min = P->press_pos_y - y;
    y_box_max = P->press_pos_y + y;

    p_counter = 0;
    sigma_1 = 0;
    sigma_3 = 0;
    mean = 0;
    differential = 0;

    runParticle = &refParticle; // start

    for (i = 0; i < numParticles; i++) // look which particles are in the box
    {
      if (runParticle->ypos > y_box_min) // larger than ymin
      {
        if (runParticle->ypos < y_box_max) // smaller than ymax
        {
          if (runParticle->xpos > x_box_min) // larger than xmin
          {
            if (runParticle->xpos < x_box_max) // smaller than xmax
            {
              p_counter = p_counter + 1; // count particles

              // set stress tensor
              sxx = runParticle->sxx;
              syy = runParticle->syy;
              sxy = runParticle->sxy;

              // eigenvalues
              smax =
                  ((sxx + syy) / 2.0) +
                  sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) + sxy * sxy);
              smin =
                  ((sxx + syy) / 2.0) -
                  sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) + sxy * sxy);

              sigma_1 = sigma_1 + smax;
              sigma_3 = sigma_3 + smin;

              mean = mean + (smax + smin) / 2.0;

              differential = differential + (smax - smin);
            }
          }
        }
      }
      runParticle = runParticle->nextP;
    }

    // and divide by number of particles
    sigma_1 = sigma_1 / p_counter;
    sigma_3 = sigma_3 / p_counter;
    mean = mean / p_counter;
    differential = differential / p_counter;

    fprintf(stat, "%le", sigma_1);
    fprintf(stat, ";%le", sigma_3);
    fprintf(stat, ";%le", mean);
    fprintf(stat, ";%le", differential);

    fprintf(stat, "\t");

    // if(j < 10.0)
    {
      z = 5.0 * j;
      j += 1.0;
    }

    x = P->delta_x * z;
    y = P->delta_y * z;
  }

  fprintf(stat, "\n");

  fclose(stat); // close file
}

void Lattice::DumpStat_Boxes_Stress_LHRes_c_d(int press_node_y,
                                              int press_node_x) {
  FILE *stat;                       // file pointer
  int i, p_counter;                 // counters
  double sxx, syy, sxy, smax, smin; // stress tensor and eigenvalues
  double sigma_1, sigma_3, mean,
      differential; // mean stress, differential stress, pressure

  float x_box_min, x_box_max, y_box_min, y_box_max;
  float j, x, y, z;

  P->press_nodes(press_node_y, press_node_x);

  // cout << P->press_pos_x  << " ::: "  << P->press_pos_y << endl;

  stat = fopen("Boxes_Stresses.txt", "a"); // open statistic output
                                           // append file
                                           // and dump the data

  j = 1.0;
  x = P->delta_x, y = P->delta_y;

  while (x <= P->press_pos_x && y <= P->press_pos_y) {
    x_box_min = P->press_pos_x - x;
    x_box_max = P->press_pos_x + x;

    y_box_min = P->press_pos_y - y;
    y_box_max = P->press_pos_y + y;

    p_counter = 0;
    sigma_1 = 0;
    sigma_3 = 0;
    mean = 0;
    differential = 0;

    runParticle = &refParticle; // start

    for (i = 0; i < numParticles; i++) // look which particles are in the box
    {
      if (runParticle->ypos > y_box_min) // larger than ymin
      {
        if (runParticle->ypos < y_box_max) // smaller than ymax
        {
          if (runParticle->xpos > x_box_min) // larger than xmin
          {
            if (runParticle->xpos < x_box_max) // smaller than xmax
            {
              p_counter = p_counter + 1; // count particles

              // set stress tensor
              sxx = runParticle->sxx;
              syy = runParticle->syy;
              sxy = runParticle->sxy;

              // eigenvalues
              smax =
                  ((sxx + syy) / 2.0) +
                  sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) + sxy * sxy);
              smin =
                  ((sxx + syy) / 2.0) -
                  sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) + sxy * sxy);

              sigma_1 = sigma_1 + smax;
              sigma_3 = sigma_3 + smin;

              mean = mean + (smax + smin) / 2.0;

              differential = differential + (smax - smin);
            }
          }
        }
      }
      runParticle = runParticle->nextP;
    }

    // and divide by number of particles
    sigma_1 = sigma_1 / p_counter;
    sigma_3 = sigma_3 / p_counter;
    mean = mean / p_counter;
    differential = differential / p_counter;

    fprintf(stat, "%le", sigma_1);
    fprintf(stat, ";%le", sigma_3);
    fprintf(stat, ";%le", mean);
    fprintf(stat, ";%le", differential);

    fprintf(stat, "\t");

    // if(j < 10.0)     // this condition together with break; command work with
    // while(1).
    {
      z = 5.0 * j;
      j += 1.0;
    }
    // else
    // break;

    x = P->delta_x * z;
    y = P->delta_y * z;
  }

  fprintf(stat, "\n");

  fclose(stat); // close file
}

void Lattice::DumpStat_UnderSealUpper_Stress_LHRes_c_d(float sealy_min,
                                                       float sealy_max) {
  FILE *stat;                       // file pointer
  int i, p_counter;                 // counters
  double sxx, syy, sxy, smax, smin; // stress tensor and eigenvalues
  double sigma_1, sigma_3, mean,
      differential; // mean stress, differential stress, pressure

  stat = fopen("UnderSealUpper_Stresses.txt", "a"); // open statistic output
                                                    // append file
                                                    // and dump the data

  p_counter = 0;
  sigma_1 = 0;
  sigma_3 = 0;
  mean = 0;
  differential = 0;

  runParticle = &refParticle;        // start
  for (i = 0; i < numParticles; i++) // look which particles are in the box
  {
    if (runParticle->ypos > 0.1) // larger than ymin
    {
      if (runParticle->ypos < sealy_min) // smaller than ymax
      {
        if (runParticle->xpos > 0.1) // larger than xmin
        {
          if (runParticle->xpos < 0.9) // smaller than xmax
          {
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

            sigma_1 = sigma_1 + smax;
            sigma_3 = sigma_3 + smin;

            mean = mean + (smax + smin) / 2.0;

            differential = differential + (smax - smin);
          }
        }
      }
    }
    runParticle = runParticle->nextP;
  }
  sigma_1 = sigma_1 / p_counter;
  sigma_3 = sigma_3 / p_counter;
  mean = mean / p_counter; // and divide by number of particles
  differential = differential / p_counter;

  fprintf(stat, "%le", sigma_1);
  fprintf(stat, ";%le", sigma_3);
  fprintf(stat, ";%le", mean);
  fprintf(stat, ";%le", differential);
  fprintf(stat, "\t");

  p_counter = 0;
  sigma_1 = 0;
  sigma_3 = 0;
  mean = 0;
  differential = 0;

  runParticle = &refParticle;        // start
  for (i = 0; i < numParticles; i++) // look which particles are in the box
  {
    if (runParticle->ypos >= sealy_min) // larger than ymin
    {
      if (runParticle->ypos <= sealy_max) // smaller than ymax
      {
        if (runParticle->xpos > 0.1) // larger than xmin
        {
          if (runParticle->xpos < 0.9) // smaller than xmax
          {
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

            sigma_1 = sigma_1 + smax;
            sigma_3 = sigma_3 + smin;

            mean = mean + (smax + smin) / 2.0;

            differential = differential + (smax - smin);
          }
        }
      }
    }
    runParticle = runParticle->nextP;
  }
  sigma_1 = sigma_1 / p_counter;
  sigma_3 = sigma_3 / p_counter;
  mean = mean / p_counter; // and divide by number of particles
  differential = differential / p_counter;

  fprintf(stat, "%le", sigma_1);
  fprintf(stat, ";%le", sigma_3);
  fprintf(stat, ";%le", mean);
  fprintf(stat, ";%le", differential);
  fprintf(stat, "\t");

  p_counter = 0;
  sigma_1 = 0;
  sigma_3 = 0;
  mean = 0;
  differential = 0;

  runParticle = &refParticle;        // start
  for (i = 0; i < numParticles; i++) // look which particles are in the box
  {
    if (runParticle->ypos > sealy_max) // larger than ymin
    {
      if (runParticle->ypos < 0.9) // smaller than ymax
      {
        if (runParticle->xpos > 0.1) // larger than xmin
        {
          if (runParticle->xpos < 0.9) // smaller than xmax
          {
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

            sigma_1 = sigma_1 + smax;
            sigma_3 = sigma_3 + smin;

            mean = mean + (smax + smin) / 2.0;

            differential = differential + (smax - smin);
          }
        }
      }
    }
    runParticle = runParticle->nextP;
  }
  sigma_1 = sigma_1 / p_counter;
  sigma_3 = sigma_3 / p_counter;
  mean = mean / p_counter; // and divide by number of particles
  differential = differential / p_counter;

  fprintf(stat, "%le", sigma_1);
  fprintf(stat, ";%le", sigma_3);
  fprintf(stat, ";%le", mean);
  fprintf(stat, ";%le", differential);
  fprintf(stat, "\t");

  fprintf(stat, "\n");
  fclose(stat); // close file
}

void Lattice::DumpStat_EdgeSeal_Stress_LHRes_c_d(float sealy_min,
                                                 float sealy_max) {
  FILE *stat;                       // file pointer
  int i, p_counter;                 // counters
  double sxx, syy, sxy, smax, smin; // stress tensor and eigenvalues
  double sigma_1, sigma_3, mean,
      differential; // mean stress, differential stress, pressure

  stat = fopen("EdgeSeal_Stresses.txt", "a"); // open statistic output
                                              // append file
                                              // and dump the data

  p_counter = 0;
  sigma_1 = 0;
  sigma_3 = 0;
  mean = 0;
  differential = 0;

  runParticle = &refParticle;        // start
  for (i = 0; i < numParticles; i++) // look which particles are in the box
  {
    if (runParticle->ypos >= sealy_min - 0.05) // larger than ymin
    {
      if (runParticle->ypos <= sealy_min + 0.05) // smaller than ymax
      {
        if (runParticle->xpos > 0.1) // larger than xmin
        {
          if (runParticle->xpos < 0.9) // smaller than xmax
          {
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

            sigma_1 = sigma_1 + smax;
            sigma_3 = sigma_3 + smin;

            mean = mean + (smax + smin) / 2.0;

            differential = differential + (smax - smin);
          }
        }
      }
    }
    runParticle = runParticle->nextP;
  }
  sigma_1 = sigma_1 / p_counter;
  sigma_3 = sigma_3 / p_counter;
  mean = mean / p_counter; // and divide by number of particles
  differential = differential / p_counter;

  fprintf(stat, "%le", sigma_1);
  fprintf(stat, ";%le", sigma_3);
  fprintf(stat, ";%le", mean);
  fprintf(stat, ";%le", differential);
  fprintf(stat, "\t");

  p_counter = 0;
  sigma_1 = 0;
  sigma_3 = 0;
  mean = 0;
  differential = 0;

  runParticle = &refParticle;        // start
  for (i = 0; i < numParticles; i++) // look which particles are in the box
  {
    if (runParticle->ypos >= sealy_max - 0.05) // larger than ymin
    {
      if (runParticle->ypos <= sealy_max + 0.05) // smaller than ymax
      {
        if (runParticle->xpos > 0.1) // larger than xmin
        {
          if (runParticle->xpos < 0.9) // smaller than xmax
          {
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

            sigma_1 = sigma_1 + smax;
            sigma_3 = sigma_3 + smin;

            mean = mean + (smax + smin) / 2.0;

            differential = differential + (smax - smin);
          }
        }
      }
    }
    runParticle = runParticle->nextP;
  }
  sigma_1 = sigma_1 / p_counter;
  sigma_3 = sigma_3 / p_counter;
  mean = mean / p_counter; // and divide by number of particles
  differential = differential / p_counter;

  fprintf(stat, "%le", sigma_1);
  fprintf(stat, ";%le", sigma_3);
  fprintf(stat, ";%le", mean);
  fprintf(stat, ";%le", differential);
  fprintf(stat, "\t");

  fprintf(stat, "\n");
  fclose(stat); // close file
}

void Lattice::DumpStat_Boxes_Stress_HRes_c_d(int press_node_y,
                                             int press_node_x) {
  FILE *stat;                       // file pointer
  int i, p_counter;                 // counters
  double sxx, syy, sxy, smax, smin; // stress tensor and eigenvalues
  double sigma_1, sigma_3, mean,
      differential; // mean stress, differential stress, pressure

  float x_box_min, x_box_max, y_box_min, y_box_max;
  float j, x, y, z;

  P->press_nodes(press_node_y, press_node_x);

  // cout << mid_node_x  << " ::: "  << mid_node_y << endl;

  stat = fopen("Boxes_Stresses.txt", "a"); // open statistic output
                                           // append file
                                           // and dump the data

  j = 1.0;
  x = P->delta_x, y = P->delta_y;

  while (x <= P->press_pos_x && y <= P->press_pos_y) {
    x_box_min = P->press_pos_x - x;
    x_box_max = P->press_pos_x + x;

    y_box_min = P->press_pos_y - y;
    y_box_max = P->press_pos_y + y;

    p_counter = 0;
    sigma_1 = 0;
    sigma_3 = 0;
    mean = 0;
    differential = 0;

    runParticle = &refParticle; // start

    for (i = 0; i < numParticles; i++) // look which particles are in the box
    {
      if (runParticle->ypos > y_box_min) // larger than ymin
      {
        if (runParticle->ypos < y_box_max) // smaller than ymax
        {
          if (runParticle->xpos > x_box_min) // larger than xmin
          {
            if (runParticle->xpos < x_box_max) // smaller than xmax
            {
              p_counter = p_counter + 1; // count particles

              // set stress tensor
              sxx = runParticle->sxx;
              syy = runParticle->syy;
              sxy = runParticle->sxy;

              // eigenvalues
              smax =
                  ((sxx + syy) / 2.0) +
                  sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) + sxy * sxy);
              smin =
                  ((sxx + syy) / 2.0) -
                  sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) + sxy * sxy);

              sigma_1 = sigma_1 + smax;
              sigma_3 = sigma_3 + smin;

              mean = mean + (smax + smin) / 2.0;

              differential = differential + (smax - smin);
            }
          }
        }
      }
      runParticle = runParticle->nextP;
    }

    // and divide by number of particles
    sigma_1 = sigma_1 / p_counter;
    sigma_3 = sigma_3 / p_counter;
    mean = mean / p_counter;
    differential = differential / p_counter;

    fprintf(stat, "%le", sigma_1);
    fprintf(stat, ";%le", sigma_3);
    fprintf(stat, ";%le", mean);
    fprintf(stat, ";%le", differential);

    fprintf(stat, "\t");

    // if(j < 10.0)
    {
      z = 10.0 * j;
      j += 1.0;
    }
    // else
    // break;

    x = P->delta_x * z;
    y = P->delta_y * z;
  }

  fprintf(stat, "\n");

  fclose(stat); // close file
}

void Lattice::DumpStat_Boxes_Pore_Perm(int press_node_y, int press_node_x) {
  FILE *stat;          // file pointer
  int i, j, p_counter; // counters
  double void_space, permeability;

  int x_box_min, x_box_max, y_box_min, y_box_max;
  int x, y, z;

  stat = fopen("Boxes_PoroPerm.txt", "a");

  x = 1, y = 1, z = 1;

  while (y <= press_node_y && x <= press_node_x) {
    // cout << x << " " << y << endl;

    y_box_min = press_node_y - y;
    if (y == press_node_y)
      y_box_max = press_node_y + y - 1;
    else
      y_box_max = press_node_y + y;

    x_box_min = press_node_x - x;
    if (x == press_node_x)
      x_box_max = press_node_x + x - 1;
    else
      x_box_max = press_node_x + x;

    p_counter = 0;
    void_space = 0;
    permeability = 0;

    for (i = P->min_size; i < P->max_size; i++) {
      for (j = P->min_size; j < P->max_size; j++) {
        if (i >= y_box_min && i <= y_box_max) {
          if (j >= x_box_min && j <= x_box_max) {
            p_counter = p_counter + 1;

            void_space = void_space + P->porosity[i][j];
            permeability = permeability + P->kappa[i][j];
          }
        }
      }
    }

    void_space =
        void_space /
        p_counter; // and divide by number of particles to get average porosity
    permeability = permeability / p_counter;

    fprintf(stat, "%lf", void_space);
    fprintf(stat, ";");
    fprintf(stat, "%le\t", permeability);

    // if(z < 6)
    {
      x = 5 * z;
      y = 5 * z;
      z += 2;
    }
    // else
    // break;
  }

  fprintf(stat, "\n");
  fclose(stat); // close file
}

void Lattice::DumpStat_Boxes_Pore_Perm_LHRes_c_d(int press_node_y,
                                                 int press_node_x) {
  FILE *stat;          // file pointer
  int i, j, p_counter; // counters
  double void_space, permeability;

  int x_box_min, x_box_max, y_box_min, y_box_max;
  int x, y, z;

  stat = fopen("Boxes_PoroPerm.txt", "a");

  x = 1, y = 1, z = 1;

  while (y <= press_node_y && x <= press_node_x) {
    // cout << x << " " << y << endl;

    y_box_min = press_node_y - y;
    y_box_max = press_node_y + y;

    x_box_min = press_node_x - x;
    x_box_max = press_node_x + x;

    p_counter = 0;
    void_space = 0;
    permeability = 0;

    for (i = P->min_size; i < P->max_size; i++) {
      for (j = P->min_size; j < P->max_size; j++) {
        if (i >= y_box_min && i <= y_box_max) {
          if (j >= x_box_min && j <= x_box_max) {
            p_counter = p_counter + 1;

            void_space = void_space + P->porosity[i][j];
            permeability = permeability + P->kappa[i][j];
          }
        }
      }
    }

    void_space =
        void_space /
        p_counter; // and divide by number of particles to get average porosity
    permeability = permeability / p_counter;

    fprintf(stat, "%le", void_space);
    fprintf(stat, " ; ");
    fprintf(stat, "%le\t", permeability);

    // if(z <= 5)
    {
      x = 5 * z;
      y = 5 * z;
      z += 1;
    }
    // else
    // break;
  }

  fprintf(stat, "\n");
  fclose(stat); // close file
}

void Lattice::DumpStat_UnderSealUpper_Pore_Perm_LHRes_c_d(float sealy_min,
                                                          float sealy_max) {
  FILE *stat;          // file pointer
  int i, j, p_counter; // counters
  double void_space, permeability;
  int y_box_min, y_box_max;

  stat = fopen("UnderSealUpper_Pore_Perm.txt", "a");

  p_counter = 0;
  void_space = 0;
  permeability = 0;
  y_box_max = sealy_min * P->max_size;

  for (i = P->min_size + 5; i < y_box_max; i++) {
    for (j = P->min_size + 5; j < P->max_size - 5; j++) {
      p_counter = p_counter + 1;

      void_space = void_space + P->porosity[i][j];
      permeability = permeability + P->kappa[i][j];
    }
  }
  void_space =
      void_space /
      p_counter; // and divide by number of particles to get average porosity
  permeability = permeability / p_counter;
  fprintf(stat, "%le", void_space);
  fprintf(stat, " ; ");
  fprintf(stat, "%le\t", permeability);

  p_counter = 0;
  y_box_min = sealy_min * P->max_size;
  y_box_max = sealy_max * P->max_size;

  for (i = y_box_min; i <= y_box_max; i++) {
    for (j = P->min_size + 5; j < P->max_size - 5; j++) {
      p_counter = p_counter + 1;

      void_space = void_space + P->porosity[i][j];
      permeability = permeability + P->kappa[i][j];
    }
  }
  void_space =
      void_space /
      p_counter; // and divide by number of particles to get average porosity
  permeability = permeability / p_counter;
  fprintf(stat, "%le", void_space);
  fprintf(stat, " ; ");
  fprintf(stat, "%le\t", permeability);

  p_counter = 0;
  for (i = P->min_size + 5; i <= y_box_max; i++) {
    for (j = P->min_size + 5; j < P->max_size - 5; j++) {
      p_counter = p_counter + 1;

      void_space = void_space + P->porosity[i][j];
      permeability = permeability + P->kappa[i][j];
    }
  }
  void_space =
      void_space /
      p_counter; // and divide by number of particles to get average porosity
  permeability = permeability / p_counter;
  fprintf(stat, "%le", void_space);
  fprintf(stat, " ; ");
  fprintf(stat, "%le\t", permeability);

  p_counter = 0;
  for (i = y_box_max + 1; i < P->max_size - 5; i++) {
    for (j = P->min_size + 5; j < P->max_size - 5; j++) {
      p_counter = p_counter + 1;

      void_space = void_space + P->porosity[i][j];
      permeability = permeability + P->kappa[i][j];
    }
  }
  void_space =
      void_space /
      p_counter; // and divide by number of particles to get average porosity
  permeability = permeability / p_counter;
  fprintf(stat, "%le", void_space);
  fprintf(stat, " ; ");
  fprintf(stat, "%le", permeability);

  fprintf(stat, "\n");
  fclose(stat); // close file
}

void Lattice::DumpStat_EdgeSeal_Pore_Perm_LHRes_c_d(float sealy_min,
                                                    float sealy_max) {
  FILE *stat;          // file pointer
  int i, j, p_counter; // counters
  double void_space, permeability;
  int y_box_min, y_box_max;

  stat = fopen("EdgeSeal_Pore_Perm.txt", "a");

  p_counter = 0;
  void_space = 0;
  permeability = 0;
  y_box_min = (sealy_min - 0.05) * P->max_size;
  y_box_max = (sealy_min + 0.05) * P->max_size;

  for (i = y_box_min; i <= y_box_max; i++) {
    for (j = P->min_size + 5; j < P->max_size - 5; j++) {
      p_counter = p_counter + 1;

      void_space = void_space + P->porosity[i][j];
      permeability = permeability + P->kappa[i][j];
    }
  }
  void_space =
      void_space /
      p_counter; // and divide by number of particles to get average porosity
  permeability = permeability / p_counter;
  fprintf(stat, "%le", void_space);
  fprintf(stat, " ; ");
  fprintf(stat, "%le\t", permeability);

  p_counter = 0;
  y_box_min = (sealy_max - 0.05) * P->max_size;
  y_box_max = (sealy_max + 0.05) * P->max_size;

  for (i = y_box_min; i <= y_box_max; i++) {
    for (j = P->min_size + 5; j < P->max_size - 5; j++) {
      p_counter = p_counter + 1;

      void_space = void_space + P->porosity[i][j];
      permeability = permeability + P->kappa[i][j];
    }
  }
  void_space =
      void_space /
      p_counter; // and divide by number of particles to get average porosity
  permeability = permeability / p_counter;
  fprintf(stat, "%le", void_space);
  fprintf(stat, " ; ");
  fprintf(stat, "%le\t", permeability);

  fprintf(stat, "\n");
  fclose(stat); // close file
}

void Lattice::DumpStat_Boxes_Pore_Perm_HRes_c_d(int press_node_y,
                                                int press_node_x) {
  FILE *stat;          // file pointer
  int i, j, p_counter; // counters
  double void_space, permeability;

  int x_box_min, x_box_max, y_box_min, y_box_max;
  int x, y, z;

  stat = fopen("Boxes_PoroPerm.txt", "a");

  x = 1, y = 1, z = 1;

  while (y <= press_node_y && x <= press_node_x) {
    y_box_min = press_node_y - y;
    y_box_max = press_node_y + y;

    x_box_min = press_node_x - x;
    x_box_max = press_node_x + x;

    // cout << x_box_min << " " << x_box_max << endl;

    p_counter = 0;
    void_space = 0;
    permeability = 0;

    for (i = P->min_size; i < P->max_size; i++) {
      for (j = P->min_size; j < P->max_size; j++) {
        if (i >= y_box_min && i <= y_box_max) {
          if (j >= x_box_min && j <= x_box_max) {
            p_counter = p_counter + 1;

            void_space = void_space + P->porosity[i][j];
            permeability = permeability + P->kappa[i][j];
          }
        }
      }
    }

    void_space =
        void_space /
        p_counter; // and divide by number of particles to get average porosity
    permeability = permeability / p_counter;

    fprintf(stat, "%le", void_space);
    fprintf(stat, " ; ");
    fprintf(stat, "%le\t", permeability);

    // if(z <= 10)
    {
      x = 10 * z;
      y = 10 * z;
      z += 1;
    }
    // else
    // break;
  }

  fprintf(stat, "\n");
  fclose(stat); // close file
}

// void Lattice::DumpStat_Press_stress_pore_perm(float y_box_min, float
// y_box_max, float x_box_min, float x_box_max)
//{
// FILE *scaledstat;	// file pointer
// int i, p_counter;	// counters
// double sxx, syy, sxy, smax, smin;	// stress tensor and eigenvalues
// double sigma_1, sigma_3, mean, differential;	// mean stress, differential
// stress, pressure double pressure;	// mean pressure double void_space,
// permeability;

// int  j, q_counter;	// counters
// int y_min, y_max, x_min, x_max;

// SymmetricMatrix sigma(2);
// Matrix V;
// DiagonalMatrix D;
// double m_sxx, m_syy, m_sxy;
// double sigma_V, sigma_H;
// double sigmaV_vec[2], sigmaH_vec[2];

// p_counter= 0;
// pressure = 0;
// sigma_1 = 0;
// sigma_3 = 0;
// mean = 0;
// differential = 0;

// runParticle = &refParticle;	// start

// for (i = 0; i < numParticles; i++)	// look which particles are in the box
//{
// if ( runParticle->ypos > y_box_min )	// larger than x_press_node
//{
// if (runParticle->ypos < y_box_max)	// smaller x_press_node
//{
// if (runParticle->xpos > x_box_min)	// larger than y_press_node
//{
// if (runParticle->xpos < x_box_max)	// smaller than y_press_node
//{
// p_counter = p_counter + 1;	// count particles

//// set stress tensor
// sxx = runParticle->sxx;
// syy = runParticle->syy;
// sxy = runParticle->sxy;

//// set tensors for eigen values and vectors
// m_sxx = m_sxx + sxx;
// m_syy = m_syy + syy;
// m_sxy = m_sxy + sxy;

//// set max and min stresses
// smax = ((sxx + syy) / 2.0) + sqrt (((sxx - syy) / 2.0) * ((sxx - syy) / 2.0)
// + sxy * sxy); smin = ((sxx + syy) / 2.0) - sqrt (((sxx - syy) / 2.0) * ((sxx
// - syy) / 2.0) + sxy * sxy);

// sigma_1 = sigma_1 + smax;
// sigma_3 = sigma_3 + smin;

// mean = mean + (smax + smin) / 2.0;

// differential = differential + (smax - smin);

// pressure = pressure +
// sqrt((runParticle->F_P_x*runParticle->F_P_x)+(runParticle->F_P_y*runParticle->F_P_y));
//}
//}
//}
//}
// runParticle = runParticle->nextP;
//}

////sigma.element(0,0) = m_sxx;
////sigma.element(1,1) = m_syy;
////sigma.element(1,0) = m_sxy;

////EigenValues(sigma, D, V);

////// get sigmax
////sigma_V = D.element(1);
////sigma_H = D.element(0);

////// read columns and put into sigmax_vec
////for (int i=0; i<2; i++)
////sigmaV_vec[i] = V.element (i,1);

////for (int i=0; i<2; i++)
////sigmaH_vec[i] = V.element (i,0);

//// and divide by number of particles to get average values of mean stress,
///diff stress and pressure
// sigma_1 = sigma_1 / p_counter;
// sigma_3 = sigma_3 / p_counter;
// mean = mean / p_counter;
// differential = differential / p_counter;
// pressure = pressure / p_counter;

// q_counter= 0;
// void_space = 0;
// permeability = 0;

// y_min = y_box_min * P->max_size;
// y_max = y_box_max * P->max_size;
// x_min = x_box_min * P->max_size;
// x_max = x_box_max * P->max_size;

// for(i = y_min; i <= y_max; i ++)
//{
// for(j = x_min; j <= x_max; j++)
//{
// q_counter = q_counter + 1;

// void_space = void_space + P->porosity[i][j];
// permeability = permeability + P->kappa[i][j];
//}
//}

// void_space = void_space / q_counter;	// and divide by number of particles to
// get average porosity permeability = permeability / q_counter;

// scaledstat = fopen ("pressbox_press_mean_diff_pore_perm.txt", "a");	// open
// statistic output
//// append file and dump the data
// fprintf (scaledstat, "%lf\t", pressure);
// fprintf (scaledstat, "%le\t", sigma_1);
// fprintf (scaledstat, "%le\t", sigma_3);
// fprintf (scaledstat, "%le\t", mean);
// fprintf (scaledstat, "%le\t", differential);
// fprintf (scaledstat, "%lf", void_space);
// fprintf (scaledstat, "\t");
// fprintf (scaledstat, "%le\n", permeability);

// fclose (scaledstat);		// close file
//}

void Lattice::DumpStat_TwoSourc_Pres_stres_por_perm_BrkBnd(float y_ModHalf) {
  FILE *stat;                     // file pointer
  int i, j, p_counter, q_counter; // counters

  double sxx, syy, sxy, smax, smin; // stress tensor and eigenvalues
  double sigma_1, sigma_3, mean,
      differential; // mean stress, differential stress, pressure
  double pressure;  // mean pressure
  double void_space, permeability;

  int y_max, x_min, x_max;

  stat = fopen("TwoSourc_Pres_stres_por_perm_BrkBnd.txt", "a");

  y_max = y_ModHalf * P->max_size;
  x_min = y_max / 2;
  x_max = x_min + (y_max / 2);

  p_counter = 0;
  pressure = 0;
  sigma_1 = 0;
  sigma_3 = 0;
  mean = 0;
  differential = 0;
  runParticle = &refParticle;        // start
  for (i = 0; i < numParticles; i++) // look which particles are in the box
  {
    if (runParticle->ypos >= 0.05) {
      if (runParticle->ypos < y_ModHalf) {
        if (runParticle->xpos > x_min / P->max_size) {
          if (runParticle->xpos < x_max / P->max_size) {
            p_counter = p_counter + 1; // count particles
            // set stress tensor
            sxx = runParticle->sxx;
            syy = runParticle->syy;
            sxy = runParticle->sxy;
            // set max and min stresses
            smax = ((sxx + syy) / 2.0) +
                   sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) + sxy * sxy);
            smin = ((sxx + syy) / 2.0) -
                   sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) + sxy * sxy);

            sigma_1 = sigma_1 + smax;
            sigma_3 = sigma_3 + smin;

            mean = mean + (smax + smin) / 2.0;
            differential = differential + (smax - smin);
            pressure =
                pressure + sqrt((runParticle->F_P_x * runParticle->F_P_x) +
                                (runParticle->F_P_y * runParticle->F_P_y));
          }
        }
      }
    }
    runParticle = runParticle->nextP;
  }
  // and divide by number of particles to get average values of mean stress,
  // diff stress and pressure
  pressure = pressure / p_counter;
  sigma_1 = sigma_1 / p_counter;
  sigma_3 = sigma_3 / p_counter;
  mean = mean / p_counter;
  differential = differential / p_counter;

  q_counter = 0;
  void_space = 0;
  permeability = 0;
  for (i = 2; i <= y_max; i++) {
    for (j = x_min; j <= x_max; j++) {
      q_counter = q_counter + 1;
      void_space = void_space + P->porosity[i][j];
      permeability = permeability + P->kappa[i][j];
    }
  }
  void_space =
      void_space /
      q_counter; // and divide by number of particles to get average porosity
  permeability = permeability / q_counter;

  fprintf(stat, "%lf\t", pressure);
  fprintf(stat, "%le\t", sigma_1);
  fprintf(stat, "%le\t", sigma_3);
  fprintf(stat, "%le\t", mean);
  fprintf(stat, "%le\t", differential);
  fprintf(stat, "%lf", void_space);
  fprintf(stat, "\t");
  fprintf(stat, "%le\n", permeability);
  fprintf(stat, "%d\t", ModLowr_bnd_nb);
  fprintf(stat, "%d\t", nbModLowr_last);
  fprintf(stat, "%d\t", nbModLowr_current);

  p_counter = 0;
  pressure = 0;
  sigma_1 = 0;
  sigma_3 = 0;
  mean = 0;
  differential = 0;
  runParticle = &refParticle;        // start
  for (i = 0; i < numParticles; i++) // look which particles are in the box
  {
    if (runParticle->ypos >= y_ModHalf) {
      if (runParticle->ypos < 0.95) {
        if (runParticle->xpos > x_min / P->max_size) {
          if (runParticle->xpos < x_max / P->max_size) {
            p_counter = p_counter + 1; // count particles
            // set stress tensor
            sxx = runParticle->sxx;
            syy = runParticle->syy;
            sxy = runParticle->sxy;
            // set max and min stresses
            smax = ((sxx + syy) / 2.0) +
                   sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) + sxy * sxy);
            smin = ((sxx + syy) / 2.0) -
                   sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) + sxy * sxy);

            sigma_1 = sigma_1 + smax;
            sigma_3 = sigma_3 + smin;

            mean = mean + (smax + smin) / 2.0;
            differential = differential + (smax - smin);
            pressure =
                pressure + sqrt((runParticle->F_P_x * runParticle->F_P_x) +
                                (runParticle->F_P_y * runParticle->F_P_y));
          }
        }
      }
    }
    runParticle = runParticle->nextP;
  }
  // and divide by number of particles to get average values of mean stress,
  // diff stress and pressure
  pressure = pressure / p_counter;
  sigma_1 = sigma_1 / p_counter;
  sigma_3 = sigma_3 / p_counter;
  mean = mean / p_counter;
  differential = differential / p_counter;

  q_counter = 0;
  void_space = 0;
  permeability = 0;
  for (i = y_max; i <= P->max_size - 2; i++) {
    for (j = x_min; j <= x_max; j++) {
      q_counter = q_counter + 1;
      void_space = void_space + P->porosity[i][j];
      permeability = permeability + P->kappa[i][j];
    }
  }
  void_space =
      void_space /
      q_counter; // and divide by number of particles to get average porosity
  permeability = permeability / q_counter;

  fprintf(stat, "%lf\t", pressure);
  fprintf(stat, "%le\t", sigma_1);
  fprintf(stat, "%le\t", sigma_3);
  fprintf(stat, "%le\t", mean);
  fprintf(stat, "%le\t", differential);
  fprintf(stat, "%lf", void_space);
  fprintf(stat, "\t");
  fprintf(stat, "%le\n", permeability);
  fprintf(stat, "%d\t", ModUppr_bnd_nb);
  fprintf(stat, "%d\t", nbModUppr_last);
  fprintf(stat, "%d\t", ModUppr_current);

  fprintf(stat, "\n");

  fclose(stat); // close file
}

void Lattice::DumpStat_PressMDMxMnStress_LHRes_d(int press_node_y,
                                                 int press_node_x) {
  FILE *statStress;                 // file pointer
  int i, j, k, l, m, p_counter;     // counters
  double pressure;                  // mean continuum pressure, elapsed time;
  double sxx, syy, sxy, smax, smin; // stress tensor and eigenvalues
  double sigma_1, sigma_3, mean,
      differential; // mean stress, differential stress, pressure

  float mid_node_x_ij, mid_node_y_ij;
  float x_ij_box_min, x_ij_box_max, y_ij_box_min, y_ij_box_max;

  statStress = fopen("Local_Stresses_d.txt", "a");

  // press_cal_time += cal_time;

  // for(i = 6; i < P->max_size-6; i +=5)
  i = 1;
  while (i < P->max_size) {
    for (j = 1; j < P->max_size; j += 5) {
      if ((i > press_node_y - 5 && i < press_node_y + 5) &&
          (j > press_node_x - 5 && j < press_node_x + 5)) {
        for (k = i - 4; k <= i + 4; k += 2) {
          for (l = j - 4; l <= j + 4; l += 2) {
            if (k != press_node_y && l != press_node_x) {
              P->press_nodes(k, l);
              mid_node_x_ij = P->press_pos_x;
              mid_node_y_ij = P->press_pos_y;

              y_ij_box_min = mid_node_y_ij - P->delta_y / 2.0;
              y_ij_box_max = mid_node_y_ij + P->delta_y / 2.0;

              x_ij_box_min = mid_node_x_ij - P->delta_x / 2.0;
              x_ij_box_max = mid_node_x_ij + P->delta_x / 2.0;

              p_counter = 0;
              pressure = 0;
              sigma_1 = 0;
              sigma_3 = 0;
              mean = 0;
              differential = 0;

              runParticle = &refParticle; // start

              for (m = 0; m < numParticles;
                   m++) // look which particles are in the box
              {
                if (runParticle->ypos >= y_ij_box_min) {
                  if (runParticle->ypos <= y_ij_box_max) {
                    if (runParticle->xpos >= x_ij_box_min) {
                      if (runParticle->xpos <= x_ij_box_max) {
                        p_counter = p_counter + 1; // count particles

                        // set local pressure
                        pressure =
                            pressure +
                            sqrt((runParticle->F_P_x * runParticle->F_P_x) +
                                 (runParticle->F_P_y * runParticle->F_P_y));

                        // set stress tensor
                        sxx = runParticle->sxx;
                        syy = runParticle->syy;
                        sxy = runParticle->sxy;

                        // set pressure tensor

                        smax = ((sxx + syy) / 2.0) +
                               sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                                    sxy * sxy);
                        smin = ((sxx + syy) / 2.0) -
                               sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                                    sxy * sxy);

                        mean = mean + (smax + smin) / 2.0;

                        differential = differential + (smax - smin);
                        sigma_1 = sigma_1 + smax;
                        sigma_3 = sigma_3 + smin;
                      }
                    }
                  }
                }
                runParticle = runParticle->nextP;
              }
              pressure =
                  pressure / p_counter; // and divide by number of particles to
                                        // get average pressure
              sigma_1 = sigma_1 / p_counter;
              sigma_3 = sigma_3 / p_counter;
              mean = mean / p_counter;
              differential = differential / p_counter;

              fprintf(statStress, "%f\t", mid_node_x_ij);
              fprintf(statStress, "%f\t", mid_node_y_ij);
              fprintf(statStress, "%le\t", pressure);
              fprintf(statStress, "%le\t", sigma_1);
              fprintf(statStress, "%le\t", sigma_3);
              fprintf(statStress, "%le\t", mean);
              fprintf(statStress, "%le", differential);

              fprintf(statStress, " \n");
            }
          }
        }
      }

      else {
        P->press_nodes(i, j);
        mid_node_x_ij = P->press_pos_x;
        mid_node_y_ij = P->press_pos_y;

        y_ij_box_min = mid_node_y_ij - P->delta_y / 2.0;
        y_ij_box_max = mid_node_y_ij + P->delta_y / 2.0;

        x_ij_box_min = mid_node_x_ij - P->delta_x / 2.0;
        x_ij_box_max = mid_node_x_ij + P->delta_x / 2.0;

        p_counter = 0;
        pressure = 0;
        sigma_1 = 0;
        sigma_3 = 0;
        mean = 0;
        differential = 0;

        runParticle = &refParticle; // start

        for (k = 0; k < numParticles;
             k++) // look which particles are in the box
        {
          if (runParticle->ypos >= y_ij_box_min) {
            if (runParticle->ypos <= y_ij_box_max) {
              if (runParticle->xpos >= x_ij_box_min) {
                if (runParticle->xpos <= x_ij_box_max) {
                  p_counter = p_counter + 1; // count particles

                  // set local pressure
                  pressure = pressure +
                             sqrt((runParticle->F_P_x * runParticle->F_P_x) +
                                  (runParticle->F_P_y * runParticle->F_P_y));

                  // set stress tensor
                  sxx = runParticle->sxx;
                  syy = runParticle->syy;
                  sxy = runParticle->sxy;

                  // set pressure tensor

                  smax = ((sxx + syy) / 2.0) +
                         sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                              sxy * sxy);
                  smin = ((sxx + syy) / 2.0) -
                         sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                              sxy * sxy);

                  mean = mean + (smax + smin) / 2.0;

                  differential = differential + (smax - smin);
                  sigma_1 = sigma_1 + smax;
                  sigma_3 = sigma_3 + smin;
                }
              }
            }
          }
          runParticle = runParticle->nextP;
        }
        pressure = pressure / p_counter; // and divide by number of particles to
                                         // get average pressure
        sigma_1 = sigma_1 / p_counter;
        sigma_3 = sigma_3 / p_counter;
        mean = mean / p_counter;
        differential = differential / p_counter;

        fprintf(statStress, "%f\t", mid_node_x_ij);
        fprintf(statStress, "%f\t", mid_node_y_ij);
        fprintf(statStress, "%le\t", pressure);
        fprintf(statStress, "%le\t", sigma_1);
        fprintf(statStress, "%le\t", sigma_3);
        fprintf(statStress, "%le\t", mean);
        fprintf(statStress, "%le", differential);

        fprintf(statStress, " \n");
      }
    }

    i += 5;
  }

  fprintf(statStress, "\n");
  fclose(statStress); // close file
}

void Lattice::DumpStat_PressMDMxMnStress_HRes_d(int press_node_y,
                                                int press_node_x) {
  FILE *statStress;                 // file pointer
  int i, j, k, l, m, p_counter;     // counters
  double pressure;                  // mean continuum pressure, elapsed time;
  double sxx, syy, sxy, smax, smin; // stress tensor and eigenvalues
  double sigma_1, sigma_3, mean,
      differential; // mean stress, differential stress, pressure

  float mid_node_x_ij, mid_node_y_ij;
  float x_ij_box_min, x_ij_box_max, y_ij_box_min, y_ij_box_max;

  statStress = fopen("Local_Stresses_d.txt", "a");

  // press_cal_time += cal_time;

  // for(i = 6; i < P->max_size-6; i +=5)
  i = 1;
  while (i < P->max_size) {
    for (j = 1; j < P->max_size; j += 10) {
      if ((i > press_node_y - 10 && i < press_node_y + 10) &&
          (j > press_node_x - 10 && j < press_node_x + 10)) {
        for (k = i - 8; k <= i + 8; k += 4) {
          for (l = j - 8; l <= j + 8; l += 4) {
            // if (k != press_node_y && l != press_node_x)
            {
              P->press_nodes(k, l);
              mid_node_x_ij = P->press_pos_x;
              mid_node_y_ij = P->press_pos_y;

              y_ij_box_min = mid_node_y_ij - P->delta_y / 2.0;
              y_ij_box_max = mid_node_y_ij + P->delta_y / 2.0;

              x_ij_box_min = mid_node_x_ij - P->delta_x / 2.0;
              x_ij_box_max = mid_node_x_ij + P->delta_x / 2.0;

              p_counter = 0;
              pressure = 0;
              sigma_1 = 0;
              sigma_3 = 0;
              mean = 0;
              differential = 0;

              runParticle = &refParticle; // start

              for (m = 0; m < numParticles;
                   m++) // look which particles are in the box
              {
                if (runParticle->ypos >= y_ij_box_min) {
                  if (runParticle->ypos <= y_ij_box_max) {
                    if (runParticle->xpos >= x_ij_box_min) {
                      if (runParticle->xpos <= x_ij_box_max) {
                        p_counter = p_counter + 1; // count particles

                        // set local pressure
                        pressure =
                            pressure +
                            sqrt((runParticle->F_P_x * runParticle->F_P_x) +
                                 (runParticle->F_P_y * runParticle->F_P_y));

                        // set stress tensor
                        sxx = runParticle->sxx;
                        syy = runParticle->syy;
                        sxy = runParticle->sxy;

                        // set pressure tensor

                        smax = ((sxx + syy) / 2.0) +
                               sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                                    sxy * sxy);
                        smin = ((sxx + syy) / 2.0) -
                               sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                                    sxy * sxy);

                        mean = mean + (smax + smin) / 2.0;

                        differential = differential + (smax - smin);
                        sigma_1 = sigma_1 + smax;
                        sigma_3 = sigma_3 + smin;
                      }
                    }
                  }
                }
                runParticle = runParticle->nextP;
              }
              pressure =
                  pressure / p_counter; // and divide by number of particles to
                                        // get average pressure
              sigma_1 = sigma_1 / p_counter;
              sigma_3 = sigma_3 / p_counter;
              mean = mean / p_counter;
              differential = differential / p_counter;

              fprintf(statStress, "%f\t", mid_node_x_ij);
              fprintf(statStress, "%f\t", mid_node_y_ij);
              fprintf(statStress, "%le\t", pressure);
              fprintf(statStress, "%le\t", sigma_1);
              fprintf(statStress, "%le\t", sigma_3);
              fprintf(statStress, "%le\t", mean);
              fprintf(statStress, "%le", differential);

              fprintf(statStress, " \n");
            }
          }
        }
      }

      else {
        P->press_nodes(i, j);
        mid_node_x_ij = P->press_pos_x;
        mid_node_y_ij = P->press_pos_y;

        y_ij_box_min = mid_node_y_ij - P->delta_y / 2.0;
        y_ij_box_max = mid_node_y_ij + P->delta_y / 2.0;

        x_ij_box_min = mid_node_x_ij - P->delta_x / 2.0;
        x_ij_box_max = mid_node_x_ij + P->delta_x / 2.0;

        p_counter = 0;
        pressure = 0;
        sigma_1 = 0;
        sigma_3 = 0;
        mean = 0;
        differential = 0;

        runParticle = &refParticle; // start

        for (k = 0; k < numParticles;
             k++) // look which particles are in the box
        {
          if (runParticle->ypos >= y_ij_box_min) {
            if (runParticle->ypos <= y_ij_box_max) {
              if (runParticle->xpos >= x_ij_box_min) {
                if (runParticle->xpos <= x_ij_box_max) {
                  p_counter = p_counter + 1; // count particles

                  // set local pressure
                  pressure = pressure +
                             sqrt((runParticle->F_P_x * runParticle->F_P_x) +
                                  (runParticle->F_P_y * runParticle->F_P_y));

                  // set stress tensor
                  sxx = runParticle->sxx;
                  syy = runParticle->syy;
                  sxy = runParticle->sxy;

                  // set pressure tensor

                  smax = ((sxx + syy) / 2.0) +
                         sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                              sxy * sxy);
                  smin = ((sxx + syy) / 2.0) -
                         sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                              sxy * sxy);

                  mean = mean + (smax + smin) / 2.0;

                  differential = differential + (smax - smin);
                  sigma_1 = sigma_1 + smax;
                  sigma_3 = sigma_3 + smin;
                }
              }
            }
          }
          runParticle = runParticle->nextP;
        }
        pressure = pressure / p_counter; // and divide by number of particles to
                                         // get average pressure
        sigma_1 = sigma_1 / p_counter;
        sigma_3 = sigma_3 / p_counter;
        mean = mean / p_counter;
        differential = differential / p_counter;

        fprintf(statStress, "%f\t", mid_node_x_ij);
        fprintf(statStress, "%f\t", mid_node_y_ij);
        fprintf(statStress, "%le\t", pressure);
        fprintf(statStress, "%le\t", sigma_1);
        fprintf(statStress, "%le\t", sigma_3);
        fprintf(statStress, "%le\t", mean);
        fprintf(statStress, "%le", differential);

        fprintf(statStress, " \n");
      }
    }

    i += 10;
  }

  fprintf(statStress, "\n");
  fclose(statStress); // close file
}

void Lattice::DumpStat_PressMDMxMnStress_LHRes_c(int press_node_y,
                                                 int press_node_x) {
  FILE *statStress;                 // file pointer
  int i, j, k, l, m, p_counter;     // counters
  double pressure;                  // mean continuum pressure, elapsed time;
  double sxx, syy, sxy, smax, smin; // stress tensor and eigenvalues
  double sigma_1, sigma_3, mean,
      differential; // mean stress, differential stress, pressure

  float mid_node_x_ij, mid_node_y_ij;
  float x_ij_box_min, x_ij_box_max, y_ij_box_min, y_ij_box_max;

  statStress = fopen("Local_Stresses_c.txt", "a");

  // press_cal_time += cal_time;

  // for(i = 5; i < P->max_size-5; i +=5)
  i = 0;
  while (i < P->max_size) {
    for (j = 0; j < P->max_size; j += 5) {
      if ((i > press_node_y - 5 && i < press_node_y + 5) &&
          (j > press_node_x - 5 && j < press_node_x + 5)) {
        for (k = i - 4; k <= i + 4; k += 2) {
          for (l = j - 4; l <= j + 4; l += 2) {
            if (k != press_node_y && l != press_node_x) {
              P->press_nodes(k, l);
              mid_node_x_ij = P->press_pos_x;
              mid_node_y_ij = P->press_pos_y;

              y_ij_box_min = mid_node_y_ij - P->delta_y / 2.0;
              y_ij_box_max = mid_node_y_ij + P->delta_y / 2.0;

              x_ij_box_min = mid_node_x_ij - P->delta_x / 2.0;
              x_ij_box_max = mid_node_x_ij + P->delta_x / 2.0;

              p_counter = 0;
              pressure = 0;
              sigma_1 = 0;
              sigma_3 = 0;
              mean = 0;
              differential = 0;

              runParticle = &refParticle; // start

              for (m = 0; m < numParticles;
                   m++) // look which particles are in the box
              {
                if (runParticle->ypos >= y_ij_box_min) {
                  if (runParticle->ypos <= y_ij_box_max) {
                    if (runParticle->xpos >= x_ij_box_min) {
                      if (runParticle->xpos <= x_ij_box_max) {
                        p_counter = p_counter + 1; // count particles

                        // set local pressure
                        pressure =
                            pressure +
                            sqrt((runParticle->F_P_x * runParticle->F_P_x) +
                                 (runParticle->F_P_y * runParticle->F_P_y));

                        // set stress tensor
                        sxx = runParticle->sxx;
                        syy = runParticle->syy;
                        sxy = runParticle->sxy;

                        // set pressure tensor

                        smax = ((sxx + syy) / 2.0) +
                               sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                                    sxy * sxy);
                        smin = ((sxx + syy) / 2.0) -
                               sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                                    sxy * sxy);

                        mean = mean + (smax + smin) / 2.0;

                        differential = differential + (smax - smin);
                        sigma_1 = sigma_1 + smax;
                        sigma_3 = sigma_3 + smin;
                      }
                    }
                  }
                }
                runParticle = runParticle->nextP;
              }
              pressure =
                  pressure / p_counter; // and divide by number of particles to
                                        // get average pressure
              sigma_1 = sigma_1 / p_counter;
              sigma_3 = sigma_3 / p_counter;
              mean = mean / p_counter;
              differential = differential / p_counter;

              fprintf(statStress, "%f\t", mid_node_x_ij);
              fprintf(statStress, "%f\t", mid_node_y_ij);
              fprintf(statStress, "%le\t", pressure);
              fprintf(statStress, "%le\t", sigma_1);
              fprintf(statStress, "%le\t", sigma_3);
              fprintf(statStress, "%le\t", mean);
              fprintf(statStress, "%le", differential);

              fprintf(statStress, " \n");
            }
          }
        }
      }

      else {
        P->press_nodes(i, j);
        mid_node_x_ij = P->press_pos_x;
        mid_node_y_ij = P->press_pos_y;

        y_ij_box_min = mid_node_y_ij - P->delta_y / 2.0;
        y_ij_box_max = mid_node_y_ij + P->delta_y / 2.0;

        x_ij_box_min = mid_node_x_ij - P->delta_x / 2.0;
        x_ij_box_max = mid_node_x_ij + P->delta_x / 2.0;

        p_counter = 0;
        pressure = 0;
        sigma_1 = 0;
        sigma_3 = 0;
        mean = 0;
        differential = 0;

        runParticle = &refParticle; // start

        for (k = 0; k < numParticles;
             k++) // look which particles are in the box
        {
          if (runParticle->ypos >= y_ij_box_min) {
            if (runParticle->ypos <= y_ij_box_max) {
              if (runParticle->xpos >= x_ij_box_min) {
                if (runParticle->xpos <= x_ij_box_max) {
                  p_counter = p_counter + 1; // count particles

                  // set local pressure
                  pressure = pressure +
                             sqrt((runParticle->F_P_x * runParticle->F_P_x) +
                                  (runParticle->F_P_y * runParticle->F_P_y));

                  // set stress tensor
                  sxx = runParticle->sxx;
                  syy = runParticle->syy;
                  sxy = runParticle->sxy;

                  // set pressure tensor

                  smax = ((sxx + syy) / 2.0) +
                         sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                              sxy * sxy);
                  smin = ((sxx + syy) / 2.0) -
                         sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                              sxy * sxy);

                  mean = mean + (smax + smin) / 2.0;

                  differential = differential + (smax - smin);
                  sigma_1 = sigma_1 + smax;
                  sigma_3 = sigma_3 + smin;
                }
              }
            }
          }
          runParticle = runParticle->nextP;
        }
        pressure = pressure / p_counter; // and divide by number of particles to
                                         // get average pressure
        sigma_1 = sigma_1 / p_counter;
        sigma_3 = sigma_3 / p_counter;
        mean = mean / p_counter;
        differential = differential / p_counter;

        fprintf(statStress, "%f\t", mid_node_x_ij);
        fprintf(statStress, "%f\t", mid_node_y_ij);
        fprintf(statStress, "%le\t", pressure);
        fprintf(statStress, "%le\t", sigma_1);
        fprintf(statStress, "%le\t", sigma_3);
        fprintf(statStress, "%le\t", mean);
        fprintf(statStress, "%le", differential);

        fprintf(statStress, " \n");
      }
    }

    i += 5;
  }

  fprintf(statStress, "\n");
  fclose(statStress); // close file
}

void Lattice::DumpStat_PressMDMxMnStress_HRes_c(int press_node_y,
                                                int press_node_x) {
  FILE *statStress;                 // file pointer
  int i, j, k, l, m, p_counter;     // counters
  double pressure;                  // mean continuum pressure, elapsed time;
  double sxx, syy, sxy, smax, smin; // stress tensor and eigenvalues
  double sigma_1, sigma_3, mean,
      differential; // mean stress, differential stress, pressure

  float mid_node_x_ij, mid_node_y_ij;
  float x_ij_box_min, x_ij_box_max, y_ij_box_min, y_ij_box_max;

  statStress = fopen("Local_Stresses_c.txt", "a");

  // press_cal_time += cal_time;

  // for(i = 10; i < P->max_size-10; i +=10)
  i = 0;
  while (i < P->max_size) {
    for (j = 0; j < P->max_size; j += 10) {
      if ((i > press_node_y - 10 && i < press_node_y + 10) &&
          (j > press_node_x - 10 && j < press_node_x + 10)) {
        for (k = i - 8; k <= i + 8; k += 4) {
          for (l = j - 8; l <= j + 8; l += 4) {
            // if (k != press_node_y && l != press_node_x)
            {
              P->press_nodes(k, l);
              mid_node_x_ij = P->press_pos_x;
              mid_node_y_ij = P->press_pos_y;

              y_ij_box_min = mid_node_y_ij - P->delta_y / 2.0;
              y_ij_box_max = mid_node_y_ij + P->delta_y / 2.0;

              x_ij_box_min = mid_node_x_ij - P->delta_x / 2.0;
              x_ij_box_max = mid_node_x_ij + P->delta_x / 2.0;

              p_counter = 0;
              pressure = 0;
              sigma_1 = 0;
              sigma_3 = 0;
              mean = 0;
              differential = 0;

              runParticle = &refParticle; // start

              for (m = 0; m < numParticles;
                   m++) // look which particles are in the box
              {
                if (runParticle->ypos >= y_ij_box_min) {
                  if (runParticle->ypos <= y_ij_box_max) {
                    if (runParticle->xpos >= x_ij_box_min) {
                      if (runParticle->xpos <= x_ij_box_max) {
                        p_counter = p_counter + 1; // count particles

                        // set local pressure
                        pressure =
                            pressure +
                            sqrt((runParticle->F_P_x * runParticle->F_P_x) +
                                 (runParticle->F_P_y * runParticle->F_P_y));

                        // set stress tensor
                        sxx = runParticle->sxx;
                        syy = runParticle->syy;
                        sxy = runParticle->sxy;

                        // set pressure tensor

                        smax = ((sxx + syy) / 2.0) +
                               sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                                    sxy * sxy);
                        smin = ((sxx + syy) / 2.0) -
                               sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                                    sxy * sxy);

                        mean = mean + (smax + smin) / 2.0;

                        differential = differential + (smax - smin);
                        sigma_1 = sigma_1 + smax;
                        sigma_3 = sigma_3 + smin;
                      }
                    }
                  }
                }
                runParticle = runParticle->nextP;
              }
              pressure =
                  pressure / p_counter; // and divide by number of particles to
                                        // get average pressure
              sigma_1 = sigma_1 / p_counter;
              sigma_3 = sigma_3 / p_counter;
              mean = mean / p_counter;
              differential = differential / p_counter;

              fprintf(statStress, "%f\t", mid_node_x_ij);
              fprintf(statStress, "%f\t", mid_node_y_ij);
              fprintf(statStress, "%le\t", pressure);
              fprintf(statStress, "%le\t", sigma_1);
              fprintf(statStress, "%le\t", sigma_3);
              fprintf(statStress, "%le\t", mean);
              fprintf(statStress, "%le", differential);

              fprintf(statStress, " \n");
            }
          }
        }
      }

      else {
        P->press_nodes(i, j);
        mid_node_x_ij = P->press_pos_x;
        mid_node_y_ij = P->press_pos_y;

        y_ij_box_min = mid_node_y_ij - P->delta_y / 2.0;
        y_ij_box_max = mid_node_y_ij + P->delta_y / 2.0;

        x_ij_box_min = mid_node_x_ij - P->delta_x / 2.0;
        x_ij_box_max = mid_node_x_ij + P->delta_x / 2.0;

        p_counter = 0;
        pressure = 0;
        sigma_1 = 0;
        sigma_3 = 0;
        mean = 0;
        differential = 0;

        runParticle = &refParticle; // start

        for (k = 0; k < numParticles;
             k++) // look which particles are in the box
        {
          if (runParticle->ypos >= y_ij_box_min) {
            if (runParticle->ypos <= y_ij_box_max) {
              if (runParticle->xpos >= x_ij_box_min) {
                if (runParticle->xpos <= x_ij_box_max) {
                  p_counter = p_counter + 1; // count particles

                  // set local pressure
                  pressure = pressure +
                             sqrt((runParticle->F_P_x * runParticle->F_P_x) +
                                  (runParticle->F_P_y * runParticle->F_P_y));

                  // set stress tensor
                  sxx = runParticle->sxx;
                  syy = runParticle->syy;
                  sxy = runParticle->sxy;

                  // set pressure tensor

                  smax = ((sxx + syy) / 2.0) +
                         sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                              sxy * sxy);
                  smin = ((sxx + syy) / 2.0) -
                         sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                              sxy * sxy);

                  mean = mean + (smax + smin) / 2.0;

                  differential = differential + (smax - smin);
                  sigma_1 = sigma_1 + smax;
                  sigma_3 = sigma_3 + smin;
                }
              }
            }
          }
          runParticle = runParticle->nextP;
        }
        pressure = pressure / p_counter; // and divide by number of particles to
                                         // get average pressure
        sigma_1 = sigma_1 / p_counter;
        sigma_3 = sigma_3 / p_counter;
        mean = mean / p_counter;
        differential = differential / p_counter;

        fprintf(statStress, "%f\t", mid_node_x_ij);
        fprintf(statStress, "%f\t", mid_node_y_ij);
        fprintf(statStress, "%le\t", pressure);
        fprintf(statStress, "%le\t", sigma_1);
        fprintf(statStress, "%le\t", sigma_3);
        fprintf(statStress, "%le\t", mean);
        fprintf(statStress, "%le", differential);

        fprintf(statStress, " \n");
      }
    }

    i += 10;
  }

  fprintf(statStress, "\n");
  fclose(statStress); // close file
}

void Lattice::DumpStat_PressMDStress_LHRes_c_d_verSeal(float min_press_limit_x,
                                                       float max_press_limit_x,
                                                       int y_nod, int x_nod) {
  FILE *statStress;                 // file pointer
  int k, p_counter;                 // counters
  double pressure;                  // mean continuum pressure, elapsed time;
  double sxx, syy, sxy, smax, smin; // stress tensor and eigenvalues
  double sigma_1, sigma_3, mean,
      differential; // mean stress, differential stress, pressure

  float source_x, source_y;
  float x_mid, y_mid, seal_x_mid;
  float x_ij_box_min, x_ij_box_max, y_ij_box_min, y_ij_box_max;

  statStress = fopen("Local_Stresses_VerSeal.txt", "a");

  P->press_nodes(y_nod, x_nod);

  // press_cal_time += cal_time;

  // for(y_mid = 0.1; y_mid <= 1.0; y_mid +=0.1)
  y_mid = 0;
  while (y_mid <= 1.0) {
    for (x_mid = 0; x_mid <= 1.0; x_mid += 0.1) {
      if (x_mid < min_press_limit_x - 2.0 * P->delta_x ||
          x_mid > max_press_limit_x + 2.0 * P->delta_x) {
        if ((y_mid >= P->press_pos_y - 0.05 &&
             y_mid <= P->press_pos_y + 0.05) &&
            (x_mid >= P->press_pos_x - 0.05 &&
             x_mid <= P->press_pos_x + 0.05)) {
          for (source_y = y_mid - 0.06; source_y < y_mid + 0.07;
               source_y += 0.03) {
            for (source_x = x_mid - 0.06; source_x < x_mid + 0.07;
                 source_x += 0.03) {
              y_ij_box_min = source_y - P->delta_y / 2.0;
              y_ij_box_max = source_y + P->delta_y / 2.0;

              x_ij_box_min = source_x - P->delta_x / 2.0;
              x_ij_box_max = source_x + P->delta_x / 2.0;

              p_counter = 0;
              pressure = 0;
              sigma_1 = 0;
              sigma_3 = 0;
              mean = 0;
              differential = 0;

              runParticle = &refParticle; // start

              for (k = 0; k < numParticles;
                   k++) // look which particles are in the box
              {
                if (runParticle->ypos >= y_ij_box_min) {
                  if (runParticle->ypos <= y_ij_box_max) {
                    if (runParticle->xpos >= x_ij_box_min) {
                      if (runParticle->xpos <= x_ij_box_max) {
                        p_counter = p_counter + 1; // count particles

                        // set local pressure
                        pressure =
                            pressure +
                            sqrt((runParticle->F_P_x * runParticle->F_P_x) +
                                 (runParticle->F_P_y * runParticle->F_P_y));

                        // set stress tensor
                        sxx = runParticle->sxx;
                        syy = runParticle->syy;
                        sxy = runParticle->sxy;

                        // set pressure tensor
                        smax = ((sxx + syy) / 2.0) +
                               sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                                    sxy * sxy);
                        smin = ((sxx + syy) / 2.0) -
                               sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                                    sxy * sxy);

                        sigma_1 = sigma_1 + smax;
                        sigma_3 = sigma_3 + smin;

                        mean = mean + (smax + smin) / 2.0;

                        differential = differential + (smax - smin);
                      }
                    }
                  }
                }
                runParticle = runParticle->nextP;
              }
              pressure =
                  pressure / p_counter; // and divide by number of particles to
                                        // get average pressure
              sigma_1 = sigma_1 / p_counter;
              sigma_3 = sigma_3 / p_counter;
              mean = mean / p_counter;
              differential = differential / p_counter;

              fprintf(statStress, "%f\t", source_x);
              fprintf(statStress, "%f\t", source_y);
              fprintf(statStress, "%le\t", pressure);
              fprintf(statStress, "%le\t", sigma_1);
              fprintf(statStress, "%le\t", sigma_3);
              fprintf(statStress, "%le\t", mean);
              fprintf(statStress, "%le\t", differential);

              fprintf(statStress, " \n");
            }
          }
        }

        else {
          y_ij_box_min = y_mid - P->delta_y / 2.0;
          y_ij_box_max = y_mid + P->delta_y / 2.0;

          x_ij_box_min = x_mid - P->delta_x / 2.0;
          x_ij_box_max = x_mid + P->delta_x / 2.0;

          p_counter = 0;
          pressure = 0;
          sigma_1 = 0;
          sigma_3 = 0;
          mean = 0;
          differential = 0;

          runParticle = &refParticle; // start

          for (k = 0; k < numParticles;
               k++) // look which particles are in the box
          {
            if (runParticle->ypos >= y_ij_box_min) {
              if (runParticle->ypos <= y_ij_box_max) {
                if (runParticle->xpos >= x_ij_box_min) {
                  if (runParticle->xpos <= x_ij_box_max) {
                    p_counter = p_counter + 1; // count particles

                    // set local pressure
                    pressure = pressure +
                               sqrt((runParticle->F_P_x * runParticle->F_P_x) +
                                    (runParticle->F_P_y * runParticle->F_P_y));

                    // set stress tensor
                    sxx = runParticle->sxx;
                    syy = runParticle->syy;
                    sxy = runParticle->sxy;

                    // set pressure tensor
                    smax = ((sxx + syy) / 2.0) +
                           sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                                sxy * sxy);
                    smin = ((sxx + syy) / 2.0) -
                           sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                                sxy * sxy);

                    sigma_1 = sigma_1 + smax;
                    sigma_3 = sigma_3 + smin;

                    mean = mean + (smax + smin) / 2.0;

                    differential = differential + (smax - smin);
                  }
                }
              }
            }
            runParticle = runParticle->nextP;
          }
          pressure = pressure / p_counter; // and divide by number of particles
                                           // to get average pressure
          sigma_1 = sigma_1 / p_counter;
          sigma_3 = sigma_3 / p_counter;
          mean = mean / p_counter;
          differential = differential / p_counter;

          fprintf(statStress, "%f\t", x_mid);
          fprintf(statStress, "%f\t", y_mid);
          fprintf(statStress, "%le\t", pressure);
          fprintf(statStress, "%le\t", sigma_1);
          fprintf(statStress, "%le\t", sigma_3);
          fprintf(statStress, "%le\t", mean);
          fprintf(statStress, "%le\t", differential);

          fprintf(statStress, " \n");
        }
      }

      else {
        for (seal_x_mid = x_mid - 1.5 * P->delta_x;
             seal_x_mid <= x_mid + 1.5 * P->delta_x;
             seal_x_mid += 1.5 * P->delta_x) {
          y_ij_box_min = y_mid - P->delta_y / 2.0;
          y_ij_box_max = y_mid + P->delta_y / 2.0;

          x_ij_box_min = seal_x_mid - P->delta_x / 2.0;
          x_ij_box_max = seal_x_mid + P->delta_x / 2.0;

          p_counter = 0;
          pressure = 0;
          sigma_1 = 0;
          sigma_3 = 0;
          mean = 0;
          differential = 0;

          runParticle = &refParticle; // start

          for (k = 0; k < numParticles;
               k++) // look which particles are in the box
          {
            if (runParticle->ypos >= y_ij_box_min) {
              if (runParticle->ypos <= y_ij_box_max) {
                if (runParticle->xpos >= x_ij_box_min) {
                  if (runParticle->xpos <= x_ij_box_max) {
                    p_counter = p_counter + 1; // count particles

                    // set local pressure
                    pressure = pressure +
                               sqrt((runParticle->F_P_x * runParticle->F_P_x) +
                                    (runParticle->F_P_y * runParticle->F_P_y));

                    // set stress tensor
                    sxx = runParticle->sxx;
                    syy = runParticle->syy;
                    sxy = runParticle->sxy;

                    // set pressure tensor
                    smax = ((sxx + syy) / 2.0) +
                           sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                                sxy * sxy);
                    smin = ((sxx + syy) / 2.0) -
                           sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                                sxy * sxy);

                    sigma_1 = sigma_1 + smax;
                    sigma_3 = sigma_3 + smin;

                    mean = mean + (smax + smin) / 2.0;

                    differential = differential + (smax - smin);
                  }
                }
              }
            }
            runParticle = runParticle->nextP;
          }
          pressure = pressure / p_counter; // and divide by number of particles
                                           // to get average pressure
          sigma_1 = sigma_1 / p_counter;
          sigma_3 = sigma_3 / p_counter;
          mean = mean / p_counter;
          differential = differential / p_counter;

          fprintf(statStress, "%f\t", seal_x_mid);
          fprintf(statStress, "%f\t", y_mid);
          fprintf(statStress, "%le\t", pressure);
          fprintf(statStress, "%le\t", sigma_1);
          fprintf(statStress, "%le\t", sigma_3);
          fprintf(statStress, "%le\t", mean);
          fprintf(statStress, "%le\t", differential);

          fprintf(statStress, " \n");
        }
      }
    }

    y_mid += 0.1;
  }

  fprintf(statStress, "\n");
  fclose(statStress); // close file
}

void Lattice::DumpStat_PressMDStress_HRes_c_d_verSeal(float min_press_limit_x,
                                                      float max_press_limit_x,
                                                      int y_nod, int x_nod) {
  FILE *statStress;                 // file pointer
  int k, p_counter;                 // counters
  double pressure;                  // mean continuum pressure, elapsed time;
  double sxx, syy, sxy, smax, smin; // stress tensor and eigenvalues
  double sigma_1, sigma_3, mean,
      differential; // mean stress, differential stress, pressure

  float source_x, source_y;
  float x_mid, y_mid, seal_x_mid;
  float x_ij_box_min, x_ij_box_max, y_ij_box_min, y_ij_box_max;

  statStress = fopen("Local_Stresses_VerSeal.txt", "a");

  P->press_nodes(y_nod, x_nod);

  // press_cal_time += cal_time;

  // for(y_mid = 0.1; y_mid <= 1.0; y_mid +=0.1)
  y_mid = 0;
  while (y_mid <= 1.0) {
    for (x_mid = 0; x_mid <= 1.0; x_mid += 0.1) {
      if (x_mid < min_press_limit_x - 4.0 * P->delta_x ||
          x_mid > max_press_limit_x + 4.0 * P->delta_x) {
        if ((y_mid >= P->press_pos_y - 0.05 &&
             y_mid <= P->press_pos_y + 0.05) &&
            (x_mid >= P->press_pos_x - 0.05 &&
             x_mid <= P->press_pos_x + 0.05)) {
          for (source_y = y_mid - 0.06; source_y < y_mid + 0.07;
               source_y += 0.03) {
            for (source_x = x_mid - 0.06; source_x < x_mid + 0.07;
                 source_x += 0.03) {
              y_ij_box_min = source_y - P->delta_y / 2.0;
              y_ij_box_max = source_y + P->delta_y / 2.0;

              x_ij_box_min = source_x - P->delta_x / 2.0;
              x_ij_box_max = source_x + P->delta_x / 2.0;

              p_counter = 0;
              pressure = 0;
              sigma_1 = 0;
              sigma_3 = 0;
              mean = 0;
              differential = 0;

              runParticle = &refParticle; // start

              for (k = 0; k < numParticles;
                   k++) // look which particles are in the box
              {
                if (runParticle->ypos >= y_ij_box_min) {
                  if (runParticle->ypos <= y_ij_box_max) {
                    if (runParticle->xpos >= x_ij_box_min) {
                      if (runParticle->xpos <= x_ij_box_max) {
                        p_counter = p_counter + 1; // count particles

                        // set local pressure
                        pressure =
                            pressure +
                            sqrt((runParticle->F_P_x * runParticle->F_P_x) +
                                 (runParticle->F_P_y * runParticle->F_P_y));

                        // set stress tensor
                        sxx = runParticle->sxx;
                        syy = runParticle->syy;
                        sxy = runParticle->sxy;

                        // set pressure tensor
                        smax = ((sxx + syy) / 2.0) +
                               sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                                    sxy * sxy);
                        smin = ((sxx + syy) / 2.0) -
                               sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                                    sxy * sxy);

                        sigma_1 = sigma_1 + smax;
                        sigma_3 = sigma_3 + smin;

                        mean = mean + (smax + smin) / 2.0;

                        differential = differential + (smax - smin);
                      }
                    }
                  }
                }
                runParticle = runParticle->nextP;
              }
              pressure =
                  pressure / p_counter; // and divide by number of particles to
                                        // get average pressure
              sigma_1 = sigma_1 / p_counter;
              sigma_3 = sigma_3 / p_counter;
              mean = mean / p_counter;
              differential = differential / p_counter;

              fprintf(statStress, "%f\t", source_x);
              fprintf(statStress, "%f\t", source_y);
              fprintf(statStress, "%le\t", pressure);
              fprintf(statStress, "%le\t", sigma_1);
              fprintf(statStress, "%le\t", sigma_3);
              fprintf(statStress, "%le\t", mean);
              fprintf(statStress, "%le\t", differential);

              fprintf(statStress, " \n");
            }
          }
        }

        else {
          y_ij_box_min = y_mid - P->delta_y / 2.0;
          y_ij_box_max = y_mid + P->delta_y / 2.0;

          x_ij_box_min = x_mid - P->delta_x / 2.0;
          x_ij_box_max = x_mid + P->delta_x / 2.0;

          p_counter = 0;
          pressure = 0;
          sigma_1 = 0;
          sigma_3 = 0;
          mean = 0;
          differential = 0;

          runParticle = &refParticle; // start

          for (k = 0; k < numParticles;
               k++) // look which particles are in the box
          {
            if (runParticle->ypos >= y_ij_box_min) {
              if (runParticle->ypos <= y_ij_box_max) {
                if (runParticle->xpos >= x_ij_box_min) {
                  if (runParticle->xpos <= x_ij_box_max) {
                    p_counter = p_counter + 1; // count particles

                    // set local pressure
                    pressure = pressure +
                               sqrt((runParticle->F_P_x * runParticle->F_P_x) +
                                    (runParticle->F_P_y * runParticle->F_P_y));

                    // set stress tensor
                    sxx = runParticle->sxx;
                    syy = runParticle->syy;
                    sxy = runParticle->sxy;

                    // set pressure tensor
                    smax = ((sxx + syy) / 2.0) +
                           sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                                sxy * sxy);
                    smin = ((sxx + syy) / 2.0) -
                           sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                                sxy * sxy);

                    sigma_1 = sigma_1 + smax;
                    sigma_3 = sigma_3 + smin;

                    mean = mean + (smax + smin) / 2.0;

                    differential = differential + (smax - smin);
                  }
                }
              }
            }
            runParticle = runParticle->nextP;
          }
          pressure = pressure / p_counter; // and divide by number of particles
                                           // to get average pressure
          sigma_1 = sigma_1 / p_counter;
          sigma_3 = sigma_3 / p_counter;
          mean = mean / p_counter;
          differential = differential / p_counter;

          fprintf(statStress, "%f\t", x_mid);
          fprintf(statStress, "%f\t", y_mid);
          fprintf(statStress, "%le\t", pressure);
          fprintf(statStress, "%le\t", sigma_1);
          fprintf(statStress, "%le\t", sigma_3);
          fprintf(statStress, "%le\t", mean);
          fprintf(statStress, "%le\t", differential);

          fprintf(statStress, " \n");
        }
      }

      else {
        for (seal_x_mid = x_mid - 3.0 * P->delta_x;
             seal_x_mid <= x_mid + 3.0 * P->delta_x;
             seal_x_mid += 3.0 * P->delta_x) {
          y_ij_box_min = y_mid - P->delta_y / 2.0;
          y_ij_box_max = y_mid + P->delta_y / 2.0;

          x_ij_box_min = seal_x_mid - P->delta_x / 2.0;
          x_ij_box_max = seal_x_mid + P->delta_x / 2.0;

          p_counter = 0;
          pressure = 0;
          sigma_1 = 0;
          sigma_3 = 0;
          mean = 0;
          differential = 0;

          runParticle = &refParticle; // start

          for (k = 0; k < numParticles;
               k++) // look which particles are in the box
          {
            if (runParticle->ypos >= y_ij_box_min) {
              if (runParticle->ypos <= y_ij_box_max) {
                if (runParticle->xpos >= x_ij_box_min) {
                  if (runParticle->xpos <= x_ij_box_max) {
                    p_counter = p_counter + 1; // count particles

                    // set local pressure
                    pressure = pressure +
                               sqrt((runParticle->F_P_x * runParticle->F_P_x) +
                                    (runParticle->F_P_y * runParticle->F_P_y));

                    // set stress tensor
                    sxx = runParticle->sxx;
                    syy = runParticle->syy;
                    sxy = runParticle->sxy;

                    // set pressure tensor
                    smax = ((sxx + syy) / 2.0) +
                           sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                                sxy * sxy);
                    smin = ((sxx + syy) / 2.0) -
                           sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                                sxy * sxy);

                    sigma_1 = sigma_1 + smax;
                    sigma_3 = sigma_3 + smin;

                    mean = mean + (smax + smin) / 2.0;

                    differential = differential + (smax - smin);
                  }
                }
              }
            }
            runParticle = runParticle->nextP;
          }
          pressure = pressure / p_counter; // and divide by number of particles
                                           // to get average pressure
          sigma_1 = sigma_1 / p_counter;
          sigma_3 = sigma_3 / p_counter;
          mean = mean / p_counter;
          differential = differential / p_counter;

          fprintf(statStress, "%f\t", seal_x_mid);
          fprintf(statStress, "%f\t", y_mid);
          fprintf(statStress, "%le\t", pressure);
          fprintf(statStress, "%le\t", sigma_1);
          fprintf(statStress, "%le\t", sigma_3);
          fprintf(statStress, "%le\t", mean);
          fprintf(statStress, "%le\t", differential);

          fprintf(statStress, " \n");
        }
      }
    }

    y_mid += 0.1;
  }

  fprintf(statStress, "\n");
  fclose(statStress); // close file
}

void Lattice::DumpStat_PressMDStress_LHRes_c_d_horSeal(float min_press_limit_y,
                                                       float max_press_limit_y,
                                                       int y_nod, int x_nod) {
  FILE *statStress;                 // file pointer
  int k, p_counter;                 // counters
  double pressure;                  // mean continuum pressure, elapsed time;
  double sxx, syy, sxy, smax, smin; // stress tensor and eigenvalues
  double sigma_1, sigma_3, mean,
      differential; // mean stress, differential stress, pressure

  float source_x, source_y;
  float x_mid, y_mid, seal_y_mid;
  float x_ij_box_min, x_ij_box_max, y_ij_box_min, y_ij_box_max;

  statStress = fopen("Local_Stresses_HorSeal.txt", "a");

  P->press_nodes(y_nod, x_nod);
  // press_cal_time += cal_time;

  // for(y_mid = 0.1; y_mid <= 1.0; y_mid +=0.1)
  y_mid = 0;
  while (y_mid <= 1.0) {
    for (x_mid = 0; x_mid <= 1.0; x_mid += 0.1) {
      if (y_mid < min_press_limit_y - 2.0 * P->delta_y ||
          y_mid > max_press_limit_y + 2.0 * P->delta_y) {
        if ((y_mid >= P->press_pos_y - 0.05 &&
             y_mid <= P->press_pos_y + 0.05) &&
            (x_mid >= P->press_pos_x - 0.05 &&
             x_mid <= P->press_pos_x + 0.05)) {
          for (source_y = y_mid - 0.06; source_y < y_mid + 0.07;
               source_y += 0.03) {
            for (source_x = x_mid - 0.06; source_x < x_mid + 0.07;
                 source_x += 0.03) {
              y_ij_box_min = source_y - P->delta_y / 2.0;
              y_ij_box_max = source_y + P->delta_y / 2.0;

              x_ij_box_min = source_x - P->delta_x / 2.0;
              x_ij_box_max = source_x + P->delta_x / 2.0;

              p_counter = 0;
              pressure = 0;
              sigma_1 = 0;
              sigma_3 = 0;
              mean = 0;
              differential = 0;

              runParticle = &refParticle; // start

              for (k = 0; k < numParticles;
                   k++) // look which particles are in the box
              {
                if (runParticle->ypos >= y_ij_box_min) {
                  if (runParticle->ypos <= y_ij_box_max) {
                    if (runParticle->xpos >= x_ij_box_min) {
                      if (runParticle->xpos <= x_ij_box_max) {
                        p_counter = p_counter + 1; // count particles

                        // set local pressure
                        pressure =
                            pressure +
                            sqrt((runParticle->F_P_x * runParticle->F_P_x) +
                                 (runParticle->F_P_y * runParticle->F_P_y));

                        // set stress tensor
                        sxx = runParticle->sxx;
                        syy = runParticle->syy;
                        sxy = runParticle->sxy;

                        // set pressure tensor
                        smax = ((sxx + syy) / 2.0) +
                               sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                                    sxy * sxy);
                        smin = ((sxx + syy) / 2.0) -
                               sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                                    sxy * sxy);

                        sigma_1 = sigma_1 + smax;
                        sigma_3 = sigma_3 + smin;

                        mean = mean + (smax + smin) / 2.0;

                        differential = differential + (smax - smin);
                      }
                    }
                  }
                }
                runParticle = runParticle->nextP;
              }
              pressure =
                  pressure / p_counter; // and divide by number of particles to
                                        // get average pressure
              sigma_1 = sigma_1 / p_counter;
              sigma_3 = sigma_3 / p_counter;
              mean = mean / p_counter;
              differential = differential / p_counter;

              fprintf(statStress, "%f\t", source_x);
              fprintf(statStress, "%f\t", source_y);
              fprintf(statStress, "%le\t", pressure);
              fprintf(statStress, "%le\t", sigma_1);
              fprintf(statStress, "%le\t", sigma_3);
              fprintf(statStress, "%le\t", mean);
              fprintf(statStress, "%le\t", differential);

              fprintf(statStress, " \n");
            }
          }
        }

        else {
          y_ij_box_min = y_mid - P->delta_y / 2.0;
          y_ij_box_max = y_mid + P->delta_y / 2.0;

          x_ij_box_min = x_mid - P->delta_x / 2.0;
          x_ij_box_max = x_mid + P->delta_x / 2.0;

          p_counter = 0;
          pressure = 0;
          sigma_1 = 0;
          sigma_3 = 0;
          mean = 0;
          differential = 0;

          runParticle = &refParticle; // start

          for (k = 0; k < numParticles;
               k++) // look which particles are in the box
          {
            if (runParticle->ypos >= y_ij_box_min) {
              if (runParticle->ypos <= y_ij_box_max) {
                if (runParticle->xpos >= x_ij_box_min) {
                  if (runParticle->xpos <= x_ij_box_max) {
                    p_counter = p_counter + 1; // count particles

                    // set local pressure
                    pressure = pressure +
                               sqrt((runParticle->F_P_x * runParticle->F_P_x) +
                                    (runParticle->F_P_y * runParticle->F_P_y));

                    // set stress tensor
                    sxx = runParticle->sxx;
                    syy = runParticle->syy;
                    sxy = runParticle->sxy;

                    // set pressure tensor
                    smax = ((sxx + syy) / 2.0) +
                           sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                                sxy * sxy);
                    smin = ((sxx + syy) / 2.0) -
                           sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                                sxy * sxy);

                    sigma_1 = sigma_1 + smax;
                    sigma_3 = sigma_3 + smin;

                    mean = mean + (smax + smin) / 2.0;

                    differential = differential + (smax - smin);
                  }
                }
              }
            }
            runParticle = runParticle->nextP;
          }
          pressure = pressure / p_counter; // and divide by number of particles
                                           // to get average pressure
          sigma_1 = sigma_1 / p_counter;
          sigma_3 = sigma_3 / p_counter;
          mean = mean / p_counter;
          differential = differential / p_counter;

          fprintf(statStress, "%f\t", x_mid);
          fprintf(statStress, "%f\t", y_mid);
          fprintf(statStress, "%le\t", pressure);
          fprintf(statStress, "%le\t", sigma_1);
          fprintf(statStress, "%le\t", sigma_3);
          fprintf(statStress, "%le\t", mean);
          fprintf(statStress, "%le\t", differential);

          fprintf(statStress, " \n");
        }
      }

      else {
        for (seal_y_mid = y_mid - 1.5 * P->delta_y;
             seal_y_mid <= y_mid + 1.5 * P->delta_y;
             seal_y_mid += 1.5 * P->delta_y) {
          y_ij_box_min = seal_y_mid - P->delta_y / 2.0;
          y_ij_box_max = seal_y_mid + P->delta_y / 2.0;

          x_ij_box_min = x_mid - P->delta_x / 2.0;
          x_ij_box_max = x_mid + P->delta_x / 2.0;

          p_counter = 0;
          pressure = 0;
          sigma_1 = 0;
          sigma_3 = 0;
          sxy = 0;
          mean = 0;
          differential = 0;

          runParticle = &refParticle; // start

          for (k = 0; k < numParticles;
               k++) // look which particles are in the box
          {
            if (runParticle->ypos >= y_ij_box_min) {
              if (runParticle->ypos <= y_ij_box_max) {
                if (runParticle->xpos >= x_ij_box_min) {
                  if (runParticle->xpos <= x_ij_box_max) {
                    p_counter = p_counter + 1; // count particles

                    // set local pressure
                    pressure = pressure +
                               sqrt((runParticle->F_P_x * runParticle->F_P_x) +
                                    (runParticle->F_P_y * runParticle->F_P_y));

                    // set stress tensor
                    sxx = runParticle->sxx;
                    syy = runParticle->syy;
                    sxy = runParticle->sxy;

                    // set pressure tensor
                    smax = ((sxx + syy) / 2.0) +
                           sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                                sxy * sxy);
                    smin = ((sxx + syy) / 2.0) -
                           sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                                sxy * sxy);

                    sigma_1 = sigma_1 + smax;
                    sigma_3 = sigma_3 + smin;

                    mean = mean + (smax + smin) / 2.0;

                    differential = differential + (smax - smin);
                  }
                }
              }
            }
            runParticle = runParticle->nextP;
          }
          pressure = pressure / p_counter; // and divide by number of particles
                                           // to get average pressure
          sigma_1 = sigma_1 / p_counter;
          sigma_3 = sigma_3 / p_counter;
          mean = mean / p_counter;
          differential = differential / p_counter;

          fprintf(statStress, "%f\t", x_mid);
          fprintf(statStress, "%f\t", seal_y_mid);
          fprintf(statStress, "%le\t", pressure);
          fprintf(statStress, "%le\t", sigma_1);
          fprintf(statStress, "%le\t", sigma_3);
          fprintf(statStress, "%le\t", mean);
          fprintf(statStress, "%le\t", differential);

          fprintf(statStress, " \n");
        }
      }
    }

    y_mid += 0.1;
  }

  fprintf(statStress, "\n");
  fclose(statStress); // close file
}

void Lattice::DumpStat_PressMDStress_LHRes_c_d_horSeal_randSource(
    float min_press_limit_y, float max_press_limit_y) {
  FILE *statStress;                 // file pointer
  int k, p_counter;                 // counters
  double pressure;                  // mean continuum pressure, elapsed time;
  double sxx, syy, sxy, smax, smin; // stress tensor and eigenvalues
  double sigma_1, sigma_3, mean,
      differential; // mean stress, differential stress, pressure

  float x_mid, y_mid, seal_x_mid, seal_y_mid;
  float x_ij_box_min, x_ij_box_max, y_ij_box_min, y_ij_box_max;

  statStress = fopen("Local_Stresses_UnderSealUpper(rand source).txt", "a");

  // press_cal_time += cal_time;

  y_mid = 0.1;
  while (y_mid < 1.0) {
    for (x_mid = 0.1; x_mid < 1.0; x_mid += 0.1) {
      if (y_mid < min_press_limit_y || y_mid > max_press_limit_y) {
        y_ij_box_min = y_mid - P->delta_y / 2.0;
        y_ij_box_max = y_mid + P->delta_y / 2.0;

        x_ij_box_min = x_mid - P->delta_x / 2.0;
        x_ij_box_max = x_mid + P->delta_x / 2.0;

        p_counter = 0;
        pressure = 0;
        sigma_1 = 0;
        sigma_3 = 0;
        mean = 0;
        differential = 0;

        runParticle = &refParticle; // start

        for (k = 0; k < numParticles;
             k++) // look which particles are in the box
        {
          if (runParticle->ypos >= y_ij_box_min) {
            if (runParticle->ypos <= y_ij_box_max) {
              if (runParticle->xpos >= x_ij_box_min) {
                if (runParticle->xpos <= x_ij_box_max) {
                  p_counter = p_counter + 1; // count particles

                  // set local pressure
                  pressure = pressure +
                             sqrt((runParticle->F_P_x * runParticle->F_P_x) +
                                  (runParticle->F_P_y * runParticle->F_P_y));

                  // set stress tensor
                  sxx = runParticle->sxx;
                  syy = runParticle->syy;
                  sxy = runParticle->sxy;

                  // set pressure tensor
                  smax = ((sxx + syy) / 2.0) +
                         sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                              sxy * sxy);
                  smin = ((sxx + syy) / 2.0) -
                         sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                              sxy * sxy);

                  sigma_1 = sigma_1 + smax;
                  sigma_3 = sigma_3 + smin;

                  mean = mean + (smax + smin) / 2.0;

                  differential = differential + (smax - smin);
                }
              }
            }
          }
          runParticle = runParticle->nextP;
        }
        pressure = pressure / p_counter; // and divide by number of particles to
                                         // get average pressure
        sigma_1 = sigma_1 / p_counter;
        sigma_3 = sigma_3 / p_counter;
        mean = mean / p_counter;
        differential = differential / p_counter;

        fprintf(statStress, "%f\t", x_mid);
        fprintf(statStress, "%f\t", y_mid);
        fprintf(statStress, "%le\t", pressure);
        fprintf(statStress, "%le\t", sigma_1);
        fprintf(statStress, "%le\t", sigma_3);
        fprintf(statStress, "%le\t", mean);
        fprintf(statStress, "%le\t", differential);
        fprintf(statStress, " \n");
      }

      else {
        for (seal_y_mid = y_mid - 2.5 * P->delta_y;
             seal_y_mid <= y_mid + 2.5 * P->delta_y;
             seal_y_mid += 1.25 * P->delta_y) {
          for (seal_x_mid = x_mid - 1.5 * P->delta_x;
               seal_x_mid <= x_mid + 1.5 * P->delta_x;
               seal_x_mid += 1.5 * P->delta_x) {
            y_ij_box_min = seal_y_mid - P->delta_y / 2.0;
            y_ij_box_max = seal_y_mid + P->delta_y / 2.0;

            x_ij_box_min = seal_x_mid - P->delta_x / 2.0;
            x_ij_box_max = seal_x_mid + P->delta_x / 2.0;

            p_counter = 0;
            pressure = 0;
            sigma_1 = 0;
            sigma_3 = 0;
            mean = 0;
            differential = 0;

            runParticle = &refParticle; // start

            for (k = 0; k < numParticles;
                 k++) // look which particles are in the box
            {
              if (runParticle->ypos >= y_ij_box_min) {
                if (runParticle->ypos <= y_ij_box_max) {
                  if (runParticle->xpos >= x_ij_box_min) {
                    if (runParticle->xpos <= x_ij_box_max) {
                      p_counter = p_counter + 1; // count particles

                      // set local pressure
                      pressure =
                          pressure +
                          sqrt((runParticle->F_P_x * runParticle->F_P_x) +
                               (runParticle->F_P_y * runParticle->F_P_y));

                      // set stress tensor
                      sxx = runParticle->sxx;
                      syy = runParticle->syy;
                      sxy = runParticle->sxy;

                      // set pressure tensor
                      smax = ((sxx + syy) / 2.0) +
                             sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                                  sxy * sxy);
                      smin = ((sxx + syy) / 2.0) -
                             sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                                  sxy * sxy);

                      sigma_1 = sigma_1 + smax;
                      sigma_3 = sigma_3 + smin;

                      mean = mean + (smax + smin) / 2.0;

                      differential = differential + (smax - smin);
                    }
                  }
                }
              }
              runParticle = runParticle->nextP;
            }
            pressure =
                pressure / p_counter; // and divide by number of particles to
                                      // get average pressure
            sigma_1 = sigma_1 / p_counter;
            sigma_3 = sigma_3 / p_counter;
            mean = mean / p_counter;
            differential = differential / p_counter;

            fprintf(statStress, "%f\t", seal_x_mid);
            fprintf(statStress, "%f\t", seal_y_mid);
            fprintf(statStress, "%le\t", pressure);
            fprintf(statStress, "%le\t", sigma_1);
            fprintf(statStress, "%le\t", sigma_3);
            fprintf(statStress, "%le\t", mean);
            fprintf(statStress, "%le\t", differential);
            fprintf(statStress, " \n");
          }
        }
      }
    }

    y_mid += 0.1;
  }

  fprintf(statStress, "\n");
  fclose(statStress); // close file
}

void Lattice::DumpStat_PressMDStress_LHRes_c_d_randSource() {
  FILE *statStress;                 // file pointer
  int k, p_counter;                 // counters
  double pressure;                  // mean continuum pressure;
  double sxx, syy, sxy, smax, smin; // stress tensor and eigenvalues
  double sigma_1, sigma_3, mean,
      differential; // mean stress, differential stress, pressure

  float x_mid, y_mid;
  float x_ij_box_min, x_ij_box_max, y_ij_box_min, y_ij_box_max;

  statStress = fopen("Local_Stresses_FullModel(rand source).txt", "a");

  // press_cal_time += cal_time;

  for (y_mid = 0; y_mid <= 1.0; y_mid += 0.1) {
    for (x_mid = 0; x_mid <= 1.0; x_mid += 0.1) {
      y_ij_box_min = y_mid - P->delta_y / 2.0;
      y_ij_box_max = y_mid + P->delta_y / 2.0;

      x_ij_box_min = x_mid - P->delta_x / 2.0;
      x_ij_box_max = x_mid + P->delta_x / 2.0;

      p_counter = 0;
      pressure = 0;
      sigma_1 = 0;
      sigma_3 = 0;
      mean = 0;
      differential = 0;

      runParticle = &refParticle; // start

      for (k = 0; k < numParticles; k++) // look which particles are in the box
      {
        if (runParticle->ypos >= y_ij_box_min) {
          if (runParticle->ypos <= y_ij_box_max) {
            if (runParticle->xpos >= x_ij_box_min) {
              if (runParticle->xpos <= x_ij_box_max) {
                p_counter = p_counter + 1; // count particles

                // set local pressure
                pressure =
                    pressure + sqrt((runParticle->F_P_x * runParticle->F_P_x) +
                                    (runParticle->F_P_y * runParticle->F_P_y));

                // set stress tensor
                sxx = runParticle->sxx;
                syy = runParticle->syy;
                sxy = runParticle->sxy;

                // set pressure tensor
                smax =
                    ((sxx + syy) / 2.0) +
                    sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) + sxy * sxy);
                smin =
                    ((sxx + syy) / 2.0) -
                    sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) + sxy * sxy);

                sigma_1 = sigma_1 + smax;
                sigma_3 = sigma_3 + smin;

                mean = mean + (smax + smin) / 2.0;

                differential = differential + (smax - smin);
              }
            }
          }
        }
        runParticle = runParticle->nextP;
      }
      pressure = pressure / p_counter; // and divide by number of particles to
                                       // get average pressure
      sigma_1 = sigma_1 / p_counter;
      sigma_3 = sigma_3 / p_counter;
      mean = mean / p_counter;
      differential = differential / p_counter;

      fprintf(statStress, "%f\t", x_mid);
      fprintf(statStress, "%f\t", y_mid);
      fprintf(statStress, "%le\t", pressure);
      fprintf(statStress, "%le\t", sigma_1);
      fprintf(statStress, "%le\t", sigma_3);
      fprintf(statStress, "%le\t", mean);
      fprintf(statStress, "%le", differential);
      fprintf(statStress, " \n");
    }
  }

  fprintf(statStress, "\n");
  fclose(statStress); // close file
}

void Lattice::DumpStat_PressMDStress_HRes_c_d_horSeal(float min_press_limit_y,
                                                      float max_press_limit_y,
                                                      int y_nod, int x_nod) {
  FILE *statStress;                 // file pointer
  int k, p_counter;                 // counters
  double pressure;                  // mean continuum pressure, elapsed time;
  double sxx, syy, sxy, smax, smin; // stress tensor and eigenvalues
  double sigma_1, sigma_3, mean,
      differential; // mean stress, differential stress, pressure

  float source_x, source_y;
  float x_mid, y_mid, seal_y_mid;
  float x_ij_box_min, x_ij_box_max, y_ij_box_min, y_ij_box_max;

  statStress = fopen("Local_Stresses_HorSeal.txt", "a");

  P->press_nodes(y_nod, x_nod);
  // press_cal_time += cal_time;

  // for(y_mid = 0.1; y_mid <= 1.0; y_mid +=0.1)
  y_mid = 0;
  while (y_mid <= 1.0) {
    for (x_mid = 0; x_mid <= 1.0; x_mid += 0.1) {
      if (y_mid < min_press_limit_y - 4.0 * P->delta_y ||
          y_mid > max_press_limit_y + 4.0 * P->delta_y) {
        if ((y_mid >= P->press_pos_y - 0.05 &&
             y_mid <= P->press_pos_y + 0.05) &&
            (x_mid >= P->press_pos_x - 0.05 &&
             x_mid <= P->press_pos_x + 0.05)) {
          for (source_y = y_mid - 0.06; source_y < y_mid + 0.07;
               source_y += 0.03) {
            for (source_x = x_mid - 0.06; source_x < x_mid + 0.07;
                 source_x += 0.03) {
              y_ij_box_min = source_y - P->delta_y / 2.0;
              y_ij_box_max = source_y + P->delta_y / 2.0;

              x_ij_box_min = source_x - P->delta_x / 2.0;
              x_ij_box_max = source_x + P->delta_x / 2.0;

              p_counter = 0;
              pressure = 0;
              sigma_1 = 0;
              sigma_3 = 0;
              mean = 0;
              differential = 0;

              runParticle = &refParticle; // start

              for (k = 0; k < numParticles;
                   k++) // look which particles are in the box
              {
                if (runParticle->ypos >= y_ij_box_min) {
                  if (runParticle->ypos <= y_ij_box_max) {
                    if (runParticle->xpos >= x_ij_box_min) {
                      if (runParticle->xpos <= x_ij_box_max) {
                        p_counter = p_counter + 1; // count particles

                        // set local pressure
                        pressure =
                            pressure +
                            sqrt((runParticle->F_P_x * runParticle->F_P_x) +
                                 (runParticle->F_P_y * runParticle->F_P_y));

                        // set stress tensor
                        sxx = runParticle->sxx;
                        syy = runParticle->syy;
                        sxy = runParticle->sxy;

                        // set pressure tensor
                        smax = ((sxx + syy) / 2.0) +
                               sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                                    sxy * sxy);
                        smin = ((sxx + syy) / 2.0) -
                               sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                                    sxy * sxy);

                        sigma_1 = sigma_1 + smax;
                        sigma_3 = sigma_3 + smin;

                        mean = mean + (smax + smin) / 2.0;

                        differential = differential + (smax - smin);
                      }
                    }
                  }
                }
                runParticle = runParticle->nextP;
              }
              pressure =
                  pressure / p_counter; // and divide by number of particles to
                                        // get average pressure
              sigma_1 = sigma_1 / p_counter;
              sigma_3 = sigma_3 / p_counter;
              mean = mean / p_counter;
              differential = differential / p_counter;

              fprintf(statStress, "%f\t", source_x);
              fprintf(statStress, "%f\t", source_y);
              fprintf(statStress, "%le\t", pressure);
              fprintf(statStress, "%le\t", sigma_1);
              fprintf(statStress, "%le\t", sigma_3);
              fprintf(statStress, "%le\t", mean);
              fprintf(statStress, "%le\t", differential);

              fprintf(statStress, " \n");
            }
          }
        }

        else {
          y_ij_box_min = y_mid - P->delta_y / 2.0;
          y_ij_box_max = y_mid + P->delta_y / 2.0;

          x_ij_box_min = x_mid - P->delta_x / 2.0;
          x_ij_box_max = x_mid + P->delta_x / 2.0;

          p_counter = 0;
          pressure = 0;
          sigma_1 = 0;
          sigma_3 = 0;
          mean = 0;
          differential = 0;

          runParticle = &refParticle; // start

          for (k = 0; k < numParticles;
               k++) // look which particles are in the box
          {
            if (runParticle->ypos >= y_ij_box_min) {
              if (runParticle->ypos <= y_ij_box_max) {
                if (runParticle->xpos >= x_ij_box_min) {
                  if (runParticle->xpos <= x_ij_box_max) {
                    p_counter = p_counter + 1; // count particles

                    // set local pressure
                    pressure = pressure +
                               sqrt((runParticle->F_P_x * runParticle->F_P_x) +
                                    (runParticle->F_P_y * runParticle->F_P_y));

                    // set stress tensor
                    sxx = runParticle->sxx;
                    syy = runParticle->syy;
                    sxy = runParticle->sxy;

                    // set pressure tensor
                    smax = ((sxx + syy) / 2.0) +
                           sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                                sxy * sxy);
                    smin = ((sxx + syy) / 2.0) -
                           sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                                sxy * sxy);

                    sigma_1 = sigma_1 + smax;
                    sigma_3 = sigma_3 + smin;

                    mean = mean + (smax + smin) / 2.0;

                    differential = differential + (smax - smin);
                  }
                }
              }
            }
            runParticle = runParticle->nextP;
          }
          pressure = pressure / p_counter; // and divide by number of particles
                                           // to get average pressure
          sigma_1 = sigma_1 / p_counter;
          sigma_3 = sigma_3 / p_counter;
          mean = mean / p_counter;
          differential = differential / p_counter;

          fprintf(statStress, "%f\t", x_mid);
          fprintf(statStress, "%f\t", y_mid);
          fprintf(statStress, "%le\t", pressure);
          fprintf(statStress, "%le\t", sigma_1);
          fprintf(statStress, "%le\t", sigma_3);
          fprintf(statStress, "%le\t", mean);
          fprintf(statStress, "%le\t", differential);

          fprintf(statStress, " \n");
        }
      }

      else {
        for (seal_y_mid = y_mid - 3.0 * P->delta_y;
             seal_y_mid <= y_mid + 3.0 * P->delta_y;
             seal_y_mid += 3.0 * P->delta_y) {
          y_ij_box_min = seal_y_mid - P->delta_y / 2.0;
          y_ij_box_max = seal_y_mid + P->delta_y / 2.0;

          x_ij_box_min = x_mid - P->delta_x / 2.0;
          x_ij_box_max = x_mid + P->delta_x / 2.0;

          p_counter = 0;
          pressure = 0;
          sigma_1 = 0;
          sigma_3 = 0;
          mean = 0;
          differential = 0;

          runParticle = &refParticle; // start

          for (k = 0; k < numParticles;
               k++) // look which particles are in the box
          {
            if (runParticle->ypos >= y_ij_box_min) {
              if (runParticle->ypos <= y_ij_box_max) {
                if (runParticle->xpos >= x_ij_box_min) {
                  if (runParticle->xpos <= x_ij_box_max) {
                    p_counter = p_counter + 1; // count particles

                    // set local pressure
                    pressure = pressure +
                               sqrt((runParticle->F_P_x * runParticle->F_P_x) +
                                    (runParticle->F_P_y * runParticle->F_P_y));

                    // set stress tensor
                    sxx = runParticle->sxx;
                    syy = runParticle->syy;
                    sxy = runParticle->sxy;

                    // set pressure tensor
                    smax = ((sxx + syy) / 2.0) +
                           sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                                sxy * sxy);
                    smin = ((sxx + syy) / 2.0) -
                           sqrt(((sxx - syy) / 2.0) * ((sxx - syy) / 2.0) +
                                sxy * sxy);

                    sigma_1 = sigma_1 + smax;
                    sigma_3 = sigma_3 + smin;

                    mean = mean + (smax + smin) / 2.0;

                    differential = differential + (smax - smin);
                  }
                }
              }
            }
            runParticle = runParticle->nextP;
          }
          pressure = pressure / p_counter; // and divide by number of particles
                                           // to get average pressure
          sigma_1 = sigma_1 / p_counter;
          sigma_3 = sigma_3 / p_counter;
          mean = mean / p_counter;
          differential = differential / p_counter;

          fprintf(statStress, "%f\t", x_mid);
          fprintf(statStress, "%f\t", seal_y_mid);
          fprintf(statStress, "%le\t", pressure);
          fprintf(statStress, "%le\t", sigma_1);
          fprintf(statStress, "%le\t", sigma_3);
          fprintf(statStress, "%le\t", mean);
          fprintf(statStress, "%le\t", differential);

          fprintf(statStress, " \n");
        }
      }
    }

    y_mid += 0.1;
  }

  fprintf(statStress, "\n");
  fclose(statStress); // close file
}

void Lattice::Real_Density_Young(float seal_y_min, float seal_y_max,
                                 double dens_seal, double dens,
                                 double yong_seal, double yong, int scal)
// scal in meters
{
  // int i;
  // dens_a = 2500 (L.st)
  // dens_b = 2600 (shale)
  // dens_c = 2500 (L.st)

  // you_a = 54e+09 (L.St)
  // you_b = 22e+09 (shale)
  // you_c = 54e+09 (L.St)

  // real_density = (1-P->porosity[i][j])*dens_fluid +
  // P->porosity[i][j]*dens_solid;

  // runParticle = &refParticle;

  // for (i=0; i<numParticles; i++)
  //{
  // if (runParticle->ypos >= y_min && runParticle->ypos <= y_max)
  //{
  // runParticle->real_density = dens;
  //}

  // runParticle = runParticle->nextP;
  //}

  int i, j;

  runParticle = &refParticle;

  for (i = 0; i < numParticles; i++) {
    // if (( runParticle->ypos >= seal_y_min ) && (runParticle->ypos <=
    // seal_y_max))
    //{
    // runParticle->real_density = dens_seal;
    // runParticle->real_young = yong_seal;
    //}
    // else
    //{
    // runParticle->real_density = dens;
    // runParticle->real_young = yong;
    //}

    runParticle->real_density = 2700;
    // runParticle->real_young = 1.0e+10;
    runParticle->real_young =
        runParticle->young * pressure_scale * pascal_scale;
    runParticle->real_radius = runParticle->radius * scal;

    // cout << runParticle->nb << " " << runParticle->real_radius << endl;

    runParticle = runParticle->nextP;
  }
}

void Lattice::Press_Initial_Gravity(int depth, int height) {
  int i;
  double stress_g, strain_g;
  float move;
  float wall_pos; // wall position

  // SetWallBoundaries(0, 1.0);

  /******************************************************
   *  F=P/A, E=stress/strain, stress=F/A, F=mg
   * strain = stress/E = F/A*E = m*g/A*E = rho*V*g/A*E = rho*g*H*A/A*E =
   * rho*g*H/E = rho*g*r/E
   * ****************************************************/

  cout << " Gravity in terms of Uniaxial Compression" << endl;

  runParticle = &refParticle; // start of list

  runParticle = runParticle->prevP; // go back one

  // wall_pos = runParticle->ypos+ runParticle->radius;	// define wall position

  /**************************************************
   * dl = (stress*L)/E = (F/A)*(L/E) = (m*g/A)*(L/E) = (rho*V*g/A)*(L/E)
   * ************************************************/
  move = 9.81 * (height + depth) * runParticle->real_density /
         runParticle->real_young;

  move = move * particlex * 0.01; //???

  move = move * 0.5;

  for (i = 0; i < numParticles; i++) {
    runParticle->ypos = runParticle->ypos - (runParticle->ypos * move);

    runParticle->ChangeBox(repBox, particlex);

    runParticle = runParticle->prevP;
  }
  // for (i=0; i < numParticles; i++) // go through all particles
  //{
  ////cout << runParticle->nb << endl;

  // if(i < particlex)
  // stress_g = runParticle->real_density * 9.8 * (runParticle->real_radius +
  // depth);
  // else
  // stress_g = runParticle->real_density * 9.8 * runParticle->real_radius;

  // strain_g = stress_g / runParticle->real_young;

  /*****************************************************
   * strain = dl/l_u = (l_d - l_u)/l_u   ->  -ve is just to cancel the sign
   * convention of strain-g
   * **************************************************/
  // runParticle->ypos = runParticle->ypos -
  //(runParticle->ypos * strain_g);

  // runParticle->ChangeBox (repBox, particlex); // always call Repbox if
  // position changes

  // runParticle = runParticle->prevP;
  //}

  // stress_g = 2700 * 9.8 * (runParticle->radius + depth);
  // stress_g = 2700 * 9.8 * runParticle->real_radius;
  // strain_g = stress_g / 3.0e+10;

  // upper_wall_pos = upper_wall_pos - wall_pos * strain_g;

  // SetRightLatticeWallPos();
}

void Lattice::SetSpringConstants() {

  double factor = sqrt(3.0) / 2.0;

  runParticle = &refParticle;

  for (int i = 0; i < numParticles; i++) {
    cout << runParticle->young * factor << endl;
    runParticle->springconst = runParticle->young * factor;

    for (int j = 0; j < 8; j++)
      runParticle->springf[j] = runParticle->young * factor;

    runParticle = runParticle->nextP;
  }
}

void Lattice::CalculateGravityOverburden(float depth, float height) {
  int i;
  double stress_overburden;

  for (i = 0; i < numParticles; i++) {
    stress_overburden = runParticle->real_density * 9.81 *
                        (depth - (runParticle->xpos * height));

    // cout << stress_overburden << endl;

    runParticle->gravforce =
        stress_overburden / (pressure_scale * pascal_scale);
    cout << runParticle->gravforce << endl;

    runParticle = runParticle->nextP;
  }
}

void Lattice::SetWeightForce() {
  runParticle = &refParticle;
  double scaled_stress;

  for (int i = 0; i < numParticles; i++) {
    //~ scaled_stress =
    //(runParticle->real_density*9.81*2.0*runParticle->real_radius) ~
    //*(runParticle->young/runParticle->real_young);
    //~ runParticle->fg =
    //scaled_stress*Pi*runParticle->radius*runParticle->radius / sqrt(3.0);

    //~ runParticle->fg =
    //Pi*pow(runParticle->radius, 2.0)*runParticle->real_density ~
    //*9.81*(1.0/runParticle->real_young)*10000 *2.0/3.0;

    // runParticle->fg =
    // Pi*pow(runParticle->radius, 2.0)*runParticle->real_density
    //			*9.81*(1.0/runParticle->real_young)*10000 *4.0/5.0;
    runParticle->fg = runParticle->real_radius * 2 * runParticle->real_density *
                      9.81 /
                      runParticle->real_young; // pascal_scale*pressure_scale
    runParticle->fg = runParticle->fg * runParticle->radius * 2.0; //?
    // cout << runParticle->fg << endl;
    //~ runParticle->fg = 2*sqrt(3)*runParticle->radius*runParticle->real_radius
    //~ *runParticle->real_density*9.81*(1.0/runParticle->real_young)*4.0/5.0;

    runParticle = runParticle->nextP;
  }
}

void Lattice::Gravity(float depth, float density) {
  int j;
  double weigth, move;

  weigth = depth * density * 9.81 / runParticle->real_young;
  // weigth = weigth * runParticle->radius * 2.0;

  cout << weigth << endl;

  for (j = 0; j < numParticles; j++) {
    runParticle->fg =
        runParticle->real_radius * 2 *
        (((1 - runParticle->average_porosity) * runParticle->real_density) +
         (runParticle->average_porosity * 1000)) *
        9.81 / runParticle->real_young; // pascal_scale*pressure_scale
    runParticle->fg = runParticle->fg * runParticle->radius * 2.0;
    if (runParticle->is_top_lattice_boundary)

      runParticle->fg = runParticle->fg + weigth;

    runParticle = runParticle->nextP;
  }

  /**************************************************
   * dl = (stress*L)/E = (F/A)*(L/E) = (m*g/A)*(L/E) = (rho*V*g/A)*(L/E)
   * ************************************************/

  move = 9.81 * (runParticle->real_radius * 2.0 * particlex) *
         runParticle->real_density / runParticle->real_young;

  move = move;

  move = move + weigth;

  cout << move << endl;

  for (i = 0; i < numParticles; i++) {
    runParticle->ypos = runParticle->ypos - (runParticle->ypos * move);

    runParticle->ChangeBox(repBox, particlex);

    runParticle = runParticle->prevP;
  }
}

void Lattice::Gravity_noload_nofluid() {
  int j;
  double weigth, move;

  for (j = 0; j < numParticles; j++) {
    runParticle->fg = runParticle->real_radius * 2 * runParticle->real_density *
                      9.81 /
                      runParticle->real_young; // pascal_scale*pressure_scale
    runParticle->fg = runParticle->fg * runParticle->radius * 2.0 * 1.5;

    runParticle = runParticle->nextP;
  }

  /**************************************************
   * dl = (stress*L)/E = (F/A)*(L/E) = (m*g/A)*(L/E) = (rho*V*g/A)*(L/E)
   * ************************************************/

  move = 9.81 * (runParticle->real_radius * 2.0 * particlex) *
         runParticle->real_density / runParticle->real_young;

  move = move * 0.5 * 1.5;

  cout << move << endl;

  for (i = 0; i < numParticles; i++) {
    runParticle->ypos = runParticle->ypos - (runParticle->ypos * move);

    runParticle->ChangeBox(repBox, particlex);

    runParticle->oldy = runParticle->ypos;

    runParticle = runParticle->prevP;
  }
}

void Lattice::Gravity_no_Move(float depth, float density) {
  int j;
  double weigth, move;

  weigth = depth * density * 9.81 * 0.5 / runParticle->real_young;

  cout << weigth << endl;

  for (j = 0; j < numParticles; j++) {
    runParticle->fg = runParticle->real_radius * 2 * runParticle->real_density *
                      9.81 /
                      runParticle->real_young; // pascal_scale*pressure_scale
    runParticle->fg = runParticle->fg * runParticle->radius * 2.0;
    if (runParticle->is_top_lattice_boundary)

      runParticle->fg = runParticle->fg + weigth;

    runParticle = runParticle->nextP;
  }
}

void Lattice::Adjust_Gravity() {
  int j;

  for (j = 0; j < numParticles; j++) {
    if (!runParticle->is_top_lattice_boundary) {
      runParticle->fg =
          runParticle->real_radius * 2 *
          (((1 - runParticle->average_porosity) * runParticle->real_density) +
           (runParticle->average_porosity * 1000)) *
          9.81 / runParticle->real_young; // pascal_scale*pressure_scale
      runParticle->fg = runParticle->fg * runParticle->radius * 2.0;
    }
    runParticle = runParticle->nextP;
  }
}

void Lattice::OutputMohoParticles(double time) {
  ofstream file;
  file.open("Moho_Particles_XY", ios::app);

  int i;
  double size;

  runParticle = firstParticle->prevP;

  for (i = 0; i < numParticles; i++) {
    if (runParticle->is_moho_particle) {
      file << runParticle->xpos << " " << time << " " << runParticle->ypos
           << endl;
    }
    runParticle = runParticle->nextP;
  }

  file << "\n";

  file.close();
}
void Lattice::MarkMohoParticles(double y) {

  int i;
  bool y_is_in_range = false;

  runParticle = &refParticle; // first particle

  for (i = 0; i < numParticles; i++) {
    if (runParticle->ypos < y)
      runParticle = runParticle->nextP;
    else {
      y_is_in_range = true;
      break;
    }
  }

  if (y_is_in_range) {
    for (i = 0; i < particlex; i++) {
      runParticle->is_moho_particle = true;
      runParticle = runParticle->nextP;
    }
  }
}

void Lattice::Healing(double heal_dist, float prob, int time_heal,
                      float springfactor, float springstrength) {
  int i, j, ii, jj, k, l, ll, kk, myPlace, neigPlace;
  int box, pos;
  int connected_neighbors;
  Particle *neig;
  double dx, dy, dd;
  float ran_nb, newprob;
  double average_springf, average_springs, average_break_Str;

  srand(time(0));

  for (i = 0; i < numParticles; i++) {
    if (runParticle->draw_break) // is bond fractured?
    {
      runParticle->done = true;

      connected_neighbors = 0;

      for (j = 0; j < 8; j++) {
        if (runParticle->neigP[j]) {
          connected_neighbors++;
          runParticle->neigP[j]->done = true;
        }
      }

      //---------------------------------------------------
      // nine positions in 2d box 3x3
      //---------------------------------------------------

      for (ii = 0; ii < 3; ii++) {
        for (jj = 0; jj < 3; jj++) {

          //------------------------------------------
          //  case  upper  or  lower  row (+-size)
          //------------------------------------------

          if (jj == 0)
            box = 2 * particlex * (-1);
          if (jj == 1)
            box = 0;
          if (jj == 2)
            box = 2 * particlex;
          //-----------------------------------------
          //  now  check  this  position
          //  all  particles  in  this  position !!
          //-----------------------------------------

          pos = runParticle->box_pos + (ii - 1) + box;

          //-----------------------------------------
          // repbox should always be positive
          //-----------------------------------------

          if (pos > 0) {
            neig = repBox[pos]; // helppointer  for  while  loop

            while (neig) // as long as there is somebody
            {

              if (!neig->done) // if not done already
              {
                dx = neig->xpos - runParticle->xpos;
                dy = neig->ypos - runParticle->ypos;

                dd = sqrt((dx * dx) + (dy * dy));

                cout << dd << endl;

                if (dd == 0)
                  cout << dd << "problem" << endl;

                if ((dd / heal_dist) < runParticle->radius * 2.0 && dd > 0)
                // if (dd < runParticle->radius*3.0)
                // if ((dd/heal_dist) < runParticle->radius*2.0)
                // if (dd < runParticle->radius*3.0)
                {
                  ran_nb = rand() / (float)RAND_MAX;

                  newprob = (runParticle->radius * 2.0 / dd) * prob;

                  // this should include the new exponential as a function of T
                  // added to this

                  cout << ran_nb << newprob << endl;

                  if (ran_nb < newprob) // now we are actually healing
                  {
                    if (connected_neighbors < 9) {
                      connected_neighbors++;

                      // get averages first

                      ll = 0;
                      average_springf = 0;
                      average_springs = 0;
                      average_break_Str = 0;

                      for (kk = 0; kk < 8; kk++) {
                        if (runParticle->neigP[kk]) {
                          ll++;
                          average_springf = +runParticle->springf[kk];
                          average_springs = +runParticle->springs[kk];
                          average_break_Str = +runParticle->break_Str[kk];
                        }
                        if (neig->neigP[kk]) {
                          ll++;
                          average_springf = +neig->springf[kk];
                          average_springs = +runParticle->springs[kk];
                          average_break_Str = +neig->break_Str[kk];
                        }
                      }

                      average_springf = average_springf / ll;
                      average_springs = average_springs / ll;
                      average_break_Str = average_break_Str / ll;

                      // find empty positions

                      for (k = 0; k < 8; k++) {
                        if (!runParticle->neigP[k]) {
                          myPlace = k;
                        }
                      }

                      for (l = 0; l < 8; l++) {
                        if (!neig->neigP[l]) {
                          neigPlace = l;
                        }
                      }

                      // cout << average_springf << endl;
                      // cout << " Young " << runParticle->young << endl;

                      runParticle->neig_spring[myPlace] = neigPlace;
                      neig->neig_spring[neigPlace] = myPlace;

                      runParticle->neigP[myPlace] = neig;
                      neig->neigP[neigPlace] = runParticle;

                      runParticle->springf[myPlace] =
                          average_springf * springfactor;
                      runParticle->springs[myPlace] =
                          average_springs * springfactor;
                      runParticle->break_Str[myPlace] =
                          average_break_Str * springstrength;

                      neig->springf[neigPlace] = average_springf * springfactor;
                      neig->springs[neigPlace] = average_springs * springfactor;
                      neig->break_Str[neigPlace] =
                          average_break_Str * springstrength;

                      runParticle->rad[myPlace] = dd / 2.0;
                      neig->rad[neigPlace] = dd / 2.0;

                      runParticle->neigpos[myPlace][0] =
                          neig->xpos - runParticle->xpos;
                      runParticle->neigpos[myPlace][1] =
                          neig->ypos - runParticle->ypos;

                      neig->neigpos[neigPlace][0] =
                          runParticle->xpos - neig->xpos;
                      neig->neigpos[neigPlace][1] =
                          runParticle->ypos - neig->ypos;

                      // cout << "initial rad_par_fluid  " <<
                      // runParticle->rad_par_fluid <<endl; cout <<
                      // runParticle->rad[k]/runParticle->radius << endl;

                      // runParticle->rad_par_fluid=runParticle->rad_par_fluid+(runParticle->rad_par_fluid*((runParticle->rad[k]/runParticle->radius)-1));
                      // runParticle->area_par_fluid =
                      // pow((runParticle->rad_par_fluid),2.0)*PI;

                      // neig->rad_par_fluid=neig->rad_par_fluid+(neig->rad_par_fluid*((neig->rad[l]/neig->radius)-1));
                      // neig->area_par_fluid =
                      // pow((neig->rad_par_fluid),2.0)*PI;

                      /*
                      cout <<  "rad_par_fluid  " << runParticle->rad_par_fluid
                      << endl;

                      if
                      (runParticle->area_par_fluid/(pow(2*runParticle->radius,2.0))>0.99)
                              {
                                      runParticle->area_par_fluid =
                      (pow(2*runParticle->radius,2.0)*0.99);
                                      runParticle->rad_par_fluid =
                      sqrt((runParticle->area_par_fluid)/PI);
                              }
                      if (neig->area_par_fluid/(pow(2*neig->radius,2.0))>0.99)
                              {
                                      neig->area_par_fluid =
                      (pow(2*neig->radius,2.0)*0.99); neig->rad_par_fluid =
                      sqrt((neig->area_par_fluid)/PI);
                              }*/

                      // runParticle->rad[myPlace] = runParticle->radius;
                      // neig->rad[neigPlace] = neig->radius;

                      runParticle->draw_break = false;
                      neig->draw_break = false;

                      runParticle->healing = time_heal;
                      neig->healing = time_heal;

                      runParticle->fracture_time = 0;
                      neig->fracture_time = 0;

                      cout << "healing" << endl;

                      nbBreak = nbBreak - 1;
                    }
                  }
                }
              }

              if (neig->next_inBox) // if there is anotherone in box
              {
                neig = neig->next_inBox; // shift the help pointer
              } else {
                break; // if not go out to next box
              }
            }
          }
        }
      }
      runParticle->done = false; // myself

      for (j = 0; j < 8; j++) // and all my neighbours
      {
        if (runParticle->neigP[j]) {
          runParticle->neigP[j]->done = false;
        }
      }
    }
    runParticle = runParticle->nextP;
  }
}

//----------------------------------------------------------------------------------------------------------
//
// new healing function 2025 based on Renaud. Uses Atomic fluctuations on the
// small scale. Probability is multiplied times time step (can be included in
// prob) and exponential of activation energy divided by Boltzmann and T in
// Kelvin. Energy is Young Modulus times volume (and can be a function if
// distance of bonds or particles relativ to each other, at the moment that is
// out of the exponential.
//
// Healing should increase with time given and with higher local temperature (T
// of particle is used)
//-------------------------------------------------------------------------------------------------------------

void Lattice::HealingTMol(double heal_dist, float prob, int time_heal,
                          float springfactor, float springstrength) {
  int i, j, ii, jj, k, l, ll, kk, myPlace, neigPlace;
  int box, pos;
  int connected_neighbors;
  Particle *neig;
  double dx, dy, dd;
  float ran_nb, newprob;
  double average_springf, average_springs, average_break_Str, factor;

  srand(time(0));

  for (i = 0; i < numParticles; i++) {
    if (runParticle->draw_break) // is bond fractured?
    {
      runParticle->done = true;

      connected_neighbors = 0;

      for (j = 0; j < 8; j++) {
        if (runParticle->neigP[j]) {
          connected_neighbors++;
          runParticle->neigP[j]->done = true;
        }
      }

      //---------------------------------------------------
      // nine positions in 2d box 3x3
      //---------------------------------------------------

      for (ii = 0; ii < 3; ii++) {
        for (jj = 0; jj < 3; jj++) {

          //------------------------------------------
          //  case  upper  or  lower  row (+-size)
          //------------------------------------------

          if (jj == 0)
            box = 2 * particlex * (-1);
          if (jj == 1)
            box = 0;
          if (jj == 2)
            box = 2 * particlex;
          //-----------------------------------------
          //  now  check  this  position
          //  all  particles  in  this  position !!
          //-----------------------------------------

          pos = runParticle->box_pos + (ii - 1) + box;

          //-----------------------------------------
          // repbox should always be positive
          //-----------------------------------------

          if (pos > 0) {
            neig = repBox[pos]; // helppointer  for  while  loop

            while (neig) // as long as there is somebody
            {

              if (!neig->done) // if not done already
              {
                dx = neig->xpos - runParticle->xpos;
                dy = neig->ypos - runParticle->ypos;

                dd = sqrt((dx * dx) + (dy * dy));

                // cout << dd << endl;

                if (dd == 0)
                  cout << dd << "problem" << endl;

                if ((dd / heal_dist) < runParticle->radius * 2.0 && dd > 0)

                {
                  ran_nb = rand() / (float)RAND_MAX;

                  newprob = (runParticle->radius * 2.0 / dd) * prob;

                  // this should include the new exponential as a function of T
                  // added to this

                  factor = exp(
                      -(2) / (runParticle->SheetT)); // 2 is Volume times Young
                                                     // divided by Boltzmann

                  newprob = newprob * factor;

                  // cout << factor << endl;

                  // cout << ran_nb << newprob << endl;

                  if (ran_nb < newprob) // now we are actually healing
                  {
                    if (connected_neighbors < 9) {
                      connected_neighbors++;

                      // get averages first

                      ll = 0;
                      average_springf = 0;
                      average_springs = 0;
                      average_break_Str = 0;

                      for (kk = 0; kk < 8; kk++) {
                        if (runParticle->neigP[kk]) {
                          ll++;
                          average_springf = +runParticle->springf[kk];
                          average_springs = +runParticle->springs[kk];
                          average_break_Str = +runParticle->break_Str[kk];
                        }
                        if (neig->neigP[kk]) {
                          ll++;
                          average_springf = +neig->springf[kk];
                          average_springs = +runParticle->springs[kk];
                          average_break_Str = +neig->break_Str[kk];
                        }
                      }

                      average_springf = average_springf / ll;
                      average_springs = average_springs / ll;
                      average_break_Str = average_break_Str / ll;

                      // find empty positions

                      for (k = 0; k < 8; k++) {
                        if (!runParticle->neigP[k]) {
                          myPlace = k;
                        }
                      }

                      for (l = 0; l < 8; l++) {
                        if (!neig->neigP[l]) {
                          neigPlace = l;
                        }
                      }

                      runParticle->neig_spring[myPlace] = neigPlace;
                      neig->neig_spring[neigPlace] = myPlace;

                      runParticle->neigP[myPlace] = neig;
                      neig->neigP[neigPlace] = runParticle;

                      runParticle->springf[myPlace] =
                          average_springf * springfactor;
                      runParticle->springs[myPlace] =
                          average_springs * springfactor;
                      runParticle->break_Str[myPlace] =
                          average_break_Str * springstrength;

                      neig->springf[neigPlace] = average_springf * springfactor;
                      neig->springs[neigPlace] = average_springs * springfactor;
                      neig->break_Str[neigPlace] =
                          average_break_Str * springstrength;

                      runParticle->rad[myPlace] = dd / 2.0;
                      neig->rad[neigPlace] = dd / 2.0;

                      runParticle->neigpos[myPlace][0] =
                          neig->xpos - runParticle->xpos;
                      runParticle->neigpos[myPlace][1] =
                          neig->ypos - runParticle->ypos;

                      neig->neigpos[neigPlace][0] =
                          runParticle->xpos - neig->xpos;
                      neig->neigpos[neigPlace][1] =
                          runParticle->ypos - neig->ypos;

                      runParticle->draw_break = false;
                      neig->draw_break = false;

                      runParticle->healing = time_heal;
                      neig->healing = time_heal;

                      runParticle->fracture_time = 0;
                      neig->fracture_time = 0;

                      cout << "healing" << endl;

                      nbBreak = nbBreak - 1;
                    }
                  }
                }
              }

              if (neig->next_inBox) // if there is anotherone in box
              {
                neig = neig->next_inBox; // shift the help pointer
              } else {
                break; // if not go out to next box
              }
            }
          }
        }
      }
      runParticle->done = false; // myself

      for (j = 0; j < 8; j++) // and all my neighbours
      {
        if (runParticle->neigP[j]) {
          runParticle->neigP[j]->done = false;
        }
      }
    }
    runParticle = runParticle->nextP;
  }
}

void Lattice::DumpStatisticPorosity(double y_box_min, double y_box_max,
                                    double x_box_min, double x_box_max,
                                    double strain) {
  FILE *stat; // file pointer

  int i, j, jj, p_counter;                    // counters
  float poro, perm, fluid_Pressure, gradient; // mean stress, differential
  // stress, pressure

  float finite_strain, breaknumber; // finite strain within the box of
  // interest

  float wall_pos, wall_Xpos;

  runParticle = &refParticle;       // start of list
  runParticle = runParticle->prevP; // go back one
  wall_pos = runParticle->ypos;     // get wall position y
  wall_Xpos = runParticle->xpos;    // get wall position x

  for (i = 0; i < particlex; i++) {
    runParticle = runParticle->prevP;
  }

  if (runParticle->xpos > wall_Xpos) {
    wall_Xpos = runParticle->xpos;
  }

  p_counter = 0;

  perm = 0.0;
  poro = 0.0;

  // finite_strain = strain * (def_time_u + def_time_p);	// calculate
  // finite_strain = def_time_p;

  runParticle = &refParticle; // start

  for (i = 0; i < numParticles; i++) // look which particles are in the
                                     // box
  {
    if (runParticle->ypos > y_box_min * wall_pos) // larger ymin
    {
      if (runParticle->ypos < y_box_max * wall_pos) // smaller ymax
      {
        if (runParticle->xpos > x_box_min * wall_Xpos) // larger xmin
        {
          if (runParticle->xpos < x_box_max * wall_Xpos) // smaller xmax
          {
            p_counter = p_counter + 1; // count particles

            // set stress tensor

            perm = perm + runParticle->average_permeability;
            poro = poro + runParticle->average_porosity;
            fluid_Pressure = fluid_Pressure + runParticle->temperature;
            gradient = gradient + runParticle->fluid_pressure_gradient;
          }
        }
      }
    }
    runParticle = runParticle->nextP;
  }

  // and divide by number of particles

  perm = perm / p_counter;
  gradient = gradient / p_counter;
  poro = poro / p_counter;
  fluid_Pressure = fluid_Pressure / p_counter;
  finite_strain = strain;
  breaknumber = nbBreak;

  stat = fopen("statisticPerm.txt", "a"); // open statistic output
  // append file

  // and dump the data

  fprintf(stat, "%f", finite_strain);
  fprintf(stat, ";%f", breaknumber);
  fprintf(stat, ";%f", fluid_Pressure);
  fprintf(stat, ";%f", gradient);
  fprintf(stat, ";%f", perm);
  fprintf(stat, ";%f\n", poro);

  fclose(stat); // close file
}

/***************************************
 * Advection Diffusion
 * **************************/

void Lattice::Initialize_Concentration_Lattice(double Press_Background) {

  P->setconcent(Press_Background);
}

void Lattice::Calculate_Concentration(double vel, double dif) { // int n;

  // P->alpha_beta_set(1.0);

  // P->dan_matrices_set_1();
  // P->dan_multiplication_1st();
  // P->test_pressure_with_source(l); // adding the source term (test please)
  // P->transpose_trans_P();
  // P->dan_invert();
  // P->dan_multiplication_inv();

  // P->dan_matrices_set_2();
  // P->dan_multiplication_1st();
  // P->transpose_trans_P();
  // P->test_pressure_with_source(l); // adding the source term (test please)
  // P->dan_invert();
  // P->dan_multiplication_inv();

  // P->SolveadvDisExplicitSteady();
  // P->Test();

  // cout<<" I am here: Advection+Diffusion!! "<<endl;

  P->SolveadvDisExplicitTransient(vel, dif);

  // P->concentration_particles(repBox);
}
void Lattice::Interpolation_Concentration() { // int n;

  // P->alpha_beta_set(1.0);

  // P->dan_matrices_set_1();
  // P->dan_multiplication_1st();
  // P->test_pressure_with_source(l); // adding the source term (test please)
  // P->transpose_trans_P();
  // P->dan_invert();
  // P->dan_multiplication_inv();

  // P->dan_matrices_set_2();
  // P->dan_multiplication_1st();
  // P->transpose_trans_P();
  // P->test_pressure_with_source(l); // adding the source term (test please)
  // P->dan_invert();
  // P->dan_multiplication_inv();

  // P->SolveadvDisExplicitSteady();
  // P->Test();

  // P->SolveadvDisExplicitTransient() ;

  P->concentration_particles(repBox);
}

void Lattice::Get_Concentration_Directly() { P->get_concentration(repBox); }

void Lattice::SetDarcyforAdv() { P->setDarcyVelocity(repBox); }

void Lattice::Set_Calcite() {
  int i;

  for (i = 0; i < numParticles; i++) // look which particles are in the
  // box
  {
    runParticle->calcite = 0.0;
    runParticle->mV = 0.00004;

    runParticle = runParticle->nextP;
  }
}

// volume input is the proposed volume change for the reaction. This should be
// adjusted to the concentration change of 0.5

void Lattice::Make_Dolomite(double reaction_time, int step, double shrink,
                            double volume, double scale) {
  int i, box_x, box_y, j;
  double ratio, dt, change, react, v_change, particle_calcite;
  float step_model;

  double rate_constant = 0.0001; // calcite see stylo paper

  for (i = 0; i < numParticles; i++) {

    if (runParticle->conc > runParticle->calcite) {
      change =
          runParticle->conc -
          runParticle->calcite; // max is 0.5, should represent volume change

      if (change > 0.0) {
        particle_calcite = 0.0;

        react = -1 * reaction_time * runParticle->mV * rate_constant *
                (1 - (runParticle->conc / 0.01));

        // calculate how much of the particle changes using the scale

        if (react > 0.0) {
          v_change = react / (2.0 * runParticle->rad_par_fluid * scale);

          particle_calcite = runParticle->calcite;

          runParticle->calcite = runParticle->calcite + (0.5 * v_change);

          if (runParticle->calcite > 0.5)
            runParticle->calcite = 0.5;

          particle_calcite =
              (runParticle->calcite - particle_calcite) /
              0.5; // ratio of particle changed for volume and fracturing
        }

        // if (react > 0.0)
        //	cout << react << endl;

        if (volume != 0.0) // porosity change in percent, negative is porosity
                           // increase, positive is porosity decrease
        {

          runParticle->rad_par_fluid =
              runParticle->rad_par_fluid +
              runParticle->rad_par_fluid * (particle_calcite * volume);

          // v_change = 1 - runParticle->area_par_fluid / (3.14 *
          // pow(runParticle->rad_par_fluid,2.0));

          // cout << v_change << endl;

          // v_change = v_change / volume;

          runParticle->area_par_fluid =
              3.14 * pow(runParticle->rad_par_fluid, 2.0);
        }
        if (shrink != 0.0) {
          runParticle->radius =
              runParticle->radius +
              runParticle->radius * (particle_calcite * shrink); // shrink

          ratio = 1 - ((particle_calcite * shrink) * -1.0);
          for (j = 0; j < 8; j++) {
            if (runParticle->neigP[j]) {
              runParticle->neigpos[j][0] *= ratio; // change
              runParticle->neigpos[j][1] *= ratio;
            }
          }
        }
      } else
        cout << "negative crystal" << endl;

      // runParticle->calcite = runParticle->calcite + (0.5*v_change);
      v_change = 0.0;
    }
    runParticle = runParticle->prevP;
  }

  // volume input is the proposed volume change for the reaction. This should be
  // adjusted to the concentration change of 0.5
}

void Lattice::Reaction_Arrhenius(double reaction_time, double scale) {
  int i, box_x, box_y, j;
  double ratio, dt, change, react, v_change, particle_calcite;
  float step_model;

  for (i = 0; i < numParticles; i++) {
    particle_calcite = 0.0;

    react =
        reaction_time * 1000 * runParticle->mV * (exp(-1 / runParticle->conc));
    // cout << react << endl;
    //  calculate how much of the particle changes using the scale

    if (react > 0.0) {
      v_change = react / (2.0 * runParticle->rad_par_fluid * scale);

      particle_calcite = runParticle->calcite;

      runParticle->calcite = runParticle->calcite + (0.5 * v_change);

      if (runParticle->calcite > 0.5)
        runParticle->calcite = 0.5;

      particle_calcite =
          (runParticle->calcite - particle_calcite) /
          0.5; // ratio of particle changed for volume and fracturing
    }

    // if (react > 0.0)
    //	cout << react << endl;
    /*
    if (volume!=0.0)   // porosity change in percent, negative is porosity
    increase, positive is porosity decrease
    {

            runParticle->rad_par_fluid = runParticle->rad_par_fluid +
    runParticle->rad_par_fluid * (particle_calcite * volume);

            //v_change = 1 - runParticle->area_par_fluid / (3.14 *
    pow(runParticle->rad_par_fluid,2.0));

            //cout << v_change << endl;

            //v_change = v_change / volume;

            runParticle->area_par_fluid = 3.14 *
    pow(runParticle->rad_par_fluid,2.0);
    }
    if (shrink!=0.0)
    {
            runParticle->radius = runParticle->radius + runParticle->radius *
    (particle_calcite * shrink);	// shrink

            ratio = 1-((particle_calcite*shrink)*-1.0);
            for (j = 0; j < 8; j++)
            {
                    if (runParticle->neigP[j])
                    {
                            runParticle->neigpos[j][0] *= ratio;	//
    change runParticle->neigpos[j][1] *= ratio;

                    }
            }
    }*/

    runParticle = runParticle->prevP;
  }

  // volume input is the proposed volume change for the reaction. This should be
  // adjusted to the concentration change of 0.5
}

void Lattice::Make_Metal(double reaction_time, int step, double shrink,
                         double volume, double scale) {
  int i, box_x, box_y, j;
  double ratio, dt, change, react, v_change, particle_calcite, conc;
  float step_model;

  double rate_constant = 0.0001; // calcite see stylo paper

  for (i = 0; i < numParticles; i++) {
    if (runParticle->conc > 0.0) {
      if (runParticle->metal_conc > 0.0) {
        if (runParticle->conc > runParticle->metal_conc)
          conc = runParticle->metal_conc;
        else
          conc = runParticle->conc;
        if ((conc) > runParticle->calcite) {
          change = (conc)-runParticle
                       ->calcite; // max is 0.5, should represent volume change

          if (change > 0.0) {
            particle_calcite = 0.0;

            react = -1 * reaction_time * runParticle->mV * rate_constant *
                    (1 - ((conc) / 0.01));

            // calculate how much of the particle changes using the scale

            if (react > 0.0) {
              v_change = react / (2.0 * runParticle->rad_par_fluid * scale);

              particle_calcite = runParticle->calcite;

              runParticle->calcite = runParticle->calcite + (0.5 * v_change);

              if (runParticle->calcite > 0.5)
                runParticle->calcite = 0.5;

              particle_calcite =
                  (runParticle->calcite - particle_calcite) /
                  0.5; // ratio of particle changed for volume and fracturing
            }

            // if (react > 0.0)
            //	cout << react << endl;

            if (volume !=
                0.0) // porosity change in percent, negative is porosity
                     // increase, positive is porosity decrease
            {

              runParticle->rad_par_fluid =
                  runParticle->rad_par_fluid +
                  runParticle->rad_par_fluid * (particle_calcite * volume);

              // v_change = 1 - runParticle->area_par_fluid / (3.14 *
              // pow(runParticle->rad_par_fluid,2.0));

              // cout << v_change << endl;

              // v_change = v_change / volume;

              runParticle->area_par_fluid =
                  3.14 * pow(runParticle->rad_par_fluid, 2.0);
            }
            if (shrink != 0.0) {
              runParticle->radius =
                  runParticle->radius +
                  runParticle->radius * (particle_calcite * shrink); // shrink

              ratio = 1 - ((particle_calcite * shrink) * -1.0);
              for (j = 0; j < 8; j++) {
                if (runParticle->neigP[j]) {
                  runParticle->neigpos[j][0] *= ratio; // change
                  runParticle->neigpos[j][1] *= ratio;
                }
              }
            }
          } else
            cout << "negative crystal" << endl;

          // runParticle->calcite = runParticle->calcite + (0.5*v_change);
          v_change = 0.0;
        }
      }
    }
    runParticle = runParticle->prevP;
  }

  /*
          dt = reaction_time;

          step_model = (step/2.0) - int(step/2.0);

          cout << step_model << endl;

          for (i = 0 ; i < 13; i++)
                  runParticle = runParticle->prevP;

          //for (j = 0; j < 10 ; j++)
          //do {
                  react = 0.0;
                          for (i = 0; i < numParticles; i++)	// look which
     particles are in the
      // box
                          {
                                  if (runParticle->calcite < 0.5)
                                  {
                                          if
     ((runParticle->conc-runParticle->calcite)*dt > react)
                                          {
                                                  react =
     (runParticle->conc-runParticle->calcite)*dt; preRunParticle = runParticle;
                                          }
                                  }
                                  if (step_model > 0.0)
                                          runParticle = runParticle->nextP;
                                  else
                                          runParticle = runParticle->prevP;
                          }
                          if (react > 0.3)
                          {
                          runParticle = preRunParticle;
                          ratio = runParticle->calcite;
                          runParticle->calcite = (runParticle->calcite +
     ((runParticle->conc-runParticle->calcite)*dt));

                          if (runParticle->react_time == 0)
                          {
                                  if (runParticle->calcite > 3.0)
                                          runParticle->react_time = step;
                          }

                          if (runParticle->calcite > 0.5)
                          {
                                  runParticle->calcite = 0.5;
                                  cout << "large number dolomite" << endl;
                          }

                          cout << "reaction" << endl;
                           //runParticle->rad_par_fluid =
     runParticle->rad_par_fluid - runParticle->rad_par_fluid *
     ((runParticle->calcite-ratio)); runParticle->rad_par_fluid =
     runParticle->rad_par_fluid - runParticle->rad_par_fluid * (0.05);
               runParticle->area_par_fluid = 3.14 *
     pow(runParticle->rad_par_fluid,2.0); runParticle->calcite = 0.5;


               box_x = runParticle->xpos*(particlex/2) + 1;
               box_y = runParticle->ypos*(particlex/2) + 1;
               change = (runParticle->calcite - ratio) / 5.0;
               //P->ChangeConcentration(box_x, box_y, change);

               //Fluid->ChangeConcentration(box_x, box_y, change);
            }

     // }while (react > 0.1);
     */
}

void Lattice::Loose_Springs(bool boundary) {
  int i, j;

  for (i = 0; i < numParticles; i++) {
    if (!runParticle->is_lattice_boundary) {
      for (j = 0; j < 8; j++) {
        if (runParticle->neigP[j]) {
          if (runParticle->is_lattice_boundary) {
            if (!runParticle->neigP[j]->is_lattice_boundary) {
              runParticle->neigP[j] = NULL;
              runParticle->springf[j] = 0.0;
              runParticle->springs[j] = 0.0;
              runParticle->break_Str[j] = 0.0;
              runParticle->neigpos[j][0] = 0.0; // change
              runParticle->neigpos[j][1] = 0.0;
            }
          } else {
            runParticle->neigP[j] = NULL;
            runParticle->springf[j] = 0.0;
            runParticle->springs[j] = 0.0;
            runParticle->break_Str[j] = 0.0;
            runParticle->neigpos[j][0] = 0.0; // change
            runParticle->neigpos[j][1] = 0.0;
          }
        }
      }
    }
    runParticle = runParticle->nextP;
  }
}

void Lattice::Shift_Particles(float variation) {
  int i, j;
  float ran_nb;

  srand(time(0));

  variation = variation / particlex;

  for (i = 0; i < numParticles; i++) {
    if (!runParticle->is_lattice_boundary) {
      ran_nb = rand() / (float)RAND_MAX;
      ran_nb = ran_nb - 0.5;
      ran_nb = ran_nb * 2.0;
      ran_nb = ran_nb * variation;

      runParticle->xpos = runParticle->xpos + ran_nb;

      ran_nb = rand() / (float)RAND_MAX;
      ran_nb = ran_nb - 0.5;
      ran_nb = ran_nb * 2.0;
      ran_nb = ran_nb * variation;

      runParticle->ypos = runParticle->ypos + ran_nb;

      runParticle->ChangeBox(repBox, particlex);
    }
    runParticle = runParticle->nextP;
  }
}

void Lattice::Adjust_Springs_Particles() {
  int i, j, springs;
  double dir_rad, newrad, small_rad;

  for (i = 0; i < numParticles; i++) {
    springs = 0;
    newrad = 0.0;
    small_rad = 1.0;

    for (j = 0; j < 8; j++) {
      if (runParticle->neigP[j]) {
        runParticle->neigpos[j][0] =
            runParticle->neigP[j]->xpos - runParticle->xpos; // change
        runParticle->neigpos[j][1] =
            runParticle->neigP[j]->ypos - runParticle->ypos; // change

        dir_rad = sqrt(runParticle->neigpos[j][0] * runParticle->neigpos[j][0] +
                       runParticle->neigpos[j][1] * runParticle->neigpos[j][1]);
        newrad = newrad + dir_rad;
        if (small_rad > dir_rad)
          small_rad = dir_rad;

        springs++;
      }
    }

    newrad = newrad / springs;
    cout << small_rad << endl;

    runParticle->radius = small_rad / 2.0;
    runParticle->average_radius = newrad / 2.0;
    runParticle->healing = springs;

    runParticle = runParticle->nextP;
  }
}

void Lattice::Connect_New_Lattice(bool changeRad) {
  int i, j, ii, jj, k, box, pos, springs, dist_count;
  double dx, dy, dd, newrad, small_rad, dir_rad;
  Particle *neig;
  Particle *neig_list[40];
  double distance[40], dist;

  for (k = 0; k < 40; k++) {
    neig_list[k] = 0;
    distance[k] = 0;
  }

  for (i = 0; i < numParticles; i++) {
    if (!runParticle->is_lattice_boundary) {
      for (ii = 0; ii < 5; ii++) {
        for (jj = 0; jj < 3; jj++) {

          //------------------------------------------
          //  case  upper  or  lower  row (+-size)
          //------------------------------------------

          if (jj == 0)
            box = 2 * particlex * (-1);
          if (jj == 1)
            box = 0;
          if (jj == 2)
            box = 2 * particlex;
          // if (jj == 3)
          // box = 4*particlex;
          //-----------------------------------------
          //  now  check  this  position
          //  all  particles  in  this  position !!
          //-----------------------------------------

          pos = runParticle->box_pos + (ii - 2) + box;

          //-----------------------------------------
          // repbox should always be positive
          //-----------------------------------------

          if (pos > 0) {
            if (repBox[pos]) {
              neig = repBox[pos]; // helppointer  for  while  loop

              while (neig) // as long as there is somebody
              {
                if (neig->nb != runParticle->nb) {
                  if (runParticle->is_lattice_boundary &&
                      neig->is_lattice_boundary) {
                    cout << " boundary spring" << endl;
                  } else // now connect again Need to fill in list with
                         // particles and distances and fill the neighbours
                         // then.

                  // also check that new equilibrium distances are given to
                  // springs from the boundaries!

                  {
                    dx = neig->xpos - runParticle->xpos;

                    dy = neig->ypos - runParticle->ypos;

                    dd = sqrt((dx * dx) + (dy * dy));

                    // here fill vector with particles and distance and use
                    // counter for it.

                    // the rest goes out of this loop into a new loop to fill
                    // from smallest distance particle.
                    // for (k = 0; k<12; k++)
                    k = 0;

                    while (neig_list[k]) {
                      k++;
                    }
                    // if (!neig_list[k])

                    cout << k << endl;

                    neig_list[k] = neig;
                    distance[k] = dd;
                  }
                }

                if (neig->next_inBox) // if there is anotherone in box
                {
                  neig = neig->next_inBox; // shift the help pointer
                } else {
                  break; // if not go out to next box
                }
              }
            }
          }
        }
      }

      for (j = 0; j < 6; j++) {
        if (runParticle->neigP[j]) {
          ; // particle already present, this should be a lattice boundary
            // spring
        } else {
          // now connect with the closest particles

          cout << "in" << endl;

          dist = 1.0;
          neig = 0;

          for (k = 0; k < 40; k++) {
            if (neig_list[k]) {
              if (distance[k] < dist) {
                neig = neig_list[k];
                dist = distance[k];
                dist_count = k;
              }
            }
          }

          // now connect with this guy

          if (dist < runParticle->radius * 3.5) {
            if (neig) {
              runParticle->neigP[j] = neig;
              runParticle->springf[j] =
                  sqrt(3.0) * runParticle->E /
                  (3.0 *
                   (1.0 - runParticle->poisson_ratio)); // set spring constant
              runParticle->springs[j] =
                  (1.0 - 3.0 * runParticle->poisson_ratio) /
                  (1.0 + runParticle->poisson_ratio) *
                  runParticle->springf[j]; // set spring to 1.0

              runParticle->springv[j] = 1.0; // viscosity

              runParticle->break_Str[j] = 0.0017;

              // runParticle->radius = runParticle->radius;

              runParticle->rad[j] = runParticle->radius; // equilibrium length

              runParticle->spring_boundary[j] = false; //??
              runParticle->neigpos[j][0] =
                  neig->xpos - runParticle->xpos; // change
              runParticle->neigpos[j][1] =
                  neig->ypos - runParticle->ypos; // change

              cout << j << endl;

              neig_list[dist_count] = 0;
            }
          }
        }
      }

      for (k = 0; k < 40; k++) {
        neig_list[k] = 0;
        distance[k] = 0;
      }
    }
    runParticle = runParticle->nextP;
  }
  for (i = 0; i < numParticles; i++) // loop
  {

    springs = 0;
    newrad = 0.0;
    small_rad = 1.0;

    for (j = 0; j < 8; j++) // run through neighbours
    {
      if (runParticle->neigP[j]) // if neighbour
      {
        if (runParticle->is_lattice_boundary &&
            !runParticle->neigP[j]->is_lattice_boundary) {
          runParticle->neigpos[j][0] =
              runParticle->neigP[j]->xpos - runParticle->xpos;
          runParticle->neigpos[j][1] =
              runParticle->neigP[j]->ypos - runParticle->ypos;
        }

        springs++;
        dir_rad = sqrt(runParticle->neigpos[j][0] * runParticle->neigpos[j][0] +
                       runParticle->neigpos[j][1] * runParticle->neigpos[j][1]);
        newrad = newrad + dir_rad;
        if (small_rad > dir_rad)
          small_rad = dir_rad;

        for (jj = 0; jj < 8; jj++) // run through neighbour neighbours
                                   // springs
        {
          if (runParticle->neigP[j]->neigP[jj]) // if
                                                // neighbour
          {
            if (runParticle->nb ==
                runParticle->neigP[j]->neigP[jj]->nb) // look
                                                      // backwards
                                                      // and
                                                      // find
                                                      // particle
                                                      // again
            {
              runParticle->neig_spring[j] = jj; // thats
              // the one
            }
          }
        }
      }
    }
    cout << springs << endl;
    newrad = newrad / springs;
    cout << small_rad << endl;
    if (changeRad) {
      runParticle->radius = small_rad / 2.0;
      runParticle->average_radius = newrad / 2.0;
      runParticle->healing = springs;
    }

    runParticle = runParticle->nextP;
  }
  cout << "out build lattice routine" << endl;
}

//****************************************************************************************************
// New fluid code written after Irfan's code in fluid lattice (see Ghani et al.
// 2013)
//
// first the fluid lattice is initialized, scale and background pressure is
// applied,
//
//
// Daniel 2018
//****************************************************************************************************

void Lattice::Initialize_Fluid_Lattice(double pressure, double scale,
                                       double conc, double time_factor) {

  Fluid = new Fluid_Lattice(particlex);
  srand(time(0));

  int i;

  runParticle = &refParticle;
  for (i = 0; i < numParticles; i++) {
    runParticle->velx = 0;
    runParticle->vely = 0;

    runParticle->oldx = runParticle->xpos;
    runParticle->oldy = runParticle->ypos;

    runParticle->rad_par_fluid = runParticle->radius;

    runParticle->real_radius = runParticle->radius * scale;

    runParticle->area_par_fluid =
        Pi * pow(runParticle->rad_par_fluid, 2.0); // area of particle radius

    runParticle->real_density = 2700;

    runParticle->real_young =
        runParticle->young * pressure_scale * pascal_scale;

    runParticle = runParticle->nextP;
  }

  Fluid->Background_Pressure(pressure, conc);

  Fluid->Set_Scale(scale, time_factor);
}

//----------------------------------------------------------------------------------------
// this is called before the pressure is calculated, reads porosity and movement
// matrix and calculated the permeability from the porosity. input is boundary
// value for applied boundary conditions and the cozeny grain size - will affect
// permeability
//
// danie 2018
//----------------------------------------------------------------------------------------

void Lattice::Fluid_Parameters(double fracture_Perm, int boundary,
                               double grain_size) {
  Fluid->Read_Porosity_Movement_Matrix(fracture_Perm, repBox, boundary);
  Fluid->Calculate_Permeability_Matrix(grain_size);
}

//----------------------------------------------------------------------------------------
// simply adds fluid at a node.
//
// danie 2018
//----------------------------------------------------------------------------------------

void Lattice::Input_Fluid_Lattice(double pressure, int x, int y) {

  Fluid->Fluid_Input(pressure, x, y);
}

//----------------------------------------------------------------------------------------
// fluid pressure calculation using ADI method, implicit solution. first set
// boundaries then copy the matrix for delta P in dt. input upper and lower
// pressure, internal interations (this could actually be internally calculated
// depending on stability. then input for boundaries and whether or not
// concentration is averaged
//
// danie 2018
//----------------------------------------------------------------------------------------

void Lattice::Calculate_Fluid_Pressure(double upper_pressure,
                                       double lower_pressure, int iteration,
                                       int sides, int average,
                                       double additional) {
  {
    int k, iter;

    // copy the old fluid pressure in an extra matrix to get the grad in t later

    // Fluid->Set_Boundaries(upper_pressure, lower_pressure, sides, 0.5,
    // additional);

    Fluid->Set_Boundaries(upper_pressure, lower_pressure, sides, 0.5,
                          additional);

    Fluid->copy_matrix();

    /************************
     * Scheme - test
     *
     * this setupt is new with corrected functions of source term
     * and pressure calculation an assignment to particles
     *
     * *********************/

    // loop to potentially make the calculation more stable

    iter = Fluid->alpha_set(1);

    cout << "iterations fluid pressure:" << iter << endl;

    for (k = 0; k < iter; k++) {

      // set the a matrix

      Fluid->alpha_set(iter);

      Fluid->matrices_set_1();
      Fluid->multiplication();
      Fluid->test_pressure_with_source(0); // adding the source term
      Fluid->transpose_trans_P();
      Fluid->invert();
      Fluid->multiplication_inv();

      Fluid->matrices_set_2();
      Fluid->multiplication();
      Fluid->transpose_trans_P();
      Fluid->test_pressure_with_source(0); // adding the source term
      Fluid->invert();
      Fluid->multiplication_inv();

      // function at the moment just simply passes the upper boundary fluid
      // pressure to the matrix, at the moment this is fixed.

      Fluid->Set_Boundaries(upper_pressure, lower_pressure, sides, 0.5,
                            additional);
    }

    // Calculate the Darcy velocity of this fluid flow step

    Fluid->Get_Fluid_Velocity();

    // function that passes back the fluid pressure values to the actual
    // particles and calculates the fluid pressure gradient and force in x and y
    // that acts on the particles. Also passes back the Darcy Velocity to the
    // particles for visualization

    Fluid->Pass_Back_Gradients(repBox, average);
  }
}

//----------------------------------------------------------------------------------------
// fluid pressure calculation using ADI method, implicit solution. first set
// boundaries then copy the matrix for delta P in dt. input upper and lower
// pressure, internal interations (this could actually be internally calculated
// depending on stability. then input for boundaries and whether or not
// concentration is averaged
//
// danie 2018
//----------------------------------------------------------------------------------------

double Lattice::Calculate_Fluid_Pressure_timeadjust(double upper_pressure,
                                                    double lower_pressure,
                                                    int iteration, int sides,
                                                    int average,
                                                    double additional) {

  double nonlinear_time;

  // copy the old fluid pressure in an extra matrix to get the grad in t later

  // Fluid->Set_Boundaries(upper_pressure, lower_pressure, sides, 0.5,
  // additional);

  Fluid->Set_Boundaries(upper_pressure, lower_pressure, sides, 0.5, additional);

  Fluid->copy_matrix();

  /************************
   * Scheme - test
   *
   * this setupt is new with corrected functions of source term
   * and pressure calculation an assignment to particles
   *
   * *********************/

  // loop to potentially make the calculation more stable

  nonlinear_time = Fluid->alpha_set_time();

  cout << "time" << nonlinear_time << endl;

  // set the a matrix

  Fluid->alpha_set_fluid(nonlinear_time);

  Fluid->matrices_set_1();
  Fluid->multiplication();
  Fluid->test_pressure_with_source(0); // adding the source term
  Fluid->transpose_trans_P();
  Fluid->invert();
  Fluid->multiplication_inv();

  Fluid->matrices_set_2();
  Fluid->multiplication();
  Fluid->transpose_trans_P();
  Fluid->test_pressure_with_source(0); // adding the source term
  Fluid->invert();
  Fluid->multiplication_inv();

  // function at the moment just simply passes the upper boundary fluid pressure
  // to the matrix, at the moment this is fixed.

  Fluid->Set_Boundaries(upper_pressure, lower_pressure, sides, 0.5, additional);

  // Calculate the Darcy velocity of this fluid flow step

  Fluid->Get_Fluid_Velocity();

  // function that passes back the fluid pressure values to the actual particles
  // and calculates the fluid pressure gradient and force in x and y that acts
  // on the particles. Also passes back the Darcy Velocity to the particles for
  // visualization

  Fluid->Pass_Back_Gradients(repBox, average);

  return nonlinear_time;
}

//----------------------------------------------------------------------------------------
// sets a hydrostatic gradient depending on the depth. Careful with this one,
// has to be scaled right otherwise the model will blow up.
//
// danie 2018
//----------------------------------------------------------------------------------------

void Lattice::Hydrostatic_Fluid(int depth) {
  Fluid->Set_hydrostatic_gradient(depth);
}

//----------------------------------------------------------------------------------------
// add random fluid pressure nodes within a layer.
//
// danie 2018
//----------------------------------------------------------------------------------------

void Lattice::Fluid_Insert_Random_Node(double ymin, double ymax,
                                       double press_increment) {
  float ran_pos_x, ran_pos_y, ran_pressure; // rand number
  double press_increment_final;
  int node_x, node_y;

  ran_pos_x = rand() / (float)RAND_MAX;
  ran_pos_y = rand() / (float)RAND_MAX;
  ran_pressure = rand() / (float)RAND_MAX;

  node_y = int((((ran_pos_y * (ymax - ymin)) + ymin) * ((particlex / 2))));
  node_x = 0 + int((ran_pos_x * ((particlex / 2))));

  press_increment_final = press_increment + (press_increment * ran_pressure);

  if (node_y == 0) {
    if (node_x == 0)
      press_increment_final = 0;
    else if (node_x == ((particlex / 2) - 1))
      press_increment_final = 0;
  }

  Fluid->Fluid_Input(press_increment_final, node_x, node_y);
}

//----------------------------------------------------------------------------------------
// add random fluid pressure nodes within a zone
//
// danie 2018
//----------------------------------------------------------------------------------------

void Lattice::Fluid_Insert_Random_Nodexy(double ymin, double ymax, double xmin,
                                         double xmax, double press_increment) {
  float ran_pos_x, ran_pos_y, ran_pressure; // rand number
  double press_increment_final;
  int node_x, node_y;

  ran_pos_x = rand() / (float)RAND_MAX;
  ran_pos_y = rand() / (float)RAND_MAX;
  ran_pressure = rand() / (float)RAND_MAX;

  node_y = int((((ran_pos_y * (ymax - ymin)) + ymin) * ((particlex / 2))));
  node_x = int((((ran_pos_x * (xmax - xmin)) + xmin) * ((particlex / 2))));

  press_increment_final = press_increment + (press_increment * ran_pressure);

  if (node_y == 0) {
    if (node_x == 0)
      press_increment_final = 0;
    else if (node_x == ((particlex / 2) - 1))
      press_increment_final = 0;
  }

  Fluid->Fluid_Input(press_increment_final, node_x, node_y);
}

//----------------------------------------------------------------------------------------
// advection not used at the moment
//
// danie 2018
//----------------------------------------------------------------------------------------

void Lattice::Calculate_Con(double vel, double dif) {

  Fluid->SolveadvDisExplicitTransient(vel, dif);
}

//----------------------------------------------------------------------------------------
// solve for advection of concentration using explicit method
//
// uses darcy velocity from the fluid pressure diffusion
//----------------------------------------------------------------------------------------

void Lattice::Calculate_Advection(double time, double space, int boundary,
                                  int method) {
  if (method == 1)
    Fluid->SolveAdvectionLax_Wenddroff(time, space / particlex, boundary);
  else if (method == 2)
    Fluid->SolveAdvectionFromm(space / particlex, time);
}

//----------------------------------------------------------------------------------------
// solve for advection of concentration using explicit method
//
// uses darcy velocity from the fluid pressure diffusion
//----------------------------------------------------------------------------------------

void Lattice::Calculate_Advection_Metal(double time, double space,
                                        int boundary) {
  Fluid->SolveAdvectionLax_Wenddroff_Metal(time, space / particlex, boundary);
}

//----------------------------------------------------------------------------------------
// solve for advection of concentration using explicit method
//
// uses darcy velocity from the fluid pressure diffusion
//----------------------------------------------------------------------------------------

void Lattice::Calculate_Advection_Temperature(double time, double space,
                                              int boundary) {
  Fluid->SolveAdvectionLax_Wenddroff_Temperature(time, space / particlex,
                                                 boundary);
}

//----------------------------------------------------------------------------------------
// solve for Metal concentration diffusion in fluid lattice. momentarily
// constant dif constant
//
// uses ADI method, implicit alternation direction
//----------------------------------------------------------------------------------------

void Lattice::Calculate_Diffusion_Metal(double constant, double time,
                                        double space) {
  Fluid->SolveMetalDiffusionImplicit(constant, time, space);
}

//----------------------------------------------------------------------------------------
// solve for Temperature diffusion in fluid lattice. momentarily constant dif
// constant
//
// uses ADI method, implicit alternation direction
//----------------------------------------------------------------------------------------

void Lattice::Calculate_Diffusion_Temperature(double constant, double time,
                                              double space) {
  Fluid->SolveTemperatureDiffusionImplicit(constant, time, space);
}

//----------------------------------------------------------------------------------------
// solve for concentration diffusion in fluid lattice. momentarily constant dif
// constant
//
// uses ADI method, implicit alternation direction
//----------------------------------------------------------------------------------------

void Lattice::Calculate_Diffusion(double constant, double time, double space) {
  Fluid->SolveDiffusionImplicit(constant, time, space);
}

// FOR SIMON !!! 2021

//----------------------------------------------------------------------------------------
// These are the new functions in lattice for the diffusion and exchange
// reaction
//
// They are all passed directly to fluid lattice at the moment
//
// This is done calling the object Fluid and point to the functions in
// fluid_lattice.cc with -> and passing the parameter further down to the fluid
// lattice function
//----------------------------------------------------------------------------------------

// first function for the diffusion of solid matter in the grains solved with
// the ADI method mixed Implicit/explicit.

void Lattice::Calculate_Diffusion_Solid(double time, double space) {
  Fluid->SolveDiffusionImplicitSolid(
      time, space); // call implicit diffusion solver in fluid_lattice for solid
  Fluid->Pass_Back_Concentration(
      repBox); // passes back new values to the particles, gets the repulsion
               // box from lattice
}

// This function is now not used ! !
// second function for the diffusion of fluid in the grain boundaries solved
// with the ADI method mixed Implicit/explicit. lost mass so not used now... may
// be reinstated later

void Lattice::Calculate_Diffusion_Fluid(double time, double space) {
  Fluid->SolveDiffusionImplicitFluid(
      time,
      space); // call implicit diffusion solver in fluid_lattice for boundaries
  Fluid->Pass_Back_Concentration(
      repBox); // passes back new values to the particles, gets the repulsion
               // box from lattice
}

// this is the function that calculates diffusion only in the grain boundaries
// it is a purely explicit method
// the method is adjusted to whatever the neighbours are and switches from one
// to 1.5 and two dimension depending on the local geometries.

void Lattice::Calculate_Diffusion_Fluid_Explicit(double time, double space) {
  Fluid->SolveDiffusionExplicitFluid(
      time, space); // call the explicit solver in fluid_lattice.cc
  Fluid->Pass_Back_Concentration(
      repBox); // passes back new values to the particles, gets the repulsion
               // box from lattice
}

// this is the Exchange function that we need to work on quite a bit.

void Lattice::Exchange_Solid_Fluid(double constant, double time, double space) {
  Fluid->Exchange(time, space,
                  repBox); // call exchange function in fluid_lattice.cc
  Fluid->Pass_Back_Concentration(
      repBox); // passes back new values to the particles, gets the repulsion
               // box from lattice
}

//-------------------------------------------------------------------------------------------------------------------------
// 	Here some initialization functions for the diffusion and exchange
//-------------------------------------------------------------------------------------------------------------------------

void Lattice::Initiate_Mineral(double background) // some background values
{
  Fluid->Init_Mineral(background); // call function in fluid_lattice.cc
}

// this function now sets the minerals, at the moment garnet and Biotite

void Lattice::Set_Mineral(int Grain,
                          int Mineral) // 1 for mineral is garnet, 2 is Biotite
{
  Fluid->Set_Mineral_Fluid(Grain, Mineral,
                           repBox); // call the function in fluid_lattice.cc
  Fluid->Pass_Back_Concentration(
      repBox); // passes back new values to the particles, gets the repulsion
               // box from lattice
}

// this function sets the parameters for the fluid diffusion
// if we want different diffusion coefficients for Mg and Fe (which we should)
// then this needs changing

void Lattice::Set_Fluid(
    double back,
    double boundary) // the second number is important for the diffusion now
{
  Fluid->Set_Boundary_Fluid(back, boundary,
                            repBox); // call things in fluid_lattice.cc
}

//---------------------------------------------------------------------------------------------------------
// end functions from diffusion and exchange Daniel and Simon 2021
//
//----------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------
// function that makes the grain boundaries more permeable for fluid flow
//
// Daniel 2018
//-----------------------------------------------------------------------------------

void Lattice::Permeable_Grain_Boundary(double factor) {
  int j;

  for (j = 0; j < numParticles; j++) {
    if (runParticle->is_boundary) {
      runParticle->rad_par_fluid =
          runParticle->rad_par_fluid - runParticle->rad_par_fluid * (factor);
      runParticle->area_par_fluid = 3.14 * pow(runParticle->rad_par_fluid, 2.0);
    }
    runParticle = runParticle->nextP;
  }
}

//----------------------------------------------------------------------------------
// function that makes the grain  more permeable for fluid flow
//
// Daniel 2018
//-----------------------------------------------------------------------------------

void Lattice::Permeable_Grain(double factor, int grain) {
  int j;

  for (j = 0; j < numParticles; j++) {
    if (runParticle->grain == grain) {
      runParticle->rad_par_fluid =
          runParticle->rad_par_fluid - runParticle->rad_par_fluid * (factor);
      runParticle->area_par_fluid = 3.14 * pow(runParticle->rad_par_fluid, 2.0);
    }
    runParticle = runParticle->nextP;
  }
}

void Lattice::Add_Concentration(int nodex, int nodey, double concentration) {
  Fluid->Set_Concentration(nodex, nodey, concentration);
}

void Lattice::Initialize_Heat(double xpos, double xpos2, double temp,
                              double critAge) {
  Fluid->SetHeatSheet(xpos, xpos2);

  Fluid->Pass_Back_SheetT(repBox, temp, critAge);
}

void Lattice::InitHeat(int T, double critAge) {
  int i, j, count;
  double breaks, springs, springf;

  count = 0;

  for (i = 0; i < numParticles; i++) {
    for (j = 0; j < 8; j++) {
      if (runParticle->neigP[j]) {
        breaks += runParticle->break_Str[j];
        springs += runParticle->springs[j];
        springf += runParticle->springf[j];

        count++;
      }
    }

    runParticle = runParticle->prevP;
  }
  average_break = breaks / count;
  average_springs = springs / count;
  average_springf = springf / count;

  Fluid->InitHeatSheet(T);

  Fluid->Pass_Back_SheetT(repBox, T, critAge);
}

void Lattice::Initialize_Heat_Box(double xpos, double xpos2, double ypos,
                                  double ypos2, int T, double critAge) {
  Fluid->SetHeatSheetbox(xpos, xpos2, ypos, ypos2, T);

  Fluid->Pass_Back_SheetT(repBox, 1200, critAge);
}

void Lattice::solve_Heat(double space, double time, double constant,
                         double temp, double critAge, int tfac) {
  int i;

  Fluid->Read_TSheet(repBox, 1);

  for (i == 0; i < tfac; i++)
    Fluid->SolveHeatSheet(space, time, constant);

  Fluid->Pass_Back_SheetT(repBox, temp, critAge);
}
void Lattice::Adjust_break() {
  int i, j;
  for (i = 0; i < numParticles; i++) {
    for (j = 0; j < 8; j++) {
      if (runParticle->neigP[j])
        runParticle->break_Str[j] =
            runParticle->break_Str[j] / (runParticle->SheetT / 80);
    }
    runParticle = runParticle->prevP;
  }
}

void Lattice::Loose_Unodes(double xmin, double xmax) {
  int i, j, ii;

  Ulistcount = 0;

  for (i = 0; i < numParticles; i++) {
    if (runParticle->xpos < xmin) {
      Ulistcount++;
      Ulist[Ulistcount] = runParticle->nb;
      ElleRemoveUnodeFromFlynn(runParticle->p_Unode->flynn(),
                               runParticle->p_Unode->id());
      runParticle->p_Unode = NULL;
      for (ii = 0; ii < 32; ii++) {
        runParticle->elle_Node[ii] = -1;
      }

    } else if (runParticle->xpos > xmax) {
      Ulistcount++;
      Ulist[Ulistcount] = runParticle->nb;
      ElleRemoveUnodeFromFlynn(runParticle->p_Unode->flynn(),
                               runParticle->p_Unode->id());
      runParticle->p_Unode = NULL;
      for (ii = 0; ii < 32; ii++) {
        runParticle->elle_Node[ii] = -1;
      }
    }

    runParticle = runParticle->prevP;
  }
}

/*
void Lattice::Grow_Particle(double posx, double posy, int boxpos, double
sheetYoung, double sheetVisc, int age)
{
        int i;
        Coords xy;

        runParticle = new Particle;
        cout << "new particle" << posx << posy << endl;
        runParticle->nextP = firstParticle;
        runParticle->prevP = firstParticle->prevP;
        firstParticle->prevP->nextP = runParticle;
        firstParticle->prevP = runParticle;
        numParticles ++;
        runParticle->nb = numParticles-1;

        runParticle->p_Unode = ElleGetParticleUnode (Ulist[Ulistcount]);
        Ulistcount = Ulistcount - 1;


        for (i = 0; i < 8; i++)
        {
                runParticle->neigP[j] = NULL;
                runParticle->springf[j] = 0.0;
                runParticle->springs[j] = 0.0;
                runParticle->break_Str[j] = 0.0;
                runParticle->neigpos[j][0] = 0.0;	// change
                runParticle->neigpos[j][1] = 0.0;
        }

        runParticle->radius = runParticle->smallest_rad = runParticle->radius *
100.0 / particlex; runParticle->area = runParticle->radius * runParticle->radius
* 3.1415927;

        runParticle->xpos = posx + runParticle->radius;
    runParticle->ypos = posy + runParticle->radius;
    runParticle->oldx = posx + runParticle->radius;
    runParticle->oldy = posy + runParticle->radius;
    runParticle->initialy = posy + runParticle->radius;

    xy.x = runParticle->xpos;	// pass x
    xy.y = runParticle->ypos;	// pass y

    runParticle->p_Unode->setPosition(&xy);	// and pass that to the
            // Unode

        runParticle->SetSprings (default_Young, default_Poisson);
        runParticle->springconst = 1.0;
        runParticle->young = 1.0;
        runParticle->next_inBox = NULL;
        runParticle->box_pos = boxpos;
        runParticle->sheet_x = posx;
        runParticle->sheet_y = posy;
        runParticle->sheet_young = sheetYoung;
        runParticle->sheet_visc = sheetVisc;
        runParticle->SheetT = 1000;
        runParticle->F_P_x = 0;
        runParticle->F_P_y = 0;
        repBox[boxpos] = runParticle;
        runParticle->react_time = age;
}
*/

void Lattice::Grow_Particle(double posx, double posy, int boxpos,
                            double sheetYoung, double sheetVisc, int age,
                            int T) {
  int i, connected_neighbors, ii, jj, box, pos, k, myPlace, neigPlace, l, kk,
      ll;
  Coords xy;
  double dx, dy, dd, average_springs, average_break_Str, average_springf,
      springfactor, springstrength;
  Particle *neig;

  springfactor = 0.2;   // 0.5
  springstrength = 0.2; // 0.1

  runParticle = new Particle;
  cout << "new particle" << posx << posy << endl;
  runParticle->nextP = firstParticle->prevP;
  runParticle->prevP = firstParticle->prevP->prevP;
  firstParticle->prevP->prevP->nextP = runParticle;
  firstParticle->prevP->prevP = runParticle;
  numParticles++;
  runParticle->nb = numParticles - 1;

  runParticle->p_Unode = ElleGetParticleUnode(Ulist[Ulistcount]);
  Ulistcount = Ulistcount - 1;

  runParticle->age = 0;

  for (i = 0; i < 8; i++) {
    runParticle->neigP[j] = NULL;
    runParticle->springf[j] = 0.0;
    runParticle->springs[j] = 0.0;
    runParticle->springv[j] = 0.0;
    runParticle->break_Str[j] = 0.0;
    runParticle->neigpos[j][0] = 0.0; // change
    runParticle->neigpos[j][1] = 0.0;
  }

  runParticle->radius = runParticle->smallest_rad =
      runParticle->radius * 100.0 / particlex;
  runParticle->area = runParticle->radius * runParticle->radius * 3.1415927;

  runParticle->xpos = posx + runParticle->radius;
  runParticle->ypos = posy + runParticle->radius;
  runParticle->oldx = posx + runParticle->radius;
  runParticle->oldy = posy + runParticle->radius;
  runParticle->initialy = posy + runParticle->radius;

  xy.x = runParticle->xpos; // pass x
  xy.y = runParticle->ypos; // pass y

  runParticle->p_Unode->setPosition(&xy); // and pass that to the
                                          // Unode

  runParticle->SetSprings(default_Young, default_Poisson);
  runParticle->springconst = 0.02; // 0.01
  runParticle->young = 0.01;
  runParticle->next_inBox = NULL;
  runParticle->box_pos = boxpos;
  runParticle->sheet_x = posx;
  runParticle->sheet_y = posy;
  runParticle->sheet_young = sheetYoung;
  runParticle->sheet_visc = sheetVisc;
  runParticle->SheetT = T;
  runParticle->Tdif = T;
  runParticle->F_P_x = 0;
  runParticle->F_P_y = 0;
  repBox[boxpos] = runParticle;
  runParticle->particle_age = age;

  driftlist[int((runParticle->ypos) * particlex)] = runParticle;

  if (runParticle->ypos < runParticle->radius * 2)
    runParticle->fix_y = true;

  if (runParticle->ypos > (1.0 - runParticle->radius * 2))
    runParticle->fix_y = true;

  connected_neighbors = 0;
  dx = dy = 0;

  for (ii = 0; ii < 3; ii++) {
    for (jj = 0; jj < 3; jj++) {

      //------------------------------------------
      //  case  upper  or  lower  row (+-size)
      //------------------------------------------

      if (jj == 0)
        box = 2 * particlex * (-1);
      if (jj == 1)
        box = 0;
      if (jj == 2)
        box = 2 * particlex;
      //-----------------------------------------
      //  now  check  this  position
      //  all  particles  in  this  position !!
      //-----------------------------------------

      pos = runParticle->box_pos + (ii - 1) + box;

      //-----------------------------------------
      // repbox should always be positive
      //-----------------------------------------

      if (pos > 0) {
        if (pos == runParticle->box_pos) {
          ;
        } else if (repBox[pos]) {

          neig = repBox[pos]; // helppointer  for  while  loop

          while (neig) // as long as there is somebody
          {
            dx += neig->xpos;
            dy += neig->ypos;
            connected_neighbors++;

            if (neig->next_inBox) // if there is anotherone in box
            {
              neig = neig->next_inBox; // shift the help pointer
            } else {
              break; // if not go out to next box
            }
          }
        }
      }
    }
  }

  if (connected_neighbors > 0) {
    dx = dx / connected_neighbors;
    dy = dy / connected_neighbors;
    cout << "new position" << dx << "/" << dy << endl;
    runParticle->xpos = dx;
    runParticle->ypos = dy;
    runParticle->oldx = dx;
    runParticle->oldy = dy;

    xy.x = runParticle->xpos; // pass x
    xy.y = runParticle->ypos; // pass y

    runParticle->p_Unode->setPosition(&xy); // and pass that to the

    runParticle->sheet_x = dx;
    runParticle->sheet_y = dy;

  } else
    cout << "no Neighbour??" << endl;

  // Relaxation();

  //---------------------------------------------------
  // nine positions in 2d box 3x3
  //---------------------------------------------------
  connected_neighbors = 0;

  for (ii = 0; ii < 3; ii++) {
    for (jj = 0; jj < 3; jj++) {

      //------------------------------------------
      //  case  upper  or  lower  row (+-size)
      //------------------------------------------

      if (jj == 0)
        box = 2 * particlex * (-1);
      if (jj == 1)
        box = 0;
      if (jj == 2)
        box = 2 * particlex;
      //-----------------------------------------
      //  now  check  this  position
      //  all  particles  in  this  position !!
      //-----------------------------------------

      pos = runParticle->box_pos + (ii - 1) + box;

      //-----------------------------------------
      // repbox should always be positive
      //-----------------------------------------

      if (pos > 0) {
        if (pos == runParticle->box_pos) {
          ;
        } else if (repBox[pos]) {

          neig = repBox[pos]; // helppointer  for  while  loop

          while (neig) // as long as there is somebody
          {

            // if (!neig->done)  // if not done already
            {
              dx = neig->xpos - runParticle->xpos;
              dy = neig->ypos - runParticle->ypos;

              dd = sqrt((dx * dx) + (dy * dy));

              // cout << dd << endl;

              // if ((dd/heal_dist) < runParticle->radius*2.0)
              if (dd < runParticle->radius * 1.7) {
                if (dd > runParticle->radius * 2.5) {
                  // ran_nb = rand () / (float) RAND_MAX;

                  // prob=(runParticle->radius*2.0/dd)*prob;

                  // cout << ran_nb << prob << endl;

                  // if (ran_nb < prob) // now we are actually healing
                  {
                    if (connected_neighbors < 8) {
                      connected_neighbors++;

                      // get averages first

                      ll = 0;
                      average_springf = 0;
                      average_springs = 0;
                      average_break_Str = 0;

                      // for (kk=0;kk<8;kk++)
                      //{
                      ////if (runParticle->neigP[kk])
                      ////{
                      ////ll ++;
                      ////average_springf =+ runParticle->springf[kk];
                      ////average_springs =+ runParticle->springs[kk];
                      ////average_break_Str =+ runParticle->break_Str[kk];
                      ////}
                      // if (neig->neigP[kk])
                      //{
                      // ll ++;
                      // average_springf =+ neig->springf[kk];
                      // average_springs =+ neig->springs[kk];
                      // average_break_Str =+ neig->break_Str[kk];
                      //}
                      //}
                      // if (ll > 0)
                      //{
                      // average_springf = average_springf / ll;
                      // average_springs = average_springs / ll;
                      // average_break_Str = average_break_Str / ll;
                      //}
                      // else
                      //{
                      // cout << "empty" << endl;
                      // average_springf = 0.05; //from 0.01 to 0.02, 0.02 is
                      // good, average_springs = 0.05; //from 0.01 to 0.02, 0.02
                      // is good, average_break_Str = 0.005;
                      //}

                      //// find empty positions

                      for (k = 0; k < 8; k++) {
                        if (!runParticle->neigP[k]) {
                          myPlace = k;
                        }
                      }

                      for (l = 0; l < 8; l++) {
                        if (!neig->neigP[l]) {
                          neigPlace = l;
                        }
                      }

                      // cout << average_springf << endl;
                      // cout << " Young " << runParticle->young << endl;

                      runParticle->neig_spring[myPlace] = neigPlace;
                      neig->neig_spring[neigPlace] = myPlace;

                      runParticle->neigP[myPlace] = neig;
                      neig->neigP[neigPlace] = runParticle;

                      runParticle->springf[myPlace] = average_springf;
                      runParticle->springs[myPlace] = average_springs;
                      runParticle->break_Str[myPlace] = average_break;

                      neig->springf[neigPlace] = average_springf;
                      neig->springs[neigPlace] = average_springs;
                      neig->break_Str[neigPlace] = average_break;

                      // runParticle->springf[myPlace] =
                      // neig->springf[neigPlace]; runParticle->springs[myPlace]
                      // = neig->springs[neigPlace];
                      // runParticle->break_Str[myPlace] =
                      // neig->break_Str[neigPlace];

                      runParticle->rad[myPlace] = dd / 2.0;
                      neig->rad[neigPlace] = dd / 2.0;

                      runParticle->neigpos[myPlace][0] =
                          neig->xpos - runParticle->xpos;
                      runParticle->neigpos[myPlace][1] =
                          neig->ypos - runParticle->ypos;

                      neig->neigpos[neigPlace][0] =
                          runParticle->xpos - neig->xpos;
                      neig->neigpos[neigPlace][1] =
                          runParticle->ypos - neig->ypos;

                      runParticle->draw_break = true;
                      // neig->draw_break = false; // may be problematic

                      // runParticle->healing = time_heal;
                      // neig->healing = time_heal;

                      // runParticle->fracture_time = 0;
                      // neig->fracture_time = 0;

                      cout << "healing" << endl;

                      nbBreak = nbBreak - 1;
                    }
                  }
                }
              }
            }

            if (neig->next_inBox) // if there is anotherone in box
            {
              neig = neig->next_inBox; // shift the help pointer
            } else {
              break; // if not go out to next box
            }
          }
        }
      }
    }
  }
}

void Lattice::Find_Empty_Box(double age, int temp) {
  int x, y, xx;
  double xpos, ypos;
  bool broken1, broken2, space;
  Particle *neig1, *neig2;

  for (x = 0; x < particlex; x++) // was -10
  {
    for (y = 0; y < particlex - 1; y++) {
      xx = (y * particlex * 2) + x;
      if (repBox[xx])
        ;
      else {
        space = true;

        if (repBox[xx + 1]) {
          neig1 = repBox[xx + 1];

          if (repBox[xx - 1]) {
            neig2 = repBox[xx - 1];

            if (neig1->xpos - neig2->xpos < (neig1->radius * 3.3)) {
              space = false;
            }

            while (neig2->next_inBox) {

              if (neig1->xpos - neig2->xpos < (neig1->radius * 3.3)) {
                space = false;
              }
              neig2 = neig2->next_inBox;
            }

          } else
            space = false;
          if (space) {
            while (neig1->next_inBox) {
              neig2 = repBox[xx - 1];

              if (neig1->xpos - neig2->xpos < (neig1->radius * 3.3)) {
                space = false;
              }

              while (neig2->next_inBox) {

                if (neig1->xpos - neig2->xpos < (neig1->radius * 3.3)) {
                  space = false;
                }
                neig2 = neig2->next_inBox;
              }
              neig1 = neig1->next_inBox;
            }
          }
        } else
          space = false;

        // if (broken1 && broken2)
        if (space) {
          xpos = double(double(x) / double(particlex));
          ypos = double(double(y) / double(particlex));
          Grow_Particle(xpos, ypos, xx, 2 * 1e9, 1 * 1e20, age, temp);
        }
      }
    }
  }
}

void Lattice::Cooling(double ref) {
  int i, kk, ll;

  for (i = 0; i < numParticles; i++) {
    runParticle->abreak = 0;
    // runParticle->E = runParticle->E + runParticle->E *
    // (ref/runParticle->Tdif); // Young MOdulus goes up with cooling
    ll = 0;
    for (kk = 0; kk < 8; kk++) {
      if (runParticle->neigP[kk]) {

        // if (runParticle->Tdif > 5)
        {
          runParticle->break_Str[kk] =
              runParticle->break_Str[kk] +
              runParticle->break_Str[kk] *
                  ((runParticle->Tdif) / (ref * ref)); // tdif/ref*ref
          runParticle->springf[kk] =
              runParticle->springf[kk] +
              runParticle->springf[kk] *
                  ((runParticle->Tdif) / (ref * ref)); // 25*ref
          runParticle->springs[kk] =
              runParticle->springs[kk] +
              runParticle->springs[kk] *
                  ((runParticle->Tdif) /
                   (ref * ref)); // Breaking Strength, goes up when cooling
        }
        // cout << runParticle->break_Str[kk]<< endl;
        runParticle->abreak += runParticle->break_Str[kk];
        ll++;
      }
    }

    if (ll > 0)
      runParticle->abreak = runParticle->abreak / ll;
    // else
    // runParticle->abreak = 1;
    // cout << runParticle->abreak << endl;

    runParticle = runParticle->prevP;
  }
  AdjustParticleConstants();
}

void Lattice::CoolingDirect(double ref, double refspring, double visc) {
  int i, kk, ll;

  for (i = 0; i < numParticles; i++) {
    if (runParticle->p_Unode) {
      runParticle->abreak = 0;
      // runParticle->E = runParticle->E + runParticle->E *
      // (ref/runParticle->Tdif); // Young MOdulus goes up with cooling
      ll = 0;
      for (kk = 0; kk < 8; kk++) {
        if (runParticle->neigP[kk]) {

          runParticle->break_Str[kk] =
              average_break *
              (200 / (runParticle->SheetT / ref)); // tdif/ref*ref
          runParticle->springf[kk] =
              average_springf *
              (200 / (runParticle->SheetT / refspring)); // 25*ref
          runParticle->springs[kk] =
              average_springs *
              (200 / (runParticle->SheetT /
                      refspring)); // Breaking Strength, goes up when cooling

          // runParticle->break_Str[kk] =  average_break*
          // (200/(runParticle->SheetT)); //tdif/ref*ref
          // runParticle->springf[kk] = average_springf*
          // (200/runParticle->SheetT)*(200/runParticle->SheetT);
          // runParticle->springs[kk] = average_springs*
          // (200/runParticle->SheetT)*(200/runParticle->SheetT);

          runParticle->abreak += runParticle->break_Str[kk];
          ll++;
        }
      }

      runParticle->young =
          average_springf * (200 / (runParticle->SheetT / refspring));

      runParticle->sheet_visc = 1e20 * (200 / (runParticle->SheetT / visc));

      if (ll > 0)
        runParticle->abreak = runParticle->abreak / ll;
      else
        // runParticle->abreak = 1;
        cout << "no springs" << endl;
    }
    runParticle = runParticle->prevP;
  }
  AdjustParticleConstants();
}

void Lattice::CoolingDirectExp(double ref, double refspring, double visc) {
  int i, kk, ll;
  double EB;

  EB = 500; // Activation Energy divided by Boltzmann constant

  for (i = 0; i < numParticles; i++) {
    if (runParticle->p_Unode) {
      runParticle->abreak = 0;
      // runParticle->E = runParticle->E + runParticle->E *
      // (ref/runParticle->Tdif); // Young MOdulus goes up with cooling
      ll = 0;
      for (kk = 0; kk < 8; kk++) {
        if (runParticle->neigP[kk]) {

          runParticle->break_Str[kk] =
              average_break * ref *
              exp(EB / runParticle->SheetT); // tdif/ref*ref
          runParticle->springf[kk] = average_springf * refspring *
                                     exp(EB / runParticle->SheetT); // 25*ref
          runParticle->springs[kk] =
              average_springs * refspring *
              exp(EB / runParticle
                           ->SheetT); // Breaking Strength, goes up when cooling
          runParticle->abreak += runParticle->break_Str[kk];
          ll++;
        }
      }

      runParticle->young =
          average_springf * refspring * exp(EB / runParticle->SheetT);

      runParticle->sheet_visc = 1e20 * visc * exp(EB / runParticle->SheetT);

      if (ll > 0)
        runParticle->abreak = runParticle->abreak / ll;
      else
        // runParticle->abreak = 1;
        cout << "no springs" << endl;
    }
    runParticle = runParticle->prevP;
  }
  AdjustParticleConstants();
}

void Lattice::Dump_Young(int start, int end) {
  FILE *stat; // file pointer
  int i;      // counters

  stat = fopen("young.txt", "a"); // open statistic output
  // append file

  runParticle = &refParticle; // start

  for (i = 0; i < numParticles; i++) // look which particles are in the
  // box
  {
    if (runParticle->nb >= start && runParticle->nb < end) // larger ymin
    {
      fprintf(stat, " %f", runParticle->particle_age);
      fprintf(stat, " %f", runParticle->SheetT);
      fprintf(stat, " %f", runParticle->abreak);
      fprintf(stat, " %f\n", runParticle->young);
    }
    runParticle = runParticle->nextP;
  }

  // fprintf (stat, "%f", finite_strain);
  // fprintf (stat, " %f", runParticle->xpos);

  fclose(stat); // close file
}

void Lattice::Get_Aperture() {
  int i, ii, jj, kk, ll, pos, box;
  double dx, dy, dd;
  Particle *neig;

  for (i = 0; i < numParticles; i++) {
    if (runParticle->draw_break) // bond is fractured
    {
      for (kk = 0; kk < 8; kk++) // try all
      {
        if (runParticle->neigP[kk]) // if neighbour is connected
        {
          runParticle->neigP[kk]->done = true; // this one was done
        }
      }

      for (ii = 0; ii < 3; ii++) {
        for (jj = 0; jj < 3; jj++) {

          //------------------------------------------
          //  case  upper  or  lower  row (+-size)
          //------------------------------------------

          if (jj == 0)
            box = 2 * particlex * (-1);
          if (jj == 1)
            box = 0;
          if (jj == 2)
            box = 2 * particlex;

          //-----------------------------------------
          //  now  check  this  position
          //  all  particles  in  this  position !!
          //-----------------------------------------

          pos = runParticle->box_pos + (ii - 1) + box;

          //-----------------------------------------
          // repbox should always be positive
          //-----------------------------------------

          if (pos > 0) {
            if (repBox[pos]) {
              neig = repBox[pos]; // helppointer  for  while  loop

              while (neig) // as long as there is somebody
              {
                if (neig->nb != runParticle->nb) {
                  if (!neig->done) {
                    dx = neig->xpos - runParticle->xpos;

                    dy = neig->ypos - runParticle->ypos;

                    dd = sqrt((dx * dx) + (dy * dy));

                    if (dd > runParticle->aperture) {
                      if (dd < (runParticle->radius * 3.0)) {
                        runParticle->aperture =
                            dd - (runParticle->radius * 2.0);
                        cout << runParticle->aperture << endl;
                      }
                    }
                  }
                }
                if (neig->next_inBox) // if there is anotherone in box
                {
                  neig = neig->next_inBox; // shift the help pointer
                } else {
                  break; // if not go out
                }
              }
            }
          }
        }
      }
      for (kk = 0; kk < 8; kk++) // try all
      {
        if (runParticle->neigP[kk]) // if neighbour is connected
        {
          runParticle->neigP[kk]->done = false; // this one was done
        }
      }
    }

    runParticle = runParticle->prevP;
  }
}

void Lattice::ShrinkTemperature(int shrink) {
  int i;
  double ratio;

  for (i = 0; i < numParticles; i++) // loop through particles
  {
    if (shrink) // if particle is in the right
    // grain
    {
      if (runParticle->Tdif > 0.0) {

        runParticle->radius =
            runParticle->radius -
            (runParticle->radius * runParticle->Tdif * 0.000008); // shrink

        ratio = 1 - (runParticle->Tdif * 0.000008);
        cout << ratio << endl;
        for (j = 0; j < 8; j++) {
          if (runParticle->neigP[j]) {
            runParticle->neigpos[j][0] *= ratio; // change
            runParticle->neigpos[j][1] *= ratio;
          }
        }

        // runParticle->radius = runParticle->radius-runParticle->radius *
        // shrink;	// shrink
        // runParticle->a = runParticle->a - runParticle->a * shrink;
        // runParticle->b = runParticle->b - runParticle->b * shrink;

        // for (j = 0; j < 8; j++)
        //{

        // if (runParticle->neigP[j])
        //  {
        // runParticle->neigpos[j][0] *= shrink;	// change
        //  runParticle->neigpos[j][1] *= shrink;
        //  individual
        //  radius
        //}
      }
    }

    runParticle = runParticle->nextP; // go to next
  }
}

void Lattice::ViscousStep(double time) {
  int i, j;
  double alen, alen2, dx, dy, dd, dif, dif1, ratio;

  for (i = 0; i < numParticles; i++) // loop through particles
  {
    if (runParticle->xpos > 0.2 & runParticle->xpos < 0.8) {
      for (j = 0; j < 8; j++) {
        if (runParticle->neigP[j]) {
          alen =
              sqrt((runParticle->neigpos[j][0] * runParticle->neigpos[j][0]) +
                   (runParticle->neigpos[j][1] * runParticle->neigpos[j][1]));
          dx = runParticle->neigP[j]->xpos - runParticle->xpos;
          dy = runParticle->neigP[j]->ypos - runParticle->ypos;
          dd = sqrt((dx * dx) + (dy * dy));

          dif = dd - alen;
          // cout << runParticle->viscosity << endl;
          if (dif > 0.0 & alen > 0.0) {
            dif1 = dif - (dif * exp((-20000000000 * time) / 1e15));
            // cout << dif1 << endl;
            alen2 = dd - dif1;
            ratio = alen2 / alen;
            // cout << ratio << endl;
            runParticle->neigpos[j][0] = runParticle->neigpos[j][0] * ratio;
            runParticle->neigpos[j][1] = runParticle->neigpos[j][1] * ratio;
          }
        }
      }
    }

    runParticle = runParticle->nextP; // go to next
  }
}
