#include "fluid_lattice.h"
#include <iostream>
#include <stdio.h>

using namespace std;

//-------------------------------------------------------------------------------
// fluid lattice February 2016, Daniel, added Advection/Diffusion 2018 Daniel
//
// based on pressure lattice by Irfan Ghani
// ADI from Till Sachau's temperature diffusion code
//
// Derivation of Pressure diffusion equation based on Renauld Toussaint's work
//
// This code deals with a compressible fluid, has a square lattice that is
// linked to the particle lattice of latte. Fluid lattice cell width is twice
// that of the repulsion box of the particle lattice. Repulsion box is used to
// link both
//
// Fluid lattice deals with variations of fluid pressure per node assuming Darcy
// flow of fluids in the background. Reads in porosity from the solid code,
// converts that into a permeability. Gives back fluid pressure gradients as
// Fluid forces in X and Y that act on the solid.
//
//-------------------------------------------------------------------------------

// ---------------------------------------------------------------
// Constructor of Fluid_Lattice class,
// is called by initialize_fluid_lattice function in lattice.cc
//
// at the moment the constructor does not do much, just defines some
// basic parameters
//
// reads in the lattice size that is passed on by the solid lattice
// = amount of particles along X, this is also the width of the
// repulsion box within the initiial configuration.
// Careful that repulsion box width in the solid is twice X in order
// to leave space on the right hand side for extension.
// ---------------------------------------------------------------

Fluid_Lattice::Fluid_Lattice(int lattice_size)
    :

      // -----------------------------------------------------------------------------
      // default variable values. Variables are declared in Pressure_lattice
      // header
      // -----------------------------------------------------------------------------

      dim(100), // the default grid width, but this is read in from the lattice
                // directly

      box(1) // ?
{
  pl_size = lattice_size; // read in the actual lattice size from the solid
                          // (amount of particles in x)

  dim = lattice_size / 2; // fluid lattice size (dim) is half of lattice size

  width = 1.0 / dim; // width of fluid lattic box, non-dimensional

  time_a = 1; // in sec, one day is 86400, default time 1 second

  mu = 1.0e-03; // viscosity default

  compressibility =
      4.5e-10; // water compressibility = 4.5e10 m.m/N at 25 C°, default

  rho_fluid = 1.0e+03; // water mass density (1000 kg/m³), default

  scale_constt = 100; // box size in meters, default, is changed later on
}

//--------------------------------------------------------------------------------------------
// Boundary and initialization conditions, functions that build up the initial
// config.
//--------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------
// this function just sets the background pressure for the fluid lattice in
// cases where there is no initial gradient.
//
// two fluid pressure matrices to calculate delta P in times
//
// called from lattice (void Lattice::Initialize_Fluid_Lattice(double pressure,
// double scale))
//--------------------------------------------------------------------------------------------

void Fluid_Lattice::Background_Pressure(double pressure, double concentration) {
  int i, j;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      oldPf[i][j] = pressure;
      Pf[i][j] = pressure;
      Con[i][j] = concentration;
    }
  }
}

// this was just for debugging, maybe remove

void Fluid_Lattice::Background_Concentration(double concentration) {
  int i, j;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      Con[i][j] = concentration;
    }
  }
}

//---------------------------------------------------------------------------------------------
// Setting the scale (and non-dimensional area of node)
//
// now also includes a change in time so that this has not to be changed in
// fluid lattice
//
// called from lattice (void Lattice::Initialize_Fluid_Lattice(double pressure,
// double scale))
//---------------------------------------------------------------------------------------------

void Fluid_Lattice::Set_Scale(double scale, double time_factor) {
  scale_constt = scale;
  area_node = width * width;
  time_a = time_a * time_factor;
}

//--------------------------------------------------------------------------------------------------
// setting a boundary condition at the start with hydrostatic conditions within
// the box
//--------------------------------------------------------------------------------------------------

void Fluid_Lattice::Set_hydrostatic_gradient(int depth) {
  int i, j;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {

      Pf[i][j] = (depth + (scale_constt * (float((dim - 1) - j) / (dim - 1)))) *
                 1000 * 9.81;
    }
  }
}

//------------------------------------------------------------------------------------------------
// Sending information to the fluid lattice, need porosity or solid fraction per
// node and the solid velocity.
//
// function receives the particle list, which is a pointer to the repulsion box,
// so that we can find the particles easily. Each of the boxes contains pointers
// to the actual particle/particles in the box, so that parameters in the
// particle can be easily accessed.
//
// fracture effect not yet used
//
// The solid fraction is determined by summing up the radii of particles, but
// the fluid radii. Particles have a radius for the elastic part and a DIFFERENT
// radius for the fluid representing the "area" of the particle
//
// called from lattice: void Lattice::Fluid_Parameters()
//-------------------------------------------------------------------------------------------------

void Fluid_Lattice::Read_Porosity_Movement_Matrix(double fracture_effect,
                                                  Particle **list,
                                                  int boundary) {
  int i, j, k, particle_in_box, count;
  double solid_average, vel_x_av, vel_y_av;

  // loop through the matrix, i is x, and j is y (note that Irfan did the
  // opposite in pressure lattice) Each fluid box contains 4 repulsion boxes.

  // calculate an average for the boundary conditions. These boundary conditions
  // turn out to be a major pain there are three different ones at the moment,
  // the average is number 3

  solid_average = 0.0;
  vel_x_av = 0.0;
  vel_y_av = 0.0;
  count = 0;

  for (i = 1; i < dim - 1; i++) // loop in x = i and y = j except for boundaries
  {
    for (j = 1; j < dim - 1; j++) {

      particle_in_box = 0; // count sum of particles in the box (FLUID box)
      rho = 0;             // solid fraction sum up

      vel_x[i][j] = 0; // velocity x matrix
      vel_y[i][j] = 0; // velocity y matrix

      solid_fraction[i][j] = 0.9; // solid fraction matrix

      //------------------------------------------------------------------------------------------------
      // In order to get the right solid fraction and to avoid grid effects
      // converting a triangular solid lattice into the square fluid lattice we
      // use a smoothing function and consider solid in the four repulsion boxes
      // of the fluid node plus the surrounding 12 repulsion boxes. The solid is
      // then weight using the smoothing function depending on how far the
      // particle is from the centre of the pressure node.
      //
      // this is also applied to the velocities.
      //-------------------------------------------------------------------------------------------------

      //-------------------------------------------------------------------------------------------------
      // Because of the above the outer rows of the matrix have to be treated
      // separately, otherwise we look into boxes that dont exist.
      //-------------------------------------------------------------------------------------------------

      // if((i > 0 && i < dim-1) && (j > 0 && j < dim-1))
      {
        //--------------------------------------------------------------------------------------------
        // Now define counters for all 16 repulsion boxes that we have to look
        // into to define the solid fraction for the pressure node as a function
        // of i and j. start in the lower left hand corner with k = 0, then go
        // to the right. So 0,1,2,3 are below the actual pressure node, 5,6,9,10
        // are in the actual pressure node.
        //
        // the repulsion box is a one-dimensional list starting in the lower
        // left hand corner of the lattice with 0 and running towards the right
        // with TWICE the lattice width (so for a 100 lattice up to 199). The
        // box in the second horizontal row above 0 is then 200 and so on.
        //
        // at the moment this routine is not using the boundaries because the
        // boxes there are either not completely full or are empty. There is
        // already the try to double up full boxes for the boundaries, but that
        // does not work yet, needs further debugging
        //--------------------------------------------------------------------------------------------

        for (k = 0; k < 16; k++) {
          if (k == 0) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 1 +
                  ((i - 1) * 2); // pl_size = resolution of particle lattice
            if (i == 0)
              box =
                  (((j * 2) - 1) * (2 * pl_size)) + 3 + ((i - 1) * 2); // take 2
            if (j == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 1 +
                    ((i - 1) * 2); // take 8
            if (j == 0 && i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 3 +
                    ((i - 1) * 2); // take 10

          } else if (k == 1) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 2 + ((i - 1) * 2);
            if (i == 0)
              box =
                  (((j * 2) - 1) * (2 * pl_size)) + 4 + ((i - 1) * 2); // take 3
            if (j == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 2 +
                    ((i - 1) * 2); // take 9
            if (j == 0 && i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 3 +
                    ((i - 1) * 2); // take 10

          } else if (k == 2) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 3 + ((i - 1) * 2);
            if (i == dim - 1)
              box =
                  (((j * 2) - 1) * (2 * pl_size)) + 1 + ((i - 1) * 2); // take 0
            if (j == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 3 +
                    ((i - 1) * 2); // take 10
            if (j == 0 && i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 2 +
                    ((i - 1) * 2); // take 9

          } else if (k == 3) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 4 + ((i - 1) * 2);
            if (i == dim - 1)
              box =
                  (((j * 2) - 1) * (2 * pl_size)) + 2 + ((i - 1) * 2); // take 1
            if (j == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 4 +
                    ((i - 1) * 2); // take 11
            if (j == 0 && i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 2 +
                    ((i - 1) * 2); // take 9

          }

          else if (k == 4) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 1 +
                  ((i - 1) * 2);
            if (i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 3 +
                    ((i - 1) * 2); // take 6
            if (j == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 6 * pl_size + 3 +
                    ((i - 1) * 2); // take 14
            if (j == 0 && i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 3 +
                    ((i - 1) * 2); // take 10

          } else if (k == 5) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 2 +
                  ((i - 1) * 2);
            if (i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 4 +
                    ((i - 1) * 2); // take 7
            if (j == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 6 * pl_size + 2 +
                    ((i - 1) * 2); // take 13
            if (j == 0 && i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 3 +
                    ((i - 1) * 2); // take 10
          } else if (k == 6) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 3 +
                  ((i - 1) * 2);
            if (i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 1 +
                    ((i - 1) * 2); // take 4
            if (j == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 6 * pl_size + 3 +
                    ((i - 1) * 2); // take 14
            if (j == 0 && i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 2 +
                    ((i - 1) * 2); // take 9

          } else if (k == 7) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 4 +
                  ((i - 1) * 2);
            if (i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 2 +
                    ((i - 1) * 2); // take 5
            if (j == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 6 * pl_size + 4 +
                    ((i - 1) * 2); // take 15
            if (j == 0 && i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 2 +
                    ((i - 1) * 2); // take 9

          } else if (k == 8) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 1 +
                  ((i - 1) * 2);
            if (i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 3 +
                    ((i - 1) * 2); // take 10
            if (j == dim - 1)
              box =
                  (((j * 2) - 1) * (2 * pl_size)) + 1 + ((i - 1) * 2); // take 0
            if (j == dim - 1 && i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 3 +
                    ((i - 1) * 2); // take 6

          } else if (k == 9) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 2 +
                  ((i - 1) * 2);
            if (i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 4 +
                    ((i - 1) * 2); // take 11
            if (j == dim - 1)
              box =
                  (((j * 2) - 1) * (2 * pl_size)) + 2 + ((i - 1) * 2); // take 1
            if (j == dim - 1 && i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 3 +
                    ((i - 1) * 2); // take 6

          } else if (k == 10) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 3 +
                  ((i - 1) * 2);
            if (i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 1 +
                    ((i - 1) * 2); // take 8
            if (j == dim - 1)
              box =
                  (((j * 2) - 1) * (2 * pl_size)) + 3 + ((i - 1) * 2); // take 2
            if (j == dim - 1 && i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 2 +
                    ((i - 1) * 2); // take 5

          } else if (k == 11) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 4 +
                  ((i - 1) * 2);
            if (i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 2 +
                    ((i - 1) * 2); // take 9
            if (j == dim - 1)
              box =
                  (((j * 2) - 1) * (2 * pl_size)) + 4 + ((i - 1) * 2); // take 3
            if (j == dim - 1 && i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 2 +
                    ((i - 1) * 2); // take 5

          } else if (k == 12) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 6 * pl_size + 1 +
                  ((i - 1) * 2);
            if (i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 6 * pl_size + 3 +
                    ((i - 1) * 2); // take 14
            if (j == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 1 +
                    ((i - 1) * 2); // take 4
            if (j == dim - 1 && i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 3 +
                    ((i - 1) * 2); // take 6

          } else if (k == 13) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 6 * pl_size + 2 +
                  ((i - 1) * 2);
            if (i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 6 * pl_size + 4 +
                    ((i - 1) * 2); // take 15
            if (j == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 2 +
                    ((i - 1) * 2); // take 5
            if (j == dim - 1 && i == 0)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 3 +
                    ((i - 1) * 2); // take 6

          } else if (k == 14) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 6 * pl_size + 3 +
                  ((i - 1) * 2);
            if (i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 6 * pl_size + 1 +
                    ((i - 1) * 2); // take 12
            if (j == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 3 +
                    ((i - 1) * 2); // take 6
            if (j == dim - 1 && i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 2 +
                    ((i - 1) * 2); // take 5

          } else if (k == 15) {
            box = (((j * 2) - 1) * (2 * pl_size)) + 6 * pl_size + 4 +
                  ((i - 1) * 2);
            if (i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 6 * pl_size + 2 +
                    ((i - 1) * 2); // take 13
            if (j == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 4 +
                    ((i - 1) * 2); // take 7
            if (j == dim - 1 && i == dim - 1)
              box = (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 2 +
                    ((i - 1) * 2); // take 5
          }

          //-----------------------------------------------------------
          // If the box is filled deal with the particle
          //-----------------------------------------------------------

          if (list[box]) {
            help = list[box]; // help now points at particle

            //------------------------------------------------------------------------------
            // first the smoothening function, get distance from centre node to
            // particle in x and y
            //------------------------------------------------------------------------------

            x_diff = abs(help->xpos - ((double(i) * width) + (width / 2.0)));
            y_diff = abs(help->ypos - ((double(j) * width) + (width / 2.0)));

            // adjust if too small

            if (y_diff < (help->radius / 2.0))
              y_diff = help->radius / 2.0;

            if (x_diff < (help->radius / 2.0))
              x_diff = help->radius / 2.0;

            // should not be larger than the width

            if (x_diff > width)
              x_diff = width;
            if (y_diff > width)
              y_diff = width;

            // final smooth function

            smooth_func = (1.0 - (x_diff / width)) * (1.0 - (y_diff / width));

            // add up solid fraction

            rho = rho + smooth_func * help->area_par_fluid;

            // add up velocities

            vel_x[i][j] =
                vel_x[i][j] +
                (help->velx * smooth_func); // velocity (flux) per node
            vel_y[i][j] = vel_y[i][j] + (help->vely * smooth_func);

            // count particles

            particle_in_box++;

            //-----------------------------------------------------------------------
            // check if there is an additional particle in the repulsion box
            // the additional particle would be attached to the help particle
            // with the next_inBox pointer.
            //-----------------------------------------------------------------------

            while (help->next_inBox) {
              // and do the same stuff as above.

              x_diff = abs(help->next_inBox->xpos -
                           ((double(i) * width) + (width / 2.0)));
              y_diff = abs(help->next_inBox->ypos -
                           ((double(j) * width) + (width / 2.0)));

              if (y_diff < (help->radius / 2.0))
                y_diff = help->radius / 2.0;
              if (x_diff < (help->radius / 2.0))
                x_diff = help->radius / 2.0;

              if (x_diff > width)
                x_diff = width;
              if (y_diff > width)
                y_diff = width;

              smooth_func = (1.0 - (x_diff / width)) * (1.0 - (y_diff / width));

              rho = rho + smooth_func * help->next_inBox->area_par_fluid;

              vel_x[i][j] =
                  vel_x[i][j] + (help->next_inBox->velx * smooth_func);
              vel_y[i][j] =
                  vel_y[i][j] + (help->next_inBox->vely * smooth_func);

              particle_in_box++;

              help = help->next_inBox;
            }
          }
        }

        // calculate solid fraction

        solid_fraction[i][j] = rho / area_node;

        solid_average =
            solid_average + solid_fraction[i][j]; // for boundary condition
        count++;

        if (solid_fraction[i][j] >= 1.0) // cannot be more than 1
          solid_fraction[i][j] = 0.999;

        if (solid_fraction[i][j] <= 0.01) // not too small, might check
          solid_fraction[i][j] = 0.01;

        vel_x[i][j] = vel_x[i][j] / rho; // average velocity(flux) per node
        vel_y[i][j] = vel_y[i][j] / rho;

        vel_x_av += vel_x[i][j]; // for boundary condition
        vel_y_av += vel_y[i][j];
      }
    }
  }

  if (boundary == 1) {
    // boundary conditions, apply after matrix filled,
    // condition 1: simply copy the next row from the matrix
    // to the boundary row on sides as well as at bottom and
    // top. Works with random background noise. Becomes
    // problematic when Elle grains are used (boundaries then
    // seem to slow things down because they are mainly oriented 90 degrees
    // to the boundary.
    // condition 2: copy neighbour row to boundary row for upper and lower
    // boundary and wrap rows for left and right boundary. Again can be
    // problematic when grain boundaries are used, but does help avoid
    // "doubling" effects during dynamic permeability development condition 3:
    // take an average of the whole lattice and apply that to boundaries. Needs
    // testing.

    for (j = 1; j < dim - 1; j++) // copy neighbour row
    {
      solid_fraction[0][j] = solid_fraction[1][j];
      vel_x[0][j] = vel_x[1][j];
      vel_y[0][j] = vel_y[1][j];

      solid_fraction[dim - 1][j] = solid_fraction[dim - 2][j];
      vel_x[dim - 1][j] = vel_x[dim - 2][j];
      vel_y[dim - 1][j] = vel_y[dim - 2][j];
    }

    for (i = 0; i < dim; i++) // copy neighbour row
    {
      solid_fraction[i][0] = solid_fraction[i][1];
      vel_x[i][0] = vel_x[i][1];
      vel_y[i][0] = vel_y[i][1];

      solid_fraction[i][dim - 1] = solid_fraction[i][dim - 2];
      vel_x[i][dim - 1] = vel_x[i][dim - 2];
      vel_y[i][dim - 1] = vel_y[i][dim - 2];
    }
  }
  if (boundary == 2) {

    for (j = 1; j < dim - 1; j++) // left and right hand side "wrap" rows
    {
      solid_fraction[0][j] = solid_fraction[dim - 2][j];
      vel_x[0][j] = vel_x[dim - 2][j];
      vel_y[0][j] = vel_y[dim - 2][j];

      solid_fraction[dim - 1][j] = solid_fraction[1][j];
      vel_x[dim - 1][j] = vel_x[1][j];
      vel_y[dim - 1][j] = vel_y[1][j];
    }

    for (i = 0; i < dim; i++) // bottom and top copy neighbour
    {
      solid_fraction[i][0] = solid_fraction[i][1];
      vel_x[i][0] = vel_x[i][1];
      vel_y[i][0] = vel_y[i][1];

      solid_fraction[i][dim - 1] = solid_fraction[i][dim - 2];
      vel_x[i][dim - 1] = vel_x[i][dim - 2];
      vel_y[i][dim - 1] = vel_y[i][dim - 2];
    }
  }
  if (boundary == 3) // use an average
  {
    // at the moment the average works best but needs to be multiplied
    // by a factor, otherwise the boundaries are too permeable
    // the main problem are the side boundaries if the fluid and reaction
    // travel upwards

    for (j = 1; j < dim - 1; j++) // check this boundary condition
    {
      solid_fraction[0][j] = 1.12 * solid_average / count;
      vel_x[0][j] = vel_x_av / count;
      vel_y[0][j] = vel_y_av / count;

      solid_fraction[dim - 1][j] = 1.12 * solid_average / count;
      vel_x[dim - 1][j] = vel_x_av / count;
      vel_y[dim - 1][j] = vel_y_av / count;
    }

    for (i = 0; i < dim; i++) // check this boundary condition
    {
      solid_fraction[i][0] = solid_average / count;
      vel_x[i][0] = vel_x_av / count;
      vel_y[i][0] = vel_y_av / count;

      solid_fraction[i][dim - 1] = solid_average / count;
      vel_x[i][dim - 1] = vel_x_av / count;
      vel_y[i][dim - 1] = vel_y_av / count;
    }
  }
}

//---------------------------------------------------------------------------------
// Caculate permeability from the porosity using Karmen Cozeny relation
// The grain size is now read in directly and changes the advection a lot,
// because of the permeability
//---------------------------------------------------------------------------------

void Fluid_Lattice::Calculate_Permeability_Matrix(double Coseny_grain_size) {
  int i, j;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      kappa[i][j] = 0;

      kappa[i][j] =
          (pow(Coseny_grain_size, 2.0) * pow((1 - solid_fraction[i][j]), 3.0));
      kappa[i][j] /= 45.0 * pow((solid_fraction[i][j]), 2.0);

      if (kappa[i][j] >
          0.00000001) // reset if too large, probably set too small (e11)
      {
        kappa[i][j] = 0.00000001;
        cout << "reset permeability";
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CALCULATION MATRIX INVERSION
// based on implicit alternating direction method, here x and y are the same
// length, so only alpha needed
//
// alpha-set
//
// matrices set 1
//
// multiplication
//
// invert
//
// matrices set 2
//
// multiplication invers
//
// Transpose
//
// source term addition
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------------------------------------------------------
// This just fills the parameter matrix for the calculation. Note that pressure
// is also in this, meaning that the calucation is easily unstable because it
// varies so much if pressure varies. at the moment the matrix size is preset to
// 200 200 200 for latte simulations with 400 particles in x we are not using
// sparse matrices, this is something that should be changed, would speed up
// things this is the prefactor for each node in the grid where pf, permeability
// and porosity vary this prefactor makes the implicit method potentially
// unstable if its too large, this should be checked and corrected for in a
// future version to make the method more stable
//--------------------------------------------------------------------------------------------------------------------

int Fluid_Lattice::alpha_set(int adjust) // need only alpha
{
  int i, j;
  double max;
  max = 0.0;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      alpha[i][j] = 0;

      {
        // alpha[i][j] = (1.0 + (compressibility *
        // Pf[i][j]))*density_factor[i][j];
        alpha[i][j] = (1.0 + (compressibility * Pf[i][j]));
        alpha[i][j] *= kappa[i][j] * time_a;
        alpha[i][j] /=
            (2.0 * mu * compressibility * (1 - solid_fraction[i][j]) *
             pow(width * scale_constt, 2.0));

        alpha[i][j] = alpha[i][j] / adjust;

        if (alpha[i][j] > max)
          max = alpha[i][j];
      }
    }
  }
  cout << max << endl;

  if (max < 10) // was 50
    return 1;
  else
    return (int)max / 10;
}

//-----------------------------------------------------------------------------------
// set the matrix for the first calculation step
// this is the matrix that could be sparse potentially
//-----------------------------------------------------------------------------------

void Fluid_Lattice::matrices_set_1() {
  // creats L.H.S (al[][][])  and R.H.S (bl[][][])coeffiecient matrix
  // at first half time step of ADI method

  int k, j, i;

  for (k = 0; k < dim; k++) {
    for (i = 0; i < dim; i++) {
      for (j = 0; j < dim; j++) {
        al[k][i][j] = 0;
        bl[k][i][j] = 0;

        if (i == j) // for diagonal elements
        {
          if (i == 0 || i == dim - 1) // set boundary conditions
          {
            al[k][i][j] = 1.0; // column set for any given column j
            bl[k][i][j] = 1.0; // row set for any given row i
          } else {
            al[k][i][j] = 1 + 2 * alpha[k][i];
            bl[k][i][j] = 1 - 2 * alpha[i][k];
          }
        } else if ((i == j - 1 && i > 0) || (i == j + 1 && i < (dim - 1))) {
          al[k][i][j] = -alpha[k][i];

          bl[k][i][j] = alpha[i][k];
        }
      }
    }
  }
}

//---------------------------------------------------------------------------
// this just copies the fluid pressure matrix in order to get the change in
// time per node for the darcy velocities that are then need for the
// advection code
//----------------------------------------------------------------------------

void Fluid_Lattice::copy_matrix() {
  int i, j;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      oldPf[i][j] = Pf[i][j];
    }
  }
}

//----------------------------------------------------------------------------
// multiply diagonal matrix with the fluid pressure matrix
//-----------------------------------------------------------------------------

void Fluid_Lattice::multiplication() {
  int i, j, k;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {

      interim_Pf[i][j] = 0;

      for (k = 0; k < dim; k++) {
        interim_Pf[i][j] += bl[j][i][k] * Pf[k][j];
      }
    }
  }
}

//----------------------------------------------------------------------------
// Matrix inversion
//----------------------------------------------------------------------------

void Fluid_Lattice::invert() {
  double factor;
  int g, i, j, k;

  for (g = 0; g < dim; g++) {
    // Creation of identity matrix
    for (i = 0; i < dim; i++) {
      for (j = 0; j < dim; j++) {
        a_invl[g][i][j] = 0;

        a_invl[g][i][j] = (i == j) ? 1 : 0;
      }
    }

    // Development of upper-triangular matrix from the input matrix
    // along with the subsequent changes in identity matrix.
    for (i = 0; i < dim - 1; i++) {
      for (k = i + 1; k < dim; k++) {
        factor = (-al[g][k][i]) / al[g][i][i];
        if (!isinf(factor) && !isnan(factor)) {
          for (j = 0; j < dim; j++) {
            al[g][k][j] += al[g][i][j] * factor;
            a_invl[g][k][j] += a_invl[g][i][j] * factor;
          }
        }
      }
    }

    for (i = 0; i < dim; i++) {
      for (j = 0; j < dim; j++) {
        if (j < i)
          al[g][i][j] = 0;
      }
    }

    // Development of diagonal matrix from developed upper-triangular
    // matrix along with subsequent changes in modified identity matrix.

    for (i = dim - 1; i > 0; i--) {
      for (k = i - 1; k >= 0; k--) {
        factor = (-al[g][k][i]) / al[g][i][i];
        if (!isinf(factor) && !isnan(factor)) {
          for (j = 0; j < dim; j++) {
            al[g][k][j] += al[g][i][j] * factor;
            a_invl[g][k][j] += a_invl[g][i][j] * factor;
          }
        }
      }
    }

    for (i = 0; i < dim; i++) {
      for (j = 0; j < dim; j++) {
        if (j > i)
          al[g][i][j] = 0;
      }
    }

    // Conversion of given matrix into identity matrix and on contrary
    // achievement of inverse matrix from the identity matrix.

    for (i = 0; i < dim; i++) {
      if (al[g][i][i] != 0) {
        factor = al[g][i][i];

        for (j = 0; j < dim; j++) {
          al[g][i][j] /= factor;
          a_invl[g][i][j] /= factor;
        }
      }
    }
  }
}

//------------------------------------------------------------------------------------------
// Set the matrix for the second step of the calculation
//------------------------------------------------------------------------------------------

void Fluid_Lattice::matrices_set_2() {
  // creats L.H.S (al[][][])  and R.H.S (bl[][][])coeffiecient matrix
  // at second half time step of ADI method

  int k, i, j;

  for (k = 0; k < dim; k++) {
    for (i = 0; i < dim; i++) {
      for (j = 0; j < dim; j++) {
        al[k][i][j] = 0;
        bl[k][i][j] = 0;

        if (i == j) // for diagonal elements
        {
          if (i == 0 || i == dim - 1) // set boundary conditions
          {
            al[k][i][j] = 1.0; // column set for any given column j
            bl[k][i][j] = 1.0; // row set for any given row i
          } else {
            al[k][i][j] = 1 + 2 * alpha[i][k];
            bl[k][i][j] = 1 - 2 * alpha[k][i];
          }
        } else if ((i == j - 1 && i > 0) || (i == j + 1 && i < (dim - 1))) {
          al[k][i][j] = -alpha[i][k];

          bl[k][i][j] = alpha[k][i];
        }
      }
    }
  }
}

//--------------------------------------------------------------------------------
// multiply interim matrix with inverted matrix
//--------------------------------------------------------------------------------

void Fluid_Lattice::multiplication_inv() {
  int i, j, k;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      Pf[i][j] = 0;

      for (k = 0; k < dim; k++)
        Pf[i][j] += a_invl[j][i][k] * interim_Pf[k][j];
    }
  }
}

//--------------------------------------------------------------------------------------
// Source term that is added to the pressure diffusion equation, uses pressure
// difference and difference in velocity of solid.
//
// check that velocity is scaled (m/second)
//--------------------------------------------------------------------------------------

void Fluid_Lattice::test_pressure_with_source() {
  int i, j;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      source_a[i][j] = 0;

      if ((i > 0 && i < dim - 1) && (j > 0 && j < dim - 1)) {
        source_a[i][j] = ((vel_x[i][j + 1] - vel_x[i][j - 1]) /
                          (2.0 * width * scale_constt));
        source_a[i][j] += ((vel_y[i + 1][j] - vel_y[i - 1][j]) /
                           (2.0 * width * scale_constt));
        source_a[i][j] *= (1.0 / compressibility + Pf[i][j]);
        source_a[i][j] *= time_a / (2.0 * (1 - solid_fraction[i][j]));
      }
    }
  }

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      source_b[i][j] = 0;

      if ((i > 0 && i < dim - 1) && (j > 0 && j < dim - 1)) {
        source_b[i][j] = vel_x[i][j] * ((Pf[i][j + 1] - Pf[i][j - 1]) /
                                        (2 * width * scale_constt));
        source_b[i][j] += vel_y[i][j] * ((Pf[i + 1][j] - Pf[i - 1][j]) /
                                         (2.0 * width * scale_constt));
        source_b[i][j] *= time_a / 2.0;
      }
    }
  }

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      source[i][j] = 0;

      source[i][j] = source_a[i][j] + source_b[i][j];
    }
  }

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++)
      interim_Pf[i][j] -= source[i][j];
  }
}

//------------------------------------------------------------------------------------
// Transpose the matrix (switch x and y)
//------------------------------------------------------------------------------------

void Fluid_Lattice::transpose_trans_P() {
  double pressure_trans[dim][dim];
  int i, j;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      pressure_trans[i][j] = 0;

      pressure_trans[i][j] = interim_Pf[j][i];
    }
  }

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      interim_Pf[i][j] = 0;

      interim_Pf[i][j] = pressure_trans[i][j];
    }
  }
}

//-----------------------------------------------------------------------------------------------
// concentration calculation Advective part
// at the moment an explicit forward calculation where Ci time n+1 is Ci time n
// minus velocity vector times delta t/delta x times Ci at time n minus Ci-1 at
// time n. Needs to determine the velocity vector first to find direction of
// Ci-1 (backward to flow) need to check stability, and is supposed to induce
// negative diffusion
//
// the prefactor is adjusted internally for the worst case scenario in the code
// to make it stable however, that means every time step is potentially
// different so need to adjust that internally
//------------------------------------------------------------------------------------------------

void Fluid_Lattice::SolveAdvectionLax_Wenddroff(double time, double space,
                                                int boundary) {
  double C_for_x, C_back_x, C_for_y, C_back_y;
  int i, j, k, v, iteration;
  double dx, dy, dt;
  double prefactor, maxprefactor;

  dx = space;
  dy = space;
  dt = time;

  // this boundary condition should probably be applied externally?

  for (v = 0; v < dim; v++) {
    Con[v][0] = 0.5; // setting the lower boundary
    Con[v][dim - 1] = 0.0;
  }

  // determine the largest prefactor

  maxprefactor = 0.0;

  for (i = 1; i < dim - 1; i++) // x loop
  {
    for (j = 1; j < dim - 1; j++) // y loop
    {
      prefactor = fabs(vel_darcyx[i][j] * dt / dx);
      if (prefactor > maxprefactor)
        maxprefactor = prefactor;
      prefactor = fabs(vel_darcyy[i][j] * dt / dy);
      if (prefactor > maxprefactor)
        maxprefactor = prefactor;
    }
  }

  cout << "advection prefactor:" << maxprefactor << endl;

  // adjust prefactor to be 0.2, this should be an internal loop
  // so that the time is constant.

  iteration = maxprefactor / 0.02;

  if (iteration < 1)
    iteration = 1;

  cout << "iteration advection:" << iteration << endl;

  // step 1 calculation of half space back and forward for x and y.

  // step 2
  for (k = 0; k < iteration; k++) {
    for (i = 0; i < dim; i++) // x loop
    {
      for (j = 1; j < dim - 1; j++) // y loop
      {
        if (i == 0) // boundary conditions are letting fluid flow in y parallel
                    // to boundary, x term constant
        {
          if (boundary == 1) {
            C_for_x =
                0.5 * (Con[i][j] + Con[i + 1][j]) -
                0.5 * vel_darcyx[i][j] * dt / dx * (Con[i + 1][j] - Con[i][j]);
            C_back_x = 0.5 * (Con[i][j] + Con[dim - 1][j]) +
                       0.5 * vel_darcyx[i][j] * dt / dx *
                           (Con[dim - 1][j] - Con[i][j]);

            // step 2 x
            dxc[i][j] =
                Con[i][j] - vel_darcyx[i][j] * dt / dx * (C_for_x - C_back_x);
          }
          if (boundary == 2) //
            dxc[i][j] = Con[i][j];
        } else if (i == dim - 1) // wrapping right boundary
        {
          if (boundary == 1) {
            C_for_x =
                0.5 * (Con[i][j] + Con[0][j]) -
                0.5 * vel_darcyx[i][j] * dt / dx * (Con[0][j] - Con[i][j]);
            C_back_x =
                0.5 * (Con[i][j] + Con[i - 1][j]) +
                0.5 * vel_darcyx[i][j] * dt / dx * (Con[i - 1][j] - Con[i][j]);

            // step 2 x
            dxc[i][j] =
                Con[i][j] - vel_darcyx[i][j] * dt / dx * (C_for_x - C_back_x);
          }
          if (boundary == 2) //
            dxc[i][j] = Con[i][j];
        } else {

          // step 1 x
          C_for_x =
              0.5 * (Con[i][j] + Con[i + 1][j]) -
              0.5 * vel_darcyx[i][j] * dt / dx * (Con[i + 1][j] - Con[i][j]);
          C_back_x =
              0.5 * (Con[i][j] + Con[i - 1][j]) +
              0.5 * vel_darcyx[i][j] * dt / dx * (Con[i - 1][j] - Con[i][j]);

          // step 2 x
          dxc[i][j] =
              Con[i][j] - vel_darcyx[i][j] * dt / dx * (C_for_x - C_back_x);
        }

        // step 1 y
        C_for_y =
            0.5 * (Con[i][j] + Con[i][j + 1]) -
            0.5 * vel_darcyy[i][j] * dt / dy * (Con[i][j + 1] - Con[i][j]);
        C_back_y =
            0.5 * (Con[i][j] + Con[i][j - 1]) +
            0.5 * vel_darcyy[i][j] * dt / dy * (Con[i][j - 1] - Con[i][j]);

        // step 2 y
        dyc[i][j] =
            Con[i][j] - vel_darcyy[i][j] * dt / dy * (C_for_y - C_back_y);
      }
    }
    for (i = 0; i < dim; i++) // x loop
    {
      for (j = 1; j < dim - 1; j++) // y loop
      {
        // add up x and y

        Con[i][j] = 0.5 * (dxc[i][j] + dyc[i][j]);

        // final adjustment, no overflow and no negative concentration to
        // improve stability

        if (Con[i][j] > 0.5)
          Con[i][j] = 0.5;
        if (Con[i][j] < 0.0)
          Con[i][j] = 0.0;
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------------
// this is not used at the moment, but can potentially be used to change
// concentration after a reaction
//
// sept 2018
//--------------------------------------------------------------------------------------------------------

void Fluid_Lattice::ChangeConcentration(int box_x, int box_y, double change) {
  Con[box_y][box_x] = Con[box_y][box_x] - change;
  if (Con[box_y][box_x] < 0.0)
    Con[box_y][box_x] = 0.0;
}

//------------------------------------------------------------------------------------------
// Passing back the values to the solid lattice including the fluid pressure
// gradients as force in x and y that is added to the force ballance in
// relaxation, but not to the calculated stress, so that the stress remains to
// be the solid stress only.
//
// This takes each pressure node and looks into its 4 repulsion boxes and gives
// the values to the particles in there.
//
// can potentially calculate average concetration between concentration nodes,
// but this does not seem to work well at the moment sept 2018
//
// The pressure gradient is taken from the particle repulsion box to the
// neighbouring node
//-------------------------------------------------------------------------------------------

void Fluid_Lattice::Pass_Back_Gradients(Particle **list, int av_conc) {
  int i, j, k;
  Particle *help;
  double smooth, grad_x, grad_y, conc_av_x, conc_av_y;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      //--------------------------------------------------------------------
      // 0 is lower left box, 1 lower right box, 2 upper left box
      // and 3 upper right box.
      //
      // we need special boundary conditions for the boundary boxes
      // left and right can be wrappinp, but up and down may be more
      // problematic
      //
      // at the moment right and left wrap, up and down copy next row
      //--------------------------------------------------------------------

      for (k = 0; k < 4; k++) // go to four quadrants of pressure node
      {
        if (k == 0)
          box =
              (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 2 + ((i - 1) * 2);
        else if (k == 1)
          box =
              (((j * 2) - 1) * (2 * pl_size)) + 2 * pl_size + 3 + ((i - 1) * 2);
        else if (k == 2)
          box =
              (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 2 + ((i - 1) * 2);
        else if (k == 3)
          box =
              (((j * 2) - 1) * (2 * pl_size)) + 4 * pl_size + 3 + ((i - 1) * 2);

        if (list[box]) {
          help = list[box];
          x_diff = abs(help->xpos - ((double(i) * width) + (width / 2.0)));
          y_diff = abs(help->ypos - ((double(j) * width) + (width / 2.0)));
          smooth = (1 - (x_diff / width)) * (1 - (y_diff / width));

          if (k == 0) // lower left box
          {
            if (i == 0) // left boundary, deal with x
            {
              grad_x = Pf[dim - 1][j] - Pf[i][j]; // wrap
              conc_av_x = (Con[dim - 1][j] + Con[i][j]) / 2.0;
            } else // normal
            {
              grad_x = Pf[i - 1][j] - Pf[i][j];
              conc_av_x = (Con[dim - 1][j] + Con[i][j]) / 2.0;
            }
            if (j == 0) // lower boundary, deal with y
            {
              grad_y =
                  Pf[i][j] - Pf[i][j + 1]; // use the next row and copy that
              conc_av_y = (Con[i][j] + Con[i][j + 1]) / 2.0;
            } else // normal
            {
              grad_y = Pf[i][j - 1] - Pf[i][j];
              conc_av_y = (Con[i][j - 1] + Con[i][j]) / 2.0;
            }
          } else if (k == 1) // lower right box
          {
            if (i == dim - 1) // wrap
            {
              grad_x = Pf[i][j] - Pf[0][j];
              conc_av_x = (Con[i][j] + Con[0][j]) / 2.0;
            } else // normal
            {
              grad_x = Pf[i][j] - Pf[i + 1][j];
              conc_av_x = (Con[i][j] + Con[i + 1][j]) / 2.0;
            }
            if (j == 0) {
              grad_y = Pf[i][j] - Pf[i][j + 1];
              conc_av_y = (Con[i][j] + Con[i][j + 1]) / 2.0;
            } else {
              grad_y = Pf[i][j - 1] - Pf[i][j];
              conc_av_y = (Con[i][j - 1] + Con[i][j]) / 2.0;
            }
          } else if (k == 2) {
            if (i == 0) // left boundary
            {
              grad_x = Pf[dim - 1][j] - Pf[i][j];
              conc_av_x = (Con[dim - 1][j] + Con[i][j]) / 2.0;
            } else {
              grad_x = Pf[i - 1][j] - Pf[i][j];
              conc_av_x = (Con[i - 1][j] + Con[i][j]) / 2.0;
            }
            if (j == dim - 1) {
              grad_y = Pf[i][j - 1] - Pf[i][j];
              conc_av_y = (Con[i][j - 1] + Con[i][j]) / 2.0;
            } else {
              grad_y = Pf[i][j] - Pf[i][j + 1];
              conc_av_y = (Con[i][j] + Con[i][j + 1]) / 2.0;
            }

          } else if (k == 3) {
            if (i == dim - 1) {
              grad_x = Pf[i][j] - Pf[0][j];
              conc_av_x = (Con[i][j] + Con[0][j]) / 2.0;
            } else {
              grad_x = Pf[i][j] - Pf[i + 1][j];
              conc_av_x = (Con[i][j] + Con[i + 1][j]) / 2.0;
            }
            if (j == dim - 1) {
              grad_y = Pf[i][j - 1] - Pf[i][j];
              conc_av_x = (Con[i][j - 1] + Con[i][j]) / 2.0;
            } else {
              grad_y = Pf[i][j] - Pf[i][j + 1];
              conc_av_x = (Con[i][j] + Con[i][j + 1]) / 2.0;
            }
          }

          // pass the pressure gradient to each particle as force in x and y

          help->F_P_y = grad_y * 3.14 * help->real_radius * help->real_radius /
                        solid_fraction[i][j];
          help->F_P_x = grad_x * 3.14 * help->real_radius * help->real_radius /
                        solid_fraction[i][j];
          help->F_P_y = help->F_P_y / (width * scale_constt);
          help->F_P_x = help->F_P_x / (width * scale_constt);

          // pass the darcy velocities in x and y. note that they are calculated
          // in a differnt function

          help->v_darcyx = vel_darcyx[i][j];
          help->v_darcyy = vel_darcyy[i][j];

          // pass aditional values for visualization purposes to the particles

          help->temperature = Pf[i][j];
          help->fluid_pressure_gradient = help->F_P_y;
          help->average_porosity = 1 - solid_fraction[i][j];

          // use or dont use the average concentration

          if (av_conc)
            help->conc = (conc_av_x + conc_av_y) / 2.0;
          else
            help->conc = Con[i][j];

          while (help->next_inBox) // and do the whole thing for other particles
                                   // in the box
          {
            help = help->next_inBox;

            x_diff = abs(help->xpos - ((double(i) * width) + (width / 2.0)));
            y_diff = abs(help->ypos - ((double(j) * width) + (width / 2.0)));
            smooth = (1 - (x_diff / width)) * (1 - (y_diff / width));

            if (k == 0) // lower left box
            {
              if (i == 0) // left boundary, deal with x
              {
                grad_x = Pf[dim - 1][j] - Pf[i][j]; // wrap
                conc_av_x = (Con[dim - 1][j] + Con[i][j]) / 2.0;
              } else // normal
              {
                grad_x = Pf[i - 1][j] - Pf[i][j];
                conc_av_x = (Con[i - 1][j] + Con[i][j]) / 2.0;
              }
              if (j == 0) // lower boundary, deal with y
              {
                grad_y =
                    Pf[i][j] - Pf[i][j + 1]; // use the next row and copy that
                conc_av_y = (Con[i][j] + Con[i][j + 1]) / 2.0;
              } else // normal
              {
                grad_y = Pf[i][j - 1] - Pf[i][j];
                conc_av_y = (Con[i][j - 1] + Con[i][j]) / 2.0;
              }
            } else if (k == 1) // lower right box
            {
              if (i == dim - 1) // wrap
              {
                grad_x = Pf[i][j] - Pf[0][j];
                conc_av_x = (Con[i][j] + Con[0][j]) / 2.0;
              } else // normal
              {
                grad_x = Pf[i][j] - Pf[i + 1][j];
                conc_av_x = (Con[i][j] + Con[i + 1][j]) / 2.0;
              }
              if (j == 0) {
                grad_y = Pf[i][j] - Pf[i][j + 1];
                conc_av_y = (Con[i][j] + Con[i][j + 1]) / 2.0;
              } else {
                grad_y = Pf[i][j - 1] - Pf[i][j];
                conc_av_y = (Con[i][j - 1] + Con[i][j]) / 2.0;
              }
            } else if (k == 2) {
              if (i == 0) // left boundary
              {
                grad_x = Pf[dim - 1][j] - Pf[i][j];
                conc_av_x = (Con[dim - 1][j] + Con[i][j]) / 2.0;
              } else {
                grad_x = Pf[i - 1][j] - Pf[i][j];
                conc_av_x = (Con[i - 1][j] + Con[i][j]) / 2.0;
              }
              if (j == dim - 1) {
                grad_y = Pf[i][j - 1] - Pf[i][j];
                conc_av_y = (Con[i][j - 1] + Con[i][j]) / 2.0;
              } else {
                grad_y = Pf[i][j] - Pf[i][j + 1];
                conc_av_y = (Con[i][j] + Con[i][j + 1]) / 2.0;
              }

            } else if (k == 3) {
              if (i == dim - 1) {
                grad_x = Pf[i][j] - Pf[0][j];
                conc_av_x = (Con[i][j] + Con[0][j]) / 2.0;
              } else {
                grad_x = Pf[i][j] - Pf[i + 1][j];
                conc_av_x = (Con[i][j] + Con[i + 1][j]) / 2.0;
              }
              if (j == dim - 1) {
                grad_y = Pf[i][j - 1] - Pf[i][j];
                conc_av_y = (Con[i][j - 1] + Con[i][j]) / 2.0;
              } else {
                grad_y = Pf[i][j] - Pf[i][j + 1];
                conc_av_y = (Con[i][j] + Con[i][j + 1]) / 2.0;
              }
            }

            // fluid forces from pressure gradients for lattice

            help->F_P_y = grad_y * 3.14 * help->real_radius *
                          help->real_radius / solid_fraction[i][j];
            help->F_P_x = grad_x * 3.14 * help->real_radius *
                          help->real_radius / solid_fraction[i][j];
            help->F_P_y = help->F_P_y / (width * scale_constt);
            help->F_P_x = help->F_P_x / (width * scale_constt);

            // output values for visualization

            help->v_darcyx = vel_darcyx[i][j];
            help->v_darcyy = vel_darcyy[i][j];
            help->temperature = Pf[i][j];
            help->fluid_pressure_gradient = help->F_P_y;
            help->average_porosity = 1 - solid_fraction[i][j];

            if (av_conc)
              help->conc = (conc_av_x + conc_av_y) / 2.0;
            else
              help->conc = Con[i][j];
          }
        }
      }
    }
  }
}

//--------------------------------------------------------------------------------------------------------------
// This just sets some boundaries for the fluid lattice plus potentially one for
// the sides, but that is not necessarily need because the IDE method using
// parallel flow along the boundaries already
//--------------------------------------------------------------------------------------------------------------

void Fluid_Lattice::Set_Boundaries(double UpperBoundaryPressure,
                                   double LowerBoundaryPressure,
                                   int double_sides, double conc,
                                   double additional) {
  int i, j;

  for (i = 0; i < dim; i++) {
    if (UpperBoundaryPressure != 0.0) {
      Pf[i][dim - 1] = UpperBoundaryPressure + additional;
      Con[i][dim - 1] = conc;
    }

    if (LowerBoundaryPressure != conc)

    {
      Pf[i][0] = LowerBoundaryPressure + additional;
      // Pf[i][0] = Pf[i][1];
      Con[i][0] = conc;
    }
    if (double_sides) {
      Pf[0][i] = UpperBoundaryPressure + additional;
      Con[0][i] = conc;
      Pf[dim - 1][i] = UpperBoundaryPressure + additional;
      Con[dim - 1][i] = conc;
      Pf[0][i] = Pf[1][i];
      Pf[dim - 1][i] = Pf[dim - 2][i];
    }
  }
}

//-------------------------------------------------------------------------------------------
// simple increase of pressure at a single node
//-------------------------------------------------------------------------------------------

void Fluid_Lattice::Fluid_Input(double pressure, int x, int y) {
  Pf[x][y] = Pf[x][y] + pressure;
  Con[x][y] += 0.5;
}

//-------------------------------------------------------------------------------------------
// this is the core for the advection code, gets the darcy velocity
// at the moment a function of changes in pressure in the pressure diffusion
// code. that means the concentration travels as a function of pressure changes
// for the advection
//
// first get delta P and then the unit vectors for old grad in x and y for the
// direction of flow.
//
// This could potentially be adjusted but works well at the moment sept 2018
//-------------------------------------------------------------------------------------------

void Fluid_Lattice::Get_Fluid_Velocity() {
  int i, j;
  double grad_xold, grad_yold, unit_vector;
  double deltaPf;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      // get delta P, Pf change per time step

      deltaPf = Pf[i][j] - oldPf[i][j];

      if (i == 0) // flow parallel to left and right hand boundaries
      {
        if (j == 0) {
          grad_xold = 0.0;
          grad_yold = 1.0;
        } else if (j == dim - 1) {
          grad_xold = 0.0;
          grad_yold = 1.0;
        } else {
          grad_xold = 0.0;
          grad_yold = oldPf[i][j - 1] - oldPf[i][j + 1];
        }
      } else if (i == dim - 1) // right boundary wrap
      {
        if (j == 0) {
          grad_xold = 0.0;
          grad_yold = 1.0;
        } else if (j == dim - 1) {
          grad_xold = 0.0;
          grad_yold = 1.0;
        } else {
          grad_xold = 0.0;
          grad_yold = oldPf[i][j - 1] - oldPf[i][j + 1];
        }
      } else if (j == 0) // bottom static only upwards
      {
        grad_xold = 0.0;
        grad_yold = 1.0;
      } else if (j == dim - 1) // top static only upwards
      {
        grad_xold = 0.0;
        grad_yold = 1.0;
      } else // normal case
      {
        grad_xold = oldPf[i - 1][j] - oldPf[i + 1][j];
        grad_yold = oldPf[i][j - 1] - oldPf[i][j + 1];
      }
      // get the unit vectors for x and y

      unit_vector = sqrt((grad_xold * grad_xold) + (grad_yold * grad_yold));

      if (unit_vector != 0.0) {
        grad_xold = grad_xold / unit_vector;
        grad_yold = grad_yold / unit_vector;
      } else {
        grad_xold = 0.0;
        grad_yold = 0.0;
      }

      // and finally calculate the darcy velocity in x and y

      vel_darcyx[i][j] =
          grad_xold * deltaPf * kappa[i][j] /
          (mu * (1 - solid_fraction[i][j]) * scale_constt * width);
      vel_darcyy[i][j] =
          grad_yold * deltaPf * kappa[i][j] /
          (mu * (1 - solid_fraction[i][j]) * scale_constt * width);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CALCULATION MATRIX INVERSION for the concentration, simpler than pressure,
// does not include source term
//
// based on implicit alternating direction method, here x and y are the same
// length, so only alpha needed
//
// alpha-set
//
// matrices set 1
//
// multiplication
//
// invert
//
// matrices set 2
//
// multiplication invers
//
// Transpose
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------------------------------------------------------
// This just fills the parameter matrix for the calculation.
// at the moment the matrix size is preset to 200 200 200 for latte simulations
// with 400 particles in x. We are not using sparse matrices, this is something
// that should be changed, would speed up things. This is the prefactor for each
// node in the grid with a constant diffusion coefficient for the moment. this
// prefactor makes the implicit method potentially unstable if its too large,
// this should be checked and corrected for in a future version to make the
// method more stable.
//--------------------------------------------------------------------------------------------------------------------

void Fluid_Lattice::SolveDiffusionImplicit(double dif_const, double time,
                                           double space) {
  alpha_set_dif(dif_const, time);

  matrices_set_1();
  multiplication_dif();
  transpose_trans_dif();
  invert();
  multiplication_inv_dif();

  matrices_set_2();
  multiplication_dif();
  transpose_trans_dif();
  invert();
  multiplication_inv_dif();
}

//--------------------------------------------------------------------------------------------------------------------
// This just fills the parameter matrix for the calculation. Note that pressure
// is also in this, meaning that the calucation is easily unstable because it
// varies so much if pressure varies.
//
// this difusion constant is constant at the moment, but could be variable
// depending on porosity/permeability of nodes
//--------------------------------------------------------------------------------------------------------------------

void Fluid_Lattice::alpha_set_dif(double dif_constant,
                                  double time_dif) // need only alpha
{
  int i, j;
  double max;

  max = 0;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      alpha[i][j] = 0;
      {
        alpha[i][j] = time_dif * dif_constant / pow(width * scale_constt, 2.0);
        if (alpha[i][j] > max)
          max = alpha[i][j];
      }
    }
  }

  cout << "dif prefactor:" << max << endl;
}

//----------------------------------------------------------------------------
// multiply diagonal matrix with the fluid pressure matrix
//-----------------------------------------------------------------------------

void Fluid_Lattice::multiplication_dif() {
  int i, j, k;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {

      ddyc[i][j] = 0; // ddyc is just used as interim here

      for (k = 0; k < dim; k++) {
        ddyc[i][j] += bl[j][i][k] * Con[k][j];
      }
    }
  }
}

//------------------------------------------------------------------------------------
// Transpose the matrix (switch x and y)
//------------------------------------------------------------------------------------

void Fluid_Lattice::transpose_trans_dif() {
  double con_trans[dim][dim];
  int i, j;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      con_trans[i][j] = 0;

      con_trans[i][j] = ddyc[j][i];
    }
  }

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      ddyc[i][j] = 0;

      ddyc[i][j] = con_trans[i][j];
    }
  }
}

//--------------------------------------------------------------------------------
// multiply interim matrix with inverted matrix
//--------------------------------------------------------------------------------

void Fluid_Lattice::multiplication_inv_dif() {
  int i, j, k;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      Con[i][j] = 0;

      for (k = 0; k < dim; k++)
        Con[i][j] += a_invl[j][i][k] * ddyc[k][j];
    }
  }
}

// this is an old version that did not work well and was not stable

void Fluid_Lattice::SolveadvDisExplicitTransient(double fluid_v,
                                                 double dif_const) {
  int i, j, v;
  double D, dt, dx, dy, dx2, dy2, conc_x, conc_y;

  D = dif_const;

  dx = width; //
  dy = width; //
  dx2 = dx * dx;
  for (v = 0; v < dim; v++) {
    Con[v][0] = 0.5; // setting the lower boundary
    Con[v][1] = 0.5;
    // Con[v][dim-1]=0;
  }
  dy2 = dy * dy;

  dt = 0.002;

  // boundary condition (background concentration is already given)

  for (i = 0; i < dim - 1; i++) // x loop
  {
    for (j = 1; j < dim - 2; j++) // y loop
    {

      ddyc[i][j] = (Con[i + 1][j] - 2 * Con[i][j] + Con[i - 1][j]) /
                   dy2; // x has to vary in y, so takes neighbours in i (=y)
      ddxc[i][j] = (Con[i][j + 1] - 2 * Con[i][j] + Con[i][j - 1]) /
                   dx2; // y has to vary in x, so takes neighbours in j (=x)

      // forward approximation takes the difference between node and node behind
      // (with respect to velocity) and adds to node

      // if (vel_darcyx[i][j]>=0.0) // x direction positive, to the right
      if (vel_darcyx[i][j] >= 0.0) {
        dxc[i][j] = (Con[i][j] - Con[i - 1][j]) / dx;
        if (dxc[i][j] > 0.0)
          dxc[i][j] = 0.0;
      } else // x direction negative, to the left
      {
        dxc[i][j] = (Con[i][j] - Con[i + 1][j]) / dx;
        if (dxc[i][j] > 0.0)
          dxc[i][j] = 0.0;
      }
      if (vel_darcyy[i][j] >= 0.0) // y direction positive, upwards
      {
        dyc[i][j] = (Con[i][j] - Con[i][j - 1]) / dy;
        if (dyc[i][j] > 0.0)
          dyc[i][j] = 0.0;
      } else // y direction positive, downwards
      {
        dyc[i][j] = (Con[i][j] - Con[i][j + 1]) / dy;
        if (dyc[i][j] > 0.0)
          dyc[i][j] = 0.0;
      }
    }
  }

  // you now have all the values to fill in the matrix and update the
  // concentration. You have to do that after the calculation, otherwise you get
  // a trend in the calculation depending on where you start because some have
  // been updated and others not. Need to update in on step.

  for (i = 0; i < dim - 1; i++) // y
  {
    for (j = 1; j < dim - 2; j++) // x
    {
      // diffusion
      // Con[i][j] = Con[i][j]+dt*D*(ddxc[i][j]+ddyc[i][j]);
      // cout << " too large y" << dt*abs(vel_darcyy[i][j])*fluid_v*dyc[i][j] <<
      // endl;
      // advection
      // if (dt*abs(vel_darcyx[i][j])*fluid_v*dxc[i][j]>0.1)
      // cout << " too large x" << dt*abs(vel_darcyx[i][j])*fluid_v*dxc[i][j] <<
      // endl; if (dt*abs(vel_darcyy[i][j])*fluid_v*dyc[i][j]<-0.0) cout << "
      // too large y" << dt*abs(vel_darcyy[i][j])*fluid_v*dyc[i][j] << endl;
      // else

      conc_x = Con[i][j] - (dt * abs(vel_darcyx[i][j]) * fluid_v * dxc[i][j]);
      conc_y = Con[i][j] - (dt * abs(vel_darcyy[i][j]) * fluid_v * dyc[i][j]);

      if (j < 4) {
        cout << "concentration x" << conc_x << endl;
        cout << "concentration y" << conc_y << endl;
      }

      Con[i][j] = conc_x + conc_y;

      // else if (vel_darcyy[i][j]< -0.00000001)
      //{
      //     Con[i][j] = Con[i][j]-(dt*vel_darcyx[i][j]*fluid_v*dxc[i][j]);
      //     Con[i][j] = Con[i][j]-(dt*vel_darcyy[i][j]*fluid_v*dyc[i][j]);
      // }

      // original

      // concent[i][j] = concent[i][j]-(dt*1000000000000*kappa[i][j]*dyc[i][j]);
      if (Con[i][j] < 0.0)
        Con[i][j] = 0.0;
      else if (Con[i][j] > 0.5)
        Con[i][j] = 0.5;
    }
  }

  for (v = 1; v < dim - 1;
       v++) // boundary condition for sides, not really wrapping... sort of
  {
    Con[0][v] = Con[1][v];

    Con[dim - 1][v] = Con[dim - 2][v];
  }
}

void Fluid_Lattice::SetFluidPhreeqc(const char *pore_fluid,
                                    const char *inf_fluid, const char *d_base) {
  iphreeqc_obj.Version();
  iphreeqc_obj.SetDatabase(d_base);
  iphreeqc_obj.SetPoreFluid(pore_fluid);
  iphreeqc_obj.SetInfFluid(inf_fluid);
  iphreeqc_obj.InitializeFluidVector();

  iphreeqc_obj.AdjustUnits(iphreeqc_obj.Conc_inf);

  // Set the pore fluid
}

void Fluid_Lattice::TestPhreeqc() {
  // double conc_Zn = iphreeqc.FindConcentration("Zn");
  // cout << "\nConcentration of Zn in the infiltrating fluid: " << conc_Zn <<
  // endl;
  /**We have two types of fluid data structures, one large container that holds
  the concentrations for the entire domain and one simple vector that contains
  the concentrations of infiltrating fluid. In fluid_lattice we have an object
  adn a pointer to the phreeqc class. However, for teh parallelization we need
  to create new objects adn load teh databse for each of them **/

  for (PhreeqC::FLUID ::iterator it = iphreeqc_obj.Conc_inf.begin();
       it != iphreeqc_obj.Conc_inf.end(); it++) {
    if (it->first != "ph") {
      const char *type = it->first; // nmae of the element or parameter
      double inf_conc = it->second; // concentration in filitrating fluid

      double conc[dimX][dimY]; //	this is the matrix we use for transport
                               //calculation
      iphreeqc_obj.GetConcentration(
          type, conc); // this is the function that fills the matrix

      // Transport

      iphreeqc_obj.SetConcentration(
          type,
          conc); // this function passes the matric back to the fluid array
    }
  }

  iphreeqc_obj.Equilibrate(iphreeqc_obj.Conc_inf, "calcite");

  // here we do the parallelization
  cout << "Reaction calculation " << endl;
  double react[dimX][dimY]; // The reaction matrix, we get the dimensions from
                            // the macros defined in fluid_lattice.h
  boost::progress_display *show_progress =
      new boost::progress_display(dimX * dimY);
  time_t start, end;
  time(&start);
#pragma omp parallel shared(react) num_threads(iphreeqc_obj.nb_threads)
  {
    // this should be private for every thread
    PhreeqC iphreeqc_o = iphreeqc_obj; // object from copy constructor
    iphreeqc_o.loader(
        iphreeqc_obj.p_database); // load the databse for this object
    PhreeqC *i_ptr;

// this is the parallel loop, we collapse the nested loop because we can :-)
#pragma omp for collapse(2)
    for (int i = 0; i < dimX; i++) {
      for (int j = 0; j < dimY; j++) {
        // printf("i = %d, j= %d, threadId = %d \n", i, j,
        // omp_get_thread_num());				//for degugging
        // we loop which thread got which entry this actually contains two
        // functions: SI which calculates the saturation index defined by the
        // second argument (mineral name) and the function GetFluid which build
        // the fluid composition at a specific point
        double si_calcite =
            iphreeqc_o.SI(iphreeqc_o.GetFluid(i, j), "calcite", false);

        // here is the horrible function that calulates reaction rates (full of
        // magic numbers). To see what they are look in phreeqc.cc
        iphreeqc_o.React(iphreeqc_o.GetFluid(i, j), -13.7, 0.67, 22.6, 0.162,
                         103, 15, 1.5768e8, "Quartz", "SiO2", false);

        react[i][j] = si_calcite; // fill the matrix
        ++(*show_progress);
      }
    }
  }
  cout << endl;

  time(&end);
  // Calculating total time taken by the program.
  double time_taken = double(end - start);
  cout << "Time taken by phreeqc: " << fixed << time_taken << setprecision(10);
  cout << " sec " << endl;
}
