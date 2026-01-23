/*****************************************************
 * Copyright: (c) L. A. Evans
 * File:      $RCSfile: erand.cc,v $
 * Revision:  $Revision: 1.1 $
 * Date:      $Date: 2006/03/05 05:29:47 $
 * Author:    $Author: levans $
 *
 ******************************************************/
/*!
        \file		erand.cc
        \brief		random number generator
        \par		Description:
                Interface for gsl random number generators
*/
#include "erand.h"

Erand::Erand() {
  _rgen = gsl_rng_alloc(gsl_rng_ran3);
  /*
   * if linking to the gsl lib, default is mt19937
   * can set generator with env var GSL_RNG
  _rgen = gsl_rng_alloc(gsl_rng_env_setup());
   */
  gsl_rng_set(_rgen, gsl_rng_default_seed);
}

Erand::Erand(unsigned long int s) {
  _rgen = gsl_rng_alloc(gsl_rng_ran3);
  /*
   * if linking to the gsl lib, default is mt19937
   * can set generator with env var GSL_RNG
  _rgen = gsl_rng_alloc(gsl_rng_env_setup());
   */
  gsl_rng_set(_rgen, s);
}

void Erand::seed(unsigned long int seed) { gsl_rng_set(_rgen, seed); }

unsigned long int Erand::ran() { return (gsl_rng_get(_rgen)); }

double Erand::randouble() { return (gsl_rng_uniform(_rgen)); }
