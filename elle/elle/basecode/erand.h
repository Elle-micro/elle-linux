/*****************************************************
 * Copyright: (c) L. A. Evans
 * File:      $RCSfile: erand.h,v $
 * Revision:  $Revision: 1.3 $
 * Date:      $Date: 2006/03/05 05:29:47 $
 * Author:    $Author: levans $
 *
 ******************************************************/
/*!
        \file		erand.h
        \brief		header for random number
        \par		Description:
                Class for accessing random number generators
*/
#ifndef E_erand_h
#define E_erand_h
/*************************************************************
 *	INCLUDE FILES
 */
#include "gsl/gsl_rng.h"
#include "stdio.h"
/*************************************************************
 *	CONSTANT DEFINITIONS
 */
/*************************************************************
 *	MACRO DEFINITIONS
 */
/************************************************************
 *	ENUMERATED DATA TYPES
 */
/*************************************************************
 *	STRUCTURE DEFINITIONS
 */
/*************************************************************
 *	IN-LINE FUNCTION DEFINITIONS
 */
/*************************************************************
 *	CLASS DECLARATIONS
 */
class Erand {
  gsl_rng *_rgen;

public:
  Erand();
  Erand(unsigned long int s);
  void seed(unsigned long int seed);
  unsigned long int ran();
  double randouble();
  ~Erand() { gsl_rng_free(_rgen); }
};
/*************************************************************
 *	EXTERNAL DATA DECLARATIONS
 */
GSL_VAR unsigned long int gsl_rng_default_seed;
/*************************************************************
 *	EXTERNAL FUNCTION PROTOTYPES
 */
#endif // E_erand_h
