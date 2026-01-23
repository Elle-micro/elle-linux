/*****************************************************
 * Copyright: (c) Dr J.K. Becker
 * File:      $RCSfile: gz_utils.h,v $
 * Revision:  $Revision: 1.1 $
 * Date:      $Date: 2005/07/12 06:53:55 $
 * Author:    $Author: levans $
 *
 ******************************************************/
/*!
        \file		gz_utils.h
        \brief		header for functions used in r/w of zip files
        \par		Description:
                                function declarations
*/
#if !defined(_E_gz_utils_h)
#define _E_gz_utils_h
/*************************************************************
 *	INCLUDE FILES
 */
#include <zlib.h>
/* modernize: make header C-safe by guarding C++ includes (analysis fast-track,
 * behavior-neutral) */
#ifdef __cplusplus
#include <string> // modernize: C++ only
#endif
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
/*************************************************************
 *	EXTERNAL DATA DECLARATIONS
 */
/*************************************************************
 *	EXTERNAL FUNCTION PROTOTYPES
 */
#ifdef __cplusplus
// modernize: hide C++ API from C translation units (analysis fast-track,
// behavior-neutral)
std::string gzReadLineSTD(gzFile in);
std::string gzReadSingleString(gzFile in);
#endif
#endif // _E_gz_utils_h
