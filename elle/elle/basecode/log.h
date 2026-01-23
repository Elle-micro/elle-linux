/*****************************************************
 * Copyright: (c) L. A. Evans
 * File:      $RCSfile: log.h,v $
 * Revision:  $Revision: 1.2 $
 * Date:      $Date: 2012/08/16 06:42:34 $
 * Author:    $Author: levans $
 *
 ******************************************************/
/*!
        \file		log.h
        \brief		header for Log fns
        \par		Description:
                Log messages are written to a file
                which can be displayed if running a GUI
*/
#if !defined(_E_log_h)
#define _E_log_h
/*************************************************************
 *	INCLUDE FILES
 */
/*************************************************************
 *	CONSTANT DEFINITIONS
 */
#define E_LOGBUFSIZ 4096
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
extern char logbuf[E_LOGBUFSIZ];
/*************************************************************
 *	EXTERNAL FUNCTION PROTOTYPES
 */
#ifdef __cplusplus
extern "C" {
#endif
void Log(int loglevel, const char *msg);
#ifdef __cplusplus
}
#endif
#endif // _E_log_h
