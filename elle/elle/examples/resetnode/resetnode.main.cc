
/*
 *  main.cc
 */

#include <stdio.h>
#include <stdlib.h>

#include "error.h"
#include "parseopts.h"
#include "init.h"
#include "runopts.h"
#include "file.h"
#include "setup.h"
#include "reset.h"

main(int argc, char **argv)
{
    int err=0;
    UserData udata;
    extern int InitReset(void);

    /*
     * initialise
     */
    ElleInit();

    /*
     * set the function to the one in your process file
     */
    ElleSetInitFunction(InitReset);

    /* initialise node id, attribute value and default attribute id */
    ElleUserData(udata);
    udata[NodeId] = (double)NO_NB;
    udata[AttribValue] = 0.0;
    udata[AttribId] = CONC_A;
    ElleSetUserData(udata);

    if (err=ParseOptions(argc,argv))
        OnError("",err);

    /*
     * this process only runs once
     */
    ElleSetStages(1);

    /*
     * set up the X window
     */
    if (ElleDisplay()) SetupApp(argc,argv);

    /*
     * set the base for naming statistics and elle files
     */
    ElleSetSaveFileRoot("resetnode");

    /*
     * run your initialisation function and start the application
     */
    StartApp();

    return(0);
} 
