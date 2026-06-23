
/*
 *  main.cc
 */

#include <stdio.h>
#include <stdlib.h>

#include "error.h"
#include "parseopts.h"
#include "init.h"
#include "runopts.h"
#include "setup.h"

main(int argc, char **argv)
{
    int err=0;
    UserData udata;
    extern int InitExpand(void);

    ElleInit();

    ElleSetInitFunction(InitExpand);

    /*ElleSetDisplay(0);*/
    /*ElleSetSaveFrequency(10);*/
    ElleUserData(udata);
    /* initialise maxarea, speedup and flag (step controlled by area) */
    udata[0] = udata[1] = udata[2] = 1.0;
    ElleSetUserData(udata);

    if (argc>1)
        if (err=ParseOptions(argc,argv))
            OnError("",err);

    if (ElleDisplay()) SetupApp(argc,argv);

    ElleSetSaveFileRoot("expand");

    /*
     * run init and run functions
     */
    StartApp();

    return(0);
} 
