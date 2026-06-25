#include <stdio.h>
#include <math.h>
#include "attrib.h"
#include "nodes.h"
#include "display.h"
#include "error.h"
#include "runopts.h"
#include "file.h"
#include "interface.h"
#include "stats.h"
#include "init.h"
#include "unodes.h"
#include "reset.h"

static int Count; /* how many iterations have been completed */
                  /* in the ProcessFunction (relate to time) */

int InitReset(), ResetNodeAttribute();

/*
 * this function will be run when the application starts,
 * when an elle file is opened or
 * if the user chooses the "Rerun" option
 */
int InitReset()
{
    char *infile;
    int err=0;
    
    /*
     * clear the data structures
     */
    ElleReinit();

    ElleSetRunFunction(ResetNodeAttribute);

    /*
     * read the data
     */
    infile = ElleFile();
    if (strlen(infile)>0) {
        if (err=ElleReadData(infile)) OnError(infile,err);
        /*
         * check for any necessary attributes which may
         * not have been in the elle file
         */
    }

    Count = 0;
}

int ResetNodeAttribute()
{
    int i, k, l;
    double val=0;
    int interval=0,st_interval=0,err=0,max;
    UserData udata;
    
    interval = EllesaveInterval();
   
    st_interval = ES_statsInterval();
    if (st_interval>0) ES_WriteStatistics(0);
   
    if (ElleDisplay()) EllePlotRegions(Count);
    
    ElleUserData(udata);

    for (i=0;i<EllemaxStages();i++) {
    /*
     * reset the value of the node attribute
     */
        if (ElleNodeIsActive((int)udata[NodeId]))
            ElleSetNodeAttribute((int)udata[NodeId],
                                      udata[AttribValue],
                                 (int)udata[AttribId]);

        Count++;
        /*
         * update the display
         */
        if (ElleDisplay()) {
            EllePlotRegions( Count );
            ElleShowStages( Count );
        }
        /*
         * check whether to write an elle file
         */
        if (interval>0) {
            if (Count>0 && Count%interval==0) {
                ElleAutoWriteFile(Count);
            }
        }
        /*
         * check whether to append the current statistics
         */
        if (st_interval>0) {
            if (Count%st_interval==0) {
                ES_WriteStatistics(Count);
            }
        }
    }
    if (!ElleDisplay()) CleanUp();
} 
