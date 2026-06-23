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
#include "update.h"
#include "reset.h"

static int Count; /* how many iterations have been completed */
                  /* in the ProcessFunction (relate to time) */

int InitReset(), ResetBndAttribute();

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

    ElleSetRunFunction(ResetBndAttribute);

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

int ResetBndAttribute()
{
    int i, k=0, j, m;
    int nb, start, curr, prev, rgns[3];
    int reset_cnt=0;
    double val=0;
    int interval=0,st_interval=0,err=0,max;
    UserData udata;
Coords xy;
    
    interval = EllesaveInterval();
   
    st_interval = ES_statsInterval();
    if (st_interval>0) ES_WriteStatistics(0);
   
    if (ElleDisplay()) EllePlotRegions(Count);
    
    ElleUserData(udata);

    for (i=0;i<EllemaxStages();i++) {
    /*
     * reset the value of the node attribute
     */
        start = (int)udata[NodeId];
/*ElleNodePosition(start,&xy);*/
/*cout << xy.y << ' ' <<endl;*/
        if (ElleNodeIsActive(start) && ElleNodeIsDouble(start)) {
            ElleRegions(start,rgns);
            ElleSetNodeAttribute(start,udata[AttribValue],
                                 (int)udata[AttribId]);
            while (k<3) {
                if (rgns[k]!=NO_NB) {
    /*vector<int> unodelist;*/
    /*ElleGetFlynnUnodeList(rgns[k],unodelist);*/
                  curr = prev = start;
                  reset_cnt=0;
                  while (reset_cnt<(int)(udata[ResetNum]) && 
                        (nb = ElleFindBndNode(curr,prev,rgns[k]))!=start
                          && ElleNodeIsDouble(nb)) {
                      ElleSetNodeAttribute(nb, udata[AttribValue],
                                       (int)udata[AttribId]);
/*ElleNodePosition(nb,&xy);*/
/*cout << xy.y << ' ' <<endl;*/
                      prev = curr;
                      curr = nb;
                      reset_cnt++;
                  }
                }
                k++;
            }
        }
        ElleUpdate();
    }
    if (!ElleDisplay()) CleanUp();
} 
