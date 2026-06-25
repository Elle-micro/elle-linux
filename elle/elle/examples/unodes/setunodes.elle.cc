#include <stdio.h>
#include <math.h>
#include "attrib.h"
#include "nodes.h"
#include "update.h"
#include "error.h"
#include "runopts.h"
#include "file.h"
#include "interface.h"
#include "init.h"
#include "triattrib.h"
#include "unodes.h"

int InitSetUnodes(), SetUnodes();
void SetUnodeAttributeFromNbs(int flynnid,int attr_id);

/*
 * this function will be run when the application starts,
 * when an elle file is opened or
 * if the user chooses the "Rerun" option
 */
int InitSetUnodes()
{
    char *infile;
    int err=0;
    
    /*
     * clear the data structures
     */
    ElleReinit();

    ElleSetRunFunction(SetUnodes);

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
}

int SetUnodes()
{
    int i, j, k;
    int max_stages, max_flynns, max_unodes;
    
    ElleCheckFiles();
    
    max_stages = EllemaxStages();
    max_flynns = ElleMaxFlynns();
    max_unodes = ElleMaxUnodes();
    for (i=0;i<max_stages;i++) {
        for (j=0;j<max_flynns;j++) {
            if (ElleFlynnIsActive(j)) {
                ElleClearTriAttributes();
                TriangulateUnodes(j,MeshData.tri);
                SetUnodeAttributeFromNbs(j,CONC_A);
            }
        }
        ElleUpdate();
    }
} 

void SetUnodeAttributeFromNbs(int flynnid,int attr_id)
{
    int i, j, count=0;
    int id, num_nbs;
    double val, total;

    vector<int> unodelist;
    ElleGetFlynnUnodeList(flynnid,unodelist);
    count = unodelist.size();
    for (i=0; i<count; i++) {
        vector<int> nbnodes;
        ElleGetTriPtNeighbours(unodelist[i],nbnodes,0); //internal pts
        num_nbs = nbnodes.size();
        for (total=0,j=0; j<num_nbs; j++)  {
            ElleGetUnodeAttribute(nbnodes[j],attr_id,&val);
            total += val;
        }
        if (num_nbs>0) {
            val = total/num_nbs;
            ElleSetUnodeAttribute(unodelist[i],attr_id,val);
        }
    }
}
