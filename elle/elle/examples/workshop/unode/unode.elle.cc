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
    
    max_stages = EllemaxStages(); // number of stages 
    max_flynns = ElleMaxFlynns(); // maximum flynn number used
    max_unodes = ElleMaxUnodes(); // maximum unode number used
    for (i=0;i<max_stages;i++)   // cycle through stages
	{
        for (j=0;j<max_flynns;j++) // cycle through flynns
		{
            if (ElleFlynnIsActive(j))  // process if flynn is active 
			{
                ElleClearTriAttributes();  // reset triangulation
                TriangulateUnodes(j,MeshData.tri); // triangulate one flynns unodes
                SetUnodeAttributeFromNbs(j,CONC_A);  // process this triangulation
            }
        }
        ElleUpdate();
    }
} 

void SetUnodeAttributeFromNbs(int flynnid,int attr_id)
{
    int i, j;
    int id, num_nbs, count;
    double val, total;

    vector<int> unodelist;  // create a vector list of unodes
    ElleGetFlynnUnodeList(flynnid,unodelist); // get the list of unodes for a flynn
    count = unodelist.size();
    for (i=0; i<count; i++) {	
				
		vector<int> nbnodes, bndflag;
        ElleGetTriPtNeighbours(unodelist[i],nbnodes,bndflag,0); //get the  list of neighbours for a unode
        num_nbs = nbnodes.size();

        for (total=0,j=0; j<num_nbs; j++)  // loop through all neighbouring unodes
		{
            ElleGetUnodeAttribute(nbnodes[j],attr_id,&val); // get a neighbouring unode attribute value
            total += val;
        }
		ElleGetUnodeAttribute(unodelist[i],attr_id,&val); // get old unode attribute value
		total += val;
		
        val = total/(num_nbs+1);
		ElleSetUnodeAttribute(unodelist[i],attr_id, val); // set new unode attribute value

    }
}
