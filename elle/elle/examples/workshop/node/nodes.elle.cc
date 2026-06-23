#include <stdio.h>
#include <math.h>
#include <iostream.h>
#include "attrib.h"
#include "nodes.h"
#include "update.h"
#include "error.h"
#include "runopts.h"
#include "file.h"
#include "init.h"
#include "interface.h"

int DoSomethingToNode(int node);
int ProcessFunction(), InitThisProcess();

/*
 * this function will be run when the application starts,
 * when an elle file is opened or
 * if the user chooses the "Rerun" option
 */
int InitThisProcess()
{
    char *infile;
    int err=0;

    /*
     * clear the data structures
     */
    ElleReinit();

    ElleSetRunFunction(ProcessFunction);

    /*
     * read the data
     */
    infile = ElleFile();
    if (strlen(infile)>0) {
        if (err=ElleReadData(infile)) OnError("",err);
        /*
         * check for necessary node attributes which
         * are not in the input file
         * if necessary call
        ElleAttributeNotInFile(infile,validAttr);
         */
    }
}

int ProcessFunction()
{
    int i, j, k;
    int err=0,max;

    ElleCheckFiles();

    for (i=0;i<EllemaxStages();i++) // loop for number of stages
	{
        max = ElleMaxNodes();	// index of maximum node used in model
        /*
         * alternate up and down the node array
         */
        if (i%2) {
            for (j=0;j<max;j++) {
                if (ElleNodeIsActive(j)) // process node if it is active
				{
                    err = DoSomethingToNode(j);
                    // check topology if new node position
                }
            }
        }
        else {
            for (j=max-1;j>=0;j--) {
                if (ElleNodeIsActive(j)) {
                    err = DoSomethingToNode(j);
                    // check topology if new node position
                }
            }
        }
        ElleUpdate();
    }
}

int DoSomethingToNode(int node)
{
    double total;
    int nbnodes[3];
	int i;
	int count;
	
    ElleNeighbourNodes(node,nbnodes); 	 // get list of neighbouring nodes

	for(i=0,total=0.0,count=0;i<3;i++) 	// loop through list of neighbouring nodes
	{	
		if (nbnodes[i]!=NO_NB) 			// if node exists
		{
			total+=ElleNodeAttribute(nbnodes[i],CONC_A);  // add neighbour node attribute to total
			count++;									  // count works out if there are 2 or 3 neighbours
		}
	}	
	
	total+=ElleNodeAttribute(node,CONC_A); 	// add node attribute to total
    total=total/(count+1);					// get average of all neighbours and central node
	
	ElleSetNodeAttribute(node,total,CONC_A); // set node attribute to average
}
