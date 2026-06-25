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
        max = ElleMaxNodes();       // index of maximum node used in model
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
        /*
         * update the count and redisplay
         * check whether to write an elle file
         */
        ElleUpdate();
    }
}

int DoSomethingToNode(int node)
{
    double total_node,total_neighbour,sep,del;
    int nbnodes[3];
	int i;
	int count,nb_count;
	Coords node_pos,neighbour_pos,rel_pos;
	
	ElleNodePosition(node, &node_pos); // get position of node
	
    ElleNeighbourNodes(node,nbnodes); // get list of neighbouring nodes

	for(i=0,nb_count=0;i<3;i++) // loop through list of neighbouring nodes	
		if (nbnodes[i]!=NO_NB) // if node exists
			nb_count++;
			
	for(i=0;i<3;i++) // loop through list of neighbouring nodes
	{	
		if (nbnodes[i]!=NO_NB) // if node exists
		{
			ElleRelPosition(&node_pos, nbnodes[i], &rel_pos, &sep); 		// get relative position of neighbour
			total_node=ElleNodeAttribute(node,CONC_A);						// get conc value of node
			total_neighbour=ElleNodeAttribute(nbnodes[i],CONC_A);			// get conc value of neighbour

			if(rel_pos.y >0)
			{
				del=rel_pos.y*total_node/(sep*2);							// get parcel of concentration to move
				ElleSetNodeAttribute(node,total_node-del,CONC_A);			// remove parcel from node
				ElleSetNodeAttribute(nbnodes[i],total_neighbour+del,CONC_A);// add parcel to neighbour
			}
			else
			{
				del=-rel_pos.y*total_neighbour/(sep*2);						// get parcel of concentration to move
				ElleSetNodeAttribute(node,total_node+del,CONC_A);			// add parcel to node
				ElleSetNodeAttribute(nbnodes[i],total_neighbour-del,CONC_A);// remove parcel from neighbour
			}
		}
	}	
	
}
