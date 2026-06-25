#include <stdio.h>
#include <math.h>
#include <vector>
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

using namespace std;

int InitSetUnodes(), SetUnodes();
void SetUnodeAttributeFromNbs(int flynnid,int attr_id);

#define SPINS 6  // number of possible orientations for any given site
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
	time_t now;
	
	time(&now); 
	
	srand((long) now); // initialise random number generator from current time
    
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
				//break; // uncomment for debugging purposes to get first flynn working
            }
        }
        ElleUpdate();
    }
} 

void SetUnodeAttributeFromNbs(int flynnid,int attr_id)
{
    int i,ii,j,k,loops;
    int id, num_nbs, count;
    double val,prob,m;
    int histo[SPINS],maxneigh,maxval,nodeval,total1,total2;
	double transprobs[6];
	int energy;
	
	for(j=0;j<6;j++) // set up transition probabilities based on number of similar neighbours
	{
		transprobs[j]=0.5*exp((double) -j);
	}
	
    vector<int> unodelist; // define vector list of unodes
    ElleGetFlynnUnodeList(flynnid,unodelist); // get the list of unodes for a flynn
    count = unodelist.size();

	for(loops=0;loops<1;loops++)
	{
	    	for (i=0; i<count; i++) {	
			ii=((int)(rand()/(double)RAND_MAX)*count)%count; 						// randomly select a unode from within flynn
			ii=(int)((rand()/(double)RAND_MAX)*count)%count; 						// randomly select a unode from within flynn
			
			histo[0]=histo[1]=histo[2]=histo[3]=0;	 		// reset histogram of site attributes to 0	
			vector<int> nbnodes,bndflag;
        	ElleGetTriPtNeighbours(unodelist[ii],nbnodes,bndflag,0); //get the  list of neighbours for a unode
        	num_nbs = nbnodes.size();

			ElleGetUnodeAttribute(unodelist[ii],attr_id, &val); //  get the unodes attribute
    	   
			nodeval=(int)val; 								// integer value for attribute

        	for (j=0,total1=0; j<num_nbs; j++)  			// loop through all neighbouring unodes
			{
            	ElleGetUnodeAttribute(nbnodes[j],attr_id,&val); // get neighbouring unode attribute value
  				histo[(int) val] ++; 						// increment histogram
				if((int) val != nodeval)
					total1++;
        	}
			
			if(histo[nodeval]==num_nbs)						// dont't bother if node is surrounded by equals
				continue;
			do												// randomly select a neighbour
			{
				k=(int)((rand()/(double)RAND_MAX)*SPINS);
			}
			while(histo[k]==0);
            
			for (j=0,total2=0; j<num_nbs; j++)  			// loop through all neighbouring unodes
			{
            	ElleGetUnodeAttribute(nbnodes[j],attr_id,&val); // get neighbouring unode attribute value
				if((int) val != k)							// increment histogram
					total2++;
        	}

			energy=(total1-total2);							// potential energy of change of node state
			if(energy >= 0)									
				prob=1.0;									// probability of change of node state
			else
				prob=transprobs[-energy];					// probability of change of node state
				
			m=rand()/(double)RAND_MAX;									// get random number
			
			if(m<prob)
				ElleSetUnodeAttribute(unodelist[ii],attr_id, (double) k ); // set new unode attribute value

    	}
	}
}
