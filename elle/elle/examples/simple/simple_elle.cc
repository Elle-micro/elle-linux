#include "simple_elle.h"


int InitElle(int argc, char **argv, char *rootname, int ProcessType)
{
    int err=0;

    /*
     * initialise
     */
    ElleInit();

    /*
     * set the function to the one in your process file
     */
    if(ProcessType==SIMPLE_FLYNN)
		ElleSetInitFunction(InitProcessFlynn);
	else if(ProcessType==SIMPLE_BNODE)
		ElleSetInitFunction(InitProcessBnode);
	else 
		ElleSetInitFunction(InitProcessUnode);
	
    if (err=ParseOptions(argc,argv))
        OnError("",err);

    /*
     * set the interval for writing to the stats file
    ES_SetstatsInterval(100);
     */

    /*
     * set the base for naming statistics and elle files
     */
    ElleSetSaveFileRoot(rootname);

    /*
     * set up the X window
     */
    if (ElleDisplay()) SetupApp(argc,argv);

    /*
     * run your initialisation function and start the application
     */
    StartApp();

    CleanUp();

    return(0);



}

/*
 * this function will be run when the application starts,
 * when an elle file is opened or
 * if the user chooses the "Rerun" option
 */

int InitProcessFlynn()
{
    char *infile;
    int err=0;
    
    /*
     * clear the data structures
     */
    ElleReinit();

    ElleSetRunFunction(ProcessFlynns);

    /*
     * read the data
     */
    infile = ElleFile();
    if (strlen(infile)>0) 
        if (err=ElleReadData(infile))
			OnError(infile,err);
    
}

int InitProcessBnode()
{
    char *infile;
    int err=0;
    
    /*
     * clear the data structures
     */
    ElleReinit();

    ElleSetRunFunction(ProcessBnodes);

    /*
     * read the data
     */
    infile = ElleFile();
    if (strlen(infile)>0) 
        if (err=ElleReadData(infile))
			OnError(infile,err);
    
}

int InitProcessUnode()
{
    char *infile;
    int err=0;
    
    /*
     * clear the data structures
     */
    ElleReinit();

    ElleSetRunFunction(SetUpUnodesByFlynn);

    /*
     * read the data
     */
    infile = ElleFile();
    if (strlen(infile)>0) 
        if (err=ElleReadData(infile))
			OnError(infile,err);
    
}

int ProcessFlynns()
{
    int i, k, l;
    int err=0,max;
    
    ElleCheckFiles();

    for (i=0;i<EllemaxStages();i++) // loop for number of stages
	{
        max = ElleMaxFlynns();		// index of maximum flynn used in model
        for (k=0;k<max;k++) 		// loop though all flynns
		{
            if (ElleFlynnIsActive(k)) // process flynn if it is active
			{
                err = MyFlynnProcess(k);  // as it suggests
            }
        }
        /*
         * update the count and redisplay
         * check whether to write an elle file
         */
        ElleUpdate();  
    }
    return(err);
} 

int ProcessBnodes()
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
                    err = MyBnodeProcess(j);
                    // check topology if new node position
                }
            }
        }
        else {
            for (j=max-1;j>=0;j--) {
                if (ElleNodeIsActive(j)) {
                    err = MyBnodeProcess(j);
                    // check topology if new node position
                }
            }
        }
        ElleUpdate();
    }
}

int SetUpUnodesByFlynn()
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
                ProcessUnodesByFlynn(j);  // process this triangulation
            }
        }
        ElleUpdate();
    }
} 

void ProcessUnodesByFlynn(int flynnid)
{
    int i, j;
    int id, num_nbs, count;
    double val, total;


    std::vector<int> unodelist;  // create a vector list of unodes
    ElleGetFlynnUnodeList(flynnid,unodelist); // get the list of unodes for a flynn
    count = unodelist.size();
    for (i=0; i<count; i++) 
	{	
		MyUnodeProcess(unodelist[i], flynnid);	

    }
}

int MyFlynnProcess(){};
int MyBnodeProcess(){};
int MyUnodeProcess(){};
