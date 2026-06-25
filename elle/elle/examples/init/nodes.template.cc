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

    for (i=0;i<EllemaxStages();i++) {
        max = ElleMaxNodes();
        /*
         * alternate up and down the node array
         */
        if (i%2) {
            for (j=0;j<max;j++) {
                if (ElleNodeIsActive(j)) {
                    err = DoSomethingToNode(j);
                    .
                    .
                }
            }
        }
        else {
            for (j=max-1;j>=0;j--) {
                if (ElleNodeIsActive(j)) {
                    err = DoSomethingToNode(j);
                    .
                    .
                }
            }
        }
        ElleUpdate();
    }
}

int DoSomethingToNode(int node)
{
    .
    .
}
