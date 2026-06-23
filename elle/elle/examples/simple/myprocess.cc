#include "simple_elle.h"

/*!
 \brief	Generic elle process template system 
 \par	Description:
 		This code allows easy processing of flynns, bnodes & unodes, 
		where the process consists of code that handles each individual object in 
		turn, once for each timestep. For example to update a bnode attribute once for 
		each timestep. For the unodes it handles them by flynn, so at the moment
		acting on all unodes regardless of parent flynn in not possible.
		All the basic loops are found in simple_elle.cc but do not have to be modified.
		
		Usage: 
		
		1) Alter the third argument to InitElle from "trial" to change to root filename for output elle files
		
		2) Alter the fourth argument to InitElle from SIMPLE_UNODE to SIMPLE_FLYNN or SIMPLE_BNODE to define which of the
		   three elle object types will be processed
		   
		3) Alter the corresponding MyProcess functions below to create a new process (Leave the other two untouched or delete them, 
		   they will be ignored in any case).
		
		According to whether SIMPLE_UNODE, SIMPLE_FLYNN, or SIMPLE_BNODE is passed to InitElle(), 
		the examples below will set:
		
			the flynn viscosity equal to the flynn id OR
			the bnode concentration equal to the node id OR 
			the unode concentration equal to the flynn id
			
		This template works with the file simple.elle found in this directory. Note that
		at the moment the unode attribute U_CONC_A must also be set.
		
 */
 
main(int argc, char **argv)
{
	InitElle( argc, argv, "trial", SIMPLE_UNODE);
} 


int MyFlynnProcess(int flynnid)
{
	ElleSetFlynnRealAttribute(flynnid,(double) flynnid,F_ATTRIB_A); 	// change this to alter flynn process
}

int MyBnodeProcess(int bnodeid)
{
	ElleSetNodeAttribute(bnodeid,(double) bnodeid,N_ATTRIB_A); 		// change this to alter bnode process
}

int MyUnodeProcess(int unodeid, int flynnid)
{
	ElleSetUnodeAttribute(unodeid,(double) flynnid, U_ATTRIB_A); 		// change this to alter unode process
}
