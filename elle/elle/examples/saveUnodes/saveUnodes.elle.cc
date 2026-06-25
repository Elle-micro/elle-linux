#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nodes.h"
#include "display.h"
#include "error.h"
#include "runopts.h"
#include "file.h"
#include "interface.h"
#include "stats.h"
#include "init.h"
#include "unodes.h"
#include "general.h"

int saveUnodes();
int InitsaveUnodes();
void display_row(int k);

int InitsaveUnodes()
{
    int err=0;
    int max;
    char *infile;
    extern int Count;

    ElleReinit();
    ElleSetRunFunction(saveUnodes);

    infile = ElleFile();
    if (strlen(infile)>0) 
    {
        if (err=ElleReadData(infile)) OnError(infile,err);        
    
        if (!ElleUnodesActive())
            OnError("No unodes in file",0);

    }
    else /* warn that there is no file open */;
    

}

int saveUnodes()
{
    int i, j, row=0, col=0, max;
    int colcnt, rowcnt;
    int write_cx = 1, write_ry=1;
    FILE *outr,*outc,*outr_x,*outc_y;
    Coords xy;
    double col_x, row_y;
    double conc;    
    UserData udata;

    ElleUserData(udata);
    if (udata[0]<1) row=ElleFindUnodeRow(udata[0]);
    else row = (int) udata[0];
    if (udata[1]<1) col=ElleFindUnodeColumn(udata[1]);
    else col = (int) udata[1];
    
    outr=fopen("row","a");
    outc=fopen("col","a");
    outr_x=fopen("row_x","a");
    outc_y=fopen("col_y","a");
    
    fprintf(outr,"%s\trow= %d\t",ElleFile(),row);
    fprintf(outc,"%s\tcol= %d\t",ElleFile(),col);
    fprintf(outr_x,"%s\trow= %d\t",ElleFile(),row);
    fprintf(outc_y,"%s\tcol= %d\t",ElleFile(),col);

    max = ElleMaxUnodes();
    ElleGetUnodePosition(0,&xy);
    col_x = xy.x;  row_y = xy.y;
    rowcnt = colcnt = 0;
    for(j=0;j<max;j++) {
        ElleGetUnodePosition(j,&xy);
        if (xy.y!=row_y) {
            row_y = xy.y;
            rowcnt++;
            colcnt=0;
            if (rowcnt%2) colcnt=1;
        }
        if (rowcnt==0) {
            if (colcnt==col) 
                col_x = xy.x;
            if (row==0) {
                ElleGetUnodeAttribute(j,CONC_A,&conc);
                fprintf(outr,"%e\t",conc);
                if (write_ry) {
                    fprintf(outr_x,"y= %e\t",xy.y);
                    write_ry=0;
                }
                fprintf(outr_x,"%e\t",xy.x);
            }
        }
        else if (rowcnt==row) {
            ElleGetUnodeAttribute(j,CONC_A,&conc);
            fprintf(outr,"%e\t",conc);
            if (write_ry) {
                fprintf(outr_x,"y= %e\t",xy.y);
                write_ry=0;
            }
            fprintf(outr_x,"%e\t",xy.x);
        }
        if (colcnt==col) {
            ElleGetUnodeAttribute(j,CONC_A,&conc);
            fprintf(outc,"%e\t",conc);
            if (write_cx) {
                fprintf(outc_y,"x= %e\t",xy.x);
                write_cx=0;
            }
            fprintf(outc_y,"%e\t",xy.y);
        }
        colcnt+=2;
    }
    fprintf(outr,"\n");
    fprintf(outc,"\n");
    fprintf(outr_x,"\n");
    fprintf(outc_y,"\n");
    
    fclose(outr);
    fclose(outc);
    fclose(outr_x);
    fclose(outc_y);

}
