#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iomanip>
#include <algorithm>
#include "display.h"
#include "nodes.h"
#include "unodes.h"
#include "error.h"
#include "string_utils.h"
#include "interface.h"
#include "check.h"
#include "crossings.h"
#include "runopts.h"
#include "polygon.h"
#include "file.h"
#include "init.h"
#include "convert.h"
#include "general.h"

int MovingNode(int node, ERegion *rgn);
int MoveNodeOut(int node1,int expflynn,bool use_area);
int CalcMove(int node, int *nbs, Coords *incr,bool use_area);
int ExpandGrain(), InitExpand();

int Count;
bool delete_attrib=0;

int InitExpand()
{
    char *infile;
    int err=0;
    int max;
    extern char InFile[];
    extern int Count;

    ElleReinit();
    ElleSetRunFunction(ExpandGrain);

    infile = ElleFile();
    if (strlen(infile)>0) {
        if (err=ElleReadData(infile)) OnError("",err);
        /*
         * check for any necessary attributes
         */
        if (!ElleFlynnAttributeActive(EXPAND)) 
            ElleAttributeNotInFile(infile,EXPAND);
        ElleAddDoubles();
        if (ElleUnodesActive() &&
                !ElleUnodeAttributeActive(U_ATTRIB_A)) {
            ElleInitUnodeAttribute(U_ATTRIB_A);
            ElleSetDefaultUnodeAttribute(1.0,U_ATTRIB_A);
            delete_attrib=1;
        }
    }

    Count = 1;
}

#ifdef DEBUG
extern void display_row(int k);
#endif

double MaxVExp = 1.0e-4;
double VExp = 0.0;
int ExpandGrain()
{
    int i, j, k, n, exp, expflynn;
    int interval=0,err=0,max,maxf, update;
    int nbs[3];
    vector<int> seq;
    bool use_area=1, reset=0;
    double speedup;
    UserData udata;
    double maxarea,currarea,eps=1e-10;
    double mass = 0.0, reset_me=0;
    double nmass = 0.0;
    double totmass = 0.0;
    double tmp = 0.0;

/*FILE *fp = fopen("area","w");*/
    ElleUserData(udata);
    maxarea = udata[0]; speedup = udata[1];  
    if (udata[2]<eps && udata[2]>-eps) use_area = 0;

    reset = ElleNodeAttributeActive(N_ATTRIB_A);
    if (reset) {
        if (!ElleNodeAttributeActive(N_ATTRIB_B))
            ElleInitNodeAttribute(N_ATTRIB_B);
    }

        max = ElleMaxNodes();
        for (j=0;j<max;j++) seq.push_back(j);
    if (ElleDisplay() && Count==1) EllePlotRegions(Count);
    for (i=0;i<EllemaxStages();i++,Count++) {
        for (j=0;j<max;j++)
            if (reset && ElleNodeIsActive(j) &&
                    ElleNodeAttribute(j,N_ATTRIB_A)!=0.0 &&
                        ElleNodeAttribute(j,CONC_A)<1)
                            ElleSetNodeAttribute(j,1,CONC_A);
        if (!(i%10)) {
        nmass = ElleTotalNodeMass(CONC_A);
        mass = ElleTotalUnodeMass(CONC_A);
        totmass = mass+nmass;
        cout << i << "umass " << mass << "  nmass " << nmass << 
             "  tot " << totmass << endl;
        }
        maxf = ElleMaxFlynns();
        for (j=0,currarea=0.0;j<maxf;j++) {
            if (ElleFlynnIsActive(j) && !isParent(j)) {
                ElleGetFlynnIntAttribute(j,&exp,EXPAND);
                if (exp!=0) currarea += ElleRegionArea(j);
            }
        }
        /*
         * rate of increase in flynn area will decrease as area
         * approaches maxarea (0 when area of expanding grains is maxarea,
         *                     1 when area of expanding grains is 0)
         */
        /*MaxVExp = (double)ElleminNodeSep()*0.01*speedup**/
                     /*(pow(2.0,(1.0-currarea/maxarea)) - 1.0);*/
        VExp = MaxVExp*speedup*(pow(2.0,(1.0-currarea/maxarea)) - 1.0);

        /*if (i%2) {*/
            /*for (j=0;j<max;j++)*/
random_shuffle(seq.begin(),seq.end());
            for (n=0;n<max;n++) {
                j=seq[n];
                if (ElleNodeIsActive(j)) {
                    if (MovingNode(j,&expflynn))
                        MoveNodeOut(j,expflynn,use_area);
                    if (ElleNodeIsActive(j))
                      if (ElleNodeIsDouble(j)) ElleCheckDoubleJ(j);
                      else if (ElleNodeIsTriple(j)) ElleCheckTripleJ(j);
                }
#if XY
                    if (reset && ElleNodeIsActive(j)) {
                        if ( ElleNodeAttribute(j,N_ATTRIB_A)!=0.0) {
                            ElleSetNodeAttribute(j,
                                ElleNodeAttribute(j,CONC_A),N_ATTRIB_B);
                            ElleSetNodeAttribute(j,1,CONC_A);
                        }
                        ElleNeighbourNodes(j,nbs);
                        for (k=0;k<3;k++)
                            if ( nbs[k]!=NO_NB &&
                               ElleNodeAttribute(nbs[k],N_ATTRIB_A)!=0.0) {
                                ElleSetNodeAttribute(nbs[k],
                                        ElleNodeAttribute(nbs[k],CONC_A),
                                                       N_ATTRIB_B);
                                ElleSetNodeAttribute(nbs[k],1,CONC_A);
                            }
                    }
#endif
                }
        /*}
        else {
            for (j=max-1;j>=0;j--)
                if (ElleNodeIsActive(j)) {
                    if (MovingNode(j,&expflynn))
                        MoveNodeOut(j,expflynn,use_area);
                    if (ElleNodeIsActive(j))
                      if (ElleNodeIsDouble(j)) ElleCheckDoubleJ(j);
                      else if (ElleNodeIsTriple(j)) ElleCheckTripleJ(j);
                }
        } */
#ifdef DEBUG
display_row(50);
#endif
        if (ElleDisplay()) {
            EllePlotRegions(Count);
            ElleShowStages( Count );
        }
        if ((interval=EllesaveInterval())>0) {
            if (Count>0 && Count%interval==0) {
                if (delete_attrib)
                    ElleRemoveDefaultUnodeAttribute(U_ATTRIB_A);
                if (err=ElleAutoWriteFile(Count))
                    OnError("",err);
            }
        }
        if (Count%20==0) {
            maxf = ElleMaxFlynns();
            for (j=0;j<maxf;j++) {
                if (ElleFlynnIsActive(j) && !isParent(j)) {
                    ElleGetFlynnIntAttribute(j,&exp,EXPAND);
                    /*if (exp!=0) printf("%lf ",ElleRegionArea(j));*/
                }
            }
            /*printf("\n");*/
        }
        if (!(i%10)) {
        nmass = ElleTotalNodeMass(CONC_A);
        mass = ElleTotalUnodeMass(CONC_A);
        totmass = mass+nmass;
        cout << i << "umass " << mass << "  nmass " << nmass << 
             "  tot " << totmass << endl;
        }
    }
}

/*
 * Allow node to move if it is on the boundary of
 * an expanding grain but not on the boundary
 * between two expanding grains
 */
int MovingNode(int node, ERegion *rgn)
{
    int i,nbnodes[3],move=0,cnt=0,val=0;
    ERegion tmp;

    ElleNeighbourNodes(node,nbnodes);
    for (i=0;i<3;i++) {
        if (nbnodes[i]!=NO_NB) {
            ElleNeighbourRegion(node,nbnodes[i],&tmp);
            ElleGetFlynnIntAttribute(tmp,&val,EXPAND);
            if (val) {
                if (!move) *rgn = tmp;
                move++;
            }
            cnt++;
        }
    }
    if (move==cnt) move=0;
    return(move);
}

int MoveNodeOut(int node1,int expflynn, bool use_area)
{
    int i, nghbr[2], nbnodes[3], rgns[3], tmp, err;
    int  cnt=0, shrinking=NO_NB, common_nb,  val;
    Coords movedist;

    /*
     * find the node numbers of the neighbours
     * and put them in the order, walking around expanding flynn,
     * nghbr[0], node1, nghbr[1]
     */
    if (err=ElleNeighbourNodes(node1,nbnodes))
        OnError("MoveNodeOut",err);
    ElleRegions(node1,rgns);
    /*
     * if it is a triple with 2 expanding regions, calculate
     * vector for each region and use average (allows for
     * expanding grains with different areas)
     */
    for (i=0,cnt=0;i<3;i++) {
        if (rgns[i]!=NO_NB) {
            ElleGetFlynnIntAttribute(rgns[i],&val,EXPAND);
            if (val) cnt++;
            else shrinking=rgns[i];
        }
    }
    if (cnt>1) {
        if ((nghbr[0] = ElleFindBndNb(node1,shrinking))==NO_NB)
            OnError("MoveNodeOut",NONB_ERR);
        i=0;
        while (i<3 && (nbnodes[i]==NO_NB || nbnodes[i]==nghbr[0])) i++;
       if (ElleNodeOnRgnBnd(nbnodes[i],shrinking))  {
            nghbr[1] = nbnodes[i];
        }
        else {
            i++;
            while (i<3 && (nbnodes[i]==NO_NB || nbnodes[i]==nghbr[0])) i++;
            if (ElleNodeOnRgnBnd(nbnodes[i],shrinking) && i<3)
                nghbr[1] = nbnodes[i];
        }
        if (i==3) OnError("MoveNodeOut - bnd nb not found",0);
    }
    else {
        if ((nghbr[1] = ElleFindBndNb(node1,expflynn))==NO_NB)
            OnError("MoveNodeOut",NONB_ERR);
        i=0;
        while (i<3 && (nbnodes[i]==NO_NB || nbnodes[i]==nghbr[1])) i++;
        if (ElleNodeOnRgnBnd(nbnodes[i],expflynn))  {
            nghbr[0] = nbnodes[i];
        }
        else {
            i++;
            while (i<3 && (nbnodes[i]==NO_NB || nbnodes[i]==nghbr[1])) i++;
            if (ElleNodeOnRgnBnd(nbnodes[i],expflynn) && i<3)
                nghbr[0] = nbnodes[i];
        }
        if (i==3) OnError("MoveNodeOut - bnd nb not found",0);
    }

    CalcMove(node1,nghbr,&movedist,use_area);
    ElleCrossingsCheck(node1,&movedist);
}

int CalcMove(int node, int *nbs, Coords *incr, bool use_area)
{
    int err=0, i, j, rgns[2];
    double a, factor, factor2, step;
    double ang, eps=1e-6;
    Coords xy[3];
    double dir_x,dir_y,newxy[2];
    double xpts[3], ypts[3], area_old, area_new;
    float dist, dist2;

    ElleNodePosition(node,&xy[1]);
    ElleNodePlotXY(nbs[0],&xy[0],&xy[1]);
    ElleNodePlotXY(nbs[1],&xy[2],&xy[1]);

    for (i=0;i<3;i++) {
        xpts[i] = xy[i].x;
        ypts[i] = xy[i].y;
    }
    area_old = polyArea(xpts,ypts,3);

    angle(xpts[1],ypts[1],
          xpts[2],ypts[2],
          xpts[0],ypts[0],&ang);

    /*
     * step is proportional to sqrt flynn area
     */
    ElleNeighbourRegion(node,nbs[1],&rgns[1]);
    step = VExp;
    if (use_area) {
        factor = sqrt(fabs(ElleRegionArea(rgns[1])));
        ElleNeighbourRegion(nbs[0],node,&rgns[0]);
        if (rgns[0]!=rgns[1]) {
            /*  
             *  triple - with 2 expanding grains so
             *    use average from nb regions
             */
            factor2 = sqrt(fabs(ElleRegionArea(rgns[0])));
            factor = (factor+factor2)*0.5;
        }
        step = VExp*(1.0-factor);
    }
    /*
     * bisect the angle
     */
    a = ang*0.5;
    if (ang<0) a += PI;
    rotate_coords(xpts[2],ypts[2],
                  xpts[1],ypts[1],
                  &dir_x, &dir_y, a);
    if (dir_y==ypts[1] && dir_x==xpts[1]) a=0;
    else 
        a = atan2(dir_y-ypts[1],dir_x-xpts[1]);
    newxy[0] = eps*cos(a) + xpts[1];
    newxy[1] = eps*sin(a) + ypts[1];
    xpts[1] = newxy[0]; ypts[1] = newxy[1];
    area_new = polyArea(xpts,ypts,3);

    /*
     * local area change should increase flynn area
     */
    if (area_old > area_new) {
        incr->x = (-step*cos(a));
        incr->y = (-step*sin(a));
    }
    else {
        incr->x = (step*cos(a));
        incr->y = (step*sin(a));
    }

    return(err);
}

#if XY
int WriteCircleNodeFile(float centrex,float centrey,float rad,
                      char *filename)
{
    char label[20];
    float max,incr,ang;
    int j, err=0, max_node, nn, num[360];
    int nn1,nn2,indx1,indx2;
    Coords current;
    FILE *fp;

    max = M_PI * 2.0;
    if ((fp=fopen(filename,"w"))==NULL) return(OPEN_ERR);
    incr = M_PI/90.0*0.45/rad;  /* 2 deg if rad=0.45 */
    max_node = (int)(max/incr);
    for (j=0;j<360;j++) num[j] = -1;
    RandomiseSequence(num,max_node);
    
    ang = 0.0;
    indx1 = 0;
    nn = num[indx1++];
    if (!id_match(FileKeys,LOCATION,label)) return(KEY_ERR);
    fprintf(fp,"%s\n",label);
    do {
        current.x = centrex + rad*cos(ang);
        current.y = centrey + rad*sin(ang);
        fprintf(fp,"%d %f %f\n",nn,current.x,current.y);
        ang += incr;
        nn = num[indx1++];
    } while( nn != -1 );

    fclose(fp);
    return( err );
}

int WriteCircleFile(char *nodefname, char *grainfname)
{
    char label[20];
    int nn, num, grain=2, nnarray[300], i=0, j;
    float x,y;
    FILE *fin, *fout;

    fin = fopen(nodefname,"r");
    fout = fopen(grainfname,"w");
    if ((num = fscanf(fin,"%19s\n", label))!=1) return(READ_ERR);
    if (!id_match(FileKeys,FLYNNS,label)) return(KEY_ERR);
    fprintf(fout,"%s\n",label);
    while (!feof(fin)) {
        if ((num = fscanf(fin,"%d %f %f\n",
                            &nn, &x, &y)) != 3) {
            printf("Read error for node file\n");
            exit(1);;
        }
        nnarray[i++] = nn;
    }
    fprintf(fout,"%d",grain);
    fprintf(fout," %d",i);
    for (j=0;j<i;j++) fprintf(fout," %d",nnarray[j]);
    fprintf(fout,"\n");
    grain = 1;
    fprintf(fout,"%d",grain);
    fprintf(fout," %d",i);
    for (j=i-1;j>=0;j--) fprintf(fout," %d",nnarray[j]);
    fprintf(fout,"\n");

    rewind(fin);
    while ((i=fgetc(fin))!=EOF) fputc(i,fout);

    fclose(fin);
    fclose(fout);
}
#endif
