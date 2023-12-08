#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "chapeau_obj.h"
#include "chapeau.h"
// Este programa lee los .bsp.retart y saca un .bsp que corresponde al restart
// todo, que lea tambien el .bsp y saque un .bsp en el tiempo 

//#define _PARANOIA_ 1
//


int main ( int argc, char * argv[] ) {
  int j;
  chapeau * ch;
  chapeau * cho;
  double rmin[2],rmax[2];
  
  fprintf(stdout,"reading %s\n",argv[1]); 
  ch=chapeau_allocloadstate(argv[1]);
  
  for (j=0;j<ch->dm;j++) {
    fprintf(stderr,"--- ch dm %d/%d\n",j,ch->dm); 
    fprintf(stderr,"--- ch N %d\n",ch->N[j]); 
    fprintf(stderr,"--- ch rmin %.5f\n",ch->rmin[j]); 
    fprintf(stderr,"--- ch rmax %.5f\n",ch->rmax[j]); 
    fprintf(stderr,"--- ch dr %.5f\n",ch->dr[j]); 
    fprintf(stderr,"--- ch idr %.5f\n",ch->idr[j]); 
  }
   
  rmin[0]=atof(argv[2]);
  rmax[0]=atof(argv[3]);
  if (ch->dm==2) {
    rmin[1]=atof(argv[4]);
    rmax[1]=atof(argv[5]);
  }

  cho=chapeau_crop(ch,rmin,rmax);
  fprintf(stderr,"-----------\n"); 
  
  // Output of the first chapeau before add in it
  chapeau_setupoutput(cho,"croped.bsp","croped",1,1);
 
  for (j=0;j<ch->dm;j++) {
    fprintf(stderr,"--- ch dm %d/%d\n",j,cho->dm); 
    fprintf(stderr,"--- ch N %d\n",cho->N[j]); 
    fprintf(stderr,"--- ch rmin %.5f\n",cho->rmin[j]); 
    fprintf(stderr,"--- ch rmax %.5f\n",cho->rmax[j]); 
    fprintf(stderr,"--- ch dr %.5f\n",cho->dr[j]); 
    fprintf(stderr,"--- ch idr %.5f\n",cho->idr[j]); 
  }

  cho->nupdate=1;
  chapeau_solve(cho);
  // chapeau_solve_secure(cho);
  chapeau_output(cho,1);
  chapeau_savestate (cho,"croped.ch");
             
  fflush(cho->ofp);
                     
  return 0;

}
