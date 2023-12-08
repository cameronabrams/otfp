#include <stdlib.h>
#include <string.h>
#include <math.h> 
#include "chapeau_obj.h"
#include "chapeau.h"
// Este programa lee un .bsp.restart para un chapeau 2D y simetriza
// la estadistica en linea y=x (sirve para cuando se ve la posicion de 2 cosas
// identicas que si se intercambian son lo mismo, nav 1.2 por ejemplo)

int main ( int argc, char * argv[] ) {
  chapeau* ch;
  int j;
  
  
  fprintf(stdout,"reading %s\n",argv[1]); 
  ch=chapeau_allocloadstate(argv[1]);
  
  for (j=0;j<ch->dm;j++) {
    fprintf(stderr,"--- ch dm %d/%d\n",j,ch->dm); 
    fprintf(stderr,"--- ch rmin %.5f\n",ch->rmin[j]); 
    fprintf(stderr,"--- ch rmax %.5f\n",ch->rmax[j]); 
    fprintf(stderr,"--- ch dr %.5f\n",ch->dr[j]); 
    fprintf(stderr,"--- ch idr %.5f\n",ch->idr[j]); 
  }

  chapeau_simetrize2D(ch);

  ch->nupdate=1;
  chapeau_solve(ch);
  chapeau_setupoutput(ch,"sim.bsp","sim",1,1);
  chapeau_output(ch,1);
  fflush(ch->ofp);
                     
  return 0;

}
