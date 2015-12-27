#include "chapeau.h"
// Este programa lee los .bsp.retart y saca un .bsp que corresponde al restart
// todo, que lea tambien el .bsp y saque un .bsp en el tiempo 

//#define _PARANOIA_ 1
//


int main ( int argc, char * argv[] ) {
  chapeau ** ch;
  int i;
  
  
  ch=(chapeau**)malloc((argc-1)*sizeof(chapeau*));
  
  for (i=1;i<argc;i++) {
    fprintf(stdout,"reading %s\n",argv[i]); 
    ch[i-1]=chapeau_allocloadstate(argv[i]);
  
    fprintf(stdout,"ch %i info:\n",i-1); 
    fprintf(stdout,"--- ch rmin %.5f\n",ch[i-1]->rmin); 
    fprintf(stdout,"--- ch rmax %.5f\n",ch[i-1]->rmax); 
    fprintf(stdout,"--- ch dr %.5f\n",ch[i-1]->dr); 
    fprintf(stdout,"--- ch idr %.5f\n",ch[i-1]->idr); 
    fprintf(stdout,"--- alpha  %.5f\n",ch[i-1]->alpha); 
 
  }
  
  // Output of the first chapeau before add in it
  chapeau_setupoutput(ch[0],"chaps.bsp","chaps",1,1);
  ch[0]->nupdate=1;
  chapeau_update_peaks(ch[0]);
  chapeau_output(ch[0],1);

  //Prepare the output of the rest to the same file
  for (i=1;i<argc-1;i++) {
    ch[i]->outputFreq =ch[0]->outputFreq;
    ch[i]->outputLevel=ch[0]->outputLevel;
    ch[i]->ofp=ch[0]->ofp;
    ch[i]->nupdate=1;
  }

  // start adding and output of each chapeau
  for (i=1;i<argc-1;i++) {
    chapeau_update_peaks(ch[i]);
    chapeau_output(ch[i],1);
    chapeau_sum(ch[0],ch[i]);
  }

  // Output of the sum
  chapeau_update_peaks(ch[0]);
  chapeau_output(ch[0],1);
                   
  fflush(ch[0]->ofp);
                     
  return 0;

}
