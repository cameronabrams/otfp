#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "chapeau_obj.h"
#include "chapeau.h"
#include "chapeau_old.h"

int main ( int argc, char * argv[] ) {
  chapeau * ch;
  old_chapeau * och;
  int i,j;
  int * periodic;

  och=old_chapeau_allocloadstate(argv[1]);

  // Allocating the object
  ch=(chapeau*)malloc(sizeof(chapeau));

  // Copy data
  ch->ofp=0;
  ch->dm=och->dm;
  ch->m=och->m;
  ch->periodic=calloc(ch->dm,sizeof(int));
  for (i=0;i<ch->dm;i++) {
    ch->periodic[i]=0;
  }
  if (ch->dm==1) {
    ch->ku=1;
  } else{
    ch->ku=och->N[0];
  }  
  ch->ldad=2*ch->ku+1;
   
  // Link to arrays
  ch->rmin=och->rmin;
  ch->rmax=och->rmax;
  ch->dr=och->dr;
  ch->idr=och->idr;
  ch->N=och->N;
  ch->r=och->r;
  ch->f=och->f;

  ch->lam=och->lam;
  ch->hits=och->hits;
  ch->A=och->A;
  ch->b=och->b;
  ch->Afull=och->Afull;
  ch->bfull=och->bfull;
 
  if (ch->dm==1) {
    ch->accumulate=accumulate_1D;
  } else {
    ch->accumulate=accumulate_2D;
  }
          
     
  // Output of the first chapeau before add in it
  chapeau_setupoutput(ch,"coverted.bsp","converted",1,1);

  //Prepare the output of the rest to the same file
  ch->outputFreq =och->outputFreq;
  ch->outputLevel=och->outputLevel;
  ch->ofp=ch->ofp;
  ch->nupdate=1;

  // Output of the sum
  chapeau_solve(ch);
  chapeau_output(ch,1);
  chapeau_savestate (ch,"converted.ch");
                   
  fflush(ch->ofp);
                     
  return 0;

}
