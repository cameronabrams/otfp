#include "chapeau.h"

//#define _PARANOIA_ 1

int main ( int argc, char * argv[] ) {
  chapeau ** ch;
  int i;
  
  
  ch=(chapeau**)malloc((argc-1)*sizeof(chapeau*));
  
  for (i=1;i<argc;i++) {
    fprintf(stdout,"reading %s\n",argv[i]); 
    ch[i-1]=chapeau_allocloadstate(argv[i]);
  }
  
  // Output first
  chapeau_setupoutput(ch[0],"chaps",1,1);
  ch[0]->updateinterval=1;
  chapeau_update_peaks(ch[0],1,1);
  chapeau_output(ch[0],1);

  // start adding
  for (i=2;i<argc;i++) {
    fprintf(stdout,"adding %s\n",argv[i]); 
    chapeau_sum(ch[0],ch[i-1]);
    chapeau_update_peaks(ch[0],1,1);
    chapeau_output(ch[0],1);
    fflush(ch[0]->ofp);
  }

}
