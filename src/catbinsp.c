/* catbinsp: reads binary-format data files generated by
   cfacv_otfp/chapeau_output(() */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main ( int argc, char * argv[] ) {
  FILE * fp=NULL;
  unsigned int fileOutputLevel;
  unsigned int desiredOutputLevel=1;
  int n,i,ln,timestep;
  int reading;
  double * data[4];
  int * intdata;
  char flag[4];

  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-f")) fp=fopen(argv[++i],"r");
    else if (!strcmp(argv[i],"-ol")) desiredOutputLevel=atoi(argv[++i]);
  }
  if (!fp) {
    fprintf(stderr,"ERROR: no input file specified.\n");
    exit(-1);
  }

  fread(flag,sizeof(char),4,fp);
  fprintf(stdout,"#flag is %c%c%c%c\n",flag[0],flag[1],flag[2],flag[3]);
  fread(&fileOutputLevel,sizeof(unsigned int),1,fp);
  fread(&n,sizeof(int),1,fp);
  fprintf(stderr,"#INFO: outputLevel %u nKnots %i\n",fileOutputLevel,n);
  if (fileOutputLevel!=desiredOutputLevel) {
    fprintf(stderr,"#WARNING: desired output level does not match file output level.  Intersecting.\n");
    desiredOutputLevel &= fileOutputLevel;
  }

  for (i=0;i<2;i++) data[i]=(double*)calloc(n,sizeof(double));
  intdata=(int*)calloc(n,sizeof(int));


  fread(&timestep,sizeof(int),1,fp);
  if (fileOutputLevel & 1) reading=fread(data[0],sizeof(double),n,fp);
  if (fileOutputLevel & 2) reading=fread(intdata,sizeof(int),n,fp);

  ln=0;
  while (reading) {

    if (desiredOutputLevel & 1) {
      fprintf(stdout,"KNOTS    %i ",timestep);
      for (i=0;i<n;i++) fprintf(stdout,"% 11.5lf",data[0][i]);
      fprintf(stdout,"\n");
    } 
    if (desiredOutputLevel & 2) {
      fprintf(stdout,"HITS     %i ",timestep);
      for (i=0;i<n;i++) fprintf(stdout,"% 1i",intdata[i]);
      fprintf(stdout,"\n");
    } 
    ln++;

    fread(&timestep,sizeof(int),1,fp);
    if (fileOutputLevel & 1) reading=fread(data[0],sizeof(double),n,fp);
    if (fileOutputLevel & 2) reading=fread(intdata,sizeof(int),n,fp);

  }
  fclose(fp);

  return 0;
}
