#pragma once
#ifndef _CHAPEAU_H_
#define _CHAPEAU_H_
 
#include <stdio.h>
   
// Chapeau type definition:  This defines the pair potential as 
// a linear expansion of chapeau functions ("peaks").  It also
// contains the variables responsible for evolving the vector
// of coefficients \lambda; these are the single-particle sums
// and the global accumulators
typedef struct CHAPEAU {
  // Conjunto de funciones base, sets de chapeaus que
  // constituyen una base para alguna funcion reconstruida por elementos
  // finitos. Los coeficientes tambien se almacenana aca.

  // number of data (not timsteps) to acumulate before update
  int nupdate;
 
  // Dimension
  int dm; 

  // Discretizacion del dominio por dimension
  int * N;       // size in mult-dimensional real space 

  int m;       // m[0] size in 1-dimensional wrap space 

  int ku;      // number of uperdiagonals 
  int ldad;    // dimension of the packed matrix
                 
               
  int * periodic; 
  double * rmin;
  double * rmax;
  double * dr;
  double * idr;  

  double * r; // Vector para guardar la posicion actual
  double * f; // Vector para guardar el gradiente de la funcion

  // Para sacar los coeficientes
  FILE * ofp; 
  int outputFreq;
  int outputLevel;

  // Since restart file is open and closed in the subrroutine, might be better
  // to store the name and not the unit? 
  char filename[255]; 

  // Variables that accumulate partial statistics to optimze FEP coeficients.
  // This is used to send between replicas and is a private copy of the
  // information acquired for the self sampling. After all the replicas
  // comuncates, and if this replica is not the center replica that add all the
  // other contributions, this variables are set to cero. If this is the
  // central replica, or non replica scheme is used, this variables contains
  // the full statistics of the sampling.
  double * b;
  double * lam;
  int * hits; 

  // This indicates the norm factor of A and b (needed to add chapeau or
  // accumulate new data). A, Afull and b are normalized to avoid long
  // numbers.... not sure if that is true. (We need something that lead
  // to a FES that does not have a huge additive constant.... ).
  double norm; 

  //matrix A will have dimensions (dm-1)*3*(ch->N[0]-1)+1 x m
  double ** A;

  // If this replica is not the center replica that add all the other
  // contributions, this variables accumulates the full statistics (including
  // all other replicas) infromation needed to optimze FEP coeficients. It this
  // replica is the center replica, ot non replica scheme is used, this
  // variables are always empty.
  double * bfull;
  double ** Afull;
   
  //pointer to procedure
  int (*accumulate)(struct CHAPEAU * self, double bias);
                
} chapeau;
 
#endif
