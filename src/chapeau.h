extern chapeau * chapeau_alloc ( int dm, double * rmin, double * rmax, int * N, int * periodic );
extern int chapeau_init ( chapeau * ch, int dm, double * rmin, double * rmax, int * N, int * periodic);
extern void chapeau_free ( chapeau * ch );
extern int chapeau_comparesize ( chapeau * ch1,  chapeau * ch2);
extern int chapeau_comparegrid ( chapeau * ch1,  chapeau * ch2);

extern void chapeau_sum ( chapeau * ch1, chapeau * ch2 );
extern chapeau * chapeau_crop (chapeau * ch, double * rmin, double * rmax);


// Output system
extern void chapeau_setupoutput ( chapeau * ch,  char * outfile, char * restartfile, int outputFreq, int outputLevel );
extern void chapeau_output ( chapeau * ch, int timestep );

// Restart system
extern void chapeau_savestate ( chapeau * ch, char * filename );
extern chapeau * chapeau_allocloadstate ( char * filename );
extern void chapeau_loadstate ( chapeau * ch, char * filename );
extern void chapeau_loadlambda ( chapeau * ch, char * filename );

extern void chapeau_solve ( chapeau * ch );
extern void chapeau_solve_secure ( chapeau * ch );

extern double chapeau_evalf_1simplex ( chapeau * ch, double z );
extern double chapeau_evalf_2simplex ( chapeau * ch, double z1, double z2 );
extern double chapeau_evalfg_2simplex ( chapeau * ch, double z1, double z2, double * g);
extern char * chapeau_serialize ( chapeau * ch );
extern void chapeau_addserialized ( chapeau * ch, char * str );
extern void chapeau_setserialized ( chapeau *ch, char * str );

extern void chapeau_set_peaks ( chapeau * ch, char * filename );
//void chapeau_baselinehits ( chapeau * ch );
//void chapeau_setmref ( chapeau * ch, double z );


extern int chapeau_simetrize2D ( chapeau * ch);
extern int accumulate_1D ( chapeau * ch, double bias);
extern int accumulate_2D ( chapeau * ch, double bias );

extern int chapeau_caneval_f1 ( chapeau * ch, double z );
extern int chapeau_caneval_f2 ( chapeau * ch, double z1, double z2 );

