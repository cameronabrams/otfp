extern void cfacvBanner ( void );

extern tamdOptStruct * New_tamdOptStruct ( double g, double kt, double dt, int riftyp);

extern smdOptStruct * New_smdOptStruct ( double target, int t0, int t1);

extern restraint * New_restraint ( double k, double z, int nCV, double * cvc, char * rftypstr, double zmin, double zmax, char * boundstr, double boundk, char * outfile, int outputFreq);
    

// Other subroutines
extern FILE * my_fopen ( char * name, char * code ) ;
extern DataSpace * NewDataSpace ( int Nctr, int Ncv, int Nrstr, long int seed );
extern int DataSpace_SetupChapeau ( DataSpace * ds, int numrep, int dm, double * min,
    int * nKnots, double * max, int * periodic, int beginaccum, int beginsolve, 
    int useTAMDforces, char * outfile, int outfreq, int outlevel, int nupdate);
extern chapeau * DataSpace_get_chapeauadress ( DataSpace * ds, int i );
extern int DataSpace_getN ( DataSpace * ds );
extern double * DataSpace_centerPos ( DataSpace * ds, int i );
extern int DataSpace_AddAtomCenter ( DataSpace * ds, int n, int * ind, double * m );
extern cv * DataSpace_add_cv ( DataSpace * ds, char * typ, int nind, int * ind,
		      double zmin, double zmax,char * boundf, double boundk, char * outfile, int outputFreq );

extern restraint * DataSpace_AddRestr  ( DataSpace * ds, double k, double targ, int nCV, double * cvc, char * rftypstr, double zmin, double zmax,char * boundf, double boundk,char * outfile, int outputFreq);
extern int restr_UpdateTamdOpt ( restraint * r, double g, double kt, double dt );
extern int restr_AddTamdOpt ( restraint * r, double g, double kt, double dt, int chid  , int chdm );
extern int restr_AddSmdOpt  ( restraint * r, double target, int t0, int t1 );
extern void restr_output  ( restraint * r );
extern double DataSpace_ComputeCVs ( DataSpace * ds);
extern int DataSpace_RestrainingForces ( DataSpace * ds, int first, int timestep, double bias);
extern double DataSpace_RestraintEnergy ( DataSpace * ds );
extern void DataSpace_ReportAll ( DataSpace * ds );
extern void DataSpace_ReportCV ( DataSpace * ds, int * active, double * res );
extern int DataSpace_checkdata ( DataSpace * ds );
extern int DataSpace_dump ( DataSpace * ds ); 
extern FILE * my_binfopen ( char * name, char * code, unsigned int outputLevel, DataSpace * ds );
extern void DataSpace_BinaryReportRestraints ( DataSpace * ds, int step, int outputlevel, FILE * fp );


// Evolve Functions
extern int cbd ( restraint * r, double f );
extern int uniformvelocity ( restraint * r, double f );


// Boundaries Functions
extern int SoftUpperWall ( restraint * r );
extern int SoftLowerWall ( restraint * r );
extern int SoftWalls ( restraint * r );
extern int pbc ( restraint * r );
extern int Periodic ( restraint * r );
extern int nada ( restraint * r );

// Potential Function
extern int HarmonicCart ( restraint * r );
extern int HarmonicCart_cutoff ( restraint * r );
extern int HarmonicCart_pbc ( restraint * r );
extern int HarmonicCart_cutoff_pbc ( restraint * r );

extern double restr_getz ( restraint * r );
extern double restr_getu ( restraint * r );
extern int restr_set_rchid ( restraint * r, DataSpace * ds, int chid);

extern void ds_saverestrains ( DataSpace * ds, char * filename );
extern void ds_loadrestrains ( DataSpace * ds, char * filename );


