
#define P_ARGS(A) ()                                                           
#define REAL float

int sm_limits__ P_ARGS(( double, double, double, double )) { return 0; }
void sm_label__ P_ARGS(( char * )) {}
void sm_draw__ P_ARGS(( double, double )){}
void sm_relocate__ P_ARGS(( double, double )) {}
int sm_connect__ P_ARGS(( REAL [], REAL [], int )) { return 0; }
void sm_lweight__ P_ARGS(( double )){}
void sm_ltype__ P_ARGS(( int )){}
void sm_ylabel__ P_ARGS(( char * )){}
void sm_box__ P_ARGS(( int, int, int, int )){}
void sm_expand__ P_ARGS(( double )){}
void sm_ticksize__ P_ARGS(( double, double, double, double )){}
void sm_window__ P_ARGS(( int, int, int, int, int, int )){}
void sm_defvar__ P_ARGS(( char *, char * )){}
int sm_histogram__ P_ARGS(( REAL [], REAL [], int )){ return 0; }
void sm_points__ P_ARGS(( REAL [], REAL [], int )){}
void sm_ptype__ P_ARGS(( REAL *, int )){}
void sm_graphics__ P_ARGS(( void )){}
void sm_alpha__ P_ARGS(( void )){}
void sm_gflush__ P_ARGS (( void )){}
void sm_erase__ P_ARGS(( void )){}
void sm_curs__ P_ARGS(( REAL *, REAL *, int *)){}
int sm_device__ P_ARGS(( char * )){ return 0; }
void sm_hardcopy__ P_ARGS(( void )){}
void sm_putlabel__ P_ARGS(( int, char * )){}
void sm_angle__ P_ARGS(( double )){}
int sm_location__ P_ARGS(( int, int, int, int )){ return 0; }
void sm_xlabel__ P_ARGS(( char * )){}
void sm_grid__ P_ARGS(( int, int )){}
void sm_ctype__ P_ARGS(( char * )){}
