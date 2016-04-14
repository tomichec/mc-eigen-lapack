/** 
 *  Copyright 2015 Vadim Biktashev, Tomas Stary
 *  
 *  Free software under GNU GPLv3.
 *  See <http://www.gnu.org/licenses/>.
 */

/* voltage range  */
#define VMIN	-100.0
#define DV	0.01
#define VMAX	70.0

/* file with values of voltage */
#define FVM	"vm_INa_%s.dat"
/* file with eigen values */
#define FEVINA	"evals_INa_%s.dat"
/* file with left eigenvectors */
#define FLEVINA	"left_evecs_INa_%s.dat"
/* file with right eigenvectors */
#define FREVINA	"right_evecs_INa_%s.dat"

#define ERROR(msg){fprintf(stderr,msg);exit(1);}
#define CALLOC(p,a,b) if(0==(p=calloc(a,b)))ERROR("not enough memory\n")
#define MAKE_ARRAY(type, name, length) type * name; CALLOC(name, length, sizeof(type));
#define OPEN_FILE(fileid, name,suffix)     \
  FILE * fileid;	\
  sprintf(filename,name,suffix);					\
  if ( ( fileid = fopen(filename,"w")) == NULL){			\
    fprintf(stderr,"Error while openning the file: %s.\n",filename);	\
    exit(1);}				
#define WRITE_EVAL(fileid, dimension, eval_Re, eval_Im) {	\
    int		ii;						\
    for( ii = 0; ii < dimension; ii++ ) {	\
      if( eval_Im[ii] == (double)0.0 ) {	\
	fprintf(fileid, " %.10e", eval_Re[ii] );	\
      } else {					\
	ERROR("The imaginary part is not zero.\n");	\
	/* fprintf(fileid, "(%.10e,%10e)", eval_Re[ii], eval_Im[ii] ); */	\
      }						\
      /* separators */				\
      (ii < (dimension-1))? fprintf(fileid, "\t"): fprintf(fileid,"\n");		\
    }						\
  }
#define WRITE_EVEC(fileid, dimension, eval_Im, evec) {		\
    int ii, jj;								\
    for( jj = 0; jj < dimension; jj++ ) {				\
      ii      = 0;								\
      while( ii < dimension ) {						\
	if( eval_Im[ii] == (double)0.0 ) {				\
	  fprintf(fileid, "%.10e", evec[jj*dimension+ii] );		\
	  ii++;								\
	} else {							\
	  ERROR("The imaginary part is not zero.\n");	\
	  /* fprintf(fileid,  "(%.10e,%.10e)", evec[jj*dimension+ii], evec[jj*dimension+(ii+1)] ); */ \
	  /* fprintf(fileid,  "(%.10e,%.10e)", evec[jj*dimension+ii], -evec[jj*dimension+(ii+1)] ); */ \
	  /* ii += 2; */							\
	}								\
	/* separators */						\
	((jj*dimension+ii) < (dimension*dimension) )? fprintf(fileid, "\t"): fprintf(fileid,"\n"); \
      }									\
    }									\
  }

/* Enumerate the markov chain states */
enum
  {
    #define _(n,i) markov_##n,
    #include "ina_markov.h"
    #undef _
    DIM		    		/* total number of Markov variables */
  };

void
ina_trans_rates_matrix(double V, double *tr)
{
  /* Updates the transition rates matrix of INa Markov chain model
     published by Clancy, Rudy (2002) */
  /* input variables: V -- membrane voltage; tr -- pointer to the
     matrix */
  int i;
  for (i=0; i<DIM*DIM; i++)
    {
      /* reset entries */
      tr[i]=0; 
    }
  /* recompute the tr matrix for new value of V */
  #define _VFUN(name,expression) double name=expression;
  #define _RATE(from,to,direct,reverse)		\
    tr[markov_##to*DIM+markov_##from]=direct;	\
    tr[markov_##from*DIM+markov_##from]-=direct;	\
    tr[markov_##from*DIM+markov_##to]=reverse;	\
    tr[markov_##to*DIM+markov_##to]-=reverse;
  #include "ina_rates.h"
  #undef _VFUN
  #undef _RATE
}
