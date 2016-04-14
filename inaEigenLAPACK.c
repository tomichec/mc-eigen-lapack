/** 
 *  Copyright 2015 Vadim Biktashev, Tomas Stary
 *  
 *  Free software under GNU GPLv3.
 *  See <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <lapacke.h>
#include <math.h>
#include "inaEigen.h"

/* dimension of the system */
#define LDA	DIM
#define LDVL	DIM
#define LDVR	DIM

int
main (int argc, const char * argv[] )
{
  /* initialize lapack variables */
  /* ld stands for leading dimension of an array -- for example it
     would be 20 for an array (20, 10) */
  const int	n = DIM, lda = LDA, ldvl = LDVL, ldvr = LDVR;
  int		lapack_exit_status;

  /****************************************/
  /* create array and allocate the memory */
  /* transition rates matrix of INa */
  MAKE_ARRAY(double, trans_rates_matrix, DIM*DIM);
  /* real part of eigenvalues */
  MAKE_ARRAY(double, eval_real, DIM);
  /* imaginary part of eigenvalues */
  MAKE_ARRAY(double, eval_imag, DIM);
  /* left eigenvectors */
  MAKE_ARRAY(double, evec_left, LDVL*DIM);
  /* right eigenvectors */
  MAKE_ARRAY(double, evec_right, LDVR*DIM);

  /**********************************************/
  /* create pointers and open the output files  */
  MAKE_ARRAY(char, filename, 64);          /* name of output file */
  OPEN_FILE(fvolt, FVM, "LAPACK");         /* file for voltage */	    
  OPEN_FILE(feval, FEVINA, "LAPACK");	   /* file for eigenvalues */       
  OPEN_FILE(fevec_left, FLEVINA, "LAPACK");  /* file for left eigenvectors */ 
  OPEN_FILE(fevec_right, FREVINA, "LAPACK"); /* file for right eigenvectors */
  free(filename);

  /* membrane potential in mV */
  double volt;
  /* number of voltage steps */
  const int volt_steps_N = (( VMAX - VMIN )/DV);	
  
  int	i;
  for (i = 0;i <= volt_steps_N;i++ )
    {/* Membrane potential loop */
      /* get membrane potential */
      volt = VMIN+i*DV;
      /* get transition rates matrix */
      ina_trans_rates_matrix(volt, trans_rates_matrix);
      /* calculate the eigenvalues and right and left eigenvectors */
      lapack_exit_status =
	LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'V', 'V', n, trans_rates_matrix,
		      lda, eval_real, eval_imag, evec_left, ldvl, evec_right, ldvr); 
      /* Check for convergence */
      if( lapack_exit_status != 0 ) ERROR( "The algorithm failed to compute eigenvalues.\n" );
      /* write results */
      fprintf(fvolt, "%.2f\n", volt);
      WRITE_EVAL(feval, n, eval_real, eval_imag);
      WRITE_EVEC(fevec_left, n, eval_imag, evec_left);
      WRITE_EVEC(fevec_right, n, eval_imag, evec_right);
    }				/* end of membrane potential loop */
  /* free memory */
  free(trans_rates_matrix);
  free(eval_imag);
  free(evec_left);
  free(evec_right);
  
  /* close files */
  fclose(fvolt);
  fclose(feval);
  fclose(fevec_left);
  fclose(fevec_right);

  return 0;
}
