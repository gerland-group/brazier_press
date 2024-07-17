#ifndef _NLCG_
#define _NLCG_

#include <stdio.h>
#include <string.h> // form memcpy

#define     DIMV             3                                      /* Number of dimensions per particles (R^3) */
#define     MEAN_CHECK      5                                      /* Number of iterations the gradient should not have changed to scape */

                              
struct NLCG_params
{
    int MAX_ITS_NLCG_FACT;                                       /* Number of CG steps is a factor of the number of particles */          
    int MAX_ITS_LINE_SEARCH;                                     /* Number of line search iterations*/
    double TOLERANCE_CRIT_LINE_SEARCH;                           /* Error tolerance for the secant method */
    double TOLERANCE_CRIT_NLCG;                                  /* Error tolerance of the force in the CG minimisation */
    double SIGMA_0;                                              /* Factor of the line search (by the average rest length of the bonds)  */   
};



                              
// Function declarations
void nlcg(struct Shell *shell, struct lids *caps, struct NLCG_params *minparams);

void movedir(double *v, double *d, double alpha, int N);
double ddotn(double *v1, double *v2, int N);
void changesign(double *v, int N);
double mean_sqforce_vertex(double *r, int N);
double dmean(double *v, int N);


#endif
