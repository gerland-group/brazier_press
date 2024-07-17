#include "nlcg.hpp"
#include "lids.hpp"
#include "grad.hpp"
#include "allvec.hpp"


/*****************************************************************************************/
/*                                       nlcg                                            */
/* Minimise an objective function with a provided gradient following the non-linear      */
/* conjugate gradient of N particles and  initial coordinates provided in a vector form  */
/* x = [x1, y1, z1, x2, y2, z2, ..., xN, yN, zN].                                        */
/* This implementation is based on the description provided in: An Introduction to the   */
/* Non-linear conjugate gradient method  without the agonizing pain. Shewchuk J.R., 1994 */
/*****************************************************************************************/
void nlcg(struct Shell *shell, struct lids *caps, struct NLCG_params *minparams) 
{
    
    double L0_CHAR;                     /* Characteristic distance between bonds */

    double delta_d;                     /* Square modulus of the search direction before line search. Used for tolerance crit. */
    double delta_new, delta_o;          /* Square norm of the gradient before and after a conjugate gradient iteration */
    double delta_mid;                   /* Dot producto between the current configuration gradient and the former */
    
    double eta, eta_prv;                /* Relations for the line search */
    double alpha;                       /* Scaling factor for the line search minimisation */
    double beta;                        /* Scaling factor for the nlcg minimisation (Polak-Riviere)*/

    double *dcg;                        /* The search direction for the conjugate gradient */
    double *r;                          /* Gradient for the set of particles */ 
    double *rpi;                        /* Gradient for the set of particles of the previus interaction */ 

    const int NDIM = shell->N*DIMV;      /* Effective length of the vector */
    
    dcg = dvector(NDIM); 
    r   = dvector(NDIM); 
    rpi = dvector(NDIM); 
    
    shell->dx  = dvector(shell->nedges); 
    shell->dy  = dvector(shell->nedges); 
    shell->dz  = dvector(shell->nedges); 
    shell->d  = dvector(shell->nedges); 
     
     
     
    // 1. INTIAL CALCULATIONS 
    // 1. a) For the potential (calculate eq. distances of springs)
    calcdr(shell->r, shell->dx, shell->dy, shell->dz,  shell->d, shell->edges, shell->nedges);
    L0_CHAR = dmean(shell->eql, shell->nedges);

    
    // 1. b) For the NLCG
    grad(shell, caps, r);

    
    changesign(r, NDIM);
    memcpy(dcg, r, NDIM*sizeof(double));
    memcpy(rpi, r, NDIM*sizeof(double));
    
    delta_new =  ddotn(r, dcg, NDIM);

    //*****************************************************************************//
    /////////////////     2. NON LINEAR CONJUGATE GRADIENT     //////////////////////
    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
        
    int k = 0;                              // Controls the reset of conjugate directions
    int itcg = 0;                           // Iterations NLCG
    int its_meanchk=0;                      // Counts the subset of iteration between 0 and MEAN_CHECK-1 (escape condition)
    double sqgrad_vertex[MEAN_CHECK];       // Squared force at a given NLCG iteration (escape condition)

    while(itcg < minparams->MAX_ITS_NLCG_FACT*NDIM) 
    {
        
        // Inexact Line search (Secant method)
        int j = 0;
        
        delta_d = ddotn(dcg, dcg, NDIM);
        alpha  = -L0_CHAR*(minparams->SIGMA_0);
        movedir(shell->r, dcg, minparams->SIGMA_0, NDIM);
        grad(shell, caps, r);
        eta_prv =  ddotn(r, dcg, NDIM);
        movedir(shell->r, dcg, -minparams->SIGMA_0, NDIM); // Reset to the original position

        do {
            grad(shell, caps, r);
            eta = ddotn(r, dcg, NDIM);

            if(ddotn(r,r,NDIM) == 0.0 || (eta_prv - eta) == 0.0) 
                break;
            
            alpha *= eta/(eta_prv - eta); 
            movedir(shell->r, dcg, alpha, NDIM);
            eta_prv = eta;
            j++;
        } while(j < minparams->MAX_ITS_LINE_SEARCH && alpha*alpha*delta_d > minparams->TOLERANCE_CRIT_LINE_SEARCH);
        

        grad(shell, caps, r);
        changesign(r ,NDIM);
        
        
        // Polak-Riviere+ search direction (see Nocedal and Wright)
        delta_o = delta_new;
        delta_new =  ddotn(r, r, NDIM);    if( delta_new  == 0.0) break;
        
        delta_mid = ddotn(r, rpi, NDIM);
        memcpy(rpi, r, NDIM*sizeof(double));
        beta = (delta_new - delta_mid)/delta_o;
        k++;
        
        // Reset condition
        if(k == NDIM || beta <= 0 || abs(delta_mid)/delta_new >= 0.1) {
            beta = 0;
            k=0;
        }
       
        //Update search direction
        for(int i=0; i<NDIM; i++)
            dcg[i] = r[i] + beta*dcg[i];
       
        // Correct for negative orientation
        if(ddotn(dcg,r,NDIM)<=0)
            memcpy(dcg, r, NDIM*sizeof(double));
        
        
        
        // Break condition if the force does not change compared with previous steps
        sqgrad_vertex[its_meanchk] = mean_sqforce_vertex(r, shell->N);
        if(itcg>2*MEAN_CHECK) {
            if( dmean(sqgrad_vertex, MEAN_CHECK) < minparams->TOLERANCE_CRIT_NLCG && 
                sqgrad_vertex[its_meanchk] < minparams->TOLERANCE_CRIT_NLCG)
                    break;
        }
        its_meanchk++; 
        if(its_meanchk==MEAN_CHECK) { its_meanchk=0; }
            
       
        itcg++;
        
        
    }
    cout << " (# Iterations: " << itcg << ")" << endl << endl;
    
    
    /*-----------------------------------------------------------------------------*/
    
    
    /****** Free memory ******/
    
    /* Free vectors */   
    free_dvector(shell->dx);
    free_dvector(shell->dy);
    free_dvector(shell->dz);
    free_dvector(shell->d);
    
    
    free_dvector(r);
    free_dvector(rpi);
    free_dvector(dcg);
    
    
}




/*****************************************************************************************/
/*                                      movedir                                          */
/* Moves the vector v in the direction given by the vector d and scaled by the factor    */
/* alpha                                                                                 */ 
/*****************************************************************************************/
void movedir(double *v, double *d, double alpha, int N)
{
    for(int i=0; i<N; i++) 
        v[i] += alpha*d[i];
}



/*****************************************************************************************/
/*                                      ddotn                                            */
/* Calculates the dot product of a vector v1 with a vector v2 with dimensions N          */
/*****************************************************************************************/
double ddotn(double *v1, double *v2, int N)
{
    double inp = 0;
    #pragma omp parallel for reduction(+:inp) num_threads(2)
    for(int i=0; i<N; i++)
        inp += v1[i]*v2[i];
    return inp;
}


/*****************************************************************************************/
/*                                      changesing                                       */
/* Change the sign of the components in the vector v with size N.                        */
/*****************************************************************************************/
void changesign(double *v, int N)
{
    for(int i=0;i<N;i++)
        v[i] = -v[i];
}


/*****************************************************************************************/
/*                                    mean_sqforce_vertex                                */
/* Calculates an average square forces per vertex                                        */
/*****************************************************************************************/
double mean_sqforce_vertex(double *r, int N)
{
    double mean_fsq=0;    
    #pragma omp parallel for reduction(+:mean_fsq) num_threads(2)
    for(int i=0; i<N; i++){
        double fx,fy,fz;
        fx = r[i*DIMV];
        fy = r[i*DIMV+1];
        fz = r[i*DIMV+2];
        mean_fsq += fx*fx + fy*fy + fz*fz;
    }
    mean_fsq /= N;
    return mean_fsq;
}

/*****************************************************************************************/
/*                                      dmean                                            */
/* Calculate the mean of the input vector v elements from 0 to the element N-1           */
/*****************************************************************************************/
double dmean(double *v, int N)
{
    double m = 0;
    for(int i=0; i<N; i++)
        m += v[i];
    return m/N;
    
}

