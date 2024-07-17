#include "grad.hpp"
#include "lids.hpp"
#include "allvec.hpp"
#include "tensor3.hpp"


void grad(struct Shell *shell, struct lids *caps, double *gf)
{
    
    int v1, v2, v3, v4;             // Indexes of the vertexes of the the triangle or the two edges and its oposed vertexes

    
    zeros_dvector(gf, shell->N*DIMV);
    
    //*******************************************************************//
    /////////               1. Pressure gradient               ////////////
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    
    double p_eff = ONESIXTH*shell->p;
    #pragma omp parallel for reduction(+:gf[:shell->N*DIMV]) private(v1,v2,v3)
    for(int t=0; t < shell->N_tri; t++) {
        
        double x1, x2, x3;              // x component of the position for every vertex in triangle
        double y1, y2, y3;              // y component of the position for every vertex in triangle
        double z1, z2, z3;              // z component of the position for every vertex in triangle
        
        v1 = shell->tri[t][0];
        v2 = shell->tri[t][1];
        v3 = shell->tri[t][2];
        
        x1 = shell->r[v1*DIMV];
        x2 = shell->r[v2*DIMV];
        x3 = shell->r[v3*DIMV];
        
        y1 = shell->r[v1*DIMV + 1];
        y2 = shell->r[v2*DIMV + 1];
        y3 = shell->r[v3*DIMV + 1];
        
        z1 = shell->r[v1*DIMV + 2];
        z2 = shell->r[v2*DIMV + 2];
        z3 = shell->r[v3*DIMV + 2];
        
        // Vertex 1
        gf[v1*DIMV]     += p_eff*(y2*z3 - y3*z2);
        gf[v1*DIMV + 1] += p_eff*(x3*z2 - x2*z3);
        gf[v1*DIMV + 2] += p_eff*(x2*y3 - y2*x3);
        
        // Vertex 2
        gf[v2*DIMV]     += p_eff*(y3*z1 - y1*z3);
        gf[v2*DIMV + 1] += p_eff*(z3*x1 - x3*z1);
        gf[v2*DIMV + 2] += p_eff*(x3*y1 - y3*x1);
        
        // Vertex 3
        gf[v3*DIMV]     += p_eff*(y1*z2 - y2*z1);
        gf[v3*DIMV + 1] += p_eff*(x2*z1 - x1*z2);
        gf[v3*DIMV + 2] += p_eff*(x1*y2 - x2*y1);
        
    }

    
    //--------------------------------------------------------------//
    
    
    //First, calculate distances between edges in the current configuration
    calcdr(shell->r, 
           shell->dx,
           shell->dy, 
           shell->dz, 
           shell->d, 
           shell->edges,
           shell->nedges);
    
    #pragma omp parallel for reduction(+:gf[:shell->N*DIMV]) private(v1,v2,v3,v4)
    for(int i=0; i<shell->nedges; i++) { 
        
        //*********************************************************************//
        /////////               2. Stretching gradient               ////////////
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
        v1 = shell->edges[i][0];
        v2 = shell->edges[i][1];

        double f, fx, fy, fz;
        double dd;
        dd = (shell->d[i] - shell->eql[i]);
        //const double A_CONST = 2.0;
        //f = (shell->k_spr_eff[i])*dd/(shell->d[i]) + 3*A_CONST*dd*dd/(shell->d[i]);
        f = (shell->k_spr_eff[i])*dd/(shell->d[i]);
        
        fx = f*shell->dx[i];
        fy = f*shell->dy[i];
        fz = f*shell->dz[i];
        
        gf[v1*DIMV]     += fx;
        gf[v1*DIMV + 1] += fy;
        gf[v1*DIMV + 2] += fz;
        
        gf[v2*DIMV]     -= fx;
        gf[v2*DIMV + 1] -= fy;
        gf[v2*DIMV + 2] -= fz;
        //-------------------------------------------------------------------//
        
        //*******************************************************************//
        /////////               3. Bending gradient               ////////////
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
        double r1[3], r2[3], r3[3], r4[3];              // Coordinates of the vectors of the edge and oposes vertexes  
        double r21[3], r31[3], r24[3], r34[3], r23[3];  // Vectors defining the edges of the triangles
        double nI[3], nJ[3];                            // Normal vectors for a couple of triangles sharing an edge
        double norm_nI, norm_nJ;                        // Norms of the normal vectors
        double gnIbnJ[3], gnJbnI[3];                    // Result of the operation (I - n*nT)*u/|n|  
        double f_bend[3];                               // Unormalized bending force
        
        v1 = shell->op_edges[i][0];
        v4 = shell->op_edges[i][1];
        v2 = shell->edges[i][0];
        v3 = shell->edges[i][1];
        
        dvec_1D_list3(shell->r, v1, r1);
        dvec_1D_list3(shell->r, v2, r2);
        dvec_1D_list3(shell->r, v3, r3);
        dvec_1D_list3(shell->r, v4, r4);
        
        dvdiff3(r3, r1, r31);
        dvdiff3(r2, r1, r21);
        dvdiff3(r2, r4, r24);
        dvdiff3(r3, r4, r34);
        dvdiff3(r2, r3, r23);
        
        // 3.2 Calculate normal vector (cross product and normalization to unity)
        dcross3(r31, r21, nI);
        dcross3(r24, r34, nJ);
        
        norm_nI = dnorm3(nI);
        norm_nJ = dnorm3(nJ);
        
        dvscal3(nI, norm_nI);
        dvscal3(nJ, norm_nJ);
        
        graduniv_byu(nI, norm_nI, nJ, gnIbnJ);
        graduniv_byu(nJ, norm_nJ, nI, gnJbnI);
        
        
        // 3.3 Gradient on each vertex
        // Vertex 1
        dcross3(r23, gnIbnJ, f_bend);
        gf[v1*DIMV]     += shell->lambda_eff[i]*f_bend[0];
        gf[v1*DIMV + 1] += shell->lambda_eff[i]*f_bend[1];
        gf[v1*DIMV + 2] += shell->lambda_eff[i]*f_bend[2];
        
        // Vertex 2
        dcross3(r31, gnIbnJ, f_bend); 
        gf[v2*DIMV]     += shell->lambda_eff[i]*f_bend[0];
        gf[v2*DIMV + 1] += shell->lambda_eff[i]*f_bend[1];
        gf[v2*DIMV + 2] += shell->lambda_eff[i]*f_bend[2];
        
        dcross3(r34, gnJbnI, f_bend); 
        gf[v2*DIMV]     -= shell->lambda_eff[i]*f_bend[0]; // (Minus sign because r34 insted of r43)
        gf[v2*DIMV + 1] -= shell->lambda_eff[i]*f_bend[1];
        gf[v2*DIMV + 2] -= shell->lambda_eff[i]*f_bend[2];
        
        
        // Vertex 3
        dcross3(r21, gnIbnJ, f_bend); 
        gf[v3*DIMV]     -= shell->lambda_eff[i]*f_bend[0]; // (Minus sign because r21 insted of r12)
        gf[v3*DIMV + 1] -= shell->lambda_eff[i]*f_bend[1];
        gf[v3*DIMV + 2] -= shell->lambda_eff[i]*f_bend[2];
        
        dcross3(r24, gnJbnI, f_bend); 
        gf[v3*DIMV]     += shell->lambda_eff[i]*f_bend[0];
        gf[v3*DIMV + 1] += shell->lambda_eff[i]*f_bend[1];
        gf[v3*DIMV + 2] += shell->lambda_eff[i]*f_bend[2];
        
        // Vertex 4
        dcross3(r23, gnJbnI, f_bend); 
        gf[v4*DIMV]     -= shell->lambda_eff[i]*f_bend[0]; // (Minus sign because r32 insted of r23)
        gf[v4*DIMV + 1] -= shell->lambda_eff[i]*f_bend[1];
        gf[v4*DIMV + 2] -= shell->lambda_eff[i]*f_bend[2];
    
    }
    
    
    //*******************************************************************//
    /////////          4. Full tube bending gradients          ////////////
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
   
    
    #pragma omp parallel for reduction(+:gf[:shell->N*DIMV]) private(v1,v2,v3) num_threads(2)
    for(int t=0; t < caps->n_tri_lid_l; t++) {
        double L;
        double f_tube_bend[3];                      // General bending gradient
        double r1[3], r2[3], r3[3];                  // Coordinates of the vectors of the edge and oposes vertexes  
        double r21[3],r31[3], r23[3];                // Vectors defining the edges of the triangles
        double nI[3];                                // Normal vectors for a couple of triangles sharing an edge
        double norm_nI;                              // Norms of the normal vectors
        double gnIbnJ[3];                            // Result of the operation (I - n*nT)*u/|n|  
        
        
        v1 = caps->tri_l[t][0];
        v2 = caps->tri_l[t][1];
        v3 = caps->tri_l[t][2]; 
         
        dvec_1D_list3(shell->r, v1, r1);
        dvec_1D_list3(shell->r, v2, r2);
        dvec_1D_list3(shell->r, v3, r3);
        
        dvdiff3(r2, r1, r21);
        dvdiff3(r2, r3, r23);
        dvdiff3(r3, r1, r31);
        
        dcross3(r21, r31, nI);
        norm_nI = dnorm3(nI);
        dvscal3(nI, norm_nI);
        graduniv_byu(nI, norm_nI, caps->nl0, gnIbnJ);
        
        L = -caps->K_BRAZIER;
    
        // Vertex 1
        dcross3(r23, gnIbnJ, f_tube_bend); 
        gf[v1*DIMV]     += L*f_tube_bend[0]; 
        gf[v1*DIMV + 1] += L*f_tube_bend[1];
        gf[v1*DIMV + 2] += L*f_tube_bend[2];
        
        // Vertex 2
        dcross3(r31, gnIbnJ, f_tube_bend); 
        gf[v2*DIMV]     += L*f_tube_bend[0]; 
        gf[v2*DIMV + 1] += L*f_tube_bend[1];
        gf[v2*DIMV + 2] += L*f_tube_bend[2];
        
        // Vertex 3
        dcross3(gnIbnJ, r21, f_tube_bend); 
        gf[v3*DIMV]     += L*f_tube_bend[0]; 
        gf[v3*DIMV + 1] += L*f_tube_bend[1];
        gf[v3*DIMV + 2] += L*f_tube_bend[2];
        
    }
    
    #pragma omp parallel for reduction(+:gf[:shell->N*DIMV]) private(v1,v2,v3) num_threads(2)
    for(int t=0; t < caps->n_tri_lid_r; t++) {
        double L;
        double f_tube_bend[3];                      // General bending gradient
        double r1[3], r2[3], r3[3];                  // Coordinates of the vectors of the edge and oposes vertexes  
        double r21[3],r31[3], r23[3];                // Vectors defining the edges of the triangles
        double nI[3];                                // Normal vectors for a couple of triangles sharing an edge
        double norm_nI;                              // Norms of the normal vectors
        double gnIbnJ[3];                            // Result of the operation (I - n*nT)*u/|n|  
        
        
        v1 = caps->tri_r[t][0];
        v2 = caps->tri_r[t][1];
        v3 = caps->tri_r[t][2]; 
         
        dvec_1D_list3(shell->r, v1, r1);
        dvec_1D_list3(shell->r, v2, r2);
        dvec_1D_list3(shell->r, v3, r3);
        
        dvdiff3(r2, r1, r21);
        dvdiff3(r2, r3, r23);
        dvdiff3(r3, r1, r31);
        
        
        dcross3(r21, r31, nI);
        norm_nI = dnorm3(nI);
        dvscal3(nI, norm_nI);
        graduniv_byu(nI, norm_nI, caps->nr0, gnIbnJ);
        
        L = -caps->K_BRAZIER;
        
        // Vertex 1
        dcross3(r23, gnIbnJ, f_tube_bend); 
        gf[v1*DIMV]     += L*f_tube_bend[0]; 
        gf[v1*DIMV + 1] += L*f_tube_bend[1];
        gf[v1*DIMV + 2] += L*f_tube_bend[2];
        
        // Vertex 2
        dcross3(r31, gnIbnJ, f_tube_bend); 
        gf[v2*DIMV]     += L*f_tube_bend[0]; 
        gf[v2*DIMV + 1] += L*f_tube_bend[1];
        gf[v2*DIMV + 2] += L*f_tube_bend[2];
        
        // Vertex 3
        dcross3(gnIbnJ, r21, f_tube_bend); 
        gf[v3*DIMV]     += L*f_tube_bend[0]; 
        gf[v3*DIMV + 1] += L*f_tube_bend[1];
        gf[v3*DIMV + 2] += L*f_tube_bend[2];
        
    }
    
     //*******************************************************************//
    /////////          5. Center of mass lids alignment         ////////////
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    
    double Rcm_r[3] = {0., 0., 0.},
           Rcm_l[3] = {0., 0., 0.};
    
    for(int i=0; i<shell->N; i++) {
        if(caps->r_r[i]==1) {
            Rcm_r[0] += shell->r[i*DIMV];
            Rcm_r[1] += shell->r[i*DIMV+1];
            Rcm_r[2] += shell->r[i*DIMV+2];
        }
        
        if(caps->r_l[i]==1) {
            Rcm_l[0] += shell->r[i*DIMV];
            Rcm_l[1] += shell->r[i*DIMV+1];
            Rcm_l[2] += shell->r[i*DIMV+2];
        }
    }
    
    Rcm_r[0] /= caps->nv_lid_r;
    Rcm_r[1] /= caps->nv_lid_r;
    Rcm_r[2] /= caps->nv_lid_r;
    
    Rcm_l[0] /= caps->nv_lid_l;
    Rcm_l[1] /= caps->nv_lid_l;
    Rcm_l[2] /= caps->nv_lid_l;
        
    double DRcm_x = Rcm_r[0] - Rcm_l[0],
           DRcm_z = Rcm_r[2] - Rcm_l[2];
    
    #pragma omp parallel for num_threads(2)
    for(int i=0; i<shell->N; i++) {
        
        if(caps->r_r[i]==1) {
            gf[i*DIMV]   += caps->K_BRAZIER_CM*(DRcm_x)/caps->nv_lid_r;
            gf[i*DIMV+2] += caps->K_BRAZIER_CM*(DRcm_z)/caps->nv_lid_r;
        }
        
        if(caps->r_l[i]==1) {
            gf[i*DIMV]   -= caps->K_BRAZIER_CM*(DRcm_x)/caps->nv_lid_l;
            gf[i*DIMV+2] -= caps->K_BRAZIER_CM*(DRcm_z)/caps->nv_lid_l;
        }
    }
    
    
    //-------------------------------------------------------------------//
    
} 


/*****************************************************************************************/
/*                                      calcdr                                           */
/* Calculates the distance components and the absolute values of the distances between   */
/* every edge of the mesh                                                                */
/*****************************************************************************************/
void calcdr(double *r, 
            double *dx,
            double *dy, 
            double *dz, 
            double *d,
            int **edges, 
            int nedges)
{
    #pragma omp parallel for
    for(int i=0; i<nedges; i++) {
        
        int v1, v2;
        v1 = edges[i][0];
        v2 = edges[i][1];
        
        dx[i] = r[v1*DIMV] - r[v2*DIMV];
        dy[i] = r[v1*DIMV + 1] - r[v2*DIMV + 1];
        dz[i] = r[v1*DIMV + 2] - r[v2*DIMV + 2];
        d[i]  = sqrt(dx[i]*dx[i]  + dy[i]*dy[i] + dz[i]*dz[i]);
    }
}



/*****************************************************************************************/
/*                                   graduniv_byu                                        */
/* This function performs the operation [(I- (v_n)*(v_n)^T)^T]*u / |v|, where v_n is the */
/* normalized vector v and u is other vector. This expresio arises after extracting the  */
/* gradient of a unitary vector and is necessary to calculate the gradient of the        */
/* bending energy                                                                        */
/*****************************************************************************************/
void graduniv_byu(double *v_uni, double normv, double *u, double *out_v)
{
    double x,y,z;
    double a,b,c;
    x = v_uni[0];
    y = v_uni[1];
    z = v_uni[2];
    
    a = u[0];
    b = u[1];
    c = u[2];

    out_v[0] = (1 - x*x)*a - x*y*b - x*z*c;
    out_v[1] = (1 - y*y)*b - x*y*a - y*z*c;
    out_v[2] = (1 - z*z)*c - x*z*a - y*z*b;
    
    out_v[0] /= normv;
    out_v[1] /= normv;
    out_v[2] /= normv;
}

