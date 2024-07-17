#include "lids.hpp"
#include "grad.hpp"
#include "meshes.hpp"
#include "tensor3.hpp"
#include "allvec.hpp"






/*****************************************************************************************/
/*                                   find_n_lids                                         */
/* Extract the total number of particles that conforms one lid NOT including the         */
/* terminal hoop of the tube                                                             */
/*****************************************************************************************/
int find_n_lids(double *r, int N)
{
    // Find the spacing and radius (simple quick-and-dirty)
    double ddx, ddy, ddz, sp, R, thrr;
    
    ddx = r[0] - r[3];
    ddy = r[1] - r[4];
    ddz = r[2] - r[5];
    sp = sqrt(ddx*ddx + ddy*ddy + ddz*ddz);
    
    R = sqrt(r[0]*r[0] + r[2]*r[2]);
    thrr = R - (sp*sqrt(3.)/2.0)*0.25;
  
    int N_lid_tot = 0;
    for(int i=0; i<N; i++) {
        double x = r[i*3+0]; 
        double z = r[i*3+2];
        double radius = sqrt(x*x + z*z);
        
        if(radius < thrr) {
            N_lid_tot++;
            //cout << i << endl;
        }
    }
    
    return N_lid_tot/2;
}



/*****************************************************************************************/
/*                                      find_tri_lids                                    */
/* Determine the triangules that collectively define the lids of the tube                */
/*****************************************************************************************/
void find_tri_lids(struct Shell *shell, struct lids *caps)
{
    std::vector<std::vector<int>> right_lid_tri;
    std::vector<std::vector<int>> left_lid_tri;
    
    
    caps->r_l = ivector(shell->nedges);
    caps->r_r = ivector(shell->nedges);
    zeros_ivector(caps->r_r, shell->nedges);
    zeros_ivector(caps->r_l, shell->nedges);
    
    double j_left[3]  = {0.0,-1.0, 0.0},
           j_rigth[3] = {0.0, 1.0, 0.0};
           
    caps->n_tri_lid_r = 0;
    caps->n_tri_lid_l = 0;
    for(int t = 0; t< shell->N_tri; t++) {
        
        int v1,v2,v3;
        v1 = shell->tri[t][0];
        v2 = shell->tri[t][1];
        v3 = shell->tri[t][2];
        
        double r1[3],r2[3],r3[3];   
        double r21[3],r31[3];
        double nI[3];          
        double norm_nI;
        
        dvec_1D_list3(shell->r, v1, r1);
        dvec_1D_list3(shell->r, v2, r2);
        dvec_1D_list3(shell->r, v3, r3);
        
        dvdiff3(r3, r1, r31);
        dvdiff3(r2, r1, r21);
        dcross3(r31, r21, nI);
        
        norm_nI = dnorm3(nI);
        dvscal3(nI, norm_nI);
        
        
        if( ddot3(nI, j_left)>LID_ALIGNMENT_AXIS)
        { 
            right_lid_tri.push_back(std::vector<int>());
            right_lid_tri[caps->n_tri_lid_r].push_back(v1);
            right_lid_tri[caps->n_tri_lid_r].push_back(v2);
            right_lid_tri[caps->n_tri_lid_r].push_back(v3);
            
            caps->r_r[v1] = 1; 
            caps->r_r[v2] = 1;
            caps->r_r[v3] = 1;
            
            caps->n_tri_lid_r++;
        }
        
        if( ddot3(nI, j_rigth)>LID_ALIGNMENT_AXIS) 
        {
            left_lid_tri.push_back(std::vector<int>());
            left_lid_tri[caps->n_tri_lid_l].push_back(v1);
            left_lid_tri[caps->n_tri_lid_l].push_back(v2);
            left_lid_tri[caps->n_tri_lid_l].push_back(v3);
            
            caps->r_l[v1] = 1; 
            caps->r_l[v2] = 1;
            caps->r_l[v3] = 1;
            
            caps->n_tri_lid_l++;
        }
            
    }
    
    
    // Number of particles in lids
    caps->nv_lid_r=0;
    caps->nv_lid_l=0;
    for(int i=0;i<shell->N; i++){
        if(caps->r_l[i]==1)
            caps->nv_lid_l++;
        if(caps->r_r[i]==1)
            caps->nv_lid_r++;
    }
    
    
    // Copy the data from vector to matrix  
    caps->tri_l = imatrix(caps->n_tri_lid_l, 3);
    caps->tri_r = imatrix(caps->n_tri_lid_r, 3);
    
    for(int t=0; t<caps->n_tri_lid_l; t++) {
        caps->tri_l[t][0] = left_lid_tri[t][0];
        caps->tri_l[t][1] = left_lid_tri[t][1];
        caps->tri_l[t][2] = left_lid_tri[t][2];
    }
    
    for(int t=0; t<caps->n_tri_lid_r; t++) {
        caps->tri_r[t][0] = right_lid_tri[t][0];
        caps->tri_r[t][1] = right_lid_tri[t][1];
        caps->tri_r[t][2] = right_lid_tri[t][2];
    }
    
    // Clear the vector 
    right_lid_tri.clear();
    left_lid_tri.clear();

}




void lids_mechanics(struct Shell *shell, struct lids *caps, double REINF_LIDS)
{
    caps->edges_str = ivector(shell->nedges);
    caps->edges_bend = ivector(shell->nedges);
    zeros_ivector(caps->edges_str, shell->nedges);
    zeros_ivector(caps->edges_bend, shell->nedges);

    for(int i=0; i<shell->nedges; i++) {
        
        // Change stiffness of the lids
        for(int t=0; t< caps->n_tri_lid_r; t++) {
            if(edge_in_tri(caps->tri_r[t], shell->edges[i][0], shell->edges[i][1])>=0) 
            {
                shell->k_spr_eff[i] *= REINF_LIDS;
                shell->lambda_eff[i] *= REINF_LIDS;
                caps->edges_str[i] = 1;
                caps->edges_bend[i] = 1;
                break;
            }
        }
        
        for(int t=0; t< caps->n_tri_lid_l; t++) {
            if(edge_in_tri(caps->tri_l[t], shell->edges[i][0], shell->edges[i][1])>=0) 
            {
                shell->k_spr_eff[i] *= REINF_LIDS;
                shell->lambda_eff[i] *= REINF_LIDS;
                caps->edges_str[i] = 1;
                caps->edges_bend[i] = 1;
                break;
            }
        }
        
        // Exclude bending for the transition lid-body
        int v1,v2,v3,v4;
        
        v1 = shell->op_edges[i][0];
        v4 = shell->op_edges[i][1];
        v2 = shell->edges[i][0];
        v3 = shell->edges[i][1];
    
        
        double r1[3],r2[3],r3[3],r4[3];   
               
        double r21[3],r31[3],         
               r24[3],r34[3];
        
        double nI[3], nJ[3];          
        
        dvec_1D_list3(shell->r, v1, r1);
        dvec_1D_list3(shell->r, v2, r2);
        dvec_1D_list3(shell->r, v3, r3);
        dvec_1D_list3(shell->r, v4, r4);
        
        dvdiff3(r3, r1, r31);
        dvdiff3(r2, r1, r21);
        dvdiff3(r2, r4, r24);
        dvdiff3(r3, r4, r34);        
        
        dcross3(r31, r21, nI);
        dcross3(r24, r34, nJ);
        
        double norm_nI = dnorm3(nI);
        double norm_nJ = dnorm3(nJ);
        
        dvscal3(nI, norm_nI);
        dvscal3(nJ, norm_nJ);
        
        
        if(ddot3(nI,nJ)<0.2) {
            shell->lambda_eff[i] = 0;
            caps->edges_bend[i] = 0;
        }
    }
}
