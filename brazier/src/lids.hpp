#include <math.h>
#include <vector>


#ifndef PI
    #define PI 3.14159265358979323846
#endif


#ifndef _LIDS_
#define _LIDS_

#define LID_ALIGNMENT_AXIS 0.99

struct lids
{
    
    // Triangles composing the lids
    int n_tri_lid_r;     // Number of particles of righth lid
    int n_tri_lid_l;     // Number of particles of left lid

    int **tri_l;         // Triangles composing the left lid
    int **tri_r;         // Triangles composing the rigth lid
    
    // Particles composing the lids
    int nv_lid_r;       // Particles right lid
    int nv_lid_l;       // Particles left lid
    int *r_l;           // Particles composing the left lids
    int *r_r;           // Particles composing teh right lids;
    
    int *edges_str;      // Edges marked with different reinforcement factor
    int *edges_bend;     // Edges marked with different bending factor
    
    // Equilibrium normal vectors
    double nl0[3];      // Normal vector left
    double nr0[3];      // Normal vector right
    
    
    // Stifness of the lids
    double K_BRAZIER;       // Stifness of the bending trap
    double K_BRAZIER_CM;    // Stifness to keep the lids aligned
};



int find_n_lids(double *r, int N);
void find_tri_lids(struct Shell *shell, struct lids *caps);
void lids_mechanics(struct Shell *shell, struct lids *caps, double REINF_LIDS);

#endif

