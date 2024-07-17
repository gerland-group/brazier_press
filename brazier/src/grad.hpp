#include <cmath>



#ifndef DIMV
#define  DIMV               3                               // R^3 for vector arrangment of particles in the matrix
#endif 

#define ONESIXTH    0.166666666666666666666667



struct Shell {
    
    double *r;                // Coordinates per vertex
    int N;                    // Number of vertices
    
    int **tri;                // Triangulation array
    int N_tri;                // Number of vertexex/particles in the mesh

    int nedges;               // Number of edges (determined from the triangulation)
    int **edges;              // Array with the pair of edges
    int **op_edges;           // Pair of vertexes defining the oposed vertexes to each edge

    double l0;                // Average distance between vertexes
    double *eql;              // Rest length of the springs
    double *dx, *dy, *dz, *d; // Distance and distance components between particles 
    
    // Physical properties
    double p;                 // Pressure
    double *k_spr_eff;        // Stretching stiffness for every edge modified for helix domains
    double *lambda_eff;       // Effective bending stiffness for a given edge (exclude lids)
};



void grad(struct Shell *shell, struct lids *caps, double *gf);


void calcdr(double *r, 
            double *dx, 
            double *dy,
            double *dz,
            double *d,
            int **edges,
            int nedges);

void graduniv_byu(double *v_uni, double normv, double *u, double *out_v);

