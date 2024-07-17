#include "meshes.hpp"
#include "cpputils.hpp"
#include "allvec.hpp"



/******************************************************************************/
/*                                                                            */
/*                               MAIN FUNCTIONS                               */
/*                                                                            */
/******************************************************************************/



/******************************************************************************/
/*                                 markedges                                  */
/* Returns a vector with either 1 or val if the edge is in the list of        */
/* indexes provided by ids.                                                   */
/******************************************************************************/
double *markedges(int *ids, int **edges, double val, int nids, int nedges)
{
    double *edgesel = dvector(nedges);
    
    
    for(int i=0; i<nedges; i++)
        edgesel[i] = 1.0;
    
    for(int j=0; j<nedges; j++) {
        for(int i=0; i<nids; i++) {
            if( edges[j][0] == ids[i] || edges[j][1] == ids[i])
                edgesel[j] = val;
        }
    }
    
    return edgesel;
    
}


/******************************************************************************/
/*                              gen_nnlist                                    */
/* This function determines the indexes of the first ring neighbourhood for   */
/* every vertex                                                               */
/* This list can also be producded from the edges list                        */
/******************************************************************************/
int **gen_nnlist(int **tri, int *cn, int N_tri, int Np)
{
    int max_q;                  // Maximum coordination number, i.e. max number nn
    int **nnlist;               // Array with the list of ids of the particles interacting                   
    int *tmp_tri;               // Temporal variable for the ids in current triangle
    int *tmp_nn;                // Temporal variable for the coordination of the current particle
    
    max_q = imax(cn, Np);
    
    nnlist = imatrix(Np, max_q);
    zeros_imatrix(nnlist, Np, max_q);
    
    tmp_tri = ivector(3);
    tmp_nn = ivector(max_q);
    
    int v1,v2;
    for(int i=0; i<Np; i++){
        
        int nq = 0;
        
        for(int t=0; t<N_tri; t++) {
            
            if( tri[t][0] == i || tri[t][1]== i || tri[t][2] == i ) {
                
                tmp_tri[0] =  tri[t][0];
                tmp_tri[1] =  tri[t][1];
                tmp_tri[2] =  tri[t][2];

                //Exchange the positions to have the particle i in firts place and check
                xchangepos(tmp_tri, i);
                
                v1 = tmp_tri[1];
                v2 = tmp_tri[2];
                
                // Check if the other two are on the list of neightbours
                if( isnn( tmp_nn, v1, nq) == 0)  {
                    nnlist[i][nq] = v1;
                    tmp_nn[nq] = v1;
                    nq++;
                }
                
                if( isnn( tmp_nn, v2, nq) == 0)  {
                    nnlist[i][nq] = v2;
                    tmp_nn[nq] = v2;
                    nq++;
                }
                
            }
        }
        
    }
    
    return nnlist;
}



/******************************************************************************/
/*                                coordnum                                    */
/* This function returns a vector with the coordination number of each        */
/* vertex. Takes as input a vector for the triangulation, the number of       */
/* triangles and the number of particles.                                     */
/******************************************************************************/
int *coordnum(int **tri, int N_tri, int Np)
{
    int *cn;
    cn = ivector(Np);
    zeros_ivector(cn, Np);
    
    for(int i = 0; i<N_tri; i++) {
        cn[tri[i][0]]++;     
        cn[tri[i][1]]++;
        cn[tri[i][2]]++;
    }
    
    return cn;
    
}



/*******************************************************************************/
/*                              find_patch                                     */
/* This function returns a patch of particles around a given particle,         */
/* including itself. THe function takes as input the index of the particle     */
/* of interest pid, the list of nearest neightbours nnlist and its coodination */
/* number. The number of concentric layers is provided as n_layers.            */
/* patch_vtx has to be a vector with size 6*(n_layers+1) and is going to be    */
/* the output together with the  total number of particles (both non-repeated) */
/*******************************************************************************/
int find_patch(int pid,  int **nnlist, int *cn, int *patch_vtx, int n_layers)
{
    
    int v;                  // Vertex under exploration
    int cnv;                // Coordination number of v
    int tv;                 // Vertex of the current particle
    int tcn_patch;          // Patch size until the current layer l
    int cn_patch;           // Number of vertexes in the patch
    
    patch_vtx[0] = pid;
    cn_patch=1;
    
    int l=0;
    while(l<n_layers) {
        
        tcn_patch=cn_patch;
        for(int i=0; i<tcn_patch; i++) {
            v = patch_vtx[i];
            cnv = cn[v];
            
            for(int j=0; j<cnv; j++) {
                tv = nnlist[v][j];
                if(isnn(patch_vtx, tv, cn_patch) == 0) {
                    patch_vtx[cn_patch] = tv;
                    cn_patch++;
                }
            }
        }
        l++;
    }

    return cn_patch;
}




/*****************************************************************************************/
/*                                  find_edges                                           */
/* Returns a matrix with nedges x 2, where each row indicates the interaction of a pair  */
/* of vertexes considering a closed surface                                              */
/*****************************************************************************************/
int **find_edges(int **tri, int N_tri)
{
    int v1, v2, v3;             /* Vertexes for notation */
    int nedges;                 /* Expected number of edges of a closed surface */
    int te;                     /* Temporal number of edges found */
    int **edges;                /* List of the interactions (v1 interacts with v2) Nx2 array */
    
    nedges = 3*N_tri/2;
    edges = imatrix(nedges, 2);
    te = 0;
    
    for(int t=0; t<N_tri; t++) {
        v1 = tri[t][0];
        v2 = tri[t][1];
        v3 = tri[t][2];
        if( pairedge(edges, v1, v2, te) == 0) { 
            edges[te][0] = v1; 
            edges[te][1] = v2; 
            te++; 
        }
        
        if( pairedge(edges, v3, v2, te) == 0) { 
            edges[te][0] = v3; 
            edges[te][1] = v2; 
            te++; 
        }
        
        if( pairedge(edges, v1, v3, te) == 0) { 
            edges[te][0] = v1; 
            edges[te][1] = v3; 
            te++; 
        }
    }

    return edges;
    
}



/*****************************************************************************************/
/*                                  find_oposed_vertexes                                 */
/* In a closed surface every edge is shared by two triangles. This function returns a    */
/* a list of nedge rows by 2 columns with the two oposed vertexes to the edge. Takes as  */ 
/* input the triangulation, the number of triangles, the edges and the number of edges   */ 
/*****************************************************************************************/
int **find_oposed_vertexes(int **tri, int N_tri, int **edges, int nedges)
{
    
    int v1,v2;
    int nef;
    int oe;
    int **op_edges;
    
    op_edges = imatrix(nedges, 2);
    
    
    for(int e=0; e<nedges; e++) {
        
        nef = 0;
        v1 = edges[e][0];
        v2 = edges[e][1];
        
        for(int t=0; t<N_tri; t++) {
            oe = oposed_vertex_tri(tri, t, v1, v2);
            if(oe > -1) {
                op_edges[e][nef] = tri[t][oe];
                nef++;
            }
            
            if(nef>1) break;
        }
    }
    
    return op_edges;
}




/******************************************************************************/
/*                                                                            */
/*                          SECONDARY FUNCTIONS                               */
/*                                                                            */
/******************************************************************************/



/******************************************************************************/
/*                                   xchangepos                               */
/* Exchange the positions in a single triangle stri to have the index with    */
/* value val in the first position of the vector.                             */
/* Called by gennlist                                                         */
/******************************************************************************/
void xchangepos(int *stri, int val)
{
    
    if(stri[1] == val) {
        
        stri[0] = stri[0] + stri[1];
        stri[1] = stri[0] - stri[1];
        stri[0] = stri[0] - stri[1];
        
    } else if (stri[2] == val) {
        
        stri[0] = stri[0] + stri[2];
        stri[2] = stri[0] - stri[2];
        stri[0] = stri[0] - stri[2];
    }
    
}


/******************************************************************************/
/*                                  isnn                                      */
/* Returns true if nn_query is in the vector nnlist with maxnn entries        */
/* Called by gennlist                                                         */
/******************************************************************************/
int isnn(int *nnlist, int nn_query, int maxnn)
{
    int b = 0;
    for(int i=0; i<maxnn; i++) {
        if(nnlist[i] == nn_query) {
           b = 1; 
           break;
        }
    }
    
    return b;
    
}


/*****************************************************************************************/
/*                                      pairedge                                         */
/* This function is used in conjuntion with find_edges. The function determines if two   */ 
/* vertexes in the triangulation are already in the list of edges                        */
/*****************************************************************************************/
int pairedge(int **edges, int v, int u, int N)
{
    int is=0;
    
    for(int i=0; i<N; i++) {
        if( (edges[i][0] == v && edges[i][1] == u)   || 
            (edges[i][0] == u && edges[i][1] == v)   ) {
            is = 1;
            break;
        }
    }
    
    return is;
}


/*****************************************************************************************/
/*                                  oposed_vertex_tri                                    */
/* Determines the position of the vertex that is not forming the edge v1, v2 in the      */
/* triangle T. If the edge is not found returns -1.                                      */
/*****************************************************************************************/
int oposed_vertex_tri(int **tri, int T, int v1, int v2)
{
    int n=-1;
    
    if(tri[T][0] == v1 && tri[T][1] == v2 )
        n=2;
    else if(tri[T][0] == v2 && tri[T][1] == v1 )
        n=2;
    else if(tri[T][0] == v1 && tri[T][2] == v2 )
        n=1;
    else if(tri[T][0] == v2 && tri[T][2] == v1 )
        n=1;
    else if(tri[T][1] == v1 && tri[T][2] == v2 )
        n=0;
    else if(tri[T][1] == v2 && tri[T][2] == v1 )
        n=0;
    
    return n;
}


/*****************************************************************************************/
/*                                  edge_in_tri                                          */
/* Same as oposed_vertex_tri but taking a single triangle as input                       */
/*****************************************************************************************/
int edge_in_tri(int *tri, int v1, int v2)
{
    int n=-1;
    
    if(tri[0] == v1 && tri[1] == v2 )
        n=2;
    else if(tri[0] == v2 && tri[1] == v1 )
        n=2;
    else if(tri[0] == v1 && tri[2] == v2 )
        n=1;
    else if(tri[0] == v2 && tri[2] == v1 )
        n=1;
    else if(tri[1] == v1 && tri[2] == v2 )
        n=0;
    else if(tri[1] == v2 && tri[2] == v1 )
        n=0;
    
    return n;
}
