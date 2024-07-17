#ifndef _MESHES_
#define _MESHES_



// Main functions
double *markdges(int *ids, int *edges, double val, int nids, int nedges);
int **gen_nnlist(int **tri, int *cn, int N_tri, int Np);
int *coordnum(int **tri, int N_tri, int Np);
int find_patch(int pid,  int **nnlist, int *cn, int *patch_vtx, int n_layers);
int **find_edges(int **tri, int N_tri);
int **find_oposed_vertexes(int **tri, int N_tri, int **edges, int nedges);

// Accesory functions
int pairedge(int **edges, int v, int u, int N);
void xchangepos(int *stri, int val);
int isnn(int *nnlist, int nn_query, int maxnn);
int oposed_vertex_tri(int **tri, int T, int v1, int v2);
int edge_in_tri(int *tri, int v1, int v2);

#endif
