#include "allvec.hpp"


/************************************************************************/
/*                                ivector                               */
/* Allocates a double vector of size M                                  */
/************************************************************************/
int *ivector(int M)
{
    int *v = new int [M];
    return v;
}


/************************************************************************/
/*                                dvector                               */
/* Allocates a double vector of size M                                  */
/************************************************************************/
double *dvector(int M)
{
    double *v = new double [M];
    return v;
}


/************************************************************************/
/*                           val_ivector                                */
/* Assign the value specified by val to the N first indexes of vecotr v */
/* allocated with ivector.                                              */
/************************************************************************/
void val_ivector(int *ivec, int val, int N)
{
    for(int i=0; i<N; i++)
        ivec[i] = val;
}


/************************************************************************/
/*                           val_dvector                                */
/* Assign the value specified by val to the N first indexes of vecotr v */
/* allocated with dvector.                                              */
/************************************************************************/
void val_dvector(double *dvec, double val, int N)
{
    for(int i=0; i<N; i++)
        dvec[i] = val;
}


/************************************************************************/
/*                           zeros_ivector                              */
/* Assign a zero value to all the positions in a M ivector              */
/************************************************************************/
void zeros_ivector(int *ivec, int M)
{
    for(int i=0; i<M; i++)
        ivec[i]=0;
}


/************************************************************************/
/*                           zeros_dvector                              */
/* Assign a zero value to all the positions in a M dvector              */
/************************************************************************/
void zeros_dvector(double *dvec, int M)
{
    for(int i=0; i<M; i++)
        dvec[i]=0;
}


/************************************************************************/
/*                           free_ivector                               */
/* Free the memory of the vector allocated with dvector                 */
/************************************************************************/
void free_ivector(int *v)
{
    delete[] v;
}


/************************************************************************/
/*                           free_dvector                               */
/* Free the memory of the vector allocated with dvector                 */
/************************************************************************/
void free_dvector(double *v)
{
    delete[] v;
}




/************************************************************************/
/************************************************************************/
/************************************************************************/



/************************************************************************/
/*                                imatrix                               */
/* Allocates a int matrix of size MxN                                   */
/************************************************************************/
int **imatrix(int M, int N)
{   

    int **mtx = new(std::nothrow) int *[M];
    if(!mtx) 
        cout << "Error in the allocation of imatrix" << endl;
    for(int i=0; i<M;i++)
        mtx[i] = new int[N];
    return mtx;
}

/************************************************************************/
/*                                dmatrix                               */
/* Allocates a double matrix of size MxN                                */
/************************************************************************/
double **dmatrix(int M, int N)
{
    double **mtx = new(std::nothrow) double *[M];
    if(!mtx) 
        cout << "Error in the allocation of dmatrix" << endl;
    for(int i=0; i<M;i++)
        mtx[i] = new double[N];
    return mtx;
}


/************************************************************************/
/*                           val_imatrix                                */
/* Assign the value specified by value to all the positions in a MxN    */
/* imatrix                                                              */ 
/************************************************************************/
void val_imatrix(int **imtx, int val, int M, int N)
{
    for(int i=0; i<M; i++) {
        for(int j=0; j<N; j++) 
            imtx[i][j] = val;
    }
}

/************************************************************************/
/*                           val_dmatrix                                */
/* Assign the value specified by value to all the positions in a MxN    */
/* imatrix                                                              */ 
/************************************************************************/
void val_dmatrix(double **dmtx, double val, int M, int N)
{
    for(int i=0; i<M; i++) {
        for(int j=0; j<N; j++) 
            dmtx[i][j] = val;
    }
}


/************************************************************************/
/*                           zeros_imatrix                              */
/* Assign a zero value to all the positions in a MxN imatrix            */
/************************************************************************/
void zeros_imatrix(int **imtx, int M, int N)
{
    for(int i=0; i<M; i++) {
        for(int j=0; j<N; j++) 
            imtx[i][j] = 0;
    }
}

/************************************************************************/
/*                           zeros_dmatrix                              */
/* Assign a zero value to all the positions in a MxN imatrix            */
/************************************************************************/
void zeros_dmatrix(double **dmtx, int M, int N)
{
    for(int i=0; i<M; i++) {
        for(int j=0; j<N; j++) 
            dmtx[i][j] = 0;
    }
}


/************************************************************************/
/*                             free_imatrix                             */
/* Free the memory for the matrix of M rows allocated with imatrix      */
/************************************************************************/
void free_imatrix(int **imtx, int M)
{
    for(int i=0; i<M;i++)
        delete[] imtx[i];        
    delete[] imtx; 
}


/************************************************************************/
/*                             free_dmatrix                             */
/* Free the memory for the matrix of M rows allocated with dmatrix      */
/************************************************************************/
void free_dmatrix(double **dmtx, int M)
{
    for(int i=0; i<M;i++)
        delete[] dmtx[i];
    delete[] dmtx; 
}




/************************************************************************/
/************************************************************************/
/************************************************************************/




/************************************************************************/
/*                            i3matrix                                  */
/* Allocated an integer 3D matrix of M rows, N columns and depth K      */
/************************************************************************/
int ***i3matrix (int M, int N, int K)
{
    int ***mtx;  
    
    mtx = new(std::nothrow) int**[M];
    if(!mtx) 
        cout << "Error in the allocation of i3matrix" << endl;
    for(int i=0; i<M; i++) {
        mtx[i] = new int*[N];
        for(int j=0; j<N; j++)
            mtx[i][j] = new int[K];
    }

    return mtx;
}


/************************************************************************/
/*                            d3matrix                                  */
/* Allocated a double 3D matrix of M rows, N columns and depth K        */
/************************************************************************/
double ***d3matrix (int M, int N, int K)
{
    double ***mtx;  
    
    mtx = new(std::nothrow) double**[M];
    if(!mtx) 
        cout << "Error in the allocation of d3matrix" << endl;

    for(int i=0; i<M; i++) {
        mtx[i] = new double*[N];
        for(int j=0; j<N; j++)
            mtx[i][j] = new double[K];
    }

    return mtx;
}


/************************************************************************/
/*                           val_i3matrix                               */
/* Assign the value specified by value to all the positions in a MxNxK  */
/* 3D matrix defined by i3matrix                                        */ 
/************************************************************************/
void val_i3matrix(int ***imtx3, int val, int M, int N, int K)
{
    for(int i=0; i<M;i++){
        for(int j=0; j<N; j++) {
            for(int k=0; k<K; k++)
                imtx3[i][j][k] = val;
        }
    }
}





/************************************************************************/
/*                           val_d3matrix                               */
/* Assign zeros to all the positions in a MxNxK 3D matrix defined by    */
/* d3matrix.                                                            */ 
/************************************************************************/
void val_d3matrix(double ***dmtx3, double val, int M, int N, int K)
{
    for(int i=0; i<M;i++){
        for(int j=0; j<N; j++) {
            for(int k=0; k<K; k++)
                dmtx3[i][j][k] = val;
        }
    }
}

/************************************************************************/
/*                           zeros_i3matrix                             */
/* Fills with 0s all the positions in a MxNxK 3D matrix defined         */
/* by i3matrix                                                          */ 
/************************************************************************/
void zeros_i3matrix(int ***imtx3, int M, int N, int K)
{
    for(int i=0; i<M;i++){
        for(int j=0; j<N; j++) {
            for(int k=0; k<K; k++)
                imtx3[i][j][k] = 0;
        }
    }
}



/************************************************************************/
/*                           zeros_d3matrix                             */
/* Fills with 0s all the positions in a MxNxK 3D matrix defined         */
/* by d3matrix                                                          */ 
/************************************************************************/
void zeros_d3matrix(double ***dmtx3, int M, int N, int K)
{
    for(int i=0; i<M;i++){
        for(int j=0; j<N; j++) {
            for(int k=0; k<K; k++)
                dmtx3[i][j][k] = 0;
        }
    }
}



/************************************************************************/
/*                             free_i3matrix                            */
/* Free the memory for the 3D matrix of M rows and N columns allocated  */
/* by using i3matrix                                                    */
/************************************************************************/
void free_i3matrix(int ***mtx, int M, int N)
{
    for(int i=0; i<M; i++) {
        for(int j=0; j<N; j++)
            delete[] mtx[i][j];
        delete[] mtx[i];
    }
    delete[] mtx;
}


/************************************************************************/
/*                             free_d3matrix                            */
/* Free the memory for the 3D matrix of M rows and N columns allocated  */
/* by using d3matrix                                                    */
/************************************************************************/
void free_d3matrix(double ***mtx, int M, int N)
{
    for(int i=0; i<M; i++) {
        for(int j=0; j<N; j++)
            delete[] mtx[i][j];
        delete[] mtx[i];
    }
    delete[] mtx;
}


