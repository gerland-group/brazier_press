#include "cpputils.hpp"




/****************************************************************************************/
/*                                                                                      */
/*                                  NUMERIC FUNCTIONS                                   */
/*                                                                                      */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/*****************************************************************************************/
/*                                      imin                                             */
/* Minimum value in the double vector input as v.                                        */
/*****************************************************************************************/
int imin(int *v, int n)
{
    int minv=v[0];
    for(int i=1; i<n; i++){
        if( v[i] < minv )
            minv=v[i];
    }
    return minv;
}

/*****************************************************************************************/
/*                                      imax                                             */
/* Maximum value in the double vector input as v.                                        */
/*****************************************************************************************/
int imax(int *v, int n)
{
    int maxv=v[0];
    for(int i=1; i<n; i++){
        if( v[i] > maxv )
            maxv=v[i];
    }
    return maxv;
}

/*****************************************************************************************/
/*                                      dmin                                             */
/* Minimum value in the double vector input as v.                                        */
/*****************************************************************************************/
double dmin(double *v, int n)
{
    double minv=v[0];
    for(int i=1; i<n; i++){
        if( v[i] < minv )
            minv=v[i];
    }
    return minv;
}

/*****************************************************************************************/
/*                                      dmax                                             */
/* Maximum value in the double vector input as v.                                        */
/*****************************************************************************************/
double dmax(double *v, int n)
{
    double maxv=v[0];
    for(int i=1; i<n; i++){
        if( v[i] > maxv )
            maxv=v[i];
    }
    return maxv;
}




/*****************************************************************************************/
/*                                      iargmin                                          */
/* Index of the maximum value in the int vector input as v.                           */
/*****************************************************************************************/
int iargmin(int *v, int n)
{
    int minv=v[0];
    int idx=0;
    for(int i=1; i<n; i++){
        if( v[i] < minv ){
            minv=v[i];
            idx=i;    
        }
    }
    return idx;
}


/*****************************************************************************************/
/*                                      iargmax                                          */
/* Index of the maximum value in the int vector input as v.                           */
/*****************************************************************************************/
int iargmax(int *v, int n)
{
    int maxv=v[0];
    int idx=0;
    for(int i=1; i<n; i++){
        if( v[i] > maxv ){
            maxv=v[i];
            idx=i;    
        }
    }
    return idx;
}


/*****************************************************************************************/
/*                                      dargmin                                          */
/* Index of the maximum value in the double vector input as v.                           */
/*****************************************************************************************/
int dargmin(double *v, int n)
{
    double minv=v[0];
    int idx=0;
    for(int i=1; i<n; i++){
        if( v[i] < minv ){
            minv=v[i];
            idx=i;    
        }
    }
    return idx;
}


/*****************************************************************************************/
/*                                      dargmax                                          */
/* Index of the maximum value in the double vector input as v.                           */
/*****************************************************************************************/
int dargmax(double *v, int n)
{
    double maxv=v[0];
    int idx=0;
    for(int i=1; i<n; i++){
        if( v[i] > maxv ){
            maxv=v[i];
            idx=i;    
        }
    }
    return idx;
}



/****************************************************************************************/
/*                                                                                      */
/*                                  STRING FUNCTIONS                                    */
/*                                                                                      */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


// Concatenate two strings str1, str2.
std::string strcat(std::string str1, std::string str2)
{
    std::stringstream ss;        
    ss << str1 << str2;
    return ss.str();
}

// Converts an int value to a string
std::string int2str(int n)
{
    std::stringstream ss;
    ss << n;
    return ss.str();
}


/*---------------------------------------------------------------------------------------*/

