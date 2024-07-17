#ifndef _CPPUTILS_
#define _CPPUTILS_

#include <string>
#include <sstream>


/*******************************************************/
/* Numeric functions                                   */
/*******************************************************/

int imin(int *v, int n);
int imax(int *v, int n);
double dmin(double *v, int n);
double dmax(double *v, int n);

int iargmin(int *v, int n);
int iargmax(int *v, int n);
int dargmin(double *v, int n);
int dargmax(double *v, int n);




/*******************************************************/
/* String Functions                                    */
/*******************************************************/
std::string strcat(std::string str1, std::string str2);
std::string int2str(int n);














#endif
