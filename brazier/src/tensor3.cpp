/****************************************************************************

    tensor3. This is a basic library to do simple operations with vectors
             (rank 1 tensors) in R3. This makes the execution faster than
             using for loops.
              
    Copyright (C) 2020  Cesar L. Pastrana

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

****************************************************************************/


#include "tensor3.hpp"







/****************************************************************/
/*                             ddot3                            */
/* Calculates the dot product between vectors v1 and v2         */
/****************************************************************/
double ddot3(double*  v1, double*  v2)
{
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}






/****************************************************************/
/*                           dvadd3                             */
/* Returns vector pointing between v1 and v2                    */
/****************************************************************/
void dvadd3(double *v1, double *v2, double *v_out)
{
    v_out[0] = v1[0] + v2[0];
    v_out[1] = v1[1] + v2[1];
    v_out[2] = v1[2] + v2[2];
}



/****************************************************************/
/*                          dmidvec                             */
/* This function returns the mid point of the edge connecting   */
/* the vectors v1 and v2                                        */
/****************************************************************/
void dmidvec3(double *v1, double *v2, double *v_out)
{
    v_out[0] = 0.5*(v1[0] + v2[0]);
    v_out[1] = 0.5*(v1[1] + v2[1]);
    v_out[2] = 0.5*(v1[2] + v2[2]);
}

/****************************************************************/
/*                            dvscal                            */
/* Returns the vector scaled by the value input in A            */
/****************************************************************/
void dvscal3(double *v, double A)
{
    v[0] = v[0]/A;
    v[1] = v[1]/A;
    v[2] = v[2]/A;
}



/****************************************************************/
/*                         dassignv3                            */
/* Assigns the values of the vectors v_in to the vector v_out   */
/* scaled by the constant s                                     */
/****************************************************************/
void dassignv3(double *v_in, double *v_out, double s)
{
    v_out[0] = v_in[0]*s;
    v_out[1] = v_in[1]*s;
    v_out[2] = v_in[2]*s;
}



/****************************************************************/
/*                          dvec_list3                          */
/* Exctract the vector n from a list off R3 vectors set         */
/* as a matrix with vectors in rows and vec components in cols  */
/****************************************************************/
void dvec_list3(double **list, int n, double *v)
{
    v[0] = list[n][0];
    v[1] = list[n][1];
    v[2] = list[n][2];
}



