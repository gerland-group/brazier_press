#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream> 
#include <algorithm>
#include <string>
#include <vector>
#include <ctime>        // Current day
#include <sys/time.h>     // Exect time
#include <iomanip>

 
    
    
#ifndef PI
    #define PI 3.14159265358979323846
#endif

#ifndef EPS
    #define EPS 1.0E-16                                     // Numerical error value
#endif

/***********************************************  PHYSICAL PROPERTIES  *************************************************/
#define  DIMV               3                               // R^3 for vector arrangment of particles in the matrix

#define  CONVATM            0.101325                        // Conversion factor between atm (input file) and pN/nm^2
#define  YOUNG_MESH_FACTOR  0.866025403784439               // Conversion factor between 2D Young modulus and spring constant
#define  POISSON            1.0/3.0                         // Poisson ratio of the mesh


/***********************************************  FILE PROPERTIES  *************************************************/

#define INIT_COORDS_FNAME        "init_coords.dat"                      //  Path of the initial file to load coardinates
#define IN_MESH_FNAME            "mesh.dat"                             //  File name to save trajectories (local directory)
#define PARAMETERS_FILE          "params.conf"                          //  Parameters file

#define PATH_PRESS_CONFS         "./press/"                             //  Minimised coordiantes after NLCG of pressurization
#define PATH_BENT_CONFS          "./bending/"                           //  Minimised coordiantes after NLCG of bendng
#define PRESS_COORDS_FNAME       "_press.dat"                           //  Minimised coordiantes after NLCG of bendng
#define BENT_COORDS_FNAME        "_brazier.dat"                         //  Minimised coordiantes after NLCG of bendng
#define LOG_FILE                 "./summary/output.log"                 //  Log file name with simulation properties (local directory)
#define MCMIN_PRESS_COORDS_FNAME "./summary/minim_coords_press.dat"     //  Final (minimised) coordiantes after NLCG
#define MCMIN_COORDS_FNAME       "./summary/minim_coords.dat"           //  Final (minimised) coordiantes after NLCG 
#define PRESSURES_FNAME          "./summary/pressures.dat"              //  Pressures employed in the ramp
#define BEND_ANGLES_FNAME         "./summary/bent_angles.dat"           //  Angles imposed in the bending ponteital (=/= that acquired by lids)

#define LIDS_EDGES_STR_FNAME      "./lids/edges_stretchin.dat"          // Edges reinforced on lids for stretching
#define LIDS_EDGES_BEND_FNAME     "./lids/edges_bending.dat"            // Edges reinforced on lids for bending (does not include boyd-lid connection)
#define LIDS_TRIS_LEFT_FNAME      "./lids/tris_left.dat"                // Tris left lid
#define LIDS_TRIS_RIGHT_FNAME     "./lids/tris_rigth.dat"               // Tris rigth lid
#define LIDS_VERTEX_LEFT_FNAME    "./lids/vertices_left.dat"            // Vertices left lid
#define LIDS_VERTEX_RIGTH_FNAME    "./lids/vertices_right.dat"          // Vertices rigth lid




/*********************************   Namespaces and type definitions   *******************************/

// Standard namespace
using std::cout;
using std::cin;
using std::endl;
using std::vector;
using std::ofstream; 
using std::ifstream; 
using std::string;
using std::getline;
using std::stringstream;
using std::nothrow; // Do not send exception but a null pointer in new

typedef vector<vector<double>> dvector_vl;
typedef vector<vector<int>> ivector_vl;



/***********************************    FUNCTION DECLARATIONS    ************************************/

double avedge_length(double *r, int **edges, int nedges);
void calc_eql(double *r, int **edges, int nedges, double *eql);
void save_coords_file(double *r, int N, std::string fname);
void print_gui(int N, int T, double P_atm, double Kappa, double E3D, double h, double K_SPR);

