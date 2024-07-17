/*******************************************************************************

    brazier:   Non-linear conjugate gradient minimisation of the energy function
               of a pressurized triangularized closed surface
             
    Copyright (C) 2022, Cesar L. Pastrana

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


*******************************************************************************/

#include "brazier.hpp"
#include "nlcg.hpp"
#include "grad.hpp"
#include "lids.hpp"
#include "meshes.hpp"
#include "tensor3.hpp"
#include "allvec.hpp"
#include "cpputils.hpp"




int main()
{
    
    
    // Main declarations +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    struct NLCG_params minparams;           // Struct for the nlcg parameters
    struct Shell shell;                     // Struct for the shell
    
    
    dvector_vl r_mtx;                       // Coordinates xyz of each particle (from file) [nm]        
    ivector_vl tri_vl;                      // Array of triangles defining the surface in vector libraru form
    
    double P;                               // Pressure parameter (value loaded from file) [atm -> pN nm]
    double DP;                              // Variation of pressure (value loaded from file)  [atm -> pN nm]
    double MAX_THETA;                       // Maximum bending angle
    double N_THETA;                         // Number of minimisation steps until reaching MAX_THETA
    double BRAZIER_FACTOR;                  // Stiffness of the potential for bending the tube
    double BRAZIER_FACTOR_CM;               // Stiffness of the potential for keeping the lids in place
    double KAPPA;                           // Bending stiffness (value loaded from file) [pN nm]
    double KMF;                             // Correction factor of the mesh (value loaded from file)
    double LAMBDA;                          // Effective bending stiffness in the triangular mesh [pN nm]
    double E_3D;                            // 3D Young modulus (value loaded from file) [MPa -> pN/nm^2]
    double h;                               // Shell thickness  (value loaded from file) [nm]
    double Y_2D;                            // 2D Young modulus [N/m -> pN/nm]
    double K_SPR;                           // Stretching stiffness [pN/nm]
    double REINF_LIDS;                      // Reinforcment factor of the lids
    
    
    
    
    //------------------------------------------------------------------------------------------------------
    
    
    
    //************************************************************************************//
    // 1. Load data                                                                       //
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

    // 1.1  Load the file with the initial coordiantes of the particles
    ifstream infile_coords(INIT_COORDS_FNAME);
    string line;
    int l=0;                                 // Line control for data loading
    if(infile_coords.is_open()) {
        while (getline(infile_coords, line))
        {
            double val;
            stringstream ss(line);
            r_mtx.push_back(vector<double>(0));

            while (ss >> val)
                r_mtx[l].push_back(val);
            
            ++l;
        }
        infile_coords.close();
    } else {
        cout << "Error while opening the coordinates file " << INIT_COORDS_FNAME << ". Exiting." << endl;
        exit(EXIT_FAILURE);
    }
        
    // 1.2  Load the file with the mesh data
    ifstream infile_tri(IN_MESH_FNAME);
    l = 0;
    if((infile_tri.is_open()))
    {
        while (getline(infile_tri, line))
        {
            int val;
            stringstream ss(line);
            tri_vl.push_back(vector<int>(0));

            while (ss >> val)
                tri_vl[l].push_back(val);
            
            ++l;
        }
        infile_tri.close();
    } else {
        cout << "Error while opening the mesh file " << IN_MESH_FNAME << ". Exiting." << endl;
        exit(EXIT_FAILURE);
    }
    
    
    // 1.3 Load the simulation parameters values ------------------------------------
    
    // Set default values
    DP = 0.0;
    P = 0.0;
    E_3D = 30.0;
    h = 5.0;
    K_SPR = 100.0;
    KAPPA = 10.0;
    KMF = 2.0/sqrt(3.0);
    LAMBDA = KAPPA*KMF; 
    N_THETA = 0;
    MAX_THETA = 0;
    BRAZIER_FACTOR = 1.0;
    BRAZIER_FACTOR_CM = 0.0;
    REINF_LIDS = 1.0;
    
    
    minparams.MAX_ITS_NLCG_FACT = 10;                                       
    minparams.MAX_ITS_LINE_SEARCH = 50;                                    
    minparams.TOLERANCE_CRIT_LINE_SEARCH = 1.0e-8;                          
    minparams.TOLERANCE_CRIT_NLCG = 1.0e-4;                                 
    minparams.SIGMA_0 = 0.05;
    
    
    // Load from file
    // std::ifstream is RAII, i.e. no need to call close
    ifstream infile_params(PARAMETERS_FILE);
    if (infile_params.is_open())
    {
        while (getline(infile_params, line)) 
        {
            line.erase(remove_if(line.begin(), line.end(), ::isspace),line.end());
            
            // Commments with the sharp symbol
            if (line[0] == '#' || line.empty()) continue; 

            auto delim_pos = line.find("=");
            auto tag = line.substr(0, delim_pos);
            auto val = line.substr(delim_pos + 1);
            
            // External factors
            if (tag == "PRESSURE") P = -CONVATM*(std::stod(val));       // Negative
            else if (tag == "DP") DP = -CONVATM*(std::stod(val));       // Negative
            else if (tag == "MAX_THETA") MAX_THETA = std::stod(val);       
            else if (tag == "N_THETA") N_THETA = std::stoi(val);       
            else if (tag == "BRAZIER_FACTOR") BRAZIER_FACTOR = std::stod(val);       
            else if (tag == "BRAZIER_FACTOR_CM") BRAZIER_FACTOR_CM = std::stod(val);       
            
            // Mechanical parameters
            else if (tag == "KAPPA") KAPPA = std::stod(val);           
            else if (tag == "KAPPA_MESH_FACTOR") KMF = std::stod(val);
            else if (tag == "3D_YOUNG_MODULUS") E_3D = std::stod(val);
            else if (tag == "SHELL_THICKNESS") h = std::stod(val);
            else if (tag == "REINF_LIDS") REINF_LIDS = std::stod(val);
            
            //Minimization parameters
            else if (tag == "MAX_ITS_NLCG_FACT")  minparams.MAX_ITS_NLCG_FACT = std::stoi(val);           
            else if (tag == "MAX_ITS_LINE_SEARCH") minparams.MAX_ITS_LINE_SEARCH = std::stoi(val);
            else if (tag == "TOLERANCE_CRIT_LINE_SEARCH")  minparams.TOLERANCE_CRIT_LINE_SEARCH = std::stod(val);
            else if (tag == "TOLERANCE_CRIT_NLCG") minparams.TOLERANCE_CRIT_NLCG = std::stod(val);
            else if (tag == "SIGMA_0") minparams.SIGMA_0 = std::stod(val);
        }
        
        // In principle not necessary beucase ifstream self cleans
        infile_params.close();
    }
    else 
    {
        cout << "Error while opening the parameters file " << PARAMETERS_FILE << ". Exiting." << endl;
        exit(EXIT_FAILURE);
    }
    
    // Calculate bending modulus from the dependency between E,h and v
    if(KAPPA<0.0)
        KAPPA = E_3D*h*h*h/(12*(1 - POISSON*POISSON));
        
 
    
    // Change the input of the brazier bending from degrees to rads
    MAX_THETA *=  PI/180.0;
    
    // 1.4  Calculates basic properties and change the structure of the allocated arrays
    shell.N = r_mtx.size(); 
    shell.N_tri = tri_vl.size();
    LAMBDA = KAPPA*KMF;
    Y_2D = E_3D*h;
    K_SPR = Y_2D*YOUNG_MESH_FACTOR;
    
    
    // Copy the set of particles to a vector form r = [x1, y1, z1, x2, y2, z2,...,  xN, yN, zN]
    shell.r = dvector(shell.N*DIMV);
    
    for(int i = 0; i<shell.N; i++) {
        shell.r[i*DIMV] = r_mtx[i][0];
        shell.r[i*DIMV + 1] = r_mtx[i][1];
        shell.r[i*DIMV + 2] = r_mtx[i][2];
    }
    r_mtx.clear();
    
    
    // The mesh is transformed from vector to a "new" object
    shell.tri = imatrix(shell.N_tri,3);
    for(int t = 0; t<shell.N_tri; t++) {
        shell.tri[t][0] = tri_vl[t][0];
        shell.tri[t][1] = tri_vl[t][1];
        shell.tri[t][2] = tri_vl[t][2];
    }
    tri_vl.clear();
    
    
    // Show GUI
    print_gui(shell.N, shell.N_tri, P, KAPPA, E_3D, h, K_SPR);
   
    
   
    
    //------------------------------------------------------------------------------------//
    
    //************************************************************************************//
    // 3. Extract topology of the mesh                                                    //
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//   
    
    cout << endl << " - 1) Analysis mesh properties... ";
    fflush(stdout);
    
   
    // Edges properties
    shell.nedges = 3*shell.N_tri/2;
    shell.edges = find_edges(shell.tri, shell.N_tri);
    shell.op_edges = find_oposed_vertexes(shell.tri, shell.N_tri, shell.edges, shell.nedges);
    
    // Equilibrium length
    shell.eql = dvector(shell.nedges);
    calc_eql(shell.r, shell.edges, shell.nedges, shell.eql);
    shell.l0 = avedge_length(shell.r, shell.edges, shell.nedges);
    

    
    
    //------------------------------------------------------------------------------------//
    
    
    
    //************************************************************************************//
    // 4. Calculate properties of the object	  	 				                      //
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    cout << "Done" << endl << " - 2) Analysis lids properties... "; 
    fflush(stdout);

    struct lids caps;
    
    // Original values of the normals
    caps.nl0[0] = 0;  caps.nl0[1] = 1;  caps.nl0[2] = 0; 
    caps.nr0[0] = 0;  caps.nr0[1] = -1;  caps.nr0[2] = 0; 
    
    
    // 4.1 Determine the number of particles of the lids 2*det_nlidtot(R, spacing) + 2
    //int N_lids = find_n_lids(shell.r, shell.N, R, sp);
    
    // 4.2 Find the triangles that compose the lids
    find_tri_lids(&shell, &caps);
    
    
    // Assignment of stretching stiffness and bending
    shell.k_spr_eff = dvector(shell.nedges);
    shell.lambda_eff = dvector(shell.nedges);
    val_dvector(shell.k_spr_eff, K_SPR, shell.nedges);
    val_dvector(shell.lambda_eff, LAMBDA, shell.nedges);
    lids_mechanics(&shell, &caps, REINF_LIDS);
    

    // Save the list of edges stretching or bending
    ofstream edges_stretching;
    ofstream edges_lambda;
    edges_stretching.open(LIDS_EDGES_STR_FNAME, ofstream::out);
    edges_lambda.open(LIDS_EDGES_BEND_FNAME, ofstream::out);
    for(int i = 0; i<shell.nedges; i++) {
        if(caps.edges_str[i] == 1)
            edges_stretching << shell.edges[i][0]  << "\t" << shell.edges[i][1] << endl;
        
        if(caps.edges_bend[i] == 1)
            edges_lambda << shell.edges[i][0]  << "\t" << shell.edges[i][1] << endl;
            
    }
    
    
    
    // Save list of particles in list
    ofstream fid_vertices_lid_r;
    ofstream fid_vertices_lid_l;    
    fid_vertices_lid_r.open(LIDS_VERTEX_RIGTH_FNAME, ofstream::out);
    fid_vertices_lid_l.open(LIDS_VERTEX_LEFT_FNAME, ofstream::out);
    for(int i = 0; i<shell.N; i++) {
        
        if(caps.r_r[i] == 1) 
            fid_vertices_lid_r << i << endl;
        
        if(caps.r_l[i] == 1)             
            fid_vertices_lid_l << i << endl;
    }
    fid_vertices_lid_r.close();
    fid_vertices_lid_l.close();
    
    
    // Save the list of triangles in lids
    ofstream fid_tri_r;
    ofstream fid_tri_l;
    fid_tri_l.open(LIDS_TRIS_LEFT_FNAME, ofstream::out);
    fid_tri_r.open(LIDS_TRIS_RIGHT_FNAME, ofstream::out);
    for(int i = 0; i<caps.n_tri_lid_r; i++) {
        fid_tri_l << caps.tri_l[i][0] << "\t" << caps.tri_l[i][1] << "\t" << caps.tri_l[i][2] << endl;
        fid_tri_r << caps.tri_r[i][0] << "\t" << caps.tri_r[i][1] << "\t" << caps.tri_r[i][2] << endl;
    }
    fid_tri_l.close();
    fid_tri_r.close();

  
    

    // Save LOG file with simulation parameters.
    ofstream fid_logfile;
    fid_logfile.open (LOG_FILE, ofstream::out);
    
    fid_logfile << "[NLCG_BRAZIER_VERSION] = " << 2.1 << endl;
    fid_logfile << "[COMPILATION_DATE] = " << __DATE__ << endl;
    fid_logfile << "[COMPILATION_TIME] = " << __TIME__ << endl;
    #if defined(_OPENMP)
        fid_logfile << "[OPENMP] = TRUE" << endl << endl;
    #else
        fid_logfile << "[OPENMP] = FALSE" <<  endl << endl;
    #endif
    
    fid_logfile << "[#PARTICLES] = " << shell.N << endl;
    fid_logfile << "[#FACES] = " << shell.N_tri << endl;
    fid_logfile << "[#EDGES] = " << shell.nedges << endl;
    fid_logfile << "[CHARACTERISTIC_REST_LENGTH_0](nm) = " << std::setprecision(10) << shell.l0 << endl;
    fid_logfile << "[PRESSURE](atm) = " << std::setprecision(10) << -P/CONVATM << endl;
    fid_logfile << "[DP](atm) = " << std::setprecision(10)<< -DP/CONVATM << endl;
    fid_logfile << "[BENDING_STIFNESS](pN nm) = " << std::setprecision(10) << KAPPA << endl;
    fid_logfile << "[BENDING_STIFNESS_MESH](pN nm) = " << std::setprecision(10) << LAMBDA << endl;
    fid_logfile << "[3D_YOUNG_MODULUS](MPa) = " << std::setprecision(10) << E_3D << endl;
    fid_logfile << "[2D_YOUNG_MODULUS](pN/nm) = " << std::setprecision(10) << Y_2D << endl;
    fid_logfile << "[SHELL_THICKNESS](pN/nm) = " << std::setprecision(10) << h << endl;
    fid_logfile << "[STRETCHING_STIFFNESS](pN/nm) = " << std::setprecision(10) << K_SPR << endl;   
    fid_logfile << "[MAX_ANGLE](rad) = " << std::setprecision(10) << MAX_THETA << endl;
    fid_logfile << "[BRAZIER_BEND_STIFFNESS_FACTOR] = " << std::setprecision(10) << BRAZIER_FACTOR << endl;
    
    
    
    //------------------------------------------------------------------------------------//
    
    
    
    //************************************************************************************//
    // 5. Non-Linear Conjugate Gradient Minimisations                                     //
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    
    cout << "Done" << endl << " - 3) MINIMIZATIONS RUNNING: " << endl; 
    fflush(stdout);
    
    double miniml0;                 // Length of the bonds after pressurisation (full pressure)
    double exec_time;               // Execution time
    time_t now;                     // Current date
    char *currt;                    // Current date in string format
    struct timeval t0, tf;          // OpenMP compatible calculation
    
    now = time(0);
    currt = ctime(&now);
    gettimeofday(&t0, NULL);
    
    
    // 5.1 -  Minimisation in ramp for pressure
    ofstream fid_pressure_list;
    std::string FULL_PATH_PRESS;                      // Path and file name as string including step
    
    fid_pressure_list.open(PRESSURES_FNAME, ofstream::out);
    caps.K_BRAZIER = 0.0;
    double tp=DP;
    int fcount=0;
    if(P < 0) {
        cout << "\r    - 3.1 Pressurization... ";
        fflush(stdout);
        while(tp > P) {
            cout << "\r    - 3.1 Pressurization... " << round(100.0*abs(tp/P)) << "%";
            fflush(stdout);
            
            shell.p = tp;
            nlcg(&shell, &caps, &minparams);  
            
            // Save the coordinates        
            FULL_PATH_PRESS = strcat(PATH_PRESS_CONFS, int2str(fcount));
            FULL_PATH_PRESS = strcat(FULL_PATH_PRESS, PRESS_COORDS_FNAME);
            save_coords_file(shell.r, shell.N, FULL_PATH_PRESS);
           
            // Saving list of pressures
            fid_pressure_list << fcount << "\t" <<  std::fixed << std::setprecision(10) << tp << "\t" << tp/CONVATM << endl;
            
            tp += DP;
            fcount++;
        }
        
        // Final pressure
        shell.p = P;
        nlcg(&shell, &caps, &minparams);  
    }
    
    // Save final configuration after pressurization
    FULL_PATH_PRESS = strcat(PATH_PRESS_CONFS, int2str(fcount));
    FULL_PATH_PRESS = strcat(FULL_PATH_PRESS, PRESS_COORDS_FNAME);
    save_coords_file(shell.r, shell.N, FULL_PATH_PRESS);
    save_coords_file(shell.r, shell.N, MCMIN_PRESS_COORDS_FNAME);
    fid_pressure_list << fcount << "\t" <<  std::fixed << std::setprecision(10) << P << "\t" << P/CONVATM << endl;
    fid_pressure_list.close();
    cout << "\r    - 3.1 Pressurization... " << "Done" << endl;
    fflush(stdout);
    
    miniml0 = avedge_length(shell.r, shell.edges, shell.nedges);
    fid_logfile << "[CHARACTERISTIC_REST_LENGTH_PRESS](nm) = " << miniml0 << endl;
    
    
    // 5.2  Minimisation in ramp for bending angle
    double dtheta = MAX_THETA/N_THETA;   
    fid_logfile << "[DELTA_ANGLE](rad) = " << dtheta << endl;
    

    ofstream fid_theta_angles;
    std::string FULL_PATH_BENT_TUBE;                            
   
    fid_theta_angles.open(BEND_ANGLES_FNAME, ofstream::out);
    
    caps.K_BRAZIER = K_SPR*shell.l0*shell.l0*BRAZIER_FACTOR;
    caps.K_BRAZIER_CM = K_SPR*BRAZIER_FACTOR_CM;
    
    for(int n=0; n<=N_THETA; n++) {
        
        cout <<  "\r    - 3.2 Tube bending... " << round(100.0*n/N_THETA) << "%"; 
        fflush(stdout);
        
        // Calculate the orientation vectors for the given imposed deformation
        double theta = n*dtheta;
        caps.nl0[0] = 0; 
        caps.nl0[1] = -cos(theta); 
        caps.nl0[2] = -sin(theta); 
        
        caps.nr0[0] = 0; 
        caps.nr0[1] = cos(theta); 
        caps.nr0[2] = -sin(theta); 
        
        nlcg(&shell, &caps, &minparams);
        
        // Save the coordinates
        FULL_PATH_BENT_TUBE = strcat(PATH_BENT_CONFS, int2str(n));
        FULL_PATH_BENT_TUBE = strcat(FULL_PATH_BENT_TUBE, BENT_COORDS_FNAME);
        save_coords_file(shell.r, shell.N, FULL_PATH_BENT_TUBE);
        
        //Save the angle used
        fid_theta_angles << n << "\t" << std::fixed << std::setprecision(10) << theta << endl;
    }    
    fid_theta_angles.close();
    
    
    gettimeofday(&tf, NULL);
    exec_time = ((tf.tv_sec  - t0.tv_sec) * 1000000u + tf.tv_usec - t0.tv_usec)/1.0e6;
    
    cout <<  "\r    - 3.2 Tube bending... Done" << endl;
    
    //-------------------------------------------------------------------------------------//
    
    
    
    //************************************************************************************//
    // 6. Save data                                                                       //
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    
    // Save the minimised coordiantes
    save_coords_file(shell.r, shell.N, MCMIN_COORDS_FNAME);
        
    fid_logfile << "[EXEC_DATE] = " << currt << endl;
    fid_logfile << "[EXEC_TIME](s) = " << exec_time << endl;
    fid_logfile.close();
    
    
    cout << "  Simulation finished" <<  endl << "  (execution time: " << exec_time << " seconds)" << endl << endl;
    
    
    //-------------------------------------------------------------------------------------//
    
    
    // Free memory
    free_dvector(shell.r);  
    free_dvector(shell.k_spr_eff);
    free_dvector(shell.eql);
    free_imatrix(shell.tri, shell.N_tri);
    free_imatrix(shell.edges, shell.nedges);
    
    // FREE STUFF IN LIDS
    
    return 0;
}




/*****************************************************************************************/
/*                                      calc_eql                                         */
/* This function calculates the rest length of every spring and sets on eql. Used to     */
/* determine the rest length of the springs                                              */
/*****************************************************************************************/
void calc_eql(double *r, int **edges, int nedges, double *eql)
{
    int v1,v2;
    double x,y,z;
    for(int i=0; i<nedges; i++) {
        v1 = edges[i][0];
        v2 = edges[i][1];
        x = r[v1*DIMV] - r[v2*DIMV];
        y = r[v1*DIMV + 1] - r[v2*DIMV + 1];
        z = r[v1*DIMV + 2] - r[v2*DIMV + 2];
        eql[i] = sqrt(x*x + y*y + z*z);
    }
}



/*****************************************************************************************/
/*                                    avedge_length                                      */
/* Calculates average rest length of the bonds in the configuration r                    */
/*****************************************************************************************/
double avedge_length(double *r, int **edges, int nedges)
{
    int v1, v2;
    double x, y, z;
    double L = 0;
    
    
    for(int i=0; i<nedges; i++) {
        v1 = edges[i][0];
        v2 = edges[i][1];
        x = r[v1*DIMV] - r[v2*DIMV];
        y = r[v1*DIMV + 1] - r[v2*DIMV + 1];
        z = r[v1*DIMV + 2] - r[v2*DIMV + 2];
        L += sqrt(x*x + y*y + z*z);
    }
    return L/nedges;
}





/*****************************************************************************************/
/*                                      print_gui                                        */
/* Prints the initial screen of indicating the parameters used                           */
/*****************************************************************************************/
void print_gui(int N, int T, double P_atm, double Kappa, double E3D, 
               double h, double K_SPR)
{
    
    cout << endl << "P R E S S U R I Z E D   M E S H : N L C G   M I N I M I S A T I O N" << endl; 
    cout << "---------------------------------------------------------------------" << endl << endl;
    cout << " - Particles: " <<  N << endl;
    cout << " - Triangles: " << T << endl;
    cout << " - Pressure (p): " <<  P_atm/CONVATM << " atm"<< endl;
    cout << " - Bending stiffness (Kappa): " <<  Kappa << " pN nm"<< endl;
    cout << " - 3D Young modulus (E): " <<  E3D << " MPa"<< endl;
    cout << " - 2D Young modulus (Y): " <<  E3D*h << " pN/nm"<< endl;
    cout << " - Shell thickness (h): " <<  h << " nm"<< endl;
    cout << " - Stretching stiffness (K_SPR): " <<  K_SPR << " pN/nm" << endl;
   
}



/*****************************************************************************************/
/*                                  save_coords_file                                     */
/* Saves the coordinates of the 2d array r with N rows in the file with file name fname  */
/*****************************************************************************************/
void save_coords_file(double *r, int N, std::string fname)
{
    ofstream fid_minim_coords;
    fid_minim_coords.open (fname, ofstream::out);
    for(int i = 0; i<N; i++) 
        fid_minim_coords <<  std::fixed << std::setprecision(10) << r[i*DIMV] << "\t" << r[i*DIMV + 1] << "\t" << r[i*DIMV + 2] << endl;
    fid_minim_coords.close();
}

