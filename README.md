# Mechanics and wrinkling patterns of pressurized bent tubes

Cesar L. Pastrana, 2023 (c)

This software was used in the manuscript "Mechanics and wrinkling patterns of pressurized bent tubes" by Pastrana C.L, Qiu L., Hutchinson J.W., Amir A. and Gerland U.
[arXiv](https://arxiv.org/abs/2402.17033) 2024.

## Description
The surface is described using a discrete bead and spring model. The total energy of the system is given by $E_{T} = E_{int} + E_{ext}$. The internal energy of the system is that of the tube and is given by:

$$
E_{int} = \frac{k_s}{2}\sum_{\langle i,j\rangle} (r_{ij} - r_{ij}^0)^2 +  k_b\sum_{\langle \alpha,\beta\rangle}(1 - \hat{n_\alpha}\cdot\hat{n}_\beta) - pV
$$

where the first sum runs over all pairs of connected notes $i,j$ constituting an edge of the mesh and $r_{ij} = |r_i - r_j|$; and the second sum is over pair of triangles $\alpha,\beta$ sharing an edge, $\hat{n}$ are their respective normal vectors and $k_b$ is a bending stiffness. For a detailed description of the relation between the discrete $k_x$ and the continum variables (Young modulus, thickness and Poisson ratio) see the original manuscript as well as the seminal work by Soung H.S. and Nelson D.R. \([Phys. Rev. A 38,1005--1018, 1988](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.38.1005)\).


The external potential is $E_{ext} = E_{\tau} + E_{CM}$. The component $E_\tau$ is included to induce the bending of the tube and is given by:

$$
E_\tau = k_{\tau}\sum_{ i\in \text{ right lids} } (1 - \hat{n}_i\cdot\hat{n}^l_0) + k_{\tau}\sum_{ i \in \text{ right lid} } (1 - \hat{n}_i\cdot\hat{n}_0^r)
$$

where the summation runs over all the triangles of the right (or left) lid and $\hat{n}_i$ is the normal vector of the triangle $i$. The term $\hat{n}^{l,r}_0$ is the preferred orientation vector for the right or left lid, $r$ and $l$ respectively.

The term $E_{CM}$ is used to enforce the alignment of the two lids in the plane normal to the long axis of the tube

$$
E_{CM} = \frac{k_{CM}}{2} |\vec{R}_{{CM},x}^l  - \vec{R}_{{CM},x}^r|^2  + \frac{k_{CM}}{2} |\vec{R}_{{CM},z}^l  - \vec{R}_{{CM},z}^r|^2,
$$

where $\vec{R}_{{CM},x}^{l,r}$ is the center of mass of the left/right lid and $x,z$ indicates the Cartesian coordinate. We found that this potential has no influence in the mechanical and geometrical response. However, in its absence a single kink is formed anywhere along the distance of the tube.

The approach described resembles the action of virtual hands to induce bending. The force exerted on the lids is adjusted with the constant $k_{\tau}$, given by $k_\tau := K_\tau k_{s}l_0^2$, where $l_0$ is the average bond length of the mesh. $K_\tau$ is a free parameter determined empirically to achieve the target curvatures of the tube. Similarly, $k_{CM} := k_s K_{CM}$. 

Each step, pressurization and bending, is followed by energy minimization. A custom non-linear conjugate gradient algorithm is used to find the minimal energy configuration.


## Compilation
A Makefile is provided in the `brazier` folder. To compile the code, simply execute:

```
make parallel
```
 This will generate the compiled code `brazier` in the `bin` folder. The `parallel` tag activates the OpenMP compatible compilation.

## Execution
The simulation is launched by simply executing `./brazier` in the `bin` folder. The number of nodes to use in the parallel computations can be configured by exporting the variable ``OMP_NUM_THREADS`, as `export OMP_NUM_THREADS = num threads to use`.

Pressurisation proceeds in a step-wise fashion from $0$ to a target pressure $p$ atm, steps of $\Delta p$ provided by the user. Once the target pressurisation is reached, the bending of the tube starts in steps, for a total number of $N_\theta$ steps. There is no control on the particular angle aquired, since it will also depend on the parameter $k_\tau$. Thus, some experimentation is needed to reach the target curvature. This is a tradeoff between computational speed and number of data points for each curvature.


### Input
The program requires as input the following files:

- `init_coords.dat`: Array of $N\times 3$ with the $N$ vertices coordinates defining the surface in the relax configuration. The units of the input coordiantes coordinates are nm.
- `mesh.dat`: Array of $T\times 3$ with the $T$ triangles defining the connectivity of the triangulated surface.
- `params.conf`: This file contains the simulation parameters to use. The meaning of each parameter is indicated in commented blocks in the file.


#### Generation of cylinders

A cylinder with lids having the approapiate configuration to be interpreted by the simulation code, can be generated with the included Python code included in the folder `gencylinder`. The code is executed as
```
python3 rodgen radius length spacing type_of_cylinder
```
where radius and length are the target radius and length, and spacing is the typical bond length of the sytsem, i.e., $r_{ij}^0$. Type of cylinder can take the values `circumf` or `long` depending on the desired orientation of the mesh. In the article, we mostly used circumferentially aligned cylinders. The execution of the script generates the files `init_coords.dat` and `mesh.dat` located in the `outfiles` folder. These files contain the position of the vertices and the ensuing triangulation and are ready to be used with `brazier`. A log file is also generated as a summary of the properties of the generated mesh.

The codes to generate the cylinders with lids are dependent on Python 3 and on the packages numpy and sys. The generation of the mesh relies on the Advancing Front Surface Reconstruction algorithm from the CGAL library (see the [meshgen](https://github.com/pastranacl/meshgen) repo).



### Output
The following output files are generated after execution of `brazier`:

- `press/X_press.dat`: Coordinates for each step of pressurization.
- `brazier/X_brazier.dat`: Coordinates for each step of bending at the given pressure. These are the main files to analyse.
- `summary/minim_coords_press.dat`:  Coordinate at the target pressure prior to bending.
- `summary/minim_coords_press.dat`:  Coordinate at the target pressure at the end of the bending routine.
- `summary/pressures.dat`: Pressures in each step.
- `summary/bent_angles.dat`: Imposed bent angles for each step.
- `summary/output.dat`: Summary of mechanical, geometrical, and simulation parameters employed.
- `lids/*.dat`: The files in this folder indicate the edges, triangles, and vertices that belong to the right or left lid of the tube.
