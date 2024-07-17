import numpy as np
import sys



def save_log(tube_props):
    """
        Save a log file for the tube with the properties 
        of the dictionary tube_props
    """
    fileid = open('outfiles/rod.log', 'w')
    fileid.write('[MODE] = ' + tube_props["orientation"] + '\n')
    fileid.write('[RADIUS] = ' +  str(tube_props["radius"]) + '\n')
    fileid.write('[SPACING] = ' +  str(tube_props["spacing"]) + '\n')
    fileid.write('[LENGTH_INPUT] = ' + str(tube_props["length_input"]) + '\n')
    fileid.write('[LENGTH_FINAL] = ' +  str(tube_props["length_final"]) + '\n')
    fileid.write('[#PARTICLES_TOT] = ' +   str(tube_props["total_particles"]) + '\n')
    fileid.write('[#PARTICLES_BODY] = ' +  str(tube_props["particles_body"]) + '\n')
    fileid.write('[#HOOPS] = ' +  str(tube_props["number_hoops"]) + '\n')
    fileid.write('[#N_PER_HOOP] = ' + str(tube_props["particles_per_hoop"]) + '\n')
    fileid.close




def gencircumlids(R, L, spacing):
    """
    Generates a cylinder of radius R and length L with
    with circumferential alignment and lids 
    (equivalent to gencyl3lids)
    
    
    Input: R = Desired radius of the tube 
           L = Desired length of the tube
           spacing = (Approximate) bond length
    
    Output: np.array with xyz coordinates
    
    
    """
    
    # Initial checks
    N_theta = round(2.*np.pi*R/spacing)
    spacing = 2.*np.pi*R/N_theta
    
    Dt = 2*np.pi/N_theta
    Dz = spacing*np.sqrt(3.)/2.0
    N_hops = round(L/Dz)+1

    # To ensure that both ends are symmetric
    if N_hops % 2 ==0:
        N_hops += 1
    
    N_hops = int(N_hops)

    
    # Main body +++++++++++++++++++++++++++++++++++++++++++++
    x = np.zeros(N_theta*N_hops)
    y = np.zeros(N_theta*N_hops)
    z = np.zeros(N_theta*N_hops)
    
    h = 0;
    theta = 0;
    for i in range(0, N_hops):
        for t in range(0,N_theta):
            tt = theta + t*Dt
            x[i*N_theta + t] = R*np.cos(tt)
            y[i*N_theta + t] = R*np.sin(tt)
            z[i*N_theta + t] = h
    
        h += Dz;
        theta += Dt/2.0;

    rxyz = np.vstack((x.T, y.T, z.T)).T

    N_body = rxyz.shape[0]

    z_end = z[-1];
    nc = (2.0*np.pi*R)/spacing
    nz = round(2.0*L/(spacing*np.sqrt(3.0))  )
    
    #--------------------------------------------------------
    
    
    # Lids ++++++++++++++++++++++++++++++++++++++++++++++++++
    Dt0 = Dt
    min_radius = spacing*np.sqrt(3.)/2.
    theta = Dt0/2.
    
    # Upper lid
    R_lid = R - Dz
    while R_lid > min_radius:
        tn_theta = round(2.*np.pi*R_lid/spacing)
        Dt = 2*np.pi/tn_theta
        
        for t in range (0, tn_theta):
            tt = theta + (t+1)*Dt
            xx = R_lid*np.cos(tt)
            yy = R_lid*np.sin(tt)
            rxyz = np.append(rxyz, [[xx, yy, z_end]], axis=0)

        theta = Dt/2.0
        R_lid = R_lid - Dz

    rxyz = np.append(rxyz, [[0.0, 0.0, z_end]], axis=0)
    
    
    # Lower lid
    theta = Dt0/2
    R_lid = R - Dz
    nnn = 0
    while R_lid > min_radius:
        tn_theta = round(2.0*np.pi*R_lid/spacing)
        Dt = 2*np.pi/tn_theta
        
        for t in range (0, tn_theta):
            tt = theta + (t+1)*Dt
            xx = R_lid*np.cos(tt)
            yy = R_lid*np.sin(tt)
            rxyz = np.append(rxyz, [[xx, yy, 0]], axis=0)
        
        R_lid = R_lid - Dz
        theta = Dt/2;
        
    rxyz = np.append(rxyz, [[0.0, 0.0, 0.0]], axis=0)
    Nt = rxyz.shape[0]
    
    
    
    
    
    # Final checks ++++++++++++++++++++++++++++++++++++++++++
    # Rotation and offset
    Rx = np.array([ [1, 0, 0], 
                    [0, 0, 1],
                    [0, -1, 0] ])
    rxyz_rot = np.matmul(Rx, rxyz.T).T    
    dy = max(rxyz_rot[:,1]) - min(rxyz_rot[:,1])
    rxyz_rot[:,1] = rxyz_rot[:,1] - dy/2.0;


    tube_props = {
        "orientation": "CIRCUMFERENTIAL",
        "spacing": spacing,
        "radius": R,
        "length_input": L,
        "length_final": N_hops*spacing*np.sqrt(3.0)/2.0,
        "total_particles": Nt,
        "particles_body": N_body,
        "number_hoops": N_hops,
        "particles_per_hoop": N_theta
    }

    return np.array(rxyz_rot), tube_props





def genlonglids(R, L, spacing):
    """
    Generates a cylinder of radius R and length L with
    with long axis alignment and lids 
    
    Input: R = Desired radius of the tube 
           L = Desired length of the tube
           spacing = (Approximate) bond length
    
    Output: np.array with xyz coordinates
    
    """
    # Initial checks
    lc = spacing*np.sqrt(3.)/2.
    N_theta = int(round(2.*np.pi*R/lc))
    if N_theta % 2 !=0:
        N_theta += 1
        
    lc = 2.*np.pi*R/N_theta
    spacing = 2.0*lc/np.sqrt(3.0)
 
    Dt = 2.0*np.pi/N_theta

    Dz = spacing
    N_hops = round(L/Dz)+1
    if N_hops % 2 ==0: # To ensure that both ends are symmetric
        N_hops += 1
    
    N_hops = int(N_hops)

   
    # Main body +++++++++++++++++++++++++++++++++++++++++++++
    x = np.zeros(N_theta*N_hops)
    y = np.zeros(N_theta*N_hops)
    z = np.zeros(N_theta*N_hops)
    
    for i in range(0, N_hops):
        for t in range(0,N_theta):
            tt = t*Dt
            x[i*N_theta + t] = R*np.cos(tt)
            y[i*N_theta + t] = R*np.sin(tt)
            z[i*N_theta + t] = Dz*(i +  (t%2.0)/2.0)
    

    rxyz = np.vstack((x.T, y.T, z.T)).T
    N_body = rxyz.shape[0]
    z_max = max(rxyz[:,2])
    z_min = min(rxyz[:,2])
     
    # Lids ++++++++++++++++++++++++++++++++++++++++++++++++++
    Dt0 = Dt
    min_radius = spacing*np.sqrt(3.)/2.
    theta = Dt0/2.
    
    # Upper lid
    R_lid = R - lc
    while R_lid > min_radius:
        tn_theta = round(2.*np.pi*R_lid/spacing)
        Dt = 2*np.pi/tn_theta;
        
        for t in range (0, tn_theta):
            tt = theta + (t+1)*Dt
            xx = R_lid*np.cos(tt)
            yy = R_lid*np.sin(tt)
            rxyz = np.append(rxyz, [[xx, yy, z_max]], axis=0)

        theta = Dt/2.0
        R_lid = R_lid - lc

    rxyz = np.append(rxyz, [[0.0, 0.0, z_max]], axis=0)
    
    
    # Lower lid
    R_lid = R - lc
    while R_lid > min_radius:
        tn_theta = round(2.*np.pi*R_lid/spacing)
        Dt = 2*np.pi/tn_theta;
        
        for t in range (0, tn_theta):
            tt = theta + (t+1)*Dt
            xx = R_lid*np.cos(tt)
            yy = R_lid*np.sin(tt)
            rxyz = np.append(rxyz, [[xx, yy, z_min]], axis=0)

        theta = Dt/2.0
        R_lid = R_lid - lc

    rxyz = np.append(rxyz, [[0.0, 0.0, z_min]], axis=0)
    
    
    
    # Incorporate the last points to close the rings
    # Upper ring
    Dt = 2*np.pi/N_theta
    for t in range(0,N_theta):
        if (t%2.0)/2.0 == 0:
            tt = t*Dt
            xx = R*np.cos(tt)
            yy = R*np.sin(tt)
            rxyz = np.append(rxyz, [[xx, yy, z_max]], axis=0)
    
    # Lower lid
    for t in range(0,N_theta):
        if (t%2.0)/2.0 != 0:
            tt = t*Dt
            xx = R*np.cos(tt)
            yy = R*np.sin(tt)
            rxyz = np.append(rxyz, [[xx, yy, z_min]], axis=0)
    
    
    Nt = rxyz.shape[0]
    
    
    
    
    
    # Final checks ++++++++++++++++++++++++++++++++++++++++++
    # Rotation and offset
    Rx = np.array([ [1, 0, 0], 
                    [0, 0, 1],
                    [0, -1, 0] ])
    rxyz_rot = np.matmul(Rx, rxyz.T).T    
    dy = max(rxyz_rot[:,1]) - min(rxyz_rot[:,1])
    rxyz_rot[:,1] = rxyz_rot[:,1] - dy/2.0;



    tube_props = {
        "orientation": "LONG",
        "spacing": spacing,
        "radius": R,
        "length_input": L,
        "length_final": dy,
        "total_particles": Nt,
        "particles_body": N_body,
        "number_hoops": N_hops,
        "particles_per_hoop": N_theta
    }

    return rxyz_rot, tube_props
    
    


if __name__ == "__main__":
    
    if len(sys.argv) != 5:
        print("Not enought arguments. Execute as rodgen radius length spacing type_of_cylinder")
        print("type_of_cylinder can be: circumf or long")
        sys.exit(0)
        
    R = np.float64(sys.argv[1])
    L = np.float64(sys.argv[2])
    spacing = np.float64(sys.argv[3])
    mode =  sys.argv[4]
    
    
    
    if mode == "circumf":
        rxyz, tube_props = gencircumlids(R, L, spacing)
    elif mode=="long":
        rxyz, tube_props =  genlonglids(R, L, spacing)
    else:
        print("The mode provided is not valid")
        sys.exit(0)
        
        
    
    
    np.savetxt('./outfiles/init_coords.dat', rxyz,fmt='%.10f', delimiter='\t')
    save_log(tube_props)
    
    """
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.scatter3D(rxyz[:,0], rxyz[:,1], rxyz[:,2])
    plt.show()
    """
    

    from meshgen import gentrimesh, plot_tri
    tri, n_tri = gentrimesh(rxyz)
    np.savetxt('./outfiles/mesh.dat', tri, fmt='%d', delimiter='\t')
    plot_tri(rxyz, tri)
    
    """
    # This was used to test the correct geometry od the bonds
    for t in range(0, n_tri):
        v1 = tri[t,0]
        v2 = tri[t,1]
        v3 = tri[t,2]
        l1 = np.linalg.norm(rxyz[v1,:] - rxyz[v2,:])
        l2 = np.linalg.norm(rxyz[v2,:] - rxyz[v3,:])
        l3 = np.linalg.norm(rxyz[v3,:] - rxyz[v1,:])
        
        print(rxyz[v1,:])
        print(rxyz[v2,:])
        print(rxyz[v3,:])
        print("  ")
        print(l1)
        print(l2)
        print(l3)
        print("===============")
        
    """
