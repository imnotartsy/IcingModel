from utils import *
from metro import *

"""# Cluster Heat Bath

A totally different way to think about solving the Ising model is to instead approach the system as a series of clusters interacting together rather than just a bunch of isolated particles. We'll attempt to identify and treat groups of lattice sites in unison rather than randomly hopping from site to site to site over and over again. Hopefully we'll be able to flip multiple lattice sites in one "sweep" of the lattice to create a system that is in a much lower energy state. The basic princi
"""
### Two Dimensions


J = 1
T = 1/k_B * 0.1 #Kelvin
N = 64

lattice = np.random.choice((-1,1),(N,N))

plt.imshow(lattice, cmap = 'gray')
plt.show()

#print(lattice)

def add_to_cluster_2D(lattice, i, j, cluster, T, J):
    """Adds the spin at position (i, j) to the cluster with probability p
    
    INPUT:
    
    lattice: two dimensional numpy array with up (1) or down (-1) spins
    
    i: row index in lattice
    
    j: column index in lattice
    
    T: temperature in Kelvin
    
    J: Magnetic Constant (set to 1 if you don't care)    
    
    OUTPUT:
    
    cluster: boolean array the same size of the lattice array designating
    which sites are a part of the cluster"""
    
    N = lattice.shape[0]
    
    #condition to see if the cluster will expand
    expand = False
    
    E_i = 1 * J * get_neighbors_2D(lattice,i,j)
    
    lattice[i,j] = -1 * lattice[i,j]
    
    E_f = 1 * J * get_neighbors_2D(lattice,i,j)
    
    delta_E = E_f - E_i
    
    lattice[i,j] = -1 * lattice[i,j]
    
    #probability condition for flipping
    p = (np.exp(-delta_E / (k_B*T)))
    
    if np.random.rand() < p:
        expand = True
    
    if cluster[i, j] == False and expand == True:
        cluster[i, j] = True
        if i > 0:
            add_to_cluster_2D(lattice, i-1, j, cluster, T, J)
        if i < N-1:
            add_to_cluster_2D(lattice, i+1, j, cluster, T, J)
        if j > 0:
            add_to_cluster_2D(lattice, i, j-1, cluster, T, J)
        if j < N-1:
            add_to_cluster_2D(lattice, i, j+1, cluster, T, J)
    
    return cluster

def cluster_heat_bath_2D(lattice, T, J):
    """Performs one update of the 2D Ising model using the Cluster Heat Bath algorithm
    
    INPUTS:
    
    lattice: Two dimensional numpy array of up (1) or down (-1) spins
    
    T: temperature in Kelvin
    
    J: Magnetic constant (set to 1 if you don't care)
    
    RETURNS:
    
    lattice: altered lattice from one iteration of CHB"""
    
    N = lattice.shape[0]

    cluster = np.zeros((N, N), dtype=bool)

    i, j = np.random.randint(0, N, 2)

    #create cluster
    add_to_cluster_2D(lattice, i, j, cluster, T, J)

    #flip cluster
    lattice[cluster] *= -1

    i = i + 1
    
    return lattice

i = 0
while i < 1e3:
    cluster_heat_bath_2D(lattice, T, J)
    i = i + 1

plt.imshow(lattice, cmap = 'gray')

## Let's compare this with metro algorithm using the same number of iterations

J = 1
T = 1/k_B * 0.1 #Kelvin
N = 64

lattice = np.random.choice((-1,1), (N,N))

Metro_2D(lattice,1e3,T,J)

"""As we can see the Cluster analysis can get more done in less iterations because it sweeps over the entire code. Of course it take more run time to do 1000 iterations of the cluster analysis than metro but we do get much much more order. A strange pattern that I fear may be a bug are these checkered "tendrils" throughout the cluster that ruin the nice clean blobby aesthetic.

### Two Dimensions - Running Until Equilibrium
"""

J = 1
T = 1/k_B * 1e-10
N = 64
tolerance = 0

lattice = np.random.choice((-1,1),(N,N), p = (0.7,0.3))

plt.imshow(lattice, cmap = 'gray')

plt.show()

def cluster_bath_eq_2D(lattice, T, J, tolerance):
    """Run Cluster Heat Bath algorithm until the lattice reaches some equilibrium
       within some tolerance.
       
       INPUTS:
       
       lattice: two dimensional numpy array of up (1) or down (-1) spins
       
       T: temperature in Kelvin
       
       J: Magnetic constant (set to 1 if you don't care)
       
       tolerance: largest change in hamiltonian that is still considered
       within equilibrium. The lower the number the more relaxed the system will
       get but run times may increase dramatically
       
       RETURNS:
       
       lattice: altered lattice from multiple iterations of CHB """
    
    Ham_initial = energy2D(lattice, J)
    
    delta_E = tolerance + 1
    
    while delta_E > tolerance:
        
        i = 0
        while i < 10:
            lattice = cluster_heat_bath_2D(lattice, T, J)
            i = i + 1
        
        Ham_final = energy2D(lattice, J)
        
        delta_E = abs(Ham_final - Ham_initial)
        
        Ham_initial = Ham_final
    
    return lattice

lattice = cluster_bath_eq_2D(lattice, T, J, tolerance)

plt.imshow(lattice, cmap = 'gray')

avg_mag_2D(lattice)

"""### Three Dimensions"""

J = 1
T = 1/k_B * 0.1 #Kelvin
N = 32

lattice = np.random.choice((-1,1),(N,N,N))

plt.imshow(lattice[1], cmap = 'gray')
plt.show()

#print(lattice)

def add_to_cluster_3D(lattice, k, i, j, cluster, T, J):
    """Adds the spin at position (i, j) to the cluster with probability p
    
    INPUT:
    
    lattice: two dimensional numpy array with up (1) or down (-1) spins
    
    k: layer index in lattice
    
    i: row index in lattice
    
    j: column index in lattice
    
    cluster: boolean array identical in size to lattice filled with False
    
    T: temperature in Kelvin
    
    J: Magnetic Constant (set to 1 if you don't care)
    
    OUTPUT:
    
    cluster: boolean array the same size of the lattice array designating
    which sites are a part of the cluster"""
    
    N = lattice.shape[0]
    
    expand = False
    
    E_i = 1 * J * get_neighbors_3D(lattice,k,i,j)
    
    lattice[k,i,j] = -1 * lattice[k,i,j]
    
    E_f = 1 * J * get_neighbors_3D(lattice,k,i,j)
    
    delta_E = E_f - E_i
    
    lattice[k,i,j] = -1 * lattice[k,i,j]
    
    p = (np.exp(-delta_E / (k_B*T)))
    
    if np.random.rand() < p:
        expand = True
    
    else:
        expand = False
    
    if cluster[k, i, j] == False and expand == True:
        cluster[k, i, j] = True
        if i > 0:
            add_to_cluster_3D(lattice, k, i-1, j, cluster, T, J)
        if i < N-1:
            add_to_cluster_3D(lattice, k, i+1, j, cluster, T, J)
        if j > 0:
            add_to_cluster_3D(lattice, k, i, j-1, cluster, T, J)
        if j < N-1:
            add_to_cluster_3D(lattice, k, i, j+1, cluster, T, J)
        if k > 0:
            add_to_cluster_3D(lattice, k-1, i, j, cluster, T, J)
        if k < N-1:
            add_to_cluster_3D(lattice, k+1, i, j, cluster, T, J)
    
    return cluster



def cluster_heat_bath_3D(lattice, T, J):
    """Performs one update of the 3D Ising model using the Cluster Heat Bath algorithm
    
    INPUTS:
    
    lattice: Three dimensional numpy array of up (1) or down (-1) spins
    
    T: temperature in Kelvin
    
    J: Magnetic constant (set to 1 if you don't care)
    
    RETURNS:
    
    lattice: altered lattice from one iteration of CHB"""
    
    N = lattice.shape[0]

    cluster = np.zeros((N, N, N), dtype=bool)

    k, i, j = np.random.randint(0, N, 3)

    add_to_cluster_3D(lattice, k, i, j, cluster, T, J)

    lattice[cluster] *= -1
    
    return lattice

i = 0
while i < 1e2:
    cluster_heat_bath_3D(lattice, T, J)
    i = i + 1
    
plt.imshow(lattice[1], cmap = 'gray')

J = 1
T = 1/k_B * 0.1 #Kelvin
N = 32

lattice = np.random.choice((-1,1), (N,N))

Metro_3D(N,1e2,T,J)

"""### Three Dimensions -  Running Until Equilibrium"""

J = 1
T = 1/k_B * 10
N = 16
tolerance = 1

lattice = np.random.choice((-1,1), (N,N,N), p = (0.75,0.25))

plt.imshow(lattice[1], cmap = 'gray')
plt.show()

def cluster_bath_eq_3D(lattice, T, J, tolerance):
    """Run Cluster Heat Bath algorithm until the lattice reaches some equilibrium
       within some tolerance.
       
       INPUTS:
       
       lattice: two dimensional numpy array of up (1) or down (-1) spins
       
       T: temperature in Kelvin
       
       J: Magnetic constant (set to 1 if you don't care)
       
       tolerance: largest change in hamiltonian that is still considered
       within equilibrium. The lower the number the more relaxed the system will
       get but run times may increase dramatically
       
       RETURNS:
       
       lattice: altered lattice from multiple iterations of CHB """
    
    Ham_initial = energy3D(lattice, J)
    
    delta_E = tolerance + 1
    
    
    while delta_E > tolerance:
        
        i = 0
        while i < 10:
            lattice = cluster_heat_bath_3D(lattice, T, J)
            i = i + 1
        
        Ham_final = energy3D(lattice, J)
        
        delta_E = abs(Ham_final - Ham_initial)
        
        Ham_initial = Ham_final
    
    return lattice

lattice = cluster_bath_eq_3D(lattice, T, J, tolerance)
    
plt.imshow(lattice[1], cmap = 'gray')

avg_mag_3D(lattice)

"""## Comparing the Runtime of Metro with CHB

We previously saw that when running the CHB and Metro algorithms for the same number of iterations the CHB function evolved much more. This is expected for an algorithm that considers multiple lattice sites in one iteration, but is it worth it? Does this extra work of building a cluster outweigh any potential benefit. Here we examine the wall time of each algorithm in two dimensions and three dimensions for square lattices of various lengths. We will run each lattice until equilibrium and check to see how much time has elapsed in the process of running each algorithm and plot the results.

### Two Dimensions
"""

D_2 = [4,8,16,32,64]

tolerance = 8

wall_time_2D_chb = []

i = 0

while i < len(D_2):
    lattice = np.random.choice((-1,1), (D_2[i],D_2[i]))

    start_time = time.time()

    cluster_bath_eq_2D(lattice, 1e20, J, tolerance)

    end_time = time.time()

    tot_time = end_time - start_time
    
    wall_time_2D_chb.append(tot_time)
    
    i = i + 1
    
    
D_2 = [4,8,16,32,64]

wall_time_2D_metro = []

i = 0

while i < len(D_2):
    lattice = np.random.choice((-1,1), (D_2[i],D_2[i]))

    start_time = time.time()

    Metro_2D_eq(lattice, 1e20, J, tolerance)[0]

    end_time = time.time()

    tot_time = end_time - start_time
    
    wall_time_2D_metro.append(tot_time)
    
    i = i + 1

plt.plot(D_2, wall_time_2D_chb, color = 'red', label = 'CHB')

plt.plot(D_2, wall_time_2D_metro, color = 'green', label = 'Metro')

plt.title('Comparison of Runtime of Metro and CHB in 2D')
plt.xlabel('Array Length (N)')
plt.ylabel('Runtime (s)')
plt.xlim(min(D_2),max(D_2))
plt.legend()

"""### Three Dimensions"""

D_3 = [2,4,8,16,32,64]

tolerance = 32

wall_time_3D_chb = []

i = 0

while i < len(D_3):
    lattice = np.random.choice((-1,1), (D_3[i],D_3[i],D_3[i]))

    start_time = time.time()

    cluster_bath_eq_3D(lattice, 1e20, J, 2)

    end_time = time.time()

    tot_time = end_time - start_time
    
    wall_time_3D_chb.append(tot_time)
    
    i = i + 1
    
    
N = [2,4,8,16,32,64]

wall_time_3D_metro = []

i = 0

while i < len(D_3):
    lattice = np.random.choice((-1,1), (D_3[i], D_3[i], D_3[i]))

    start_time = time.time()

    Metro_3D_eq(lattice, 1e20, J, 2)[0]

    end_time = time.time()

    tot_time = end_time - start_time
    
    wall_time_3D_metro.append(tot_time)
    
    i = i + 1

plt.plot(D_3, wall_time_3D_chb, color = 'red', label = 'CHB')

plt.plot(D_3, wall_time_3D_metro, color = 'green', label = 'Metro')

plt.title('Comparison of Runtime of Metro and CHB in 3D')
plt.xlabel('Array Length (N)')
plt.ylabel('Runtime (s)')
plt.xlim(min(D_3),max(D_3))
plt.legend()

"""Finally we see that CHB is a superior algorithm that gives more advantage in comparison to Metro the larger the lattice you begin with. One can image the cluster algorithm being more advantagous to flip tens or hundreds of spins at once to massively drop the energy of the lattice at once while the Metro algorithm has to laborously and blindly switch individual spins.
