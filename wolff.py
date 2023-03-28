from utils import *
from metro import *

"""
# Wolff Algorithm
Similarly to the heat bath algorithm, Wolff is also a cluster-based Monte Carlo algorithm. 

The main difference is that the Heat Bath algorithm updates the spin configuration by considering each spin in the lattice one at a time, and calculating the probability of that spin being in the up or down state besed on its neighbors. Whereas Wolff algorithm updates the spin configuration by growing clusters of aligned spins. It starts with a random spin, and adds to the cluster all neighbors with the same spin. It then adds neighboring spins to the cluster with a probability determined by the temperature.

In theory, the Wolff algorithm is more efficient and faster at simulating critical phenomena, such as phase transitions, because it can update large clusters of spins at once, leading to faster equilibration times.
"""

# Define a function to add neighbors to the cluster
def add_neighbors(i, j, visited, cluster):
    """
    Recursively add neighboring spins to the cluster
    """
    # Mark the current spin as visited
    visited[i, j] = True
    cluster.append((i, j))

    # Check the neighboring spins
    neighbors = [(i+1,j), (i-1,j), (i,j+1), (i,j-1)]
    for ni, nj in neighbors:
        if ni < 0 or nj < 0 or ni >= N or nj >= N:
            continue  # skip neighbors that are outside the lattice
        if visited[ni, nj]:
            continue  # skip neighbors that are already in the cluster
        if lattice[ni, nj] != lattice[i, j]:
            continue  # skip neighbors with a different spin orientation
        if random.random() < 1 - np.exp(-2*J/T): #  add it to the cluster with probability e^(-2J/T)
            cluster = add_neighbors(ni, nj, visited, cluster)

    return cluster

def wolff_2d_eq(lattice, T, J, tolerance): #runs until it finds equilibrium
    print("Running Wolff in a 2D lattice.")
    
    #Plot initial alignment
    plt.imshow(lattice, cmap = 'gray')
    plt.show()
    
    #Print Hamiltonian
    print('Initial Hamiltonian is ' + str(energy2D(lattice,J)))
    
    #Print Average Magnetation
    print('Initial Avg Mag is ' + str(avg_mag_2d(lattice)))
    
    #intialize loop
    rep = 0 #Repetitions
    it = [] #iterations, used for plotting later on
    ham = [] #hamiltonian, used to find the equilibrium point
    
    count = 0 #Counter, used to determine if hamiltonian will be calculated or not
              #If we were to calc the Ham every iteration that'd be a huge time waste
    
    print_freq = 4 #frequency of which the ham will be printed (1 calc/N iterations)
                   #frequency will increase as the code runs exponentially
                   #Until it caps at 10,000. This is also to save run time
    
    delta_ham = tolerance + 1 #change in energy between two states of lattice
    
    ##We will say the code has reached equilibrium when the change
    ##in energy of the system is below a given tolerance
    
    

    #while delta_ham > tolerance:
    for i in range(5000):
        #randomly select the index of a lattice site
        # Choose a random spin to start the cluster
        i, j = random.randint(0, N-1), random.randint(0, N-1)

        # Initialize the cluster with the chosen spin
        cluster = [(i, j)]
        visited = np.zeros((N, N), dtype=bool)
        cluster = add_neighbors(i, j, visited, cluster)

        # Flip all spins in the cluster
        for i, j in cluster:
            lattice[i, j] *= -1

        # Compute the energy and magnetization of the lattice
        # E = energy(lattice)
        E = energy2D(lattice, J)
        M = np.sum(lattice)

        # Print out the results every 100 steps
        # if step % 100 == 0:
        #     print(f"Step {step}: Energy = {E}, Magnetization = {M}")
            
        # Increase counter 
        rep = rep + 1
        count = count + 1
            
        if count == print_freq:
            
            it.append(rep-1)

            ham.append(energy2D(lattice,J))
            
            if len(ham) > 1:
                delta_ham = abs(ham[len(ham)-1] - ham[len(ham)-2])
            
            #print(delta_ham)
            
            count = 0
            
            if print_freq < 1e5:
                print_freq = print_freq * 4

    plt.imshow(lattice, cmap = 'gray')
    plt.show()

    print('Final Hamiltonian is ' + str(energy2D(lattice,J)))

    print(f"Magnetization = {M}")
    
    print('The number of random number utilized is ' + str(3 * (rep-1)))
        
    return it, ham, avg_mag_2d(lattice)

# Running Wolff
J = 1
N = 18 # size of the lattice
rep = 0
T = 50 #Kelvin
tolerance = 1
k_B = 1.380649e-23 #Boltzman constant

lattice = np.random.choice((-1,1),(N,N))

[x,y,z] = wolff_2d_eq(lattice, T, J, tolerance)

plt.plot(x,y)

plt.xlabel('Iterations')
plt.ylabel('Hamiltonian')