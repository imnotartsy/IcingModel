from utils import *
from metro import *

"""
# Glauber Algorithm
The Glauber algorithm is another Monte Carlo algorithm that can be used to simulate the behavior of a system of interacting spins, such as the Ising model. The algorithm used is as follow:

1. Choose a random particle.
2. Sum its four neighboring spins
3. Compute the change in energy if the spin x, y were to flip using the Hamiltonian.
4. Flip the spin with probability $e^{- \frac{\Delta E}{T}} / (1 + e^{- \frac{\Delta E}{T}})$
"""

def Glauber_2d_eq(lattice, T, J, tolerance): #runs until it finds equilibrium
    print("Running Glauber in a 2D lattice.")
    
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
    
    while delta_ham > tolerance:
        #randomly select the index of a lattice site
        i = np.random.randint(0,N)
        j = np.random.randint(0,N)
        
        #Boundary conditions to find edges
        
        #2D arrays have neighbors to the left,
        #right, above, and below
        
        #We define 'above/up' as the lattice site in
        #the same column but preceding row
        
        #below/down is same column but following row
        
        #Left is same row but preceding column
        
        #Right is same row but following column

        if j == 0:
            left = 0
        else:
            left = lattice[i,j]*lattice[i,j-1]

        if j == (N-1):
            right = 0
        elif j < (N-1):
            right = lattice[i,j]*lattice[i,j+1]

        if i == 0:
            up = 0
        else:
            up = lattice[i,j]*lattice[i-1,j]

        if i == (N-1):
            down = 0
        elif i < (N-1):
            down = lattice[i,j]*lattice[i+1,j]

        neighbors = up + down + left + right
        
        #calculate intial energy

        E_i = -1*J * neighbors
        
        #flip spin

        lattice[i,j] = -1 * lattice[i,j]
        
        #same conditions as before

        if j == 0:
            left = 0
        else:
            left = lattice[i,j]*lattice[i,j-1]

        if j == (N-1):
            right = 0
        elif j < (N-1):
            right = lattice[i,j]*lattice[i,j+1]

        if i == 0:
            up = 0
        else:
            up = lattice[i,j]*lattice[i-1,j]

        if i == (N-1):
            down = 0
        elif i < (N-1):
            down = lattice[i,j]*lattice[i+1,j]

        neighbors = up + down + left + right
        
        #calculate new energy energy_final

        E_f = -1* J * neighbors
        
        #condition, the higher the temp the more likely
        #the particle is to flip even if it's a higher
        #enegy state
        
        expon =  (-1*(E_f - E_i)) / T # -delta E / T  
        value = np.exp(expon) / (1+ np.exp(expon))

        #Glauber
        if E_f < E_i or np.random.rand(1)[0] < value:
            # Not flipped
            pass
        else: 
            #Flip
            lattice[i,j] = -1 * lattice[i,j]
            
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

    print('Final Avg Mag is ' + str(avg_mag_2d(lattice)))
    
    print('The number of random number utilized is ' + str(3 * (rep-1)))
        
    return it, ham, avg_mag_2d(lattice)

"""## Comparison with Metropolis

Glauber and Metropolis are both types of Markov Chain Monte Carlo algorithms and, as such, should give identical final results. 
In our tests, however, Metropolis was faster. Especially if Metropolis selects a spin deterministically (such as following a specific order) while Glauber selects each spin randomly with equal probability.
"""

# Running Metropolis
J = 1
N = 16
rep = 0
T = 1 #Kelvin
tolerance = 8
k_B = 1.380649e-23 #Boltzmann constant

lattice = np.random.choice((-1,1),(N,N))

[x,y,z] = Metro_2d_eq(lattice, T, J, tolerance)

plt.plot(x,y)

plt.xlabel('Iterations')
plt.ylabel('Hamiltonian')

# Running Glauber
J = 1
N = 16
rep = 0
T = 1 #Kelvin
tolerance = 8

lattice = np.random.choice((-1,1),(N,N))

[x,y,z] = Glauber_2d_eq(lattice, T, J, tolerance)

plt.plot(x,y)

plt.xlabel('Iterations')
plt.ylabel('Hamiltonian')

"""Over our tests, the average iterations per algorithm are:

- Glauber: 349,523
- Metropolis: 8,7379

The main difference, besides choosing a spin deterministically vs randomly, is the acceptance criterion for flipping a spin. Metropolis includes the Boltzmann constant lowering the impact that temperature has in the simulation. 

Because of this difference in the acception criterion, at thermal equilibrium, these two algorithms give identical results. However, at higher temperatures Glauber takes much longer to run.
"""