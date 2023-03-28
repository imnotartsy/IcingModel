from utils import *



"""# Metropolis Implementations and Analysis

## Metro: Barebones; Time + Space Conscious 1D & 2D
"""

def Metro1D(lat, J, T, iterations=200):
  for _ in range(iterations):
    site = np.random.randint(0, len(lat))

    if site == 0:
      delta = -1 * J * (-1* lat[site] * lat[site+1])
    elif site == len(lat)-1:
      delta = -1 * J * (lat[site-1] * -1* lat[site])
    else:
      delta = -1 * J * (lat[site-1] * -1* lat[site] + -1* lat[site] * lat[site+1])

    if delta < 0 or np.random.randint(0, 1) < np.e**(-delta/(sp.constants.k*T)):
      lat[site] *= -1

  return lat


def Metro2D(lat, J, T, iterations=200):
  for _ in range(iterations):
    site_x = np.random.randint(0, len(lat))
    site_y = np.random.randint(0, len(lat[0]))

    if site_x == 0:
      delta_x = -1 * J * (-1* lat[site_x][site_y] * lat[site_x+1][site_y])
    elif site_x == len(lat)-1:
      delta_x = -1 * J * (lat[site_x-1][site_y] * -1* lat[site_x][site_y])
    else:
      delta_x = -1 * J * (lat[site_x-1][site_y] * -1* lat[site_x][site_y] + -1* lat[site_x][site_y] * lat[site_x+1][site_y])

    if site_y == 0:
      delta_y = -1 * J * (-1* lat[site_x][site_y] * lat[site_x][site_y+1])
    elif site_y == len(lat[0])-1:
      delta_y = -1 * J * (lat[site_x][site_y-1] * -1* lat[site_x][site_y])
    else:
      delta_y = -1 * J * (lat[site_x][site_y-1] * -1* lat[site_x][site_y] + -1* lat[site_x][site_y] * lat[site_x][site_y+1])


    if delta_x + delta_y < 0 or np.random.randint(0, 1) < np.e**(-(delta_x + delta_y )/(sp.constants.k*T)):
      lat[site_x][site_y] *= -1

  return lat

"""### Metro Examples"""

# Ajacency Matrix example
#2D Adjacency Matrix Example

lat = gen2DLat(5)
adjacency = gen2DAdjacency(lat)
print2DAdjacency(lat, adjacency)
plt.imshow(adjacency, cmap='gray')

# 1D Example
N = 15
J = 1
T = 1

lat = gen1DLat(N)

print1D(lat)
init_e = energy1D(lat, J)
print("\tInitial Energy:", init_e)

lat_out = Metro1D(lat, J, T, iterations=200)
print1D(lat_out)
print("\tOutput Energy", energy1D(lat_out, J))

# 2D example
N = 4
J = 1
T = 1

lat = gen2DLat(N)
print2D(lat)


init_e = energy2D(lat, J)
print("\tInitial Energy:",init_e)



lat_out = Metro2D(lat, J, T)
print2D(lat_out)
print("\tOutput Energy", energy2D(lat_out, J))

"""## Metro: Running N iterations + Animation

### One Dimension
"""

N = 64
J = 1
T = 1e-10 #Kelvin #kelvin

    #Create Lattice
lattice = np.random.choice((-1,1),N)

#Print Initial Hamiltonian
print('Initial Hamiltonian is ' + str(energy1D(lattice, J)))
#Print Initial Average Mag
print('Initial Avg Mag is ' + str(avg_mag(lattice)))


def Metro_1D(lattice, iterations, T, J):
    """Performs N iterations on one dimensional lattice using Metropolitan Algorithm
    
    INPUTS:
    
    lattice: one dimensional numpy array of up (1) or down (-1) spins
    
    iterations: number of desired repitions of Metro function.
    One repition = one flip at one site
    
    T = temperature in kelvin
    
    J = constant, set to 1 if you don't care
    
    RETURNS:
    
    Numpy array of the altered lattice"""
    
    #intialize loop   
    rep = 0 #Repetitions
    N = len(lattice)
    while rep < iterations:
        #Select random lattice site
        i = np.random.randint(0,N)
        #Boundary condtions: In one dim
        #There are neighbors to the right
        #and left
        
        #if i == 0 there is only a one to the right
        
        #if i == N-1 there's only one to the left
        
        if i == 0:
            neighbors = lattice[i]*lattice[i+1]
        elif i == (N-1):
            neighbors = lattice[i]*lattice[i-1]
        elif i < (N-1):
            neighbors = lattice[i]*lattice[i-1] + lattice[i]*lattice[i+1]

        E_i = -1*J * neighbors
        
        #Flip the spin

        lattice[i] = -1 * lattice[i]
        
        #Redo boundary conditons, just copy-paste

        if i == 0:
            neighbors = lattice[i]*lattice[i+1]
        elif i == (N-1):
            neighbors = lattice[i]*lattice[i-1]
        elif i < (N-1):
            neighbors = lattice[i]*lattice[i-1] + lattice[i]*lattice[i+1]

        E_f = -1*J * neighbors
        
        #choose to maintain flip or not
        #the higher the temp the more likely the particle is going
        #to choose the higher energy state

        if E_f < E_i or np.random.rand(1)[0] < np.exp(-1*(E_f-E_i)/(k_B*T)):
            rep = rep + 1

        else:
            lattice[i] = -1 * lattice[i]
            rep = rep + 1
    
    return lattice

i = 0
it = []
ham = []

lattice = np.random.choice((-1,1),N)

while i < 1e3:
    it.append(i*1)
    ham.append(energy1D(lattice, J))
    lattice = Metro_1D(lattice, 1, T, J)
    i = i + 1

[x,y] = [it,ham]

plt.plot(x, y)

plt.xlabel('Iterations')
plt.ylabel('Hamiltonian')

"""### Two Dimensions"""

J = 1
N = 16
rep = 0
T = 3e2 #Kelvin
    
#Create Lattice
lattice = np.random.choice((-1,1),(N,N))

#Plot initial alignment
plt.imshow(lattice, cmap = 'gray')
plt.show()
#Print Hamiltonian
print('Initial Hamiltonian is ' + str(energy2D(lattice,J)))
#Print Average Magnetation
print('Initial Avg Mag is ' + str(avg_mag_2D(lattice)))

def Metro_2D(lattice, iterations, T, J):
    """Performs N iterations on two dimensional lattice using Metropolitan Algorithm
    
    INPUTS:
    
    lattice: two dimensional numpy array of up (1) or down (-1) spins
    
    iterations: number of desired repitions of Metro function.
    One repition = one flip at one site
    
    T = temperature in kelvin
    
    J = constant, set to 1 if you don't care
    
    RETURNS:
    
    Numpy array of the altered lattice"""
    
    #Initialize loop
    rep = 0
    
    N = len(lattice)
    while rep < iterations:
        #randomly select the index of a lattice site
        i = np.random.randint(0,N)
        j = np.random.randint(0,N)

        E_i = -1*J * get_neighbors_2D(lattice,i,j)
        
        #flip spin

        lattice[i,j] = -1 * lattice[i,j]
        
        #calculate new energy energy_final

        E_f = -1* J * get_neighbors_2D(lattice,i,j)
        
        #condition, the higher the temp the more likely
        #the particle is to flip even if it's a higher
        #energy state

        if E_f < E_i or np.random.rand(1)[0] < np.exp(-1*(E_f-E_i)/(k_B*T)):
            rep = rep + 1

        else:
            lattice[i,j] = -1 * lattice[i,j]
            rep = rep + 1
        
    return lattice

plt.imshow(Metro_2D(lattice, 1e4, T, J), cmap = 'gray')
plt.show()

i = 0
it = []
ham = []

lattice = np.random.choice((-1,1),(N,N))

while i < 1e3:
    it.append(i*100)
    ham.append(energy2D(lattice, J))
    count = 0
    lattice = Metro_2D(lattice, 100, T, J)
    i = i + 1

[x,y] = [it,ham]

plt.plot(x, y)

plt.xlabel('Iterations')
plt.ylabel('Hamiltonian')

"""### Animation"""

J = 1
N = 200
T = 3e2 #Kelvin

lattice = np.random.choice((-1,1),(N,N))

fig = plt.figure(figsize=(8,6))
ax = plt.axes()
plt.close()

def frame(w):
    ax.clear()
    global x,y
    t = w + 10000 #t increases by 1 each frame
    lat = Metro_2D(lattice, 1e4, T, J)
    ax.set_title("Quenching")
    ax.imshow(lattice, cmap = 'gray')

anim = animation.FuncAnimation(fig, frame, frames=300,
                               blit=False, repeat=True)

anim

"""### Three Dimensions"""

rep = 0
N = 16
J =  1
T = 1e-20 #Kelvin

#create lattice
lattice = np.random.choice((-1,1),(N,N,N))
#Show initial energy
print('Initial Hamiltonian is ' + str(energy3D(lattice,J)))
#Show initial average mag
print('Initial Avg Mag is ' + str(avg_mag_3D(lattice)))
#Show the second layer of the solid in initial state
plt.imshow(lattice[1], cmap = 'gray')
plt.show()

def Metro_3D(lattice,iterations,T,J):
    
    """Performs N iterations on three dimensional lattice using Metropolitan Algorithm
    
    INPUTS:
    
    lattice: three dimensional numpy array of up (1) or down (-1) spins
    
    iterations: number of desired repitions of Metro function.
    One repition = one flip at one site
    
    T = temperature in kelvin
    
    J = constant, set to 1 if you don't care
    
    RETURNS:
    
    Numpy array of the altered lattice"""
    
    #Initialize loop
    rep = 0
    
    N = len(lattice)
    while rep < iterations:
        #Select lattice site at random
        k = np.random.randint(0,N) #Layer number
        i = np.random.randint(0,N) #row number
        j = np.random.randint(0,N) #column number
        
        #initial energy

        E_i = -1* J * get_neighbors_3D(lattice,k,i,j)
        
        #Flip the spin

        lattice[k,i,j] = -1 * lattice[k,i,j]
        
        #FINAL ENERGY

        E_f = -1* J * get_neighbors_3D(lattice,k,i,j)
        
        #Condtion for keeping the flip
        #the higher the temp the more likely
        #the particle is to flip even if it's a higher
        #enegy state

        if E_f < E_i or np.random.rand(1)[0] < np.exp(-1*(E_f-E_i)/(k_B*T)):
            rep = rep + 1

        else:
            lattice[k,i,j] = -1 * lattice[k,i,j]
            rep = rep + 1
    
    return lattice

plt.imshow(Metro_3D(lattice, 1e5, T, J)[1], cmap = 'gray')
plt.show()

i = 0
it = []
ham = []

lattice = np.random.choice((-1,1),(N,N,N))

while i < 1e2:
    it.append(i*1000)
    ham.append(energy3D(lattice, J))
    count = 0
    lattice = Metro_3D(lattice, 1000, T, J)
    i = i + 1

[x,y] = [it,ham]

plt.plot(x, y)

plt.xlabel('Iterations')
plt.ylabel('Hamiltonian')

"""

### Fourth Dimensional Metro Function"""

rep = 0
N = 10
J =  1
T = 1e-20 #Kelvin

#create lattice
lattice = np.random.choice((-1,1),(N,N,N,N))
#Show initial energy
print('Initial Hamiltonian is ' + str(hamiltonian_4D(lattice, J)))
#Show initial average mag
print('Initial Avg Mag is ' + str(avg_mag_4D(lattice)))

def Metro_4D(lattice,iterations,T,J):
    
    """Performs N iterations on four dimensional lattice using Metropolitan Algorithm
    
    INPUTS:
    
    lattice: four dimensional numpy array of up (1) or down (-1) spins
    
    iterations: number of desired repitions of Metro function.
    One repition = one flip at one site
    
    T = temperature in kelvin
    
    J = constant, set to 1 if you don't care
    
    RETURNS:
    
    Numpy array of the altered lattice"""
    
    #Initialize loop
    front = 0
    rep = 0
    N = len(lattice)
    while rep < iterations:
        #Select lattice site at random
        z = np.random.randint(0,N) #Solid Number
        k = np.random.randint(0,N) #Layer number
        i = np.random.randint(0,N) #row number
        j = np.random.randint(0,N) #column number
        
        #INITIAL ENERGY

        E_i = -1* J * get_neighbors_4D(lattice,z,k,i,j)
        
        #Flip the spin

        lattice[z,k,i,j] = -1 * lattice[z,k,i,j]
        
        #FINAL ENERGY
        
        E_f = -1* J * get_neighbors_4D(lattice,z,k,i,j)
        
        #Condtion for keeping the flip
        #the higher the temp the more likely
        #the particle is to flip even if it's a higher
        #enegy state

        if E_f < E_i or np.random.rand(1)[0] < np.exp(-1*(E_f-E_i)/(k_B*T)):
            rep = rep + 1

        else:
            lattice[z,k,i,j] = -1 * lattice[z,k,i,j]
            rep = rep + 1
    
    #Print final results
    
    return lattice

i = 0
it = []
ham = []

lattice = np.random.choice((-1,1),(N,N,N,N))

while i < 1e2:
    it.append(i*1000)
    ham.append(hamiltonian_4D(lattice, J))
    lattice = Metro_4D(lattice, 1000, T, J)
    i = i + 1

[x,y] = [it,ham]

plt.plot(x, y)

plt.xlabel('Iterations')
plt.ylabel('Hamiltonian')

"""## Metro: Variable Lattice Shapes
In this section, the Metropolis algorithm is used for any arbitrary adjacency matrix and associated charges array (that represents the charges of each node in the adjacency matrix). Specific examples in this section used are the standard 2D square lattice as well as 2D triangular lattices.
"""

def MetroAdjacency(adj, charges, J, T, iterations=200):
  for _ in range(iterations):
    site = np.random.randint(0, len(adj))

    edges = [i for i in range(len(adj[site])) if adj[site][i] == 1]
    print(site, edges)

    delta = 0
    for edge in edges:
      delta += -1 * J * (charges[edge] * -1* charges[site])

    if delta < 0 or np.random.randint(0, 1) < np.e**(-(delta)/(sp.constants.k*T)):
      charges[site] *= -1

  return charges

# Adjacency Matrix test with 2D Square Latice

n = 4
T = 1
J = 1

# adjacency = genTriangualrAdjacency(10)
lat = gen2DLat(n)
adjacency = gen2DAdjacency(lat)
show(adjacency)

n = len(adjacency)
G = nx.DiGraph() 

colors = []
# charges = np.random.choice([1, -1], size=(n))
charges = [element for sublist in lat for element in sublist]
for i in range(n):
  if charges[i]== 1:
    colors.append('orange')
  elif charges[i]==-1:
    colors.append('blue')

for i in range(n): 
 for j in range(n): 
   if adjacency[i][j] == 1: 
      G.add_edge(i,j) 

nx.draw(G, node_color=colors) 
plt.show() 

iterations = 20
for i in range(iterations):

  charges = MetroAdjacency(adjacency, charges, J, T, 1)
  colors = []
  for i in range(n):
    if charges[i]== 1:
      colors.append('orange')
    elif charges[i]==-1:
      colors.append('blue')

  nx.draw(G, node_color=colors) 
  plt.show()

# Adjacency Matrix test with 2D Triangular Latice


adjacency = genTriangularAdjacency(10)
show(adjacency)

n = len(adjacency)
G = nx.DiGraph() 

colors = []
charges = np.random.choice([1, -1], size=(n))
for i in range(n):
  if charges[i]== 1:
    colors.append('orange')
  elif charges[i]==-1:
    colors.append('blue')

for i in range(n): 
 for j in range(n): 
   if adjacency[i][j] == 1: 
      G.add_edge(i,j) 

nx.draw(G, node_color=colors) 
plt.show() 

iterations = 20
for i in range(iterations):

  charges = MetroAdjacency(adjacency, charges, 1, 1, 1)
  colors = []
  for i in range(n):
    if charges[i]== 1:
      colors.append('orange')
    elif charges[i]==-1:
      colors.append('blue')

  nx.draw(G, node_color=colors) 
  plt.show()

"""## Metro: Finding equilibrium & Critical Temperatures

In this code we'll be running the Metro algorithm until we reach the desired state of equilibrium. We will say the code has reached equilibrium when the change in energy of the system is below a given tolerance. This means we'll have to keep track of the energy as the code runs to compare the energy in one state with the energy in another. We'll call this periodicy of checking the energy print_freq. As the code runs we'll increase this number so the code isn't bogged down by having to calculate the hamiltonian every single time something changes. But we'll also cap it at some point so that if the code is close to reaching equilibrium we aren't running exponentially more iterations.

### One Dimension
"""

J = 1
N = 128
rep = 0
T = 3e2 #Kelvin
tolerance = 0 #Joules

#Create Lattice
lattice = np.random.choice((-1,1),N)

#Print Initial Hamiltonian
print('Initial Hamiltonian is ' + str(energy1D(lattice, J)))

#Print Initial Average Mag
print('Initial Avg Mag is ' + str(avg_mag(lattice)))

def Metro_1D_eq(lattice, T, J, tolerance):
    """Repeats the Metropolitan Algorithm in a one dimensional lattice until equilibrum
    is reached within a desired tolerance
    
    INPUTS:
    
    lattice: one dimensional numpy array containing up (1) or down (-1) spins
    
    T: temperature in Kelvin
    
    J: Magnetic constant for material (put at 1 if you don't care)
    
    tolerance: largest change in hamiltonian that is still considered
    within equilibrium. The lower the number the more relaxed the system will
    get but run times may increase dramatically
    
    RETURNS:
    
    lattice: numpy array of lattice after reaching equilibrium
    
    it: list of numbers representing the number of iterations at which the
    Hamiltonian was calculated
    
    ham: list of Hamiltonians calculated throughout the running of the code
    
    avg_mag(lattice): The average magnetic alignment of the lattice after reaching
    equilibrium"""
    
    #intialize loop
    it = [0] #iterations, used for plotting later on
    ham = [energy1D(lattice, J)] #hamiltonian, used to find the equilibrium point
    
    rep = 0 #Counter, used to determine if hamiltonian will be calculated or not
              #If we were to calc the Ham every iteration that'd be a huge time waste
    
    r = 1e3 #frequency of which the ham will be printed (1 calc/N iterations)
                   #frequency will increase as the code runs exponentially
                   #Until it caps at 10,000. This is also to save run time
    
    delta_ham = tolerance + 1 #change in energy between two states of lattice
    
    ##We will say the code has reached equilibrium when the change
    ##in energy of the system is below a given tolerance
    
    while delta_ham > tolerance:
        
        lattice = Metro_1D(lattice, r, T, J)
            
        it.append(it[len(it) - 1] + r)

        ham.append(energy1D(lattice, J))
        
        delta_ham = abs(ham[len(ham)-1] - ham[len(ham)-2])
        
        if r < 1e5:
            r = r + 1
    
    return lattice, it, ham, avg_mag(lattice)

[lattice, x, y, z] = Metro_1D_eq(lattice, T, J, 0)

plt.plot(x, y)

plt.xlabel('Iterations')
plt.ylabel('Hamiltonian')

plt.show()

i = 0
it = []
ham = []

lattice = np.random.choice((-1,1),N)

while i < 1e3:
    it.append(i*1)
    ham.append(energy1D(lattice, J))
    lattice = Metro_1D(lattice, 1, T, J)
    i = i + 1

[x,y] = [it,ham]

plt.plot(x, y)

plt.xlabel('Iterations')
plt.ylabel('Hamiltonian')

"""### Two Dimensions"""

J = 1
N = 16
rep = 0
T = 3e2 #Kelvin
tolerance = 2


#Create Lattice
lattice = np.random.choice((-1,1),(N,N))

#Plot initial alignment
plt.imshow(lattice, cmap = 'gray')
plt.show()

#Print Hamiltonian
print('Initial Hamiltonian is ' + str(energy2D(lattice,J)))

#Print Average Magnetation
print('Initial Avg Mag is ' + str(avg_mag_2D(lattice)))

def Metro_2D_eq(lattice, T, J, tolerance):
    """Repeats the Metropolitan Algorithm in a square lattice until equilibrum
    is reached within a desired tolerance
    
    INPUTS:
    
    lattice: two dimensional numpy array containing up (1) or down (-1) spins
    
    T: temperature in Kelvin
    
    J: Magnetic constant for material (put at 1 if you don't care)
    
    tolerance: largest change in hamiltonian that is still considered
    within equilibrium. The lower the number the more relaxed the system will
    get but run times may increase dramatically
    
    RETURNS:
    
    lattice: numpy array of lattice after reaching equilibrium
    
    it: list of numbers representing the number of iterations at which the
    Hamiltonian was calculated
    
    ham: list of Hamiltonians calculated throughout the running of the code
    
    avg_mag(lattice): The average magnetic alignment of the lattice after reaching
    equilibrium"""
    
    #intialize loop
    rep = 0 #Repetitions
    it = [0] #iterations, used for plotting later on
    ham = [energy2D(lattice, J)] #hamiltonian, used to find the equilibrium point
    
    r = 64 #frequency of which the ham will be printed (1 calc/N iterations)
                   #frequency will increase as the code runs exponentially
                   #Until it caps at 10,000. This is also to save run time
    
    delta_ham = tolerance + 1 #change in energy between two states of lattice
    
    ##We will say the code has reached equilibrium when the change
    ##in energy of the system is below a given tolerance
    
    while delta_ham > tolerance:
        
        lattice = Metro_2D(lattice, r, T, J)
            
        it.append(it[len(it) - 1] + r)

        ham.append(energy2D(lattice,J))
        
        delta_ham = abs(ham[len(ham) - 1] - ham[len(ham) - 2])
        
        if r < 1e5:
            r = r * 2
        
    return lattice, it, ham, avg_mag_2D(lattice)

[lattice, x , y, z] = Metro_2D_eq(lattice, T, J, tolerance)

plt.plot(x,y)

plt.xlabel('Iterations')
plt.ylabel('Hamiltonian')
plt.show()

plt.imshow(lattice, cmap = 'gray')
plt.ylabel('Hamiltonian')

"""### Critical Temperature in 2D

To generate the critical temeprature the code will run four different lattices starting with either a 75% up or 75% down initial orientation and averaging their results to find the average spin at various temperatures. We'll consider the error to be Poison error which is correlated to the size of the arrays and number of particles used.
"""

J = 1
T = np.arange(0.5,10,0.25)*(1/k_B) #Kelvin
tolerance = 8
N = 8

i = 0
mag_1 = []

mag_2 = []

mag_3 = []

mag_4 = []

lattice_1 = np.random.choice((-1,1), size = (N,N), p = (0.25,0.75))
lattice_2 = np.random.choice((-1,1), size = (N,N), p = (0.25,0.75))
lattice_3 = np.random.choice((-1,1), size = (N,N), p = (0.25,0.75))
lattice_4 = np.random.choice((-1,1), size = (N,N), p = (0.25,0.75))

while i < len(T):
    mag_1.append(Metro_2D_eq(lattice_1, T[i], J, tolerance)[3])
    mag_2.append(Metro_2D_eq(lattice_2, T[i], J, tolerance)[3])
    mag_3.append(Metro_2D_eq(lattice_3, T[i], J, tolerance)[3])
    mag_4.append(Metro_2D_eq(lattice_4, T[i], J, tolerance)[3])
    i = i + 1
    
mag_1 = np.asarray(mag_1)

mag_2 = np.asarray(mag_2)

mag_3 = np.asarray(mag_3)

mag_4 = np.asarray(mag_4)

mag_alignment = 0.25 * (mag_1 + mag_2 + mag_3 + mag_4)

plt.plot(T,mag_alignment, label = '75% Up Lattice', color = 'green')
plt.errorbar(T, mag_alignment, yerr = (4*N**2)**-0.5, color = 'green', alpha = 0.2)

i = 0
mag_1 = []

mag_2 = []

mag_3 = []

mag_4 = []

lattice_1 = np.random.choice((-1,1), size = (N,N), p = (0.75,0.25))
lattice_2 = np.random.choice((-1,1), size = (N,N), p = (0.75,0.25))
lattice_3 = np.random.choice((-1,1), size = (N,N), p = (0.75,0.25))
lattice_4 = np.random.choice((-1,1), size = (N,N), p = (0.75,0.25))

while i < len(T):
    mag_1.append(Metro_2D_eq(lattice_1, T[i], J, tolerance)[3])
    mag_2.append(Metro_2D_eq(lattice_2, T[i], J, tolerance)[3])
    mag_3.append(Metro_2D_eq(lattice_3, T[i], J, tolerance)[3])
    mag_4.append(Metro_2D_eq(lattice_4, T[i], J, tolerance)[3])
    i = i + 1
    
mag_1 = np.asarray(mag_1)

mag_2 = np.asarray(mag_2)

mag_3 = np.asarray(mag_3)

mag_4 = np.asarray(mag_4)

mag_alignment = 0.25 * (mag_1 + mag_2 + mag_3 + mag_4)

    
plt.plot(T,mag_alignment, label = '75% Down Lattice', color = 'red')
plt.errorbar(T, mag_alignment, yerr = (4*N**2)**-0.5, color = 'red', alpha = 0.2)

plt.xlabel('Temperature (Kelvin)')
plt.ylabel('Average Magnetic Alignment')
plt.title('2D Average Magnetic Alignment vs Temperature')
plt.legend()

"""### Three Dimensions"""

J = 1
N = 16
rep = 0
T = 3e2 #Kelvin
tolerance = 128

#create lattice
lattice = np.random.choice((-1,1),(N,N,N))

#Show initial energy
print('Initial Hamiltonian is ' + str(energy3D(lattice,J)))

#Show initial average mag
print('Initial Avg Mag is ' + str(avg_mag_3D(lattice)))

#Show the second layer of the solid in initial state
plt.imshow(lattice[1], cmap = 'gray')
plt.show()


def Metro_3D_eq(lattice,T,J,tolerance):
    
    """Repeats the Metropolitan Algorithm in a cube lattice until equilibrum
    is reached within a desired tolerance
    
    INPUTS:
    
    lattice: three dimensional numpy array containing up (1) or down (-1) spins
    
    T: temperature in Kelvin
    
    J: Magnetic constant for material (put at 1 if you don't care)
    
    tolerance: largest change in hamiltonian that is still considered
    within equilibrium. The lower the number the more relaxed the system will
    get but run times may increase dramatically
    
    RETURNS:
    
    lattice: numpy array of lattice after reaching equilibrium
    
    it: list of numbers representing the number of iterations at which the
    Hamiltonian was calculated
    
    ham: list of Hamiltonians calculated throughout the running of the code
    
    avg_mag(lattice): The average magnetic alignment of the lattice after reaching
    equilibrium"""

    front = 0
    
    #intialize loop
    rep = 0 #Repetitions
    it = [0] #iterations, used for plotting later on
    ham = [energy3D(lattice, J)] #hamiltonian, used to find the equilibrium point
    
    r = 256 #frequency of which the ham will be printed (1 calc/N iterations)
                   #frequency will increase as the code runs exponentially
                   #Until it caps at 10,000. This is also to save run time
    
    delta_ham = tolerance + 1 #change in energy between two states of lattice
    
    ##We will say the code has reached equilibrium when the change
    ##in energy of the system is below a given tolerance
    
    while delta_ham > tolerance:
        lattice = Metro_3D(lattice, r, T, J)
            
        it.append(it[len(it) - 1] + r)

        ham.append(energy3D(lattice,J))
        
        delta_ham = abs(ham[len(ham) - 1] - ham[len(ham) - 2])
        
        if r < 1e5:
            r = r * 2
    
    return lattice, it, ham, avg_mag_3D(lattice)

[lattice, x, y, z] = Metro_3D_eq(lattice,T,J,tolerance)

plt.plot(x,y)

plt.xlabel('Iterations')
plt.ylabel('Hamiltonian')
plt.show()

plt.imshow(lattice[1], cmap = 'gray')

"""### Critical Temperature in 3D"""

J = 1
T = np.arange(0.5,10,0.25)*(1/k_B) #Kelvin
tolerance = 64
N = 8

i = 0
mag_1 = []

mag_2 = []

mag_3 = []

mag_4 = []

lattice_1 = np.random.choice((-1,1), size = (N,N,N), p = (0.25,0.75))
lattice_2 = np.random.choice((-1,1), size = (N,N,N), p = (0.25,0.75))
lattice_3 = np.random.choice((-1,1), size = (N,N,N), p = (0.25,0.75))
lattice_4 = np.random.choice((-1,1), size = (N,N,N), p = (0.25,0.75))

while i < len(T):
    mag_1.append(Metro_3D_eq(lattice_1, T[i], J, tolerance)[3])
    mag_2.append(Metro_3D_eq(lattice_2, T[i], J, tolerance)[3])
    mag_3.append(Metro_3D_eq(lattice_3, T[i], J, tolerance)[3])
    mag_4.append(Metro_3D_eq(lattice_4, T[i], J, tolerance)[3])
    i = i + 1
    
mag_1 = np.asarray(mag_1)

mag_2 = np.asarray(mag_2)

mag_3 = np.asarray(mag_3)

mag_4 = np.asarray(mag_4)

mag_alignment = 0.25 * (mag_1 + mag_2 + mag_3 + mag_4)

plt.plot(T,mag_alignment, label = '75% Up Lattice', color = 'green')
plt.errorbar(T, mag_alignment, yerr = (4*N)**-0.5, color = 'green', alpha = 0.2)

i = 0

mag_1 = []

mag_2 = []

mag_3 = []

mag_4 = []

lattice_1 = np.random.choice((-1,1), size = (N,N,N), p = (0.75,0.25))
lattice_2 = np.random.choice((-1,1), size = (N,N,N), p = (0.75,0.25))
lattice_3 = np.random.choice((-1,1), size = (N,N,N), p = (0.75,0.25))
lattice_4 = np.random.choice((-1,1), size = (N,N,N), p = (0.75,0.25))

while i < len(T):
    mag_1.append(Metro_3D_eq(lattice_1, T[i], J, tolerance)[3])
    mag_2.append(Metro_3D_eq(lattice_2, T[i], J, tolerance)[3])
    mag_3.append(Metro_3D_eq(lattice_3, T[i], J, tolerance)[3])
    mag_4.append(Metro_3D_eq(lattice_4, T[i], J, tolerance)[3])
    i = i + 1
    
mag_1 = np.asarray(mag_1)

mag_2 = np.asarray(mag_2)

mag_3 = np.asarray(mag_3)

mag_4 = np.asarray(mag_4)

mag_alignment = 0.25 * (mag_1 + mag_2 + mag_3 + mag_4)

    
plt.plot(T,mag_alignment, label = '75% Down Lattice', color = 'red')
plt.errorbar(T, mag_alignment, yerr = (4*N)**-0.5, color = 'red', alpha = 0.2)

plt.xlabel('Temperature (Kelvin)')
plt.ylabel('Average Magnetic Alignment')
plt.title('3D Average Magnetic Alignment vs Temperature')
plt.legend()