# Helper Functions
"""
# Generation, printing/image generation, energy calctions
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy import constants
import networkx as nx 

import random

from matplotlib import rc
rc('animation', html='jshtml')
from math import *
import matplotlib.animation as animation
import matplotlib

matplotlib.rcParams['animation.embed_limit'] = 2**128

import sys
sys.setrecursionlimit(10000)



"""## Latice Generation"""

# Generation
def gen1DLat(n):
  """Generates a 1D lattice
    Returns an np array where each element represents an element in a lattice"""
  return np.random.choice([1, -1], size=(n))

def gen2DLat(n):
  """Generates a 2D lattice
    Returns an np 2D array where each element represents an element in a lattice"""
  return np.random.choice([1, -1], size=(n, n))

def gen3DLat(n):
  """Generates a D lattice
  Returns an np 3D array where each element represents an element in a lattice"""
  return np.random.choice([1, -1], size=(n, n, n))

"""## Latice Printing"""

# Printing
def print1D(lat, tabs=0):
  """Prints a 1D lattice in matrix form"""
  for i in range(len(lat)):
    print("\t"*tabs + str(lat[i]), end=" \t")
  print()

def print2D(lat, tabs=0):
  """Prints a 2D lattice in matrix form"""
  for i in range(len(lat)):
    for j in range(len(lat[0])):
      print("\t"*tabs +str(lat[i][j]), end=" \t")
    print()
  print()

def print3D(lat, tabs=0):
  """Prints a 3D lattice in matrix form"""
  for i in range(len(lat)):
    print('layer', i)
    for j in range(len(lat[0])):
      for z in range(len(lat[0][0])):
        print("\t"*tabs +str(lat[i][j][z]), end=" \t")
      print()
    # print()
  print()

def show(lat):
  plt.imshow(lat, cmap = 'gray')
  plt.show()

"""## Lattice to Adjacency Matrix"""

# Get Elements
def get1DElements(lat):
  return list(lat) # a list of the elements is just the lat

def get2DElements(lat):
  elems = []
  coords = []
  for i in range(len(lat[0])):
    for j in range(len(lat)):
      elems.append(lat[i][j])
      coords.append([i, j])

  return elems, coords

# Getting Adjacency Matrices 

def gen2DAdjacency(lat):
  elems, coords = get2DElements(lat)

  max_x = len(lat)
  max_y = len(lat[0])

  adjacency = np.zeros((max_x * max_y, max_x * max_y))

  for idx in range(len(coords)):
    x = coords[idx][0]
    y = coords[idx][1]

    curr = x*max_x + y

    # left neighbor
    if x > 0:
      adjacency[curr][(x-1)*max_x + y] = 1
      adjacency[(x-1)*max_x + y][curr] = 1

    # right neighbor
    if x < max_x-1:
      adjacency[curr][(x+1)*max_x + y] = 1
      adjacency[(x+1)*max_x + y][curr] = 1

    # top neighbor
    if y > 0:
      adjacency[curr][x*max_x + (y-1)] = 1
      adjacency[x*max_x + (y-1)][curr] = 1

    # bottom neighbor
    if y < max_y-1:
      adjacency[curr][x*max_x + (y+1)] = 1
      adjacency[x*max_x + (y+1)][curr] = 1
    
    # self
    adjacency[curr][curr] = -1

  return adjacency

# Printing Adjacency Matrices
def print2DAdjacency(lat, arr):
  elems, coords = get2DElements(lat)

  print("       ", end="")
  for i in range(len(arr)):
    print(coords[i], end="\t")
  print()

  for i in range(len(arr)):
    print(coords[i], end="\t")
    for j in range(len(arr[0])):
      print(round(arr[i][j]), end="\t")
    print()

def gen2DAdjacency(lat):
  elems, coords = get2DElements(lat)

  max_x = len(lat)
  max_y = len(lat[0])

  adjacency = np.zeros((max_x * max_y, max_x * max_y))

  for idx in range(len(coords)):
    x = coords[idx][0]
    y = coords[idx][1]

    curr = x*max_x + y

    # left neighbor
    if x > 0:
      adjacency[curr][(x-1)*max_x + y] = 1
      adjacency[(x-1)*max_x + y][curr] = 1

    # right neighbor
    if x < max_x-1:
      adjacency[curr][(x+1)*max_x + y] = 1
      adjacency[(x+1)*max_x + y][curr] = 1

    # top neighbor
    if y > 0:
      adjacency[curr][x*max_x + (y-1)] = 1
      adjacency[x*max_x + (y-1)][curr] = 1

    # bottom neighbor
    if y < max_y-1:
      adjacency[curr][x*max_x + (y+1)] = 1
      adjacency[x*max_x + (y+1)][curr] = 1
    
    # self
    adjacency[curr][curr] = -1

  return adjacency
# Energy
def energy1D(lat, J):
  sum = 0
  for i in range(len(lat)-1):
    sum += lat[i] * lat[i+1]

  return -1 * J * sum


# 2D energy (BFS)
def find2DDownwardsNeighbors(curr_node, x_max, y_max):
  i = curr_node[0]
  j = curr_node[1]

  neighbors = []
  if i < x_max - 1:
    neighbors.append([i+1, j])
  if j < y_max - 1:
    neighbors.append([i, j+1])

  return neighbors

"""## Triangular Adjacency Matrices"""

# Triangular matrix 

def genTriangularAdjacency(nodes):
  if nodes == 3:
    # 3 nodes; 1 triangle (2nd triangular number)
    adjacency = [[-1,  1,  1], 
                [ 1, -1,  1],
                [ 1,  1, -1]]
  
  elif nodes == 6:
    # 6 nodes; 4 triangles
    adjacency =[[-1,  1,  1,  0,  0,  0], # node 1
                [ 1, -1,  1,  1,  1,  0], # node 2
                [ 1,  1, -1,  0,  1,  1], # node 3
                [ 0,  1,  0, -1,  1,  0], # node 4
                [ 0,  1,  1,  1, -1,  1], # node 5
                [ 0,  0,  1,  0,  1, -1]] # node 6
            
  elif nodes == 10:
    # 10 nodes; 9 triangles
    adjacency = [[-1,  1,  1,  0,  0,  0,  0,  0,  0,  0], # a
                [ 1, -1,  1,  1,  1,  0,  0,  0,  0,  0], # b
                [ 1,  1, -1,  0,  1,  1,  0,  0,  0,  0], # c
                [ 0,  1,  0, -1,  1,  0,  1,  1,  0,  0], # d
                [ 0,  1,  1,  1, -1,  1,  0,  1,  1,  0], # e
                [ 0,  0,  1,  0,  1, -1,  0,  0,  1,  1], # f
                [ 0,  0,  0,  1,  0,  0, -1,  1,  0,  0], # g
                [ 0,  0,  0,  1,  1,  0,  1, -1,  1,  0], # h
                [ 0,  0,  0,  0,  1,  1,  0,  1, -1,  1], # i
                [ 0,  0,  0,  0,  0,  1,  0,  0,  1, -1]] # j

  elif nodes == 15:
  # 15 nodes; 16 triangles
    adjacency = [[-1,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0], # 0, a
              [ 1, -1,  1,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0], # 1, b
              [ 1,  1, -1,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0], # 2, c
              [ 0,  1,  0, -1,  1,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0], # 3, d
              [ 0,  1,  1,  1, -1,  1,  0,  1,  1,  0,  0,  0,  0,  0,  0], # 4, e
              [ 0,  0,  1,  0,  1, -1,  0,  0,  1,  1,  0,  0,  0,  0,  0], # 5, f
              [ 0,  0,  0,  1,  0,  0, -1,  1,  0,  1,  1,  0,  0,  0,  0], # 6, g
              [ 0,  0,  0,  1,  1,  0,  1, -1,  1,  0,  1,  1,  0,  0,  0], # 7, h
              [ 0,  0,  0,  0,  1,  1,  0,  1, -1,  1,  0,  1,  1,  0,  0], # 8, i
              [ 0,  0,  0,  0,  0,  1,  1,  0,  1, -1,  1,  0,  1,  1,  0], # 9, j
              [ 0,  0,  0,  0,  0,  0,  1,  1,  0,  1, -1,  0,  0,  1,  1], #10, k
              [ 0,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0, -1,  1,  0,  1], #11, l
              [ 0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  0,  1, -1,  1,  0], #12, m
              [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  0,  1, -1,  1], #13, n
              [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  1, -1]] #14, o

  else:
    print("Node count is not a valid triangular number or the matrix has not been implemented")
    adjacency = []

  return adjacency

"""## Energy Calculations"""

def energy2D(lat, J):
  nodes = []
  nodes.append([0, 0])
  seen = []

  x_max = len(lat)
  y_max = len(lat[0])

  sum = 0
  while len(nodes) > 0:
    curr_node = nodes.pop()
    seen.append(curr_node)

    neighbors = find2DDownwardsNeighbors(curr_node, x_max, y_max)
    curr_charge = lat[curr_node[0]][curr_node[1]]
    for neighbor in neighbors:
      sum += curr_charge * lat[neighbor[0]][neighbor[1]]

      if neighbor not in seen:
        nodes.append(neighbor)
  return -1 * J * sum


# 3D energy (BFS)
def find3DDownwardsNeighbors(curr_node, x_max, y_max, z_max):
  i = curr_node[0]
  j = curr_node[1]
  k = curr_node[2]

  neighbors = []
  if i < x_max - 1:
    neighbors.append([i+1, j, k])
  if j < y_max - 1:
    neighbors.append([i, j+1, k])
  if k < z_max - 1:
    neighbors.append([i, j, k+1])

  return neighbors



def energy3D(lat, J):
  nodes = []
  nodes.append([0, 0, 0])
  
  seen = []

  x_max = len(lat)
  y_max = len(lat[0])
  z_max = len(lat[0])

  sum = 0
  while len(nodes) > 0:
    curr_node = nodes.pop()
    seen.append(curr_node)

    neighbors = find3DDownwardsNeighbors(curr_node, x_max, y_max, z_max)
    curr_charge = lat[curr_node[0]][curr_node[1]][curr_node[2]]
    for neighbor in neighbors:
      sum += curr_charge * lat[neighbor[0]][neighbor[1]][neighbor[2]]

      if neighbor not in seen:
        nodes.append(neighbor)

  return -1 * J * sum

"""## Calculating the Hamiltonian and Average Magnetic Alignment"""

##Constants and Variables
J = 1
j = 0
T = 2.0457e2 #kelvin
k_B = 1.380649e-23 #Boltzman constant

##Helper Functions

def get_neighbors_1D(lattice,i):
    """Helper function that finds the iteractions between lattice
    site i and its neighbors in a one dimensional lattice
    
    INPUTS:
    
    lattice: one dimensional numpy array with up (1) or down (-1) spins
    
    i: lattice index
    
    RETURNS
    
    neighbors: interger sum of the products of the spins of lattice site
    i with it's neighbor(s)"""
    
    N = len(lattice)

    #Boundary condtions: In one dim
    #There are neighbors to the right
    #and left

    #if i == 0 there is only a one to the right

    #if i == N-1 there's only one to the left

    neighbors = 0

    if i == 0:
        left = 0
    else:
        left = lattice[i]*lattice[i-1]

    if i == (N-1):
        right = 0
    else:
        right = lattice[i]*lattice[i+1]

        neighbors = left + right

    return neighbors
    
def get_neighbors_2D(lattice,i,j):
    """Helper function that finds the iteractions between lattice
    site [i,j] and its neighbors in a two dimensional lattice
    
    INPUTS:
    
    lattice: two dimensional numpy array with up (1) or down (-1) spins
    
    i: lattice row index
    
    j: lattice column index
    
    RETURNS
    
    neighbors: interger sum of the products of the spins of lattice site
    [i,j] with it's neighbors"""
        
    N = len(lattice)

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
    else:
        right = lattice[i,j]*lattice[i,j+1]

    if i == 0:
        up = 0
    else:
        up = lattice[i,j]*lattice[i-1,j]

    if i == (N-1):
        down = 0
    else:
        down = lattice[i,j]*lattice[i+1,j]

    neighbors = up + down + left + right

    return neighbors
    
def get_neighbors_3D(lattice,k,i,j):
    """Helper function that finds the iteractions between lattice
    site [k,i,j] and its neighbors in a three dimensional lattice
    
    INPUTS:
    
    lattice: three dimensional numpy array with up (1) or down (-1) spins
    
    k: lattice layer index
    
    i: lattice row index
    
    j: lattice column index
    
    RETURNS
    
    neighbors: interger sum of the products of the spins of lattice site
    [k,i,j] with it's neighbors"""

    N = len(lattice)

    #Boundary conditions to find edges

    ##D arrays have neighbors to the left,
    #right, above, below, front, and back

    #Back is same row and column but preceding
    #layer

    #Front is same row and column but following
    #layer

    #Above/up is same column but preceding row

    #Below/down is same column but following row

    #Left is same row but preceding column

    #Right is same row but following column

    #add up each of these contributions then multiply
    #by -J

    if k == 0:
        back = 0
    else:
        back = lattice[k,i,j]*lattice[k-1,i,j]

    if k == (N-1):
        front = 0
    else:
        front = lattice[k,i,j]*lattice[k+1,i,j]

    if j == 0:
        left = 0
    else:
        left = lattice[k,i,j]*lattice[k,i,j-1]

    if j == (N-1):
        right = 0
    else:
        right = lattice[k,i,j]*lattice[k,i,j+1]

    if i == 0:
        up = 0
    else:
        up = lattice[k,i,j]*lattice[k,i-1,j]

    if i == (N-1):
        down = 0
    else:
        down = lattice[k,i,j]*lattice[k,i+1,j]

    neighbors = up + down + left + right + front + back

    return neighbors
    
def get_neighbors_4D(lattice,z,k,i,j):
    """Helper function that finds the iteractions between lattice
    site [z,k,i,j] and its neighbors in a four dimensional lattice
    
    INPUTS:
    
    lattice: four dimensional numpy array with up (1) or down (-1) spins
    
    z: lattice solid index
    
    k: lattice layer index
    
    i: lattice row index
    
    j: lattice column index
    
    RETURNS
    
    neighbors: interger sum of the products of the spins of lattice site
    [z,k,i,j] with it's neighbors"""
        
    N = len(lattice)

            #Boundary conditions to find edges

    ##4D arrays have neighbors to the left,
    #right, above, below, front, and back. They
    #Also have two more neighbors we'll call alpha
    #And beta

    #Alpha is in the negative z direction

    #Beta is in the positive z direction

    #Back is same row and column but preceding
    #layer

    #Front is same row and column but following
    #layer

    #Above/up is same column but preceding row

    #Below/down is same column but following row

    #Left is same row but preceding column

    #Right is same row but following column

    #add up each of these contributions then multiply
    #by -J

    if z == 0:
        alpha = 0
    else:
        alpha = lattice[z,k,i,j]*lattice[z-1,k,i,j]

    if z == (N-1):
        beta = 0
    else:
        beta = lattice[z,k,i,j]*lattice[z+1,k,i,j]

    if k == 0:
        back = 0
    else:
        back = lattice[z,k,i,j]*lattice[z,k-1,i,j]

    if k == (N-1):
        front = 0
    else:
        front = lattice[z,k,i,j]*lattice[z,k+1,i,j]

    if j == 0:
        left = 0
    else:
        left = lattice[z,k,i,j]*lattice[z,k,i,j-1]

    if j == (N-1):
        right = 0
    else:
        right = lattice[z,k,i,j]*lattice[z,k,i,j+1]

    if i == 0:
        up = 0
    else:
        up = lattice[z,k,i,j]*lattice[z,k,i-1,j]

    if i == (N-1):
        down = 0
    else:
        down = lattice[z,k,i,j]*lattice[z,k,i+1,j]

    neighbors = up + down + left + right + front + back + alpha + beta

    return neighbors

##Functions


##ONE DIMENSION

def hamiltonian(lattice):
    
    i = 0
    H = 0
    N = len(lattice)
    while i < N:
        E = -J * get_neighbors_1D(lattice,i)
        H = H + E
        i = i + 1
    H = H / 2
    return H

def avg_mag(lattice):
    """Function that finds average spin in one dimensional lattice
    
    INPUT:
    
    lattice: one dimensional numpy array with up (1) or down (-1) spins
    
    RETURN:
    
    M: integer, average spin in the lattice"""
    
    N = len(lattice)
    spin_sum = np.sum(lattice)
    M = spin_sum/N
    return M

##TWO DIMENSIONS

def hamiltonian_2D(lattice):
    i = 0 #Row number
    j = 0 #Column number
    H = 0
    N = len(lattice)
    while i < N:
        while j < N:            
            E = -J * get_neighbors_2D(lattice,i,J)
            H = H + E
            j = j + 1
        i = i + 1
        j = 0
    H = H / 2
    return H

def avg_mag_2D(lattice):
    """Function that finds average spin in two dimensional lattice
    
    INPUT:
    
    lattice: two dimensional numpy array with up (1) or down (-1) spins
    
    RETURN:
    
    M: integer, average spin in the lattice"""
    
    spin_sum = np.sum(lattice)
    M = spin_sum/(N**2)
    return M

##THREE DIMENSIONS

def hamiltonian_3D(lattice):
    k = 0 #Layer number
    i = 0 #Row number
    j = 0 #Column number
    H = 0
    N = len(lattice)
    while i < N:
        while j < N:
            while k < N:
                E = -J * get_neighbors_3D(lattice,k,i,j)
                H = H + E
                k = k + 1
            j = j + 1
            k = 0
        i = i + 1
        j = 0
        k = 0
    H = H / 2
    return H

def avg_mag_3D(lattice):
    """Function that finds average spin in three dimensional lattice
    
    INPUT:
    
    lattice: three dimensional numpy array with up (1) or down (-1) spins
    
    RETURN:
    
    M: integer, average spin in the lattice"""
    
    spin_sum = np.sum(lattice)
    M = spin_sum/(N**3)
    return M


##FOUR DIMENSIONS?!?!?!?!

def hamiltonian_4D(lattice, J):
    """Function that finds the Hamiltonian of a four dimensional lattice
    
    INPUTS:
    
    lattice: four dimensional numpy array with up (1) or down (-1) spins
    
    J: magnetic constant for the material (set to 1 if you don't care)
    
    RETURNS:
    
    H: Integer, Hamiltonian of the system"""
    
    z = 0 #Solid number
    k = 0 #Layer number
    i = 0 #Row number
    j = 0 #Column number
    H = 0
    N = len(lattice)
    while i < N:
        while j < N:
            while k < N:
                while z < N:
                    E = -J * get_neighbors_4D(lattice,z,k,i,j)
                    H = H + E
                    z = z + 1
                k = k + 1
                z = 0
            j = j + 1
            k = 0
            z = 0
        i = i + 1
        j = 0
        k = 0
        z = 0
    H = H / 2
    return H

def avg_mag_4D(lattice):
    """Function that finds average spin in four dimensional lattice
    
    INPUT:
    
    lattice: four dimensional numpy array with up (1) or down (-1) spins
    
    RETURN:
    
    M: integer, average spin in the lattice"""
    
    spin_sum = np.sum(lattice)
    M = spin_sum/(N**4)
    return M