<<<<<<< HEAD
# IsPy Package
The Ising Python package (IsPy) provides implementation of four different Monte Carlo algorithms that simulates the Ising model of solids: Metropolis, Glauber, Wolff, and Cluster Heat Bath algorithms. It simulates how solids relax in low temperatures to aligned/anti-aligned spins. 
This code is flexible enough to be used for square and non-square lattices. This code relies on numpy and SciPy for array generation/random number generation and the Boltzmann constant. 

Installation
git clone https://github.com/TODO.

Recommendations
We highly recommend importing sys to increase the recursion limit of functions in order to utilize the CHB code to its best potential. We also highly recommend using matplotlib.pyplot to visualize these lattices, in particular plt.imshow() to show how these lattices change over time. We’ve also checked and the code plays well with animation in case you want to use matplotlib.animation to animate the relaxation of the system in a cool (and very satisfying) figure.

This code uses lattices of various dimensions. Glauber and Metro run well in one, two, and three dimensions. We even have a function to run Metro in four dimensions, if that’s your thing. CHB runs in two and three dimensions well.

The algorithms can be run in two separate manners: running for N iterations and running until equilibrium. As you can expect the N iterations versions asks you to specify the number of times you’d like the code to be run before returning the altered lattice while the equilibrium functions will run until reaching a certain tolerance where the change of energy between two iterations is so small it’s considered insignificant. Obviously this second version can be very taxing on runtime so make sure to choose lattice size and tolerance wisely. Our best advice is to start small and with a large tolerance and then work your way down.

Notes

Most of this code depends on the use of square matrices to run. This is a limit and an update to include rectangular matrices could be implemented with one day's work. But as of right now this is a limit the user must work within.
We have a function that generates 1D, 2D, and 3D lattices in gen1DLat(), gen2DLat(), gen3DLat(), respectively.

Cluster Heat Bath utilizes a recursive function. For larger arrays (particularly in 3D) this means that the recursion limit will be reached very easily. To avoid a recursion error we suggest raising the recursion limit to 100,000 or even 1,000,000, the limit we used in the making of this module.

License
This package is released under the MIT License. See the LICENSE file for more information.

Contributing
If you find a bug or would like to contribute to the package, please open an issue or pull request on the GitHub repository. We’ll make sure to respond immediately.
=======
# IsingModel
>>>>>>> aa2bb4b01be803328a54fbcece4e476aa93657a3
