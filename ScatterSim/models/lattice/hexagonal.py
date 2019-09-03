import numpy as np

from ScatterSim.LatticeObjects import Lattice


class HexagonalLattice(Lattice):
    def __init__(self, objects, lattice_spacing_a=1.0, lattice_spacing_b=None,
                 lattice_spacing_c=None, sigma_D=0.01, lattice_types=None,
                 lattice_coordinates=None, lattice_positions=None,
                 symmetry=None, alpha=90, beta=90, gamma=60):
        ''' Hexagonal lattice

            lattice_spacing_a : the lattice spacing in the x plane
            lattice_spacing_b : the lattice spacing in the x-y plane
                if None, set to a
            lattice_spacing_c : the lattice spacing along the z axis

            Calculations:
                lattice_spacing_a -> a
                lattice_spacing_c -> c
                a1 = (1, 0, 0)*a
                a2 = (1/2., sqrt(3)/2., 0)*b
                a3 = (0, 0, 1)*c

                Reciprocal vectors:
                b1 = 2*np.pi/a*(1, -1/np.sqrt(3), 0)
                b2 = 2*np.pi/b*(0, 2/np.sqrt(3), 0)
                b3 = 2*np.pi/c*(0, 0, 1)
        '''
        # Define the lattice
        if symmetry is None:
            symmetry = {}
            symmetry['crystal family'] = 'cubic'
            symmetry['crystal system'] = 'cubic'
            symmetry['Bravais lattice'] = 'F'
            symmetry['crystal class'] = 'hexoctahedral'
            symmetry['point group'] = 'm3m'
            symmetry['space group'] = 'Fm3m'

        if lattice_positions is None:
            lattice_positions = ['all']

        if lattice_coordinates is None:
            lattice_coordinates = [(0.0, 0.0, 0.0),
                               ]
        if lattice_types is None:
            lattice_types = [1]
        if lattice_spacing_b is None:
            lattice_spacing_b = lattice_spacing_a
        if lattice_spacing_c is None:
            lattice_spacing_c = lattice_spacing_a

        # TODO : gamma should be 60 deg I think (but alpha, beta, gamma are not
        # really used)
        super(
            HexagonalLattice,
            self).__init__(
            objects,
            lattice_spacing_a=lattice_spacing_a,
            lattice_spacing_b=lattice_spacing_b,
            lattice_spacing_c=lattice_spacing_c,
            alpha=alpha,
            beta=beta,
            gamma=gamma,
            sigma_D=sigma_D,
            symmetry=symmetry,
            lattice_positions=lattice_positions,
            lattice_coordinates=lattice_coordinates,
            lattice_types=lattice_types)

        # create the basis vectors and save
        # TODO : make more similar to other implementations (should this go in
        # super()?)
        a = lattice_spacing_a
        b = lattice_spacing_b
        c = lattice_spacing_c

        # the basis vectors for hexagonal
        # TODO : should put this in init?
        self.b1 = 2*np.pi/a*np.array([1, 1/np.sqrt(3), 0])
        self.b2 = 2*np.pi/b*np.array([0, 2/np.sqrt(3), 0])
        self.b3 = 2*np.pi/c*np.array([0, 0, 1])

    def symmetry_factor(self, h, k, l):
        """Returns the symmetry factor (0 for forbidden)."""
        return 1
        # if (h%2==0) and (k%2==0) and (l%2==0):
        # All even
        # return 4
        # elif (h%2==1) and (k%2==1) and (l%2==1):
        # All odd
        # return 4
        # else:
        # return 0

    def q_hkl(self, h, k, l):
        """Determines the position in reciprocal space for the given reflection."""

        qhkl_vector = h*self.b1 + k*self.b2 + l*self.b3

        # NOTE: Valid for ideal hexagonal only
        # NOTE : This was Kevin's old code, but should it not have been
        # "-h + 2*k" instead of "h + 2*k"?
        #qhkl_vector = ( 2*np.pi*h/(self.lattice_spacing_a), \
                        #2*np.pi*(h+2*k)/(np.sqrt(3)*self.lattice_spacing_b), \
                        #2*np.pi*l/(self.lattice_spacing_c) )
        # I moved the 2*np.pi/lattice_spacing etc to the init in the variables
        # self.b1, b2, b3 etc
        qhkl = np.sqrt( qhkl_vector[0]**2 + qhkl_vector[1]**2 + qhkl_vector[2]**2 )

        return (qhkl, qhkl_vector)



    def q_hkl_length(self, h, k, l):
        raise NotImplementedError("Add if needed")

    def unit_cell_volume(self):
        raise NotImplementedError("Add if needed")