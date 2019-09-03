import numpy as np

from ScatterSim.LatticeObjects import Lattice


class FCCLattice(Lattice):
    def __init__(self, objects, lattice_spacing_a=1.0, sigma_D=0.01):
        # Define the lattice
        symmetry = {}
        symmetry['crystal family'] = 'cubic'
        symmetry['crystal system'] = 'cubic'
        symmetry['Bravais lattice'] = 'F'
        symmetry['crystal class'] = 'hexoctahedral'
        symmetry['point group'] = 'm3m'
        symmetry['space group'] = 'Fm3m'

        # positions = ['all']
        lattice_positions = ['corner', 'faceXY', 'faceYZ', 'faceXZ']

        lattice_coordinates = [(0.0, 0.0, 0.0),
                               (0.5, 0.5, 0.0),
                               (0.0, 0.5, 0.5),
                               (0.5, 0.0, 0.5),
                               ]
        lattice_types = [1, 1, 1, 1]

        super(
            FCCLattice,
            self).__init__(
            objects,
            lattice_spacing_a=lattice_spacing_a,
            sigma_D=sigma_D,
            symmetry=symmetry,
            lattice_positions=lattice_positions,
            lattice_coordinates=lattice_coordinates,
            lattice_types=lattice_types)

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
        """Determines the position in reciprocal space for the given
        reflection."""

        prefactor = (2 * np.pi / self.lattice_spacing_a)
        qhkl_vector = np.array([prefactor * h,
                                prefactor * k,
                                prefactor * l])
        qhkl = np.sqrt(
            qhkl_vector[0]**2 +
            qhkl_vector[1]**2 +
            qhkl_vector[2]**2)

        return (qhkl, qhkl_vector)

    def q_hkl_length(self, h, k, l):
        prefactor = (2 * np.pi / self.lattice_spacing_a)
        qhkl = prefactor * np.sqrt(h**2 + k**2 + l**2)

        return qhkl

    def unit_cell_volume(self):

        return self.lattice_spacing_a**3