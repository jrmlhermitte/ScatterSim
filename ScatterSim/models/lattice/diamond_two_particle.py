import numpy as np

from ScatterSim.LatticeObjects import Lattice


class DiamondTwoParticleLattice(Lattice):
    # a.k.a. zincblende
    '''
        These are implied:
            lattice_spacing_b = lattice_spacing_a
            lattice_spacing_c = lattice_spacing_a
            alpha = radians(90)
            beta = radians(90)
            gamma = radians(90)

    '''

    def __init__(self, objects, lattice_spacing_a=1.0, sigma_D=0.01):
        lattice_spacing_a = lattice_spacing_a

        sigma_D = np.array(sigma_D)          # Lattice disorder

        # Define the lattice
        symmetry = {
            'crystal family': 'cubic',
            'crystal system': 'cubic',
            'Bravais lattice': 'P',
            'crystal class': 'hexoctahedral',
            'point group': 'm3m',
            'space group': 'Pm3m',
        }

        # positions = ['cornerFCC', 'tetrahedralFCC']
        lattice_positions = [
            'corner',
            'faceXY',
            'faceYZ',
            'faceXZ',
            'tetra1',
            'tetra2',
            'tetra3',
            'tetra4']

        lattice_coordinates = [(0.0, 0.0, 0.0),
                               (0.5, 0.5, 0.0),
                               (0.0, 0.5, 0.5),
                               (0.5, 0.0, 0.5),
                               (0.25, 0.25, 0.25),
                               (0.25, 0.75, 0.75),
                               (0.75, 0.25, 0.75),
                               (0.75, 0.75, 0.25),
                               ]

        lattice_types = [1, 1, 1, 1, 2, 2, 2, 2]

        super(
            DiamondTwoParticleLattice,
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

    def unit_cell_volume(self):
        return self.lattice_spacing_a**3