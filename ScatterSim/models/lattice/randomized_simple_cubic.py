import numpy as np

from ScatterSim.LatticeObjects import Lattice


class RandomizedSimpleCubicLattice(Lattice):
    def __init__(self, objects, lattice_spacing_a=1.0, sigma_D=0.01,
                 filling_probs=None, lattice_coordinates=None,
                 lattice_types=None, n_repeat=3):
        ''' This is an extended lattice where we randomize the filling.
            objects : the NanoObjects (or composite) to place in lattice
            lattice_spacing_a/b/c : the a/b/c lattice spacings
        '''
        if filling_probs is None:
            self.filling_probs = np.ones(len(objects))

        # prepare variables of lattice
        symmetry = {
            'crystal family': 'cubic',
            'crystal system': 'cubic',
            'Bravais lattice': 'P',
            'crystal class': 'hexoctahedral',
            'point group': 'm3m',
            'space group': 'Pm3m',
        }

        lattice_positions = ['corner'] * len(objects)  # not used

        if lattice_coordinates is None:
            lattice_coordinates = [(0.0, 0.0, 0.0)] * len(objects)
        sub_lattice_coordinates = lattice_coordinates

        lattice_coordinates = list()
        lattice_types = list()
        lattice_positions = list()
        # now prepare the elements
        for i in range(n_repeat):
            x = i / n_repeat
            for j in range(n_repeat):
                y = j / n_repeat
                for k in range(n_repeat):
                    z = k / n_repeat
                    for l, coord in enumerate(sub_lattice_coordinates):
                        xsub, ysub, zsub = x + \
                            coord[0], y + coord[1], z + coord[2]
                        if np.random.uniform() <= filling_probs[l]:
                            lattice_coordinates.append((xsub, ysub, zsub))
                            lattice_types.append(l)
                            lattice_positions.append("N/A")

        lattice_spacing_a = n_repeat * lattice_spacing_a

        # now call parent function to initialize
        super(
            RandomizedSimpleCubicLattice,
            self).__init__(
            objects,
            lattice_spacing_a=lattice_spacing_a,
            sigma_D=sigma_D,
            alpha=90,
            beta=90,
            gamma=90,
            symmetry=symmetry,
            lattice_positions=lattice_positions,
            lattice_coordinates=lattice_coordinates,
            lattice_types=lattice_types)

    def symmetry_factor(self, h, k, l):
        """Returns the symmetry factor (0 for forbidden)."""
        return 1

    def unit_cell_volume(self):
        return self.lattice_spacing_a**3