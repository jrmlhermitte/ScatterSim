import numpy as np

from ScatterSim.LatticeObjects import Lattice


class ExtendedLattice(Lattice):
    ''' An extended lattice of objects, can choose random fillings.
        fill_probs : list of probabilities
    '''

    def __init__(
            self,
            objects,
            lattice_spacing_a=1.0,
            fill_probs=None,
            lattice_spacing_b=None,
            lattice_spacing_c=None,
            alpha=90,
            beta=None,
            gamma=None,
            sigma_D=0.01,
            lattice_positions=None,
            lattice_coordinates=None,
            symmetry=None,
            lattice_types=None,
            repeat=3):
        ''' For now only supports one type.'''
        if fill_probs is None:
            fill_probs = [1]

        # sub_lattice_positions = lattice_positions
        sub_lattice_coordinates = lattice_coordinates
        # sub_lattice_types = lattice_types
        if sub_lattice_coordinates is None:
            sub_lattice_coordinates = [[0, 0, 0]]
        lattice_positions = []
        lattice_coordinates = []
        lattice_objects = []
        lattice_types = []

        # fill_prob = fill_probs[0]

        subunit_x, subunit_y, subunit_z = 1., 1., 1.
        # TODO also add positions and types
        for i in range(len(sub_lattice_coordinates)):
            xp, yp, zp = sub_lattice_coordinates[i]

            for ix in range(repeat):
                xi = ix + xp * subunit_x
                for iy in range(repeat):
                    yi = iy + yp * subunit_y
                    for iz in range(repeat):
                        zi = iz + zp * subunit_z

                        if(np.random.uniform(0.0, 1.0) <= fill_probs[i]):
                            lattice_positions.append('N/A')
                            lattice_types.append(i)
                            lattice_coordinates.append((xi, yi, zi))
                            lattice_objects.append(objects[i])

        super(
            ExtendedLattice,
            self).__init__(
            objects,
            lattice_coordinates=lattice_coordinates,
            lattice_types=lattice_types,
            lattice_positions=lattice_positions,
            lattice_spacing_a=lattice_spacing_a,
            lattice_spacing_b=lattice_spacing_b,
            lattice_spacing_c=lattice_spacing_c,
            alpha=alpha,
            beta=beta,
            gamma=gamma,
            sigma_D=sigma_D)

    def symmetry_factor(self, h, k, l):
        """Returns the symmetry factor (0 for forbidden)."""
        return 1

    def unit_cell_volume(self):
        return self.lattice_spacing_a**3