from ScatterSim.LatticeObjects import Lattice


class SimpleCubic(Lattice):
    def __init__(
            self,
            objects,
            lattice_spacing_a=1.0,
            sigma_D=0.01,
            lattice_coordinates=None,
            lattice_types=None):
        # prepare variables of lattice
        symmetry = {
            'crystal family': 'cubic',
            'crystal system': 'cubic',
            'Bravais lattice': 'P',
            'crystal class': 'hexoctahedral',
            'point group': 'm3m',
            'space group': 'Pm3m',
        }
        if not isinstance(objects, list):
            objects = [objects]

        lattice_positions = ['placement'] * len(objects)

        if lattice_coordinates is None:
            lattice_coordinates = [(0.0, 0.0, 0.0)] * len(objects)
        lattice_types = list(range(len(objects)))

        # now call parent function to initialize
        super(
            SimpleCubic,
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