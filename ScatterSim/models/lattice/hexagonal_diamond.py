import numpy as np

from ScatterSim.models.lattice.hexagonal import HexagonalLattice


# Original author : Kevin G. Yager (most of these lattices are Kevin G. Yager)
# HexagonalDiamondLattice
###################################################################
class HexagonalDiamondLattice(HexagonalLattice):
    # Initialization
    ########################################
    def __init__(self, objects, lattice_spacing_a=None, lattice_spacing_b=None,
                 lattice_spacing_c=None, alpha=90, beta=90, gamma=60,
                 sigma_D=0.01):
        ''' This is for a hexagonal diamond lattice.

            Some useful relations:
                lattice_spacing_a = bond_length*np.sqrt(2/3.)
                where bond_length is the bond length. You must specify the
                positionsi in terms of lattice spacing, not bond length so use
                this conversion if needed.
        '''
        # write the lattice spacings in terms of the bond length
        # this is different from the previous conventions but stresses
        # that hexagonal diamond really depends on bond length

        if lattice_spacing_a is None:
            raise ValueError("Error, must specify at least the"
                             " a lattice spacing.")
        if lattice_spacing_b is None:
            lattice_spacing_b = lattice_spacing_a
        if lattice_spacing_c is None:
            lattice_spacing_c = lattice_spacing_a*4.0/np.sqrt(6.)

        # Define the lattice
        symmetry = {}
        symmetry['crystal family'] = 'triclinic'
        symmetry['crystal system'] = 'triclinic'
        symmetry['Bravais lattice'] = '?'
        symmetry['crystal class'] = '?'
        symmetry['point group'] = '?'
        symmetry['space group'] = '?'

        positions = ['network1', 'network2']
        lattice_positions = ['corner bottom layer', \
                                    'strut lower', \
                                    'strut higher', \
                                    'corner midlayer', \
                                ]
        # note: the vector components in x and y are 1/3 and
        # not 1/2. This is because these vectors are non orthogonal
        lattice_coordinates = [ (0.0, 0.0, 0.0), \
                                 (1.0/3.0, 1.0/3.0, 1.0/8.0), \
                                 (1.0/3.0, 1.0/3.0, 4.0/8.0), \
                                 (0.0, 0.0, 5.0/8.0), \
                               ]
        lattice_objects = [objects[0], \
                          ]
                                    #objects[0], \
                                    #objects[0], \
                                    #objects[0], \
                                #]
        lattice_types = [1,1,1,1]

        super(
            HexagonalDiamondLattice,
            self).__init__(
            objects,
            lattice_spacing_a=lattice_spacing_a,
            lattice_spacing_b=lattice_spacing_b,
            lattice_spacing_c=lattice_spacing_c,
            sigma_D=sigma_D,
            symmetry=symmetry,
            lattice_positions=lattice_positions,
            lattice_coordinates=lattice_coordinates,
            lattice_types=lattice_types,
            alpha=alpha,
            beta=beta,
            gamma=gamma)