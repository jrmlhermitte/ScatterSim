import numpy as np

from ScatterSim.CompositeNanoObjects import CompositeNanoObject
from ScatterSim.models.nano_objects.cylinder import CylinderNanoObject


class OctahedronCylindersNanoObject(CompositeNanoObject):
    """An octahedron object made of cylinders or like object. The canonical
        (unrotated) version has the square cross-section in the x-y plane, with
        corners pointing along +z and -z.  The corners are on the x-axis and
        y-axis. The edges are 45 degrees to the x and y axes.
    The canonical (unrotated) version of the cylinders should be aligned along
    the z axis. Replace cylinders with spheres, cylindrical shells, prolate
    ellipsoids etc as you wish.

            It is best to just brute force define all the terms in one shot,
                which I chose to do here.
            notes about rotations:
                - i need to dot the rotation matrix of the octahedron with each
                    individual cylinder's rotationo element
            edgelength : length of edge
            edgespread : how much to expand the rods by (not a good name)
                positive is expansion, negative is compression
            edgesep : separation of element from edge
            rest of pargs are used for the cylinder object
                ex : radius, height
            linkerlength : if specified, adds linkers of specified length,
                centered in between the octahedra
            linkerradius : if linkerlength specified, will add linkers of this
            radius
            rho_linker : if linkerlength specified, adds linkers of this
            density (defaults to same density as cylinders in octahedra)
            linkerobject : the object to use for the linkers (defaults to
            baseObject)
    """

    def __init__(
            self,
            baseObject=None,
            linkerObject=None,
            pargs={},
            seed=None):
        if baseObject is None:
            baseObject = CylinderNanoObject
        if linkerObject is None:
            linkerObject = baseObject

        # Set defaults
        if 'edgeshift' not in pargs:
            pargs['edgeshift'] = 0.0
        if 'edgespread' not in pargs:
            pargs['edgespread'] = 0.0
        if 'linkerlength' in pargs:
            addlinkers = True
        else:
            addlinkers = False

        # raise errors for undefined parameters
        if 'edgelength' not in pargs:
            raise ValueError("Need to specify an edgelength for this object")

        # these are slight shifts per cyl along the axis
        # positive is away from COM and negative towards
        shiftlabels = [
            # these correspond to the poslist
            'CYZ1', 'CXZ1', 'CYZ2', 'CXZ2',
            'CXY1', 'CXY4', 'CXY3', 'CXY2',
            'CYZ3', 'CXZ3', 'CYZ4', 'CXZ4',
            'linker1', 'linker2', 'linker3', 'linker4',
            'linker5', 'linker6',
        ]

        # you flip x or y from original shifts to move along edge axis
        # not a good explanation but some sort of personal bookkeeping for
        # now...
        shiftfacs = [
            # top
            [0, -1, 1],
            [-1, 0, 1],
            [0, 1, 1],
            [1, 0, 1],
            # middle
            [-1, 1, 0],
            [-1, -1, 0],
            [1, -1, 0],
            [1, 1, 0],
            # bottom
            [0, 1, -1],
            [1, 0, -1],
            [0, -1, -1],
            [-1, 0, -1]
        ]

        for lbl1 in shiftlabels:
            if lbl1 not in pargs:
                pargs[lbl1] = 0.

        # calculate shift of COM from edgelength and edgespread
        fac1 = np.sqrt(2) / 2. * \
            ((.5 * pargs['edgelength']) + pargs['edgespread'])
        eL = pargs['edgelength']
        if addlinkers:
            sL = pargs['linkerlength']

        poslist = [
            # eta, theta, phi, x0, y0, z0
            # top part
            [0, 45, -90, 0, fac1, fac1],
            [0, 45, 0, fac1, 0, fac1],
            [0, 45, 90, 0, -fac1, fac1],
            [0, -45, 0, -fac1, 0, fac1],
            # now the flat part
            [0, 90, 45, fac1, fac1, 0],
            [0, 90, -45, fac1, -fac1, 0],
            [0, 90, 45, -fac1, -fac1, 0],
            [0, 90, -45, -fac1, fac1, 0],
            # finally bottom part
            [0, 45, -90, 0, -fac1, -fac1],
            [0, 45, 0, -fac1, 0, -fac1],
            [0, 45, 90, 0, fac1, -fac1],
            [0, -45, 0, fac1, 0, -fac1],
        ]

        if addlinkers:
            poslist_linker = [
                # linkers
                [0, 0, 0, 0, 0, eL / np.sqrt(2) + sL / 2.],
                [0, 0, 0, 0, 0, -eL / np.sqrt(2) - sL / 2.],
                [0, 90, 0, eL / np.sqrt(2) + sL / 2., 0, 0],
                [0, 90, 0, -eL / np.sqrt(2) - sL / 2., 0, 0],
                [0, 90, 90, 0, eL / np.sqrt(2) + sL / 2., 0],
                [0, 90, 90, 0, -eL / np.sqrt(2) - sL / 2., 0],
            ]
            for row in poslist_linker:
                poslist.append(row)
            shiftfacs_linker = [
                [0, 0, 1],
                [0, 0, -1],
                [1, 0, 0],
                [-1, 0, 0],
                [0, 1, 0],
                [0, -1, 0],
            ]
            for row in shiftfacs_linker:
                shiftfacs.append(row)

        poslist = np.array(poslist)
        shiftfacs = np.array(shiftfacs)

        # now add the shift factors
        for i in range(len(poslist)):
            poslist[i, 3:] += np.sqrt(2) / 2. * \
                shiftfacs[i] * pargs[shiftlabels[i]]

        # need to create objslist and pargslist
        objlist = list()
        pargslist = list()
        for i, pos in enumerate(poslist):
            objlist.append(baseObject)

            eta, phi, theta, x0, y0, z0 = pos
            pargstmp = dict()
            pargstmp.update(pargs)
            pargstmp['eta'] = eta
            pargstmp['theta'] = theta
            pargstmp['phi'] = phi
            pargstmp['x0'] = x0
            pargstmp['y0'] = y0
            pargstmp['z0'] = z0

            labeltmp = shiftlabels[i]
            if 'linker' in labeltmp:
                if 'rho_linker' in pargs:
                    pargstmp['rho_object'] = pargs['rho_linker']
                if 'linkerlength' in pargs:
                    pargstmp['height'] = pargs['linkerlength']
                if 'linkerradius' in pargs:
                    pargstmp['radius'] = pargs['linkerradius']

            pargslist.append(pargstmp)

        super(
            OctahedronCylindersNanoObject,
            self).__init__(
            objlist,
            pargslist,
            pargs=pargs)