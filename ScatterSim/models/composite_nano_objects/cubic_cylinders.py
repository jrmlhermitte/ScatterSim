from ScatterSim.CompositeNanoObjects import CompositeNanoObject
from ScatterSim.models.nano_objects.cylinder import CylinderNanoObject


class CubicCylindersNanoObject(CompositeNanoObject):
    """ A cylinder nano object aligned along axes, meant for cubic structure.
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

        # raise errors for undefined parameters
        if 'height' not in pargs:
            raise ValueError("Need to specify a height for this object")
        if 'radius' not in pargs:
            raise ValueError("Need to specify a radius for this object")

        eL = pargs['height']

        poslist = [
            # eta, theta, phi, x0, y0, z0
            # top part
            [0, 0, 0, 0, 0, eL / 2.],
            [0, 90, 0, eL / 2., 0, 0],
            [0, 90, 90, 0, eL / 2., 0],
        ]

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

            pargslist.append(pargstmp)

        super(
            CubicCylindersNanoObject,
            self).__init__(
            objlist,
            pargslist,
            pargs=pargs)