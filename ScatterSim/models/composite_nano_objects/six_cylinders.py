import numpy as np

from ScatterSim.CompositeNanoObjects import CompositeNanoObject
from ScatterSim.models.nano_objects.cylinder import CylinderNanoObject


class SixCylinderNanoObject(CompositeNanoObject):
    ''' A nano object comprised of six cylinders.'''

    def __init__(self, pargs={}):
        ''' Initalized the object.
        '''
        objslist = list()
        pargslist = list()

        # use default packing if not specified
        if 'packradius' not in pargs:
            pargs['packradius'] = 2 * pargs['radius'] / (2 * np.pi / 6.)

        packradius = pargs['packradius']

        dphi = 2 * np.pi / 6.
        for i in range(6):
            x0, y0, z0 = packradius * \
                np.cos(i * dphi), packradius * np.sin(i * dphi), 0
            pargs_tmp = {
                'radius': pargs['radius'],
                'height': pargs['height'],
                'x0': x0,
                'y0': y0,
                'z0': z0,
            }
            pargslist.append(pargs_tmp)
            objslist.append(CylinderNanoObject)

        super(
            SixCylinderNanoObject,
            self).__init__(
            objslist,
            pargslist,
            pargs=pargs)