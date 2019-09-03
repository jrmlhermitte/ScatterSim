import numpy as np

from ScatterSim.CompositeNanoObjects import CompositeNanoObject
from ScatterSim.models.nano_objects.sphere import SphereNanoObject


class CoreShellNanoObject(CompositeNanoObject):
    ''' A core shell nano object.
        This is a composite nanoobject.

        pargs:
            radius_inner : inner radius
            rho_object_inner : inner density
            radius_outer : outer radius
            rho_object_outer : outer density
            rho_ambient : the ambient density (solvent)

            Note : all rho's must be *positive*

        Notes
        -----
        Use of non-zero ambient has not been tested.
        This has been thought through but should be tested.

    '''

    def __init__(self, pargs={}):
        objslist = list()
        pargslist = list()

        objslist = [SphereNanoObject, SphereNanoObject]

        # check for arguments
        if 'rho_object_inner' not in pargs or\
            'radius_inner' not in pargs or\
            'rho_object_outer' not in pargs or\
            'radius_outer' not in pargs or\
                'rho_ambient' not in pargs:
            raise ValueError("Missing args, please check correct syntax")

        deltarho_inner = pargs['rho_object_inner'] - pargs['rho_object_outer']
        deltarho_outer = pargs['rho_object_outer'] - pargs['rho_ambient']

        # get the signs and take absolute value
        sign_inner = 1
        if deltarho_inner < 0:
            sign_inner = -1
            deltarho_inner = np.abs(deltarho_inner)

        sign_outer = 1
        if deltarho_outer < 0:
            sign_outer = -1
            deltarho_outer = np.abs(deltarho_outer)

        pargs_inner = {
            'radius': pargs['radius_inner'],
            'rho_object': deltarho_inner,
            'rho_ambient': 0,
            'sign': sign_inner
        }

        pargs_outer = {
            'radius': pargs['radius_outer'],
            'rho_object': deltarho_outer,
            'rho_ambient': 0,
            'sign': sign_outer
        }

        pargslist = [pargs_inner, pargs_outer]

        super(
            CoreShellNanoObject,
            self).__init__(
            objslist,
            pargslist,
            pargs=pargs)