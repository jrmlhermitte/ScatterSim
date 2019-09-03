import numpy as np

from ScatterSim.CompositeNanoObjects import CompositeNanoObject
from ScatterSim.models.composite_nano_objects.tetrahedral_cylinders import TetrahedraCylindersNanoObject


class HourglassCylindersNanoObject(CompositeNanoObject):
    ''' A nano object made up of cylinders in a tetrahedral shape.'''
    # TODO : Add density as parameter

    def __init__(self, pargs={}):
        ''' Initalized the object.
        '''

        # use default packing if not specified
        L = pargs["edgelength"]
        Rc = L*np.sqrt(3/8.)

        pargs_bot = pargs.copy()
        pargs_bot['z0'] = -Rc
        pargs_top = pargs.copy()
        pargs_top['phi'] = 180
        pargs_top['theta'] = 180
        pargs_top['z0'] = Rc

        objslist = [TetrahedraCylindersNanoObject,
                    TetrahedraCylindersNanoObject]
        pargslist = [pargs_bot, pargs_top]

        super(
            HourglassCylindersNanoObject,
            self).__init__(
            objslist,
            pargslist,
            pargs=pargs)