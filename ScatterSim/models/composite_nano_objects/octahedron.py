from ScatterSim.CompositeNanoObjects import CompositeNanoObject
from ScatterSim.models.nano_objects.pyramid import PyramidNanoObject


class OctahedronNanoObject(CompositeNanoObject):
    ''' An octahedron made of two pyramids.
        This is a composite nanoobject.
    '''

    def __init__(self, pargs=None):
        if pargs is None:
            pargs = {}
        objslist = list()
        pargslist = list()

        # need to correctly set up positions and angles
        # parameters : obj, x0, y0, z0, eta, theta, phi, sign (set to 1)
        parameters = [
            [PyramidNanoObject, 0, 0, 0, 0, 0, 0, 1],
            [PyramidNanoObject, 0, 0, 0, 0, 180, 0, 1],
        ]
        for i in range(len(parameters)):
            objslist.append(parameters[i][0])

            pargslist.append(pargs.copy())
            pargslist[i]['x0'] = parameters[i][1]
            pargslist[i]['y0'] = parameters[i][2]
            pargslist[i]['z0'] = parameters[i][3]
            pargslist[i]['eta'] = parameters[i][4]
            pargslist[i]['phi'] = parameters[i][5]
            pargslist[i]['theta'] = parameters[i][6]
            pargslist[i]['sign'] = parameters[i][7]

        super(OctahedronNanoObject, self).__init__(objslist, pargslist, pargs)