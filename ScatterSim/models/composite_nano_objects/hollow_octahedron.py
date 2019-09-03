from ScatterSim.CompositeNanoObjects import CompositeNanoObject
from ScatterSim.models.composite_nano_objects.octahedron import OctahedronNanoObject


class HollowOctahedronNanoObject(CompositeNanoObject):
    ''' An octahedron made of two pyramids.
        This is a composite nanoobject. Composite nanoobjects are just lists
            of nano objects. They're convenient in that they don't need any
            extra routines added. Everything should just work with the
            initializer and inheriting the rest.

        outer_radius : outer radius of octahedron
        radius_ratio : ratio of inner to outer radius (inner/outer). I chose it
            this way so that this could be supplied to a 1 parameter
            polydisperse object. Note this should be less than 1.

        inner_radius : inner radius of octahedron (not a used parg)

    '''

    def __init__(self, pargs={}):
        objslist = list()
        pargslist = list()

        if 'radius' not in pargs:
            raise ValueError("Need to specify the outer radius")
        if 'radius_ratio' not in pargs:
            raise ValueError("Need to specify ratio of inner to outer radius")

        outer_radius = pargs['radius']
        inner_radius = outer_radius * pargs['radius_ratio']
        # parameters : obj, x0, y0, z0, eta, theta, phi, sign (1 adds, -1
        # subtracts)
        parameters = [
            [OctahedronNanoObject, 0, 0, 0, 0, 0, 0, inner_radius, -1],
            [OctahedronNanoObject, 0, 0, 0, 0, 0, 0, outer_radius, 1],
        ]
        for i in range(len(parameters)):
            objslist.append(parameters[i][0])

            pargslist.append(pargs.copy())
            pargslist[i]['x0'] = parameters[i][1]
            pargslist[i]['y0'] = parameters[i][2]
            pargslist[i]['z0'] = parameters[i][3]
            pargslist[i]['eta'] = parameters[i][4]
            pargslist[i]['theta'] = parameters[i][5]
            pargslist[i]['phi'] = parameters[i][6]
            pargslist[i]['radius'] = parameters[i][7]
            pargslist[i]['sign'] = parameters[i][8]

        super(
            HollowOctahedronNanoObject,
            self).__init__(
            objslist,
            pargslist,
            pargs)