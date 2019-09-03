from ScatterSim.CompositeNanoObjects import CompositeNanoObject
from ScatterSim.models.nano_objects.cylinder import CylinderNanoObject
from ScatterSim.models.nano_objects.sphere import SphereNanoObject


class CylinderDumbbellNanoObject(CompositeNanoObject):
    ''' A nano object comprised of a cylinder with two spheres on edges (like a
    dummbell).'''

    def __init__(self, pargs={}):
        ''' Initalized the object.
        '''
        objslist = list()
        pargslist = list()

        # use default packing if not specified

        cyl_rad = pargs.pop('cylinder_radius')
        cyl_den = pargs.pop('cylinder_density')
        cyl_height = pargs.pop('cylinder_height')
        sph_rad = pargs.pop('sphere_radius')
        # sph_den = pargs.pop('sphere_density')

        pargs_cyl = {
            'radius': cyl_rad,
            'height': cyl_height,
            'x0': 0,
            'y0': 0,
            'z0': 0,
            'rho_ambient': cyl_den
        }

        el = (cyl_height) * .5

        pargs_sphere1 = {
            'radius': sph_rad,
            'x0': 0,
            'y0': 0,
            'z0': -el,
        }

        pargs_sphere2 = {
            'radius': sph_rad,
            'x0': 0,
            'y0': 0,
            'z0': el,
        }

        pargslist = [pargs_cyl, pargs_sphere1, pargs_sphere2]
        objslist = [CylinderNanoObject, SphereNanoObject, SphereNanoObject]

        super(
            CylinderDumbbellNanoObject,
            self).__init__(
            objslist,
            pargslist,
            pargs=pargs)