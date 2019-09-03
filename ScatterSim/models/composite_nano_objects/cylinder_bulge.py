from ScatterSim.CompositeNanoObjects import CompositeNanoObject
from ScatterSim.models.nano_objects.cylinder import CylinderNanoObject
from ScatterSim.models.nano_objects.sphere import SphereNanoObject


class CylinderBulgeNanoObject(CompositeNanoObject):
    ''' A nano object comprised of a cylinder with a sphere bulge in center.'''

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
        sph_den = pargs.pop('sphere_density')

        pargs_cyl = {
            'radius': cyl_rad,
            'height': cyl_height,
            'x0': 0,
            'y0': 0,
            'z0': 0,
            'rho_object': cyl_den
        }

        pargs_sphere = {
            'radius': sph_rad,
            'x0': 0,
            'y0': 0,
            'z0': 0,
            'rho_object': sph_den,
        }

        pargslist = [pargs_cyl, pargs_sphere]
        objslist = [CylinderNanoObject, SphereNanoObject]

        super(
            CylinderBulgeNanoObject,
            self).__init__(
            objslist,
            pargslist,
            pargs=pargs)