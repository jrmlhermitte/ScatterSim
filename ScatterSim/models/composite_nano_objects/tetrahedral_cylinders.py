import numpy as np

from ScatterSim.CompositeNanoObjects import CompositeNanoObject
from ScatterSim.models.nano_objects.cylinder import CylinderNanoObject


class TetrahedraCylindersNanoObject(CompositeNanoObject):
    ''' A nano object made up of cylinders in a tetrahedral shape.
        NOTE : Needs to be verified...
        The canonical unrotated version is aligned along the z-axis.
        The top inverted tetrahedra has one edge aligned with the x axis and
        one vertex pointing to the +y direction.
        The bottom tetrahedra is a mirror image of this.
    '''
    # TODO : Add density as parameter

    def __init__(self, pargs={}):
        ''' Initalized the object.
        '''
        objslist = list()
        pargslist = list()

        # use default packing if not specified
        r = pargs["radius"]
        L = pargs["edgelength"]
        H = pargs["height"]
        # the triangle angles (degrees)
        beta = 30
        # the triangle angles (radians)
        beta_r = np.radians(beta)
        # the face angle (radians)
        alpha_r = np.arccos(.5/np.cos(beta_r))
        # the face angle (degrees)
        alpha = np.degrees(alpha_r)
        # midway point along height
        D = L/2./np.cos(beta_r)

        # midway to intersection points
        A = D/2. # also it is: L/2.*np.tan(beta_r)
        # note this is NOT the height parameter

        # in order of x, y, z, eta, phi, theta
        posanglelist = [
                   # the base
                   [0, -A, 0, 0, 90, 0],
                   [L/2*(1-np.sin(beta_r)), L/2*(np.cos(beta_r) - np.tan(beta_r)), 0, 0, 90, 60],
                   [-L/2*(1-np.sin(beta_r)), L/2*(np.cos(beta_r) - np.tan(beta_r)), 0, 0, 90, -60],
                   # the height
                   [0,D - L/2*np.cos(alpha_r), L/2*np.sin(alpha_r), 0, 90-alpha, -90],
                   [L*np.cos(alpha_r)*np.cos(beta_r)-L/2*np.cos(alpha_r)*np.cos(beta_r),
                      -L*np.cos(alpha_r)*np.sin(beta_r)+L/2*np.cos(alpha_r)*np.sin(beta_r),
                      L/2*np.sin(alpha_r), 0, 90-alpha, 30],
                   [-L*np.cos(alpha_r)*np.cos(beta_r)+L/2*np.cos(alpha_r)*np.cos(beta_r),
                      -L*np.cos(alpha_r)*np.sin(beta_r)+L/2*np.cos(alpha_r)*np.sin(beta_r),
                      L/2*np.sin(alpha_r), 0, -(90-alpha), -30],
        ]


        eta, theta, phi = 0, np.pi/4., 0


        # the vector pointing to center of rotation for this construction
        Rc = L*np.sqrt(3/8.)
        tetra_height = L*np.sin(alpha_r)
        z_center = tetra_height - Rc
        COM = np.array([0, 0, -z_center])


        pargs_cyl = {
            'radius': r,
            'height': H,
            'x0': 0,
            'y0': 0,
            'z0': 0,
        }

        pargslist = list()
        objslist = list()
        for x0, y0, z0, eta, phi, theta in posanglelist:
            pargstmp = pargs_cyl.copy()
            pargstmp['x0'] = x0 + COM[0]
            pargstmp['y0'] = y0 + COM[1]
            pargstmp['z0'] = z0 + COM[2]
            pargstmp['eta'] = eta
            pargstmp['phi'] = phi
            pargstmp['theta'] = theta
            pargslist.append(pargstmp)
            objslist.append(CylinderNanoObject)

        super(
            TetrahedraCylindersNanoObject,
            self).__init__(
            objslist,
            pargslist,
            pargs=pargs)