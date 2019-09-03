import numpy as np

from ScatterSim.NanoObjects import NanoObject


class PyramidNanoObject(NanoObject):
    """ A square-based truncated pyramid nano-object. The canonical (unrotated)
    version has the square-base in the x-y plane, with the peak pointing along
    +z.  The base-edges are parallel to the x-axis and y-axis (i.e. the corners
    point at 45 degrees to axes.
    Edge length of base is 2*R

    A regular pyramid (for octahedra) will have equilateral triangles as faces
        and not be truncated.
        Thus H = sqrt(2)*R and face_angle = arctan(sqrt(2)) = 54.7356 degrees
        Only set "R" to get this regular pyramid.

    From [Kevin's website] (http://gisaxs.com/index.php/Form_Factor:Pyramid):
        For pyramid of base edge-length 2R, and height H. The angle of the
        pyramid walls is \alpha. If H < R/ \tan\alpha then the pyramid is
        truncated (flat top).

    Originally from:
    Lazzari, Rémi. "IsGISAXS: a program for grazing-incidence small-angle X-ray
    scattering analysis of supported islands." Journal of Applied
    Crystallography 35.4 (2002): 406-421.
    doi:10.1107/S0021889802006088
    """

    def __init__(self, pargs={}):
        super(PyramidNanoObject, self).__init__(pargs=pargs)

        # defaults
        if 'radius' not in self.pargs:
            self.pargs['radius'] = 1.
        if 'height' not in self.pargs:
            # Assume user wants a 'regular' pyramid
            self.pargs['height'] = np.sqrt(2.0) * self.pargs['radius']
        if 'pyramid_face_angle' not in self.pargs:
            self.pargs['pyramid_face_angle'] = np.degrees(
                np.arctan(np.sqrt(2)))

    def V(self, rvec):
        """Returns the intensity of the real-space potential at the
        given real-space coordinates."""

        # first rotate and shift
        rvec = self.map_rcoord(rvec)
        x, y, z = rvec

        R = self.pargs['radius']
        H = min(
            self.pargs['height'],
            R *
            np.tan(
                np.radians(
                    self.pargs['pyramid_face_angle'])))
        R_z = R - z / np.tan(np.radians(self.pargs['pyramid_face_angle']))
        V = ((z < H) * (z > 0) * (np.abs(x) < np.abs(R_z))
             * (np.abs(y) < np.abs(R_z))).astype(float)

        return V

    def volume(self):
        ''' Volume of a pyramid.
        '''
        # height = self.pargs['height']
        # base_len = self.pargs['radius']/np.sqrt(2)*2
        # base_area = base_len**2
        # return base_area*height/3.

        R = self.pargs['radius']
        H = self.pargs['height']
        tan_alpha = np.tan(np.radians(self.pargs['pyramid_face_angle']))
        volume = (4.0 / 3.0) * tan_alpha * (R**3 - (R - H / tan_alpha)**3)

        return volume

    def form_factor(self, qvec):
        """Returns the complex-amplitude of the form factor at the given
        q-coordinates.
        From Kevin's website

        Notes : I have checked that the q=0 scaling scales as volume and it
        does.
            So composite objects (sums and differences) should scale properly.
            For example, for a regular pyramid (equilateral triangles), the
            volume is 4*sqrt(2)/3*R^3 where 2*R is the edge length at the base.
            I checked and the q=0 scattering is exactly volume^2 as it should
            be (assuming other prefactors are 1).
        """

        # Phase must be retrieved *before* mapping in q
        phase = self.get_phase(qvec)

        qvec = self.map_qcoord(qvec)
        # fix divide by zero errors, threshold values (abs val)
        self.thresh_array(qvec, 1e-6)

        qx, qy, qz = qvec

        R = self.pargs['radius']
        H = self.pargs['height']
        tan_alpha = np.tan(np.radians(self.pargs['pyramid_face_angle']))
        amod = 1.0 / tan_alpha
        # volume = self.volume()

        q1 = 0.5 * ((qx - qy) * amod + qz)
        q2 = 0.5 * ((qx - qy) * amod - qz)
        q3 = 0.5 * ((qx + qy) * amod + qz)
        q4 = 0.5 * ((qx + qy) * amod - qz)

        K1 = np.sinc(q1 * H / np.pi) * np.exp(+1.0j * q1 * H) + \
            np.sinc(q2 * H / np.pi) * np.exp(-1.0j * q2 * H)
        K2 = -1.0j * np.sinc(q1 * H / np.pi) * np.exp(+1.0j * q1 * H) + \
            1.0j * np.sinc(q2 * H / np.pi) * np.exp(-1.0j * q2 * H)
        K3 = np.sinc(q3 * H / np.pi) * np.exp(+1.0j * q3 * H) + \
            np.sinc(q4 * H / np.pi) * np.exp(-1.0j * q4 * H)
        K4 = -1.0j * np.sinc(q3 * H / np.pi) * np.exp(+1.0j * q3 * H) + \
            1.0j * np.sinc(q4 * H / np.pi) * np.exp(-1.0j * q4 * H)

        F = (H / (qx * qy)) *\
            (K1 * np.cos((qx - qy) * R)
             + K2 * np.sin((qx - qy) * R)
             - K3 * np.cos((qx + qy) * R)
             - K4 * np.sin((qx + qy) * R))
        F *= phase*self.pargs['delta_rho']

        return F