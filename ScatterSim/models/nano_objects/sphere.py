import numpy as np

from ScatterSim.NanoObjects import NanoObject


class SphereNanoObject(NanoObject):
    ''' This is as the name of object describes, a sphere.'''

    def __init__(self, pargs={}):
        super(SphereNanoObject, self).__init__(pargs=pargs)

        if 'radius' not in self.pargs:
            # Set a default size
            self.pargs['radius'] = 1.0

    def V(self, rvec):
        """Returns the intensity of the real-space potential at the
        given real-space coordinates."""

        rvec = self.map_rcoord(np.array(rvec))
        R = self.pargs['radius']
        r = np.sqrt(rvec[0]**2 + rvec[1]**2 + rvec[2]**2)
        return (r < R).astype(float) * self.pargs['delta_rho']

    def volume(self):
        return 4 / 3. * np.pi * self.pargs['radius']**3

    def form_factor(self, qvec):
        ''' Compute the form factor of a sphere. '''
        phase = self.get_phase(qvec)

        qx, qy, qz = qvec
        R = self.pargs['radius']

        volume = self.volume()

        q = np.sqrt(qx**2 + qy**2 + qz**2)
        qR = q * R

        # threshold to avoid zero values
        qR = np.maximum(1e-8, qR)

        # use: 3* j1(qR)/qR (let scipy handle it)
        # numerically more stable than its equivalent:
        # 3*( np.sin(qR) - qR*np.cos(qR) )/( qR**3 )

        F = self.pargs['delta_rho'] * volume * \
            3 * spherical_jn(1, qR) / qR * phase

        return F

    # override some complex functions to save time
    def form_factor_isotropic(self, q, num_phi=50, num_theta=50):
        """Returns the particle form factor, averaged over every possible orientation.
        """
        return self.form_factor(np.array([q, 0, 0]))

    def form_factor_squared_isotropic(self, q, num_phi=50, num_theta=50):
        """Returns the particle form factor squared, averaged over every
        possible orientation.
        """
        # numpy should broadcast 0 to same length as q
        return self.form_factor_squared(np.array([q, 0, 0]))

    def to_string(self):
        """Returns a string describing the object."""
        s = "SphereNanoObject: A sphere of radius = %.3f nm" % \
            self.pargs['radius']
        return s

    def to_short_string(self):
        """Returns a short string describing the object's variables.
        (Useful for distinguishing objects of the same class.)"""
        s = "R = %.3f nm" % self.pargs['radius']
        return s