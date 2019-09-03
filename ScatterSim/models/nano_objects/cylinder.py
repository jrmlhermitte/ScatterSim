import numpy as np

from ScatterSim.NanoObjects import NanoObject


class CylinderNanoObject(NanoObject):
    """A cylinder nano-object. The canonical (unrotated) version
    has the circular-base in the x-y plane, with the length along z.

    self.pargs contains parameters:
        rho_ambient : the cylinder density
        rho_object : the solvent density
        radius : (default 1.0) the cylinder radius
        length : (default 1.0) the cylinder length

        eta,phi,eta: Euler angles
        x0, y0, z0 : the position of cylinder COM relative to origin
        The object is rotated first about origin, then translated to
            where x0, y0, and z0 define it to be.

    these are calculated after the fact:
        delta_rho : rho_ambient - rho1
    """

    def __init__(self, pargs={}):
        super(CylinderNanoObject, self).__init__(pargs=pargs)

        if 'radius' not in self.pargs:
            self.pargs['radius'] = 1.0

        if 'height' not in self.pargs:
            self.pargs['height'] = 1.0

    def V(self, rvec):
        """Returns the intensity of the real-space potential at the
        given real-space coordinates.
        Returns 1 if in the space, 0 otherwise.
        Can be arrays.
        Rotate then translate.

            rotation_matrix is an extra rotation to add on top of the built
            in rotation (from eta, phi, theta elements in object)
        """
        R = self.pargs['radius']
        L = self.pargs['height']

        rvec = self.map_rcoord(rvec)
        # could be in one step, but making explicit
        x, y, z = rvec

        r = np.hypot(x, y)
        result = np.zeros(x.shape)
        w = np.where((z <= L / 2.) * (z >= -L / 2.) * (np.abs(r) <= R))
        if len(w[0]) > 0:
            result[w] = self.pargs['delta_rho']

        return result

    def volume(self):
        return np.pi * self.pargs['radius']**2 * self.pargs['height']

    def form_factor(self, qvec):
        """Returns the complex-amplitude of the form factor at the given
        q-coordinates.

            The first set of q are assumed to be that of the root qx,qy,qz
            The next sets are the parent coordinates. These are only supplied
            if the parent is not the root The root q is only necessary for the
            form factor, where the phase factor needs to be calculated with
            respect to the origin of that coordinate system.
        """
        # Phase must be retrieved *before* mapping in q
        phase = self.get_phase(qvec)

        # first rotate just as for V
        qvec = self.map_qcoord(qvec)
        self.thresh_array(qvec, 1e-4)
        qx, qy, qz = qvec

        R = self.pargs['radius']
        H = self.pargs['height']
        volume = self.volume()

        qr = np.hypot(qx, qy)

        # NOTE : Numpy's sinc function adds
        # a factor of pi in we need to remove.
        # Why numpy... why??? ><
        F = 2 * np.sinc(qz * H / 2. / np.pi) * j1(qr * R) / qr / R + 1j * 0
        F *= phase
        F *= self.pargs['delta_rho'] * volume

        return F