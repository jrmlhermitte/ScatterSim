import numpy as np
# cylndrical and spherical Bessel functions


# NanoObject
####################################################################
class NanoObject:
    """Defines a nano-object, which can then be placed within a lattice
    for computing scattering data. A nano-object can be anisotropic.

        This is the base class all nano-objects should inherit
    """
    conversion_factor = 1E-4        # Converts units from 1E-6 A^-2 into nm^-2

    def __init__(self, pargs=None):
        '''
            pargs are the potential arguments, a dictionary of arguments
        '''
        if pargs is None:
            pargs = {}
        # look for and set defaults
        eta = self.check_arg('eta', pargs, default=0)
        theta = self.check_arg('theta', pargs, default=0)
        phi = self.check_arg('phi', pargs, default=0)
        x0 = self.check_arg('x0', pargs, default=0)
        y0 = self.check_arg('y0', pargs, default=0)
        z0 = self.check_arg('z0', pargs, default=0)

        # now update object's pargs with pargs
        self.check_arg('rho_ambient', pargs, default=0.0)
        self.check_arg('rho_object', pargs, default=15.0)

        self.pargs = dict()
        self.pargs.update(pargs)
        # delta rho is difference between sample and solvent density
        self.pargs['delta_rho'] = abs(self.pargs['rho_ambient'] -
                                      self.pargs['rho_object'])

        self.set_origin(x0, y0, z0)
        # this will set the rotation matrix
        self.set_angles(eta=eta, phi=phi, theta=theta)

    def __add__(self, other):
        ''' The add operator should now
            just create a new composite NanoObject.
            There are four cases:
                NanoObject + NanoObject
                NanoObject + CompositeNanoObject
                CompositeNanoObject + CompositeNanoObject
                CompositeNanoObject + NanoObject

            # Note : the pargs of the new object is overridden by original
            # object
        '''
        from .CompositeNanoObjects import CompositeNanoObject  # noqa
        if isinstance(self, CompositeNanoObject):
            nano_objects = self.nano_objects
            pargs = self.pargs
        else:
            nano_objects = [self]
            # no pargs for parent composite object set yet
            pargs = {}

        nano_objects.append(other)

        # TODO : should move to composite nano object code
        obj = CompositeNanoObject(nano_objects, pargs=pargs)
        obj.pargs['sign'] = 1.
        return obj

    def check_arg(self, name, pargs, default=0):
        ''' Check dictionary for a parameter. If not there, set the parg
            to the default.

            Returns the result in both cases
        '''
        if name in pargs:
            return pargs[name]
        else:
            pargs[name] = default
            return default

    def rebuild(self, pargs={}):
        """Allows the object to have its potential
        arguments (pargs) updated. Note that this doesn't
        replace the old pargs entirely. It only modifies
        (or adds) the key/values provided by the new pargs."""
        self.pargs.update(pargs)
        # also recompute the delta_rho
        self.pargs['delta_rho'] = abs(self.pargs['rho_ambient'] -
                                      self.pargs['rho_object'])

    def set_angles(self, eta=None, phi=None, theta=None):
        """Update one or multiple orientation angles (degrees).
            These are the typical Euler coordinates.

            See Also
            -------
                rotation_elements
        """
        if eta is not None:
            self.pargs['eta'] = np.copy(eta)
        if phi is not None:
            self.pargs['phi'] = np.copy(phi)
        if theta is not None:
            self.pargs['theta'] = np.copy(theta)

        self.rotation_matrix = self.rotation_elements(self.pargs['eta'],
                                                      self.pargs['phi'],
                                                      self.pargs['theta'])

    def set_origin(self, x0=None, y0=None, z0=None):
        ''' Set the origin of the sample.
        '''
        if x0 is not None:
            self.pargs['x0'] = np.copy(x0)
        else:
            x0 = self.pargs['x0']
        if y0 is not None:
            self.pargs['y0'] = np.copy(y0)
            y0 = self.pargs['y0']
        else:
            y0 = self.pargs['y0']
        if z0 is not None:
            self.pargs['z0'] = np.copy(z0)
            z0 = self.pargs['z0']
        else:
            z0 = self.pargs['z0']

        self.origin = np.array([x0, y0, z0])

    def thresh_array(self, r, val):
        ''' threshold array to have minimum value.'''
        w = np.where(np.abs(r) < val)
        if len(w[0]) > 0:
            r[w] = val

    def rotation_elements(self, eta, phi, theta):
        """Converts angles into an appropriate rotation matrix.

        Three-axis rotation:
            1. Rotate about +z by eta (counter-clockwise in x-y plane)
            2. Tilt by phi with respect to +z (rotation about y-axis,
            clockwise in x-z plane) then
            3. rotate by theta in-place (rotation about z-axis,
            counter-clockwise in x-y plane)
        """

        eta = np.radians(eta)
        phi = np.radians(phi)
        theta = np.radians(theta)

        c1 = np.cos(eta)
        c2 = np.cos(phi)
        c3 = np.cos(theta)
        s1 = np.sin(eta)
        s2 = np.sin(phi)
        s3 = np.sin(theta)

        rotation_elements = np.array([
            [c1 * c2 * c3 - s1 * s3, -c3 * s1 - c1 * c2 * s3, c1 * s2],
            [c1 * s3 + c2 * c3 * s1, c1 * c3 - c2 * s1 * s3, s1 * s2],
            [-c3 * s2, s2 * s3, c2],
        ])

        return rotation_elements

    def get_phase(self, qvec):
        ''' Get the phase factor from the shift, for the q position.'''
        phase = np.exp(1j * qvec[0] * self.pargs['x0'])
        phase *= np.exp(1j * qvec[1] * self.pargs['y0'])
        phase *= np.exp(1j * qvec[2] * self.pargs['z0'])

        return phase

    def map_qcoord(self, qcoord):
        ''' Map the reciprocal space coordinates from the parent object to
        child object (this one).  The origin shift is not needed here, and so
        this function basically just rotates the coordinates given, using the
        internal rotation matrix.
        Translation is a phase which is computed separately.

            Parameters
            ----------
            qcoord : float array, the q coordinates to rotate.
                Leftermost index are the q components: [qx, qy, qz] where qx,
                qy, qz are any dimension

            Returns
            -------
            qcoord : the rotated coordinates

            See Also
            --------
            map_coord (to specify arbitrary rotation matrix)
            set_origin
            set_angles
        '''
        return self.map_coord(qcoord, self.rotation_matrix)

    def map_rcoord(self, rcoord):
        ''' Map the real space coordinates from the parent object to child
        object's coordinates (this one), using the internal rotation matrix and
        origin (set by set_angles, set_origin). This sets the 0 and orientation
        of sample to align with internal orientation. This function is
        generally recursively called for every traversal into child
        coordinates.

            Parameters
            ----------
            rcoord : float array, the r coordinates to map
                Leftermost index are the r components: [x, y, z] where
                x, y, z are any dimension

            Returns
            -------
            qcoord : the rotated coordinates

            See Also
            --------
            map_coord (to specify your own rotations/translations)
            set_origin
            set_angles
            '''
        # tmp fix for objects that dont  set origin
        if not hasattr(self, 'origin'):
            raise ValueError(
                "Origin has not been set (please set origin to something"
                " like (0,0,0) in the object initialization")
        # update the origin from the x0, y0, z0 values
        self.set_origin()

        return self.map_coord(rcoord, self.rotation_matrix, origin=self.origin)

    def map_coord(self, coord, rotation_matrix, origin=None):
        ''' Map the real space coordinates from the parent object to child
        object's coordinates (this one), using the specified rotation matrix
        and origin. This sets the 0 and orientation of sample to align with
        internal orientation. This function is generally recursively called for
        every traversal into child coordinates.

            Parameters
            ----------
            coord : float array, the r coordinates to map
            rotation_matrix : 3x3 ndarray, the rotation matrix
            origin : length 3 array (optional), the origin (default [0,0,0])

            Returns
            -------
            coord : the rotated coordinates

            Notes
            -----
            Most objects will use map_rcoord and map_qcoord
            instead, which know what to do with the internal rotation
            matrices of the objects.

            The slowest varying index of r is the coordinate

            See Also
            --------
            map_coord (to specify your own rotations/translations)
            set_origin
            set_angles
            '''
        if coord.shape[0] != 3:
            raise ValueError("Error slowest varying dimension is not coord"
                             "(3)")

        if origin is None:
            x0, y0, z0 = 0, 0, 0
        else:
            x0, y0, z0 = origin

        # first subtract origin
        # save time, don't do it no origin specified
        if origin is not None:
            coord = np.copy(coord)
            coord[0] = coord[0] - x0
            coord[1] = coord[1] - y0
            coord[2] = coord[2] - z0

        # next dot product
        coord = np.tensordot(rotation_matrix, coord, axes=(1, 0))

        return coord

    def form_factor_numerical(self, qvector, num_points=100, size_scale=None):
        ''' This is a brute-force calculation of the form-factor, using the
        realspace potential. This is computationally intensive and should be
        avoided in preference to analytical functions which are put into the
        "form_factor(qx,qy,qz)" function.

            Parameters
            ----------
            qvector : float arrays the reciprocal space coordinates
            num_points : int, optional, the number of points to sample
            rotation_elements : rotation matrix

            Returns
            -------
            coord : the complex form factor

            Note : NOT TESTED YET
            '''
        # TODO : TEST THIS FUNCTION
        qvector = self.map_qcoord(qvector)

        if size_scale is None:
            if 'radius' in self.pargs:
                size_scale = 2.0 * self.pargs['radius']
            else:
                size_scale = 2.0

        x_vals, dx = np.linspace(-size_scale, size_scale,
                                 num_points, endpoint=True, retstep=True)
        y_vals, dy = np.linspace(-size_scale, size_scale,
                                 num_points, endpoint=True, retstep=True)
        z_vals, dz = np.linspace(-size_scale, size_scale,
                                 num_points, endpoint=True, retstep=True)

        dVolume = dx * dy * dz

        f = 0.0 + 0.0j

        # Triple-integral over 3D space
        for x in x_vals:
            for y in y_vals:
                for z in z_vals:
                    r = (x, y, z)
                    V = self.V(x, y, z, rotation_elements=self.rotation_matrix)

                    # not sure which axes need to be dotted
                    f += V*np.exp(1j*np.tensordot(qvector, r,
                                                  axes=(0, 0)))*dVolume

        return self.pargs['delta_rho'] * f

    def form_factor_squared_numerical(
            self,
            qvector,
            num_points=100,
            size_scale=None):
        """Returns the square of the form factor."""
        f = self.form_factor_numerical(
            qvector, num_points=num_points, size_scale=size_scale)
        g = f * f.conjugate()
        return g.real

    def form_factor_squared(self, qvector):
        """Returns the square of the form factor.
            Value returned is real.
            Note : Need to implement form_factor.
        """

        f = self.form_factor(qvector)
        g = f * np.conj(f)
        return g.real

    def func_orientation_spread(
            self,
            q,
            func,
            dtype,
            num_phi=50,
            num_theta=50,
            orientation_spread=None):
        """ Compute an orientation average of function func from a sphere of
        points whose radius is q

            Parameters
            ----------
            q : q points, can be an array
            func :  should take a qvector as sole argument ([qx,qy,qz])
                This is a general function used by form_factors calcs
                (isotropic or orient spread)

            dtype :  the data type (complex or float)

            orientation_spread: If set, overrides the one in pargs. If not set,
                searches in pargs.  The default distribution is uniform.
                Endpoints are given by pargs['orientation_spread'] which are
                phi_start, phi_end, theta_start, theta_end

            Returns
            -------
                F : result of function func, averaged over specified
                orientations

            Notes
            -----
            This function is intended to be used to create a effective form
            factor when particle orientations have some distribution.
            Math:
                solid angle is dS = r^2 sin(th) dth dphi
                ignore the r part (or q)

        """
        if orientation_spread is None:
            if 'orientation_spread' not in self.pargs:
                raise ValueError(
                    "form_factor_orientation_spread : Sorry, "
                    "orientation_spread not "
                    "set in defining potential arguments"
                    " (pargs). Please define this parameter and run again.")
            else:
                orientation_spread = self.pargs['orientation_spread']

        phi_start, phi_end, theta_start, theta_end = orientation_spread

        # Phi is orientation around z-axis (in x-y plane)
        phi_vals, dphi = np.linspace(
            phi_start, phi_end, num_phi, endpoint=False, retstep=True)
        # Theta is tilt with respect to +z axis
        theta_vals, dtheta = np.linspace(
            theta_start, theta_end, num_theta, endpoint=False, retstep=True)

        F = np.zeros_like(q, dtype=dtype)
        dStot = 0.

        for theta in theta_vals:
            qz = q * np.cos(theta)
            dS = np.sin(theta) * dtheta * dphi
            dStot += dS * num_phi

            for phi in phi_vals:
                qx = -q * np.sin(theta) * np.cos(phi)
                qy = q * np.sin(theta) * np.sin(phi)
                qvector = np.array([qx, qy, qz])

                F += func(qvector) * dS

        # when the orientation spread is full solid angle, still I think it's
        # better to use dStot than 4*np.pi
        F /= dStot

        return F

    def form_factor_isotropic(self, q, num_phi=50, num_theta=50):
        """Returns the particle form factor, averaged over every possible orientation.
                Math:
                    solid angle is dS = r^2 sin(th) dth dphi
                    ignore the r part (or q)
        """
        # phi_start, phi_end, theta_start, theta_end
        orientation_spread = (0, 2 * np.pi, 0, np.pi)
        return self.func_orientation_spread(
            q,
            self.form_factor,
            complex,
            num_phi=num_phi,
            num_theta=num_theta,
            orientation_spread=orientation_spread)

    def form_factor_squared_isotropic(self, q, num_phi=50, num_theta=50):
        """Returns the square of the form factor, under the assumption of
        random orientation. In other words, we average over every possible
        orientation.
        This value is denoted by P(q)

        Note this is different from form_factor_isotropic, where the difference
        is:
            |<F>|^2 versis <|F|^2>. See "Periodic lattices of arbitrary
            nano-objects: modeling and applications for self-assembled systems"
            (DOI: 10.1107/S160057671302832X)
        """
        # phi_start, phi_end, theta_start, theta_end
        orientation_spread = (0, 2 * np.pi, 0, np.pi)
        return self.func_orientation_spread(
            q,
            self.form_factor_squared,
            float,
            num_phi=num_phi,
            num_theta=num_theta,
            orientation_spread=orientation_spread)

    def form_factor_orientation_spread(
            self,
            q,
            num_phi=50,
            num_theta=50,
            orientation_spread=None):
        """Returns the particle form factor, averaged over some orientations.
        This function is intended to be used to create a effective form factor
        when particle orientations have some distribution.

        orientation_spread: If set, overrides the one in pargs. If not set,
        searches in pargs.

        The default distribution is uniform. Endpoints are given by
        pargs['orientation_spread'] which are phi_start, phi_end, theta_start,
        theta_end
        """
        return self.func_orientation_spread(
            q,
            self.form_factor,
            complex,
            num_phi=num_phi,
            num_theta=num_theta,
            orientation_spread=orientation_spread)

    def form_factor_squared_orientation_spread(
            self, q, num_phi=50, num_theta=50, orientation_spread=None):
        """Returns the particle form factor, averaged over some orientations.
        This function is intended to be used to create a effective form factor
        when particle orientations have some distribution.

        orientation_spread: If set, overrides the one in pargs. If not set,
        searches in pargs.

        The default distribution is uniform. Endpoints are given by
        pargs['orientation_spread'] which are phi_start, phi_end, theta_start,
        theta_end
        """
        return self.func_orientation_spread(
            q,
            self.form_factor_squared,
            float,
            num_phi=num_phi,
            num_theta=num_theta,
            orientation_spread=orientation_spread)

    def beta_ratio(self, q, num_phi=50, num_theta=50, approx=False):
        """Returns the beta ratio: |<F(q)>|^2 / <|F(q)|^2>
        This ratio depends on polydispersity: for a monodisperse system, beta =
        1 for all q."""
        numerator = np.abs(
            self.form_factor_isotropic(
                q,
                num_phi=num_phi,
                num_theta=num_theta))**2
        denominator = self.form_factor_squared_isotropic(
            q, num_phi=num_phi, num_theta=num_theta)
        return numerator / denominator

    def P_beta(self, q, num_phi=50, num_theta=50, approx=False):
        """Returns P (isotropic_form_factor_squared) and beta_ratio.
        This function can be highly optimized in derived classes."""

        P = self.form_factor_squared_isotropic(
            q, num_phi=num_phi, num_theta=num_theta)
        beta = self.beta_ratio(
            q,
            num_phi=num_phi,
            num_theta=num_theta,
            approx=approx)

        return P, beta

    def to_string(self):
        """Returns a string describing the object."""
        s = "Base NanoObject (zero potential everywhere)."
        return s

    def to_short_string(self):
        """Returns a short string describing the object's variables.
        (Useful for distinguishing objects of the same class.)"""
        s = "(0)"
        return s

    def form_factor(self, qvector):
        """Returns the complex-amplitude of the form factor at the given
        q-coordinates.
            qvector is an array as such: [qx, qy, qz] where qx, qy, qz need
            matching dimensions
        """

        # example:
        # qvector = self.map_qcoord(qvector)
        raise NotImplementedError(
            "Needs to be implemented by inheriting object")

    def V(self, rvec, rotation_elements=None):
        """Returns the intensity of the real-space potential at the
        given real-space coordinates.
        rvec : [x,y,z]

        This method should be overwritten.
        """
        raise NotImplementedError(
            "This needs to be implemented by the inherting object")

    def projections(self, length, npoints=100):
        ''' Compute the xy, yz, and xz projections (in that order).

            This is a convenience routine to allow one to see approximately
            what the nano object looks like. Useful when creating new composite
            nano objects.

            Parameters
            ----------

            length : length of the box to compute the projections.
                Will compute a 3D box of [-length, +length] in x, y and z

            npoints : the number of points to calculate per dimension
                default 100
                WARNING : This creates a npoints x npoints x npoints array

            Returns
            -------

            V_xy : the xy projection
            V_xz : the xz projection
            V_yz : the xz projection

            Notes
            -----
            To compute the projection, this function must first compute a 3D
            density field of the object, and add various projections.  This
            code can become very slow and memory intensive if npoints is too
            large.

            It also is strongly reliant on how well Obj.V() was coded.  When
            adding new NanoObjects, be careful to program in the Obj.V()
            function properly.
            This is a great tool for CompositeNanoObjects
        '''
        x = np.linspace(-length, length, npoints)
        # ij indexing means that we index in V[x,y,z]
        # Note that rightermost index is fastest varying index
        x, y, z = np.meshgrid(x, x, x, indexing='ij')
        V = self.V(np.array([x, y, z]))
        V_xy = np.sum(V, axis=2).T
        V_xz = np.sum(V, axis=1).T
        V_yz = np.sum(V, axis=0).T

        return V_xy, V_xz, V_yz
