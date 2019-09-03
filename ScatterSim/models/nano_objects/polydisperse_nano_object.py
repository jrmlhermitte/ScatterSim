import numpy as np

from ScatterSim.NanoObjects import NanoObject


class PolydisperseNanoObject(NanoObject):
    """Defines a polydisperse nano-object, which has a distribution in
        argname of width argstdname
        if not specified, argname defaults to the parameter 'radius'
            and argstdname defaults to 'sigma_R'
        the width is in absolute units (not a percentage or ratio of the
            argument it's varying)

        Note : this is slow, but more general. If slow, it may be worth hard
        coding the form factor (if  not too complex).
    """

    def __init__(
            self,
            baseNanoObjectClass,
            pargs={},
            argname=None,
            argstdname=None):

        # this should set the orientation and origin defaults
        NanoObject.__init__(self, pargs=pargs)

        if argname is None:
            argname = 'radius'
        if argstdname is None:
            argstdname = 'sigma_R'

        self.argname = argname
        self.argstdname = argstdname

        # More defaults specific to the polydisperse object
        if argstdname not in self.pargs:
            raise ValueError(
                "Error : did not specify a {} in pargs.".format(argstdname) +
                " Please specify, or else do not use the "
                "polydisperse modifier.")
        if argname not in self.pargs:
            self.pargs[argname] = 1.

        if 'distribution_type' not in self.pargs:
            self.pargs['distribution_type'] = 'gaussian'
        if 'distribution_num_points' not in self.pargs:
            self.pargs['distribution_num_points'] = 5

        self.baseNanoObjectClass = baseNanoObjectClass

        self.distribution_list = []

    def rebuild(self, pargs={}):
        """Allows the object to have its potential
        arguments (pargs) updated. Note that this doesn't
        replace the old pargs entirely. It only modifies
        (or adds) the key/values provided by the new pargs.

        For a polydisperse object, need to update pargs for elements in
        distribution, omitting variable that's being modified.
        """
        self.pargs.update(pargs)
        self.distribution_list = self.distribution(force=True)

    def distribution(self, spread=2.5, force=False):
        ''' Compute the probability distribution.
            Returns the previously calculated distribution
                or recomputes if doesn't exist.

            Run rebuild to reset.
        '''
        if self.distribution_list == [] or force:
            argname = self.argname
            argstdname = self.argstdname
            # Build the distribution
            mean_val = self.pargs[argname]
            rms_val = self.pargs[argstdname]
            n = self.pargs['distribution_num_points']
            if self.pargs['distribution_type'] == 'gaussian':
                self.distribution_list = \
                    self.distribution_gaussian(mean=mean_val, rms=rms_val,
                                               num_points=n, spread=spread)
            elif self.pargs['distribution_type'] == 'lognormal':
                self.distribution_list = \
                    self.distribution_lognormal(mean=mean_val, rms=rms_val,
                                                num_points=n, spread=spread)
            else:
                print("Unknown distribution type in distribution().")

        # Return the existing distribution
        return self.distribution_list

    def distribution_gaussian(
            self,
            mean=1.0,
            rms=0.01,
            num_points=11,
            spread=2.5):
        '''
            Gaussian distribution of parameters.

            mean : the mean value

            rms : the rms value

        '''

        distribution_list = []

        step = 2 * spread * rms / (num_points - 1)
        # the sampled value
        sample = mean - step * (num_points - 1) / 2.0

        prefactor = 1 / (rms * np.sqrt(2 * np.pi))

        for i in range(num_points):
            delta = mean - sample
            wt = prefactor * np.exp(- (delta**2) / (2 * (rms**2)))

            curNanoObject = self.baseNanoObjectClass(pargs=self.pargs)
            curNanoObject.rebuild(pargs={self.argname: sample})

            distribution_list.append([sample, step, wt, curNanoObject])

            sample += step

        return distribution_list

    def distribution_lognormal(
            self,
            mean=1.0,
            rms=0.01,
            num_points=91,
            spread=10):
        '''
            Lognormal distribution of parameters.

            mean : this here will mean the scale of the lognorm (the higher the
            further the peak)

            rms : this will mean the shape of the lognorm distribution

        '''
        from scipy.stats import lognorm

        distribution_list = []
        actual_mean = lognorm.mean(rms, loc=0, scale=mean)
        actual_std = lognorm.std(rms, loc=0, scale=mean)

        # use the actual mean and std from lognormal distribution, dont
        # use the parameters, which are for the underlying Gaussian
        # distribution
        step = 2 * spread * actual_std / (num_points - 1)
        # the sampled value
        sample = actual_mean - step * (num_points - 1) / 2.0

        for i in range(num_points):
            wt = lognorm.pdf(sample, rms, loc=0, scale=mean)

            curNanoObject = self.baseNanoObjectClass(pargs=self.pargs)
            curNanoObject.rebuild(pargs={self.argname: sample})

            distribution_list.append([sample, step, wt, curNanoObject])

            sample += step

        return distribution_list

    def V(self, rvec):
        """Returns the average potential"""
        return self.dist_sum('V', rvec[0].shape, float, rvec)

    def volume(self):
        ''' ret avg volume'''
        return self.dist_sum('volume', 1, float)[0]

    def dist_sum(self, funcname, shape, dtype, *args, **kwargs):
        ''' Sum the function with name 'funcname' over variable 'vec' over the current
        distribution.
        Forwards other keyword arguments to function

        Parameters
        ---------
        funcname : the function name
        shape : the shape of the result
        dtype : the data type
        args : arguments to the function
        components : specifies if vec is of form [qx,qy,qz] (True)
            or just q (False)
        kwargs : keyword arguments to function

        '''
        res = np.zeros(shape, dtype=dtype)
        cts = 0.

        for R, dR, wt, curNanoObject in self.distribution():
            res_R = getattr(curNanoObject, funcname)(*args, **kwargs)
            res += wt * res_R * dR
            cts += wt * dR

        if cts == 0.:
            raise ValueError(
                "Nothing was added to distribution? \n"
                "Distribution list is: {}".format(
                    self.distribution()))

        return res / cts

    def form_factor(self, qvec):
        """Returns the complex-amplitude of the form factor at the given
            <F>_d
        q-coordinates."""
        return self.dist_sum('form_factor', qvec[0].shape, complex, qvec)

    def form_factor_distavg_squared(self, qvec):
        '''
            |<F>_d|^2
        '''
        return np.abs(
            self.dist_sum(
                'form_factor',
                qvec[0].shape,
                complex,
                qvec))**2

    def form_factor_squared(self, qvec):
        """Returns the square of the form factor.

            <|F|^2>_d
        """
        return self.dist_sum('form_factor_squared', qvec[0].shape, float, qvec)

    def form_factor_isotropic(self, q, num_phi=50, num_theta=50):
        ''' Returns the isotropic form factor
            < <F>_iso >_d

        '''
        return self.dist_sum(
            'form_factor_isotropic',
            q.shape,
            complex,
            q,
            num_phi=num_phi,
            num_theta=num_theta)

    def form_factor_squared_isotropic(self, q, num_phi=50, num_theta=50):
        ''' Returns the isotropic form factor
            < <|F|^2>_iso >_d
        '''
        return self.dist_sum(
            'form_factor_squared_isotropic',
            q.shape,
            float,
            q,
            num_phi=num_phi,
            num_theta=num_theta)

    def beta_numerator(self, q, num_phi=50, num_theta=50):
        """Returns the numerator of the beta ratio: |<<F(q)>_d>_iso|^2"""
        return np.abs(
            self.form_factor_isotropic(
                q,
                num_phi=num_phi,
                num_theta=num_theta))**2

    def beta_numerator_iso_external(self, q, num_phi=50, num_theta=50):
        """Calculates the beta numerator under the assumption that the
        orientational averaging is done last. That is, instead of calculating
        |<<F>>_iso|^2, we calculate <|<F>|^2>_iso
        """

        return self.func_orientation_spread(
            q,
            self.form_factor_distavg_squared,
            num_phi=num_phi,
            num_theta=num_theta)

    def beta_ratio(self, q, num_phi=50, num_theta=50, approx=False):
        """Returns the beta ratio: |<<F(q)>_iso>_d|^2 / <<|F(q)|^2>_iso>_d This
        ratio depends on polydispersity: for a monodisperse system, beta = 1
        for all q. """

        if approx:
            radius = self.pargs['radius']
            sigma_R = self.pargs['sigma_R']
            beta = np.exp(-((radius * sigma_R * q)**2))
            return beta
        else:
            # numerator and denominator
            beta_num = self.beta_numerator(
                q, num_phi=num_phi, num_theta=num_theta)
            beta_den = self.form_factor_squared_isotropic(
                q, num_phi=num_phi, num_theta=num_theta)
            return beta_num / beta_den