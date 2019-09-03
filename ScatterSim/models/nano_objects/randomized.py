import numpy as np

from ScatterSim.NanoObjects import NanoObject


class RandomizedNanoObject(NanoObject):
    """Defines a nano-object of which certain parameters are randomized

        This is a little more useful than a polydisperse nano object
            in that it allows you to randomize a few parameters at once.
        The more parameters you randomize, the higher you may want to set
        nsamples for accuracy.

        Parameters
        ----------

        pargs : the base potential arguments (parameters) of the system
            In pargs, you need to specify the parameters necessary to build the
            base NanoObject Class supplied. On top of this, you also need to
            specify one more parameter:
                'distribution_num_points' : this is the number of points to
                randomly sample. Basically, this class will spawn this number
                of nano-objects, with each nano-object having parameters that
                are randomly varied by argdict (described next).
                This defaults to 21. You should set this larger for better
                accuracy.

        argdict:
            This will be randomized parameters. These will override the pargs.
            It is a dictionary, with the key of the dictionary entries being
            the parameter names. The entries themselves are dictionaries
            containing more information. For example, for a polydisperse sphere
            whose radius follows the lognormal distribution, you would put:
                {
                 'radius' : {
                             'distribution_type' : 'lognormal',
                             'mean' : 1, # average val
                             'sigma' : .1,  #std deviation
                            },
                 # add more arguments here etc...
                }

        Notes
        -----
        Call 'build_objects' to build a list of objects with parameters
        randomly varied. The form factor etc will then be computed by averaging
        over these quanties. Note that form_factor_squared != form_factor^2 in
        this case! See Beta(q) in [1] for more information.


        References
        ---------
        [1] Yager, Kevin G., et al. "Periodic lattices of arbitrary
        nano-objects: modeling and applications for self-assembled systems."
        Journal of Applied Crystallography 47.1 (2014): 118-129.

    """

    def __init__(self, baseNanoObjectClass, pargs={}, argdict=None):

        # this should set the orientation and origin defaults
        NanoObject.__init__(self, pargs=pargs)

        if argdict is None:
            argdict = dict(radius={'distribution': 'gaussian',
                                   'sigma': .1, 'mean': 1})
        if 'distribution_num_points' not in pargs:
            pargs['distribution_num_points'] = 21

        self.argdict = argdict
        self.pargs = pargs

        # check that parameters in args dict are there
        for key, val in self.argdict.items():
            if key not in self.pargs:
                raise ValueError("Error {} not present in pargs".format(key))

            if 'distribution_type' not in val:
                val['distribution_type'] = 'gaussian'

        self.baseNanoObjectClass = baseNanoObjectClass
        self.object_list = []
        self.rebuild()

    def rebuild(self, pargs={}, argdict={}):
        """Allows the object to have its potential
        arguments (pargs) updated. Note that this doesn't
        replace the old pargs entirely. It only modifies
        (or adds) the key/values provided by the new pargs.

        For a polydisperse object, need to update pargs for elements in
        distribution, omitting variable that's being modified.
        """
        self.pargs.update(pargs)
        self.argdict.update(argdict)
        self.build_objects()

    def build_objects(self):
        ''' Build a list of objects whose parameters are randomly
                sampled according to a probability distribution.

            Compute the probability distribution.
            Returns the previously calculated distribution
                or recomputes if doesn't exist.

            Since we can have multipe parameters with different
                distributions, this is a distribution list.

            Run rebuild to reset.
            # TODO : Make distribution from numbers
                when summing run same object, rebuild and re-calculate
                since it's random, you don't want to cache either
                'distribution_num_points' : 21,
                {'radius' : {'distribution_type' : 'gaussian',
                    # parameters specific to distribution
                    'avg' : 1, # average val
                    'std' : .1,  #std deviation
                    'spread' : 2.5}}
        '''
        self.object_list = list()
        for i in range(self.pargs['distribution_num_points']):
            pargs_tmp = self.pargs.copy()
            for key, entry in self.argdict.items():
                pargs_tmp[key] = self.sample_point(entry)
            self.object_list.append(self.baseNanoObjectClass(pargs=pargs_tmp))

    def sample_point(self, dist_dict):
        ''' Sample a point from a distribution specified by
            dist_dict
            Currently supported:
                'distribution_type' (case insensitive)
                    'Gaussian' or 'normal'
                        parameters :
                            'mean' : average of Gaussian (normal) distribution
                            'sigma' : standard deviation of Gaussian (normal)
                            distribution
                    'uniform' :
                        parameters :
                            'low' : lower bound of distribution
                            'high' : upper bound of distribution
                    'lognormal'
                        parameters:
                            'mean' : mean value of underlying normal
                            distribution
                            'sigma' : standard deviation of underlying normal
                            distribution


        '''
        _supported_distributions = [
            'gaussian', 'normal', 'uniform', 'lognormal']
        distribution_type = dist_dict['distribution_type'].lower()
        if distribution_type == 'gaussian' or distribution_type == 'normal':
            mean = dist_dict['mean']
            sigma = dist_dict['sigma']
            return np.random.normal(loc=mean, scale=sigma)
        elif distribution_type == 'uniform':
            low = dist_dict['low']
            high = dist_dict['high']
            return np.random.uniform(low=low, high=high)
        elif distribution_type == 'lognormal':
            mean = dist_dict['mean']
            sigma = dist_dict['sigma']
            return np.random.lognormal(mean=mean, sigma=sigma)
        else:
            errorstr = "Error, distribution {} not supported".format(
                distribution_type)
            errorstr = errorstr + \
                "\nSupported are: {}".format(_supported_distributions)
            raise ValueError(errorstr)

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

        for curNanoObject in self.object_list:
            res_R = getattr(curNanoObject, funcname)(*args, **kwargs)
            res += res_R
            cts += 1

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