from .NanoObjects import NanoObject
import numpy as np
from copy import deepcopy
# This file is where more complex nano objects can be stored
# These interface the same way as NanoObjects but since they're a little more
# complex it makes sense to separate them from NanoObjects
# TODO : Better simpler API for making composite object
# like make_composite( objectslist, pargslist)


class CompositeNanoObject(NanoObject):
    ''' This is a nano object made up of a collection of nano objects.
        You specify them with a list of objects and pargs (dictionaries).
        The pargs can contain positions and rotations as well.
        Need to redefine all the form factors etc since now they're computed
        for each element.

        objlist : either list of object identifier classes or objects
        themselves
        parglist : parameters for the object classes. If set to None, then
        objlist is assume a list of objects, and not classes

        There is a new parameter in the pargs, it is 'sign'. If it's one, the
        sample adds. If it's negative, it subtracts. Useful for things like
        core shell. Playing with 'rho' is too risky so I reserved a parameter
        for it.

        This is new compared to Kevin's old code. This is redundant with his
        Lattice class but allows to make more complex objects without
        worrying about the lattice.


        Notes
        -----

        Be very careful about subtracting objects. Make sure to think about the
        scenario where there exists an ambient solution (rho_ambient).
    '''

    def __init__(self, objlist, parglist=None, pargs={}):
        '''
            Composite object

            objlist : either list of classes or instantiated objects
                    (pargslist must be None for the latter)
        '''
        super(CompositeNanoObject, self).__init__(pargs=pargs)

        # now define the objects in a list
        self.nano_objects = list()

        if parglist is None:
            for obj in objlist:
                # just in case list has repeating object
                # make a deepcopy
                self.nano_objects.append(deepcopy(obj))
        else:
            for nobj, pargobj in zip(objlist, parglist):

                self.nano_objects.append(nobj(pargs=pargobj))

        # set defaults
        for obj in self.nano_objects:
            if 'sign' not in obj.pargs:
                # defaults to additive
                obj.pargs['sign'] = 1.

    def form_factor(self, qvector):
        """Returns the complex-amplitude of the form factor at the given
            q-coordinates.
            qvector is an array as such: [qx, qy, qz] where qx, qy, qz need
            matching dimensions
        """
        # Phase must be retrieved *before* mapping in q
        phase = self.get_phase(qvector)

        # first rotate just as for V
        qvector = self.map_qcoord(qvector)

        F = np.zeros_like(qvector[0], dtype=complex)
        for nobj in self.nano_objects:
            F += nobj.pargs['sign'] * nobj.form_factor(qvector) * phase

        return F

    def volume(self):
        ''' Return the sum of the volume of all objects.
        '''
        volume = 0.
        for nobj in self.nano_objects:
            volume += nobj.pargs['sign'] * nobj.volume()
        return volume

    def V(self, rvec):
        """Returns the intensity of the real-space potential at the
        given real-space coordinates.
        rvec : [x,y,z]
        """
        # always transform first
        rvec = self.map_rcoord(rvec)

        Vtot = np.zeros_like(rvec[0])
        for nobj in self.nano_objects:
            Vtot += nobj.pargs['sign'] * nobj.V(rvec)

        return Vtot
