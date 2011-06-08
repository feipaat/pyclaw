r"""
Module containing all Pyclaw solution objects

:Authors:
    David I. Ketcheson -- Initial version
"""

import numpy as np

class State(object):
    r"""
    Contains the current state on a particular grid, including q, t, and aux.
    Usually aux is time-independent.

    Both q and aux are initialized to None.  They cannot be accessed until
    meqn and maux (respectively) are set.

    :State Data:
    
        The arrays :attr:`q` and :attr:`aux` have variable 
        extents based on the grid dimensions and the values of 
        :attr:`meqn` and :attr:`maux`.  Note that these are initialy set to 
        None and later set to appropriately sized empty numpy arrays when
        :attr:`meqn` and :attr:`maux` are set.

    Typical usage, assuming a 1D grid::

        >>> state = State(grid)
        >>> meqn = 2
        >>> state.q[0,:] = sin(grid.x.center)
        >>> state.q[1,:] = cos(grid.x.center)
        >>> maux = 1
        >>> aux[0,:] = exp(grid.x.center)
    """

    # ========== Property Definitions ========================================
    def meqn():
        doc = r"""(int) - Number of unknowns (components of q)"""
        def fget(self):
            if self.q is None:
                raise Exception('state.meqn has not been set.')
            else: return self.q.shape[0]
        def fset(self,meqn):
            if self.q is None:
                q_shape = [meqn]
                q_shape.extend(self.grid.n)
                self.q = np.empty(q_shape,order='F')
            else:
                raise Exception('You cannot change state.meqn after q is initialized.')
        return locals()
    meqn = property(**meqn())
    def maux():
        doc = r"""(int) - Number of auxiliary fields"""
        def fset(self,maux):
            if self.aux is not None:
                raise Exception('You cannot change state.maux after aux is initialized.')
            else:
                aux_shape = [maux]
                aux_shape.extend(self.grid.n)
                self.aux = np.empty(aux_shape,order='F')
        def fget(self):
            if self.aux is not None: return self.aux.shape[0]
            else: return 0
        return locals()
    maux = property(**maux())

    # ========== Class Methods ===============================================
    def __init__(self,grid):
        import pyclaw.grid
        if not isinstance(grid,pyclaw.grid.Grid):
            raise Exception("""A PyClaw State object must be initialized with
                             a PyClaw Grid object.""")

        # ========== Attribute Definitions ===================================
        self.grid = grid
        r"""pyclaw.Grid.grid - The grid this state lives on"""
        self.aux = None
        r"""(ndarray(...,maux)) - Auxiliary array for this grid containing per 
            cell information"""
        self.q   = None
        r"""(ndarray(...,meqn)) - Cell averaged quantity being evolved."""
        self.aux_global = {}
        r"""(dict) - Dictionary of global values for this grid, 
            ``default = {}``"""
        self.t=0.
        r"""(float) - Current time represented on this grid, 
            ``default = 0.0``"""
        self.stateno = 1
        r"""(int) - State number of current state, ``default = 1``"""

    def __str__(self):
        output = "State %s:\n" % self.stateno
        output += "  t=%s mbc=%s meqn=%s\n  " % (self.t,self.mbc,self.meqn)
        output += '\n'
        if self.q is not None:
            output += "  q.shape=%s" % str(self.q.shape)
        if self.aux is not None:
            output += " aux.shape=%s" % str(self.aux.shape)
        return output

    def is_valid(self):
        r"""
        Checks to see if this state is valid
        
        The state is declared valid based on the following criteria:
            - :attr:`q` is not None
            
        A debug logger message will be sent documenting exactly what was not 
        valid.
            
        :Output:
         - (bool) - True if valid, false otherwise.
        
        """
        import logging
        valid = True
        logger = logging.getLogger('solution')
        if self.q is None:
            logger.debug('The array q has not been initialized.')
            valid = False
        return valid
 
        
    def set_cparam(self,fortran_module):
        """
        Set the variables in fortran_module.cparam to the corresponding values in
        grid.aux_global.  This is the mechanism for passing scalar variables to the
        Fortran Riemann solvers; cparam must be defined as a common block in the
        Riemann solver.

        This function should be called from solver.setup().  This seems like a fragile
        interdependency between solver and grid; perhaps aux_global should belong
        to solver instead of grid.

        This function also checks that the set of variables defined in cparam 
        all appear in aux_global.

        """
        if hasattr(fortran_module,'cparam'):
            if not set(dir(fortran_module.cparam)) <= set(self.aux_global.keys()):
                raise Exception('Some required value(s) in the cparam common block in the Riemann solver have not been set in aux_global.')
            for global_var_name,global_var_value in self.aux_global.iteritems(): 
                setattr(fortran_module.cparam,global_var_name,global_var_value)

    # ========== Copy functionality ==========================================
    def __copy__(self):
        return self.__class__(self)
        
        
    def __deepcopy__(self,memo={}):
        import copy
        result = self.__class__(copy.deepcopy(self.grid))
        result.__init__(copy.deepcopy(self.grid))
        
        for attr in ('stateno','t','meqn'):
            setattr(result,attr,copy.deepcopy(getattr(self,attr)))
        
        if self.q is not None:
            result.q = copy.deepcopy(self.q)
        if self.aux is not None:
            result.aux = copy.deepcopy(self.aux)
        result.aux_global = copy.deepcopy(self.aux_global)
        
        return result
 