from __future__ import print_function, absolute_import, division
import _my_little_magfie
import f90wrap.runtime
import logging
import numpy

class My_Little_Magfie(f90wrap.runtime.FortranModule):
    """
    Module my_little_magfie
    
    
    Defined at my_little_magfie.f90 lines 1-23
    
    """
    pass
    @staticmethod
    def _eval_field_b(x, b):
        """
        _eval_field_b(x, b)
        
        
        Defined at my_little_magfie.f90 lines 8-13
        
        Parameters
        ----------
        x : float array
        b : float array
        
        """
        _my_little_magfie.f90wrap_eval_field_b(x=x, b=b)
    
    @staticmethod
    def _eval_field_b_and_a(x, y, z, bx, by, bz, ax, ay, az):
        """
        _eval_field_b_and_a(x, y, z, bx, by, bz, ax, ay, az)
        
        
        Defined at my_little_magfie.f90 lines 15-23
        
        Parameters
        ----------
        x : float
        y : float
        z : float
        bx : float
        by : float
        bz : float
        ax : float
        ay : float
        az : float
        
        """
        _my_little_magfie.f90wrap_eval_field_b_and_a(x=x, y=y, z=z, bx=bx, by=by, bz=bz, \
            ax=ax, ay=ay, az=az)
    
    @staticmethod
    def eval_field(*args, **kwargs):
        """
        eval_field(*args, **kwargs)
        
        
        Defined at my_little_magfie.f90 lines 3-5
        
        Overloaded interface containing the following procedures:
          _eval_field_b
          _eval_field_b_and_a
        
        """
        for proc in [My_Little_Magfie._eval_field_b, \
            My_Little_Magfie._eval_field_b_and_a]:
            try:
                return proc(*args, **kwargs)
            except TypeError:
                continue
        
    
    _dt_array_initialisers = []
    

my_little_magfie = My_Little_Magfie()

