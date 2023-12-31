'''
Some convenient unit type classes for conversions
'''
import numpy as np

from astropy.units import Quantity, Unit
import astropy.units as u
import astropy.constants as const
from collections.abc import Iterable

class Flux(Quantity):

    def __init__(self, quantity):
        '''
        This is a Flux Unit Type and should have units of Energy / Distance^2 / Time
        '''
        super(Flux, self).__init__()

        self._name = 'energy flux'
        
    @staticmethod
    def isflux(quantity):
        '''
        Tests if this is a flux
        '''
        _name = 'energy flux'
        
        if not isinstance(quantity, list):
            quantity = [quantity]
            
        return all(_name in q.unit.physical_type for q in quantity)
        
    def tofluxdensity(self, wave_eff:Quantity=None, freq_eff:Quantity=None,
                      out_units:Unit='mag(AB)'):
        '''
        Converts to flux density by dividing by the given effective wave/freq

        Args:
            freq_eff [astropy.units.Quantity]: An astropy Quantity with the effective frequency
            wave_eff [astropy.units.Quantity]: An astropy Quantity with the effective wavelength
        
        Returns:
            A FluxDensity Quantity Object with the adjusted data
        '''

        # check that the input unit is a flux density
        if not FluxDensity.isfluxdensity(1*Unit(out_units)):
            raise ValueError('The output units that were given are not a flux density!!')

        # perform the conversion to out_units
        if freq_eff is not None and wave_eff is None:
            if not isinstance(freq_eff, np.ndarray):
                freq_eff = np.array(freq_eff)

            out = (self/freq_eff).to(out_units)

        elif freq_eff is None and wave_eff is not None:
            if not isinstance(wave_eff, np.ndarray):
                wave_eff = np.array(wave_eff)
    
            out = (self/wave_eff).to(out_units)
            
        else:
            raise ValueError('Either freq_eff or wave_eff must be provided (NOT both!)')

        # double check this is a fluxdensity now
        if not FluxDensity.isfluxdensity(out):
            raise ValueError('Something went wrong! This is not a flux density still!')

        return out
        
    def toflux(self, out_units=u.erg/u.cm**2/u.s, **kwargs):
        '''
        Just a wrapper that returns itself. This will make my other code cleaner
        '''
        return self.to(out_units)
        
class FluxDensity(Quantity):

    def __init__(self, quantity):
        '''
        This is a Flux Unit Type and should have units of a Energy / Distance^2 / Time / Frequency 
        
        This also handles AB/ST Magnitudes depending on if it is in frequency or 
        wavelength space
        '''
        super(FluxDensity, self).__init__()

        self._name_wave = 'spectral flux density wav'
        self._name_freq = 'spectral flux density'
        
    @staticmethod
    def isfluxdensity(quantity):
        '''
        Tests if this is a flux density
        '''
        _name_wave = 'spectral flux density wav'
        _name_freq = 'spectral flux density'
        
        if not isinstance(quantity, list):
            quantity = [quantity]
       
        test1 = all(_name_freq in q.unit.physical_type for q in quantity)
        test2 = all(_name_wave in q.unit.physical_type for q in quantity)
        return test1 or test2
        
    def toflux(self, freq_eff:Quantity=None, wave_eff:Quantity=None,
               out_units:Unit=u.erg/u.cm**2/u.s) -> Flux:
        '''
        Converts to flux by multiplying by the fiven effective wave/freq

        Args:
            freq_eff [astropy.units.Quantity]: An astropy Quantity with the effective frequency
            wave_eff [astropy.units.Quantity]: An astropy Quantity with the effective wavelength
            out_units [astropy.units.Unit]: An astropy Unit with units of flux

        Returns:
            A Flux Quantity Object with the adjusted data
        '''

        # check that the input unit is a flux
        if not Flux.isflux(1*Unit(out_units)):
            raise ValueError('The output units that were given are not a flux!!')
        
        # perform the conversion
        if freq_eff is not None and wave_eff is None:
            if not isinstance(freq_eff, np.ndarray):
                freq_eff = np.array(freq_eff)

            if 'wav' in self.unit.physical_type:
                eff = const.c / freq_eff # convert to wavelengths
            else:
                eff = freq_eff
            out = (self*eff).to(out_units)

        elif freq_eff is None and wave_eff is not None:
            if not isinstance(wave_eff, np.ndarray):
                wave_eff = np.array(wave_eff)

            if 'wav' in self.unit.physical_type:
                eff = const.c / wave_eff # convert to frequency
            else:
                eff = wave_eff
                
            out = (self*eff).to(out_units)
            
        else:
            raise ValueError('Either freq_eff or wave_eff must be provided (NOT both!)')

        # double check this is a flux now
        if not Flux.isflux(out):
            raise ValueError('Something went wrong! This is not a flux still!')

        return Flux(out)
        
    def tofluxdensity(self, out_units='mag(AB)', freq_eff:Quantity=None,
                      wave_eff:Quantity=None):
        '''
        Just a wrapper that returns itself. This will make my other code cleaner
        '''

        wav = self._name_wave
        nowav = self._name_freq
        
        # check if out_units is in frequency or wavelength space
        if wav in Unit(out_units).physical_type and nowav in self.unit.physical_type:
            # There is a discrepancy
            # we have to put self in wavelength space
            if wave_eff is not None:
                # this is ideal, nothing needs to be done
                pass
            elif freq_eff is not None:
                # we can work with this
                wave_eff = const.c / freq_eff
            else:
                raise ValueError('We need a wave_eff or freq_eff to make this conversion!')
            
            eq = u.spectral_density(wave_eff)

        elif nowav in Unit(out_units).physical_type and wav in self.unit.physical_type:
            if freq_eff is not None:
                pass
            elif wave_eff is not None: 
                freq_eff = const.c/wave_eff
            else:
                raise ValueError('We need a wave_eff to make this conversion!')
            eq = u.spectral_density(freq_eff)

        else:
            eq = None

        return self.to(out_units, equivalencies=eq)
        
def get_type(quantity):
    '''
    Finds the type of the input quantity and returns the quantity as that class
    '''

    if FluxDensity.isfluxdensity(quantity):
        # we have to be a little careful cause some of these are mags
        if str(quantity.unit) == 'mag(AB)':
            print('fouund magab')
            return FluxDensity(quantity.to(u.erg/u.cm**2/u.s/u.Hz))
        elif str(quantity.unit) == 'mag(ST)':
            return FluxDensity(quantity.to(u.erg/u.cm**2/u.s/u.nm))
        else:
            return FluxDensity(quantity)
    elif Flux.isflux(quantity):
        return Flux(quantity)
    else:
        raise ValueError('This quantity does not have valild units!')
