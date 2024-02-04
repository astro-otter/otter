'''
Some global helper functions for various processing tasks
'''
import astropy.units as u
from .constants import FILTER_MAP_WAVE

def filter_to_obstype(band_name):
    '''
    Converts a band name to either 'radio', 'uvoir', 'xray'
    '''

    try:
        wave_eff = FILTER_MAP_WAVE[band_name]*u.nm
    except KeyError as exc:
        raise Exception('No Effective Wavelength Known for this band, please add it to constants'
                        ) from exc

    if wave_eff > 1*u.mm:
        return 'radio'
    elif wave_eff <= 1*u.mm and wave_eff >= 10*u.nm:
        return 'uvoir'
    else:
        return 'xray'
    

def clean_schema(schema):
    '''
    Clean out Nones and empty lists from the given subschema
    '''
    for key, val in list(schema.items()):
        if val is None or (isinstance(val, (list, dict)) and len(val) == 0):
            del schema[key]
    return schema
