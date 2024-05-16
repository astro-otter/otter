"""
Some constants, mappings, and functions to be used across the software
"""

from __future__ import annotations
import os
import ads
import astropy.units as u

"""
Helper functions first that just don't belong anywhere else
"""


def filter_to_obstype(band_name):
    """
    Converts a band name to either 'radio', 'uvoir', 'xray'
    """

    try:
        wave_eff = FILTER_MAP_WAVE[band_name] * u.nm
    except KeyError as exc:
        raise Exception(
            f"No Effective Wavelength Known for {band_name}, please add it to constants"
        ) from exc

    if wave_eff > 1 * u.mm:
        return "radio"
    elif wave_eff <= 1 * u.mm and wave_eff >= 10 * u.nm:
        return "uvoir"
    else:
        return "xray"


def clean_schema(schema):
    """
    Clean out Nones and empty lists from the given subschema
    """
    for key, val in list(schema.items()):
        if val is None or (isinstance(val, (list, dict)) and len(val) == 0):
            del schema[key]
    return schema


def bibcode_to_hrn(bibcode):
    """
    Converts a bibcode to a human_readable_name (hrn) using ADSQuery
    """

    try:
        adsquery = list(ads.SearchQuery(bibcode=bibcode))[0]
    except IndexError:
        raise ValueError(f"Could not find {bibcode} on ADS!")

    authors = adsquery.author
    year = adsquery.year

    if len(authors) == 0:
        raise ValueError("This ADS bibcode does not exist!")
    elif len(authors) == 1:
        author = authors[0]
    elif len(authors) == 2:
        author = authors[0] + " & " + authors[1]
    else:  # longer than 2
        author = authors[0] + " et al."

    # generate the human readable name
    hrn = author + " (" + year + ")"
    return hrn


"""
Then the constants and dictionary mappings used throughout
"""

# gives the effective wavelength for each filter given
# these are all in nanometers!
FILTER_MAP_WAVE = {
    "FUV": 153.8620701901866,
    "NUV": 231.56631043707714,
    "UVW2": 207.98754332676123,
    "uvw2": 207.98754332676123,
    "W2": 207.98754332676123,
    "2": 207.98754332676123,
    "uw2": 207.98754332676123,
    "UVM2": 225.47802478793594,
    "uvm2": 225.47802478793594,
    "M2": 225.47802478793594,
    "M": 225.47802478793594,
    "um2": 225.47802478793594,
    "UVW1": 261.3713060531025,
    "uvw1": 261.3713060531025,
    "W1": 261.3713060531025,
    "1": 261.3713060531025,
    "uw1": 261.3713060531025,
    "u": 356.17887353001856,
    "u'": 356.17887353001856,
    "up": 356.17887353001856,
    "uprime": 356.17887353001856,
    "U_S": 347.06360491031495,
    "s": 347.06360491031495,
    "us": 347.06360491031495,
    "U": 353.10502283105023,
    "B": 443.0539845758355,
    "B_S": 435.912081730874,
    "b": 435.912081730874,
    "bs": 435.912081730874,
    "g": 471.8872246248687,
    "g'": 471.8872246248687,
    "gp": 471.8872246248687,
    "gprime": 471.8872246248687,
    "F475W": 471.8872246248687,
    "g-DECam": 482.6787274749997,
    "c": 540.8724658332794,
    "cyan": 540.8724658332794,
    "V": 553.7155963302753,
    "V_S": 543.0131091205997,
    "v": 543.0131091205997,
    "vs": 543.0131091205997,
    "Itagaki": 651.0535687558726,
    "white": 752.0,
    "unfilt.": 616.690135,
    "0": 616.690135,
    "C": 616.690135,
    "clear": 616.690135,
    "pseudobolometric": 616.690135,
    "griz": 616.690135,
    "RGB": 616.690135,
    "LGRB": 616.690135,
    "G": 673.5412573108297,
    "Kepler": 641.6835660569259,
    "TESS": 797.2360657697333,
    "DLT40": 615.8130149792426,
    "Open": 615.8130149792426,
    "Clear": 615.8130149792426,
    "w": 638.9300625093831,
    "o": 686.6260690394873,
    "orange": 686.6260690394873,
    "r": 618.5194476741524,
    "r'": 618.5194476741524,
    "rp": 618.5194476741524,
    "rprime": 618.5194476741524,
    "F625W": 618.5194476741524,
    "r-DECam": 643.2062638192127,
    "R": 646.9439215118385,
    "Rc": 646.9439215118385,
    "R_s": 646.9439215118385,
    "i": 749.9704174464691,
    "i'": 749.9704174464691,
    "ip": 749.9704174464691,
    "iprime": 749.9704174464691,
    "F775W": 749.9704174464691,
    "i-DECam": 782.6680306208917,
    "I": 788.558706467662,
    "Ic": 788.558706467662,
    "z_s": 867.9495480864285,
    "zs": 867.9495480864285,
    "z": 896.1488333992431,
    "z'": 896.1488333992431,
    "zp": 896.1488333992431,
    "zprime": 896.1488333992431,
    "z-DECam": 917.8949537472383,
    "y": 963.3308299506817,
    "y-DECam": 989.965087304703,
    "J": 1255.0918319447906,
    "H": 1630.5155019191195,
    "K": 2157.3600605745955,
    "Ks": 2157.3600605745955,
    "F070W": 705.5727879998312,
    "F090W": 904.2281265089156,
    "F115W": 1157.001589027877,
    "F150W": 1503.9880463410511,
    "F200W": 1993.3922957570885,
    "F225W": 2372.81,
    "F277W": 2769.332372846113,
    "F300M": 2990.7606605760484,
    "F335M": 3363.887076210947,
    "F356W": 3576.787256375927,
    "F360M": 3626.0578695682693,
    "F444W": 4415.974447587756,
    "F560W": 5645.279496731566,
    "F770W": 7663.455798629626,
    "F1000W": 9968.161727011531,
    "F1130W": 11310.984595876938,
    "F1280W": 12831.396996921212,
    "F1500W": 15091.367399905488,
    "F1800W": 18006.083119653664,
    "F2100W": 20842.526633138932,
    "F2550W": 25408.228367890282,
}
"""
Mapping for the effective wavelength in nanometers for all filters used in the dataset.
"""


# gives the effective frequency for all filters
# These are all in THz
FILTER_MAP_FREQ = {
    "FUV": 1975.086895569116,
    "NUV": 1346.3548820463916,
    "UVW2": 1531.4976984760474,
    "uvw2": 1531.4976984760474,
    "W2": 1531.4976984760474,
    "2": 1531.4976984760474,
    "uw2": 1531.4976984760474,
    "UVM2": 1360.083095675749,
    "uvm2": 1360.083095675749,
    "M2": 1360.083095675749,
    "M": 1360.083095675749,
    "um2": 1360.083095675749,
    "UVW1": 1236.8527545450256,
    "uvw1": 1236.8527545450256,
    "W1": 1236.8527545450256,
    "1": 1236.8527545450256,
    "uw1": 1236.8527545450256,
    "u": 849.2871562331687,
    "u'": 849.2871562331687,
    "up": 849.2871562331687,
    "uprime": 849.2871562331687,
    "U_S": 875.611103788721,
    "s": 875.611103788721,
    "us": 875.611103788721,
    "U": 858.321721875779,
    "B": 688.8500955332158,
    "B_S": 696.7876979144597,
    "b": 696.7876979144597,
    "bs": 696.7876979144597,
    "g": 648.9823425403824,
    "g'": 648.9823425403824,
    "gp": 648.9823425403824,
    "gprime": 648.9823425403824,
    "F475W": 648.9823425403824,
    "g-DECam": 635.8015668464043,
    "c": 580.1132515050684,
    "cyan": 580.1132515050684,
    "V": 548.3068934496129,
    "V_S": 554.9815375506427,
    "v": 554.9815375506427,
    "vs": 554.9815375506427,
    "Itagaki": 577.0861573682259,
    "white": 30079.243284322874,
    "unfilt.": 601.5655810567023,
    "0": 601.5655810567023,
    "C": 601.5655810567023,
    "clear": 601.5655810567023,
    "pseudobolometric": 601.5655810567023,
    "griz": 601.5655810567023,
    "RGB": 601.5655810567023,
    "LGRB": 601.5655810567023,
    "G": 518.6766845466752,
    "Kepler": 519.5058608954615,
    "TESS": 403.1881955125893,
    "DLT40": 629.637672549936,
    "Open": 629.637672549936,
    "Clear": 629.637672549936,
    "w": 520.8387777057242,
    "o": 451.71177203298663,
    "orange": 451.71177203298663,
    "r": 489.2629992899134,
    "r'": 489.2629992899134,
    "rp": 489.2629992899134,
    "rprime": 489.2629992899134,
    "F625W": 489.2629992899134,
    "r-DECam": 472.4459671948087,
    "R": 471.26223689126897,
    "Rc": 471.26223689126897,
    "R_s": 471.26223689126897,
    "i": 402.8409598867557,
    "i'": 402.8409598867557,
    "ip": 402.8409598867557,
    "iprime": 402.8409598867557,
    "F775W": 402.8409598867557,
    "i-DECam": 386.62233825433924,
    "I": 382.7915178046724,
    "Ic": 382.7915178046724,
    "z_s": 346.66628641927826,
    "zs": 346.66628641927826,
    "z": 337.7343708777923,
    "z'": 337.7343708777923,
    "zp": 337.7343708777923,
    "zprime": 337.7343708777923,
    "z-DECam": 328.753462451287,
    "y": 312.24818210606065,
    "y-DECam": 303.4727730182509,
    "J": 239.862442505934,
    "H": 185.33613196897403,
    "K": 139.79431978859097,
    "Ks": 139.79431978859097,
    "F070W": 431.04176743403116,
    "F090W": 336.17431986268366,
    "F115W": 262.87628654288676,
    "F150W": 201.94374815011136,
    "F200W": 152.56522352568953,
    "F277W": 110.05136786468209,
    "F300M": 100.56915203596012,
    "F335M": 89.41072625742719,
    "F356W": 85.01984846997881,
    "F360M": 82.9357933095218,
    "F444W": 68.96667222373961,
    "F560W": 53.67852315133938,
    "F770W": 39.87175477126777,
    "F1000W": 30.349460503852665,
    "F1130W": 26.53952983680919,
    "F1280W": 23.59741975845449,
    "F1500W": 20.08679352819493,
    "F1800W": 16.773842151606242,
    "F2100W": 14.581938602646188,
    "F2550W": 11.919267708332558,
}
"""
Mapping for the effective frequencies in THz for all the filters used in OTTER
"""


# x-ray telescope areas for converting
# NOTE: these are estimates from the links provided
# Since this is inherently instrument dependent they are not entirely reliable
# All are for 1-2 keV
XRAY_AREAS = {
    # https://swift.gsfc.nasa.gov/about_swift/Sci_Fact_Sheet.pdf
    "swift": 135 * u.cm**2,
    # https://heasarc.gsfc.nasa.gov/docs/rosat/ruh/handbook/node39.html#SECTION00634000000000000000
    "rosat": 400 * u.cm**2,
    # https://www.cosmos.esa.int/web/xmm-newton/technical-details-mirrors
    "xmm": 1500 * u.cm**2,
    "xmm slew": 1500 * u.cm**2,
    "xmm pointed": 1500 * u.cm**2,
    # https://cxc.harvard.edu/cdo/about_chandra
    "chandra": 600 * u.cm**2,
}
"""
X-Ray telescope areas that are used for converting from counts to other units.

NOTE: These are estimates from the following links
* https://swift.gsfc.nasa.gov/about_swift/Sci_Fact_Sheet.pdf
* https://heasarc.gsfc.nasa.gov/docs/rosat/ruh/handbook/node39.html#SECTION00634000000000000000
* https://www.cosmos.esa.int/web/xmm-newton/technical-details-mirrors
* https://cxc.harvard.edu/cdo/about_chandra
"""


# define a working base directory constant
BASEDIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
"""
Base directory for the OTTER API software package
"""

DATADIR = os.path.join(BASEDIR, "data", "base")
"""
Deprecated database directory that IS NOT always constant anymore
"""

# Overarching schema that stops once we get down to a string or list
schema = {
    "schema_version": {"value": "0", "comment": "Copied from tde.space"},
    "name": {"default_name": None, "alias": []},
    "coordinate": [],
    "distance": [],
    "classification": [],
    "reference_alias": [],
    "date_reference": [],
    "photometry": [],
    "spectra": [],
    "filter_alias": [],
}
"""
Schema dictionary to be filled with values from the subschemas
"""

# sub schemas that get filled into lists
name_alias_schema = {"value": None, "reference": None}
"""
Subschema for the name and alias dictionary
"""

coordinate_schema = {
    "ra": None,
    "dec": None,
    "l": None,
    "b": None,
    "lon": None,
    "lat": None,
    "ra_units": None,
    "dec_units": None,
    "l_units": None,
    "b_units": None,
    "lon_units": None,
    "lat_units": None,
    "ra_error": None,
    "dec_error": None,
    "l_error": None,
    "b_error": None,
    "lon_error": None,
    "lat_error": None,
    "epoch": None,
    "frame": None,
    "coord_type": None,
    "computed": None,
    "reference": None,
    "default": None,
}
"""
Subschema to describe the possible attributes for the coordinate dictionary
"""


distance_schema = {
    "value": None,
    "unit": None,
    "error": None,
    "cosmology": None,
    "reference": None,
    "computed": None,
    "uuid": None,
    "default": None,
    "distance_type": None,
}
"""
Subschema to describe the possible attributes for the distance dictionary
"""


classification_schema = {
    "object_class": None,
    "confidence": None,
    "class_type": None,
    "reference": None,
    "default": None,
}
"""
Subschema to describe the attributes for the classification dictionary
"""

reference_alias_schema = {"name": None, "human_readable_name": None}
"""
Subschema to describe the attributes for the reference alias dictionary
"""

date_reference_schema = {
    "value": None,
    "date_format": None,
    "date_type": None,
    "reference": None,
    "computed": None,
}
"""
Subschema to describe the date_reference dictionary attributes
"""

photometry_schema = {
    "raw": None,
    "raw_err": None,
    "raw_units": None,
    "value": None,
    "value_err": None,
    "value_units": None,
    "epoch_zeropoint": None,
    "epoch_redshift": None,
    "filter": None,
    "filter_key": None,
    "obs_type": None,
    "telescope_area": None,
    "date": None,
    "date_format": None,
    "date_err": None,
    "ignore": None,
    "upperlimit": None,
    "sigma": None,
    "sky": None,
    "telescope": None,
    "instrument": None,
    "phot_type": None,
    "exptime": None,
    "aperature": None,
    "observer": None,
    "reducer": None,
    "pipeline": None,
    "corr_k": None,
    "corr_av": None,
    "corr_host": None,
    "corr_hostav": None,
    "val_k": None,
    "val_s": None,
    "val_av": None,
    "val_host": None,
    "val_hostav": None,
}
"""
Subschema to describe all of the possible attributes that can be used in the photometry
dictionary
"""

spectra_schema = {
    "wavelength": None,
    "wavelength_units": None,
    "flux": None,
    "fluxerr": None,
    "raw": None,
    "raw_err": None,
    "sky": None,
    "lamp": None,
    "flux_units": None,
    "telescope": None,
    "instrument": None,
    "date": None,
    "date_format": None,
    "date_err": None,
    "exptime": None,
    "slit": None,
    "airmass": None,
    "disperser": None,
    "resolution": None,
    "resolution_units": None,
    "min_wave": None,
    "max_wave": None,
    "filter": None,
    "filter_key": None,
    "standard_name": None,
    "ignore": None,
    "spec_type": None,
    "aperture": None,
    "observer": None,
    "reducer": None,
    "pipeline": None,
    "corr_k": None,
    "corr_av": None,
    "corr_host": None,
    "corr_hostav": None,
    "corr_flux": None,
    "corr_phot": None,
    "val_k": None,
    "val_av": None,
    "val_host": None,
    "val_hostav": None,
}

filter_alias_schema = {
    "filter_key": None,
    "wave_eff": None,
    "wave_min": None,
    "wave_max": None,
    "freq_eff": None,
    "freq_min": None,
    "freq_max": None,
    "zp": None,
    "wave_units": None,
    "freq_units": None,
    "zp_units": None,
    "zp_system": None,
}
"""
Subschema to describe the attributes in the filter_alias dictionary
"""

# package the subschemas by the key used for that location in the Transient object
subschema = {
    "name/alias": name_alias_schema,
    "coordinate": coordinate_schema,
    "distance": distance_schema,
    "classification": classification_schema,
    "reference_alias": reference_alias_schema,
    "date_reference": date_reference_schema,
    "photometry": photometry_schema,
    "spectra": spectra_schema,
    "filter_alias": filter_alias_schema,
}
"""
A useful variable to describe all of the subschemas that are available and can be used
"""
