"""
Class for a transient,
basically just inherits the dict properties with some overwriting
"""

from __future__ import annotations
import warnings
from copy import deepcopy
import re
from collections.abc import MutableMapping
from typing import Callable
from typing_extensions import Self
import logging

import numpy as np
import pandas as pd

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord
from dust_extinction.parameter_averages import BaseExtRvModel, G23
from dustmaps.sfd import SFDQuery

from ..exceptions import (
    FailedQueryError,
    IOError,
    OtterLimitationError,
    TransientMergeError,
)
from ..util import XRAY_AREAS, _KNOWN_CLASS_ROOTS, _DuplicateFilter
from .host import Host

np.seterr(divide="ignore")
logger = logging.getLogger(__name__)


class Transient(MutableMapping):
    def __init__(self, d={}, name=None):
        """
        Overwrite the dictionary init

        Args:
            d (dict): A transient dictionary
            name (str): The default name of the transient, default is None and it will
                        be inferred from the input dictionary.
        """
        self.data = d

        if "reference_alias" in self:
            self.srcmap = {
                ref["name"]: ref["human_readable_name"]
                for ref in self["reference_alias"]
            }
            self.srcmap["TNS"] = "TNS"
        else:
            self.srcmap = {}

        if "name" in self:
            if "default_name" in self["name"]:
                self.default_name = self["name"]["default_name"]
            else:
                raise AttributeError("Missing the default name!!")
        elif name is not None:
            self.default_name = name
        else:
            self.default_name = "Missing Default Name"

        # Make it so all coordinates are astropy skycoords

    def __getitem__(self, keys):
        """
        Override getitem to recursively access Transient elements
        """

        if isinstance(keys, (list, tuple)):
            return Transient({key: self[key] for key in keys if key in self})
        elif isinstance(keys, str) and "/" in keys:  # this is for a path
            s = "']['".join(keys.split("/"))
            s = "['" + s
            s += "']"
            return eval(f"self{s}")
        elif (
            isinstance(keys, int)
            or keys.isdigit()
            or (keys[0] == "-" and keys[1:].isdigit())
        ):
            # this is for indexing a sublist
            return self[int(keys)]
        else:
            return self.data[keys]

    def __setitem__(self, key, value):
        """
        Override set item to work with the '/' syntax
        """

        if isinstance(key, str) and "/" in key:  # this is for a path
            s = "']['".join(key.split("/"))
            s = "['" + s
            s += "']"
            exec(f"self{s} = value")
        else:
            self.data[key] = value

    def __delitem__(self, keys):
        if "/" in keys:
            raise OtterLimitationError(
                "For security, we can not delete with the / syntax!"
            )
        else:
            del self.data[keys]

    def __iter__(self):
        return iter(self.data)

    def __len__(self):
        return len(self.data)

    def __repr__(self):
        return f"Transient(\n\tName: {self.default_name},\n\tKeys: {self.keys()}\n)"

    def keys(self):
        return self.data.keys()

    def __add__(self, other, strict_merge=True):
        """
        Merge this transient object with another transient object

        Args:
            other [Transient]: A Transient object to merge with
            strict_merge [bool]: If True it won't let you merge objects that
                                 intuitively shouldn't be merged (ie. different
                                 transient events).
        """

        # first check that this object is within a good distance of the other object
        if (
            strict_merge
            and self.get_skycoord().separation(other.get_skycoord()) > 10 * u.arcsec
        ):
            raise TransientMergeError(
                "These two transients are not within 10 arcseconds!"
                + " They probably do not belong together! If they do"
                + " You can set strict_merge=False to override the check"
            )

        # create set of the allowed keywords
        allowed_keywords = {
            "name",
            "date_reference",
            "coordinate",
            "distance",
            "filter_alias",
            "schema_version",
            "photometry",
            "classification",
            "host",
        }

        merge_subkeys_map = {
            "name": None,
            "date_reference": ["value", "date_format", "date_type"],
            "coordinate": None,  # may need to update this if we run into problems
            "distance": ["value", "distance_type", "unit"],
            "filter_alias": None,
            "schema_version": None,
            "photometry": None,
            "classification": None,
            "host": [
                "host_ra",
                "host_dec",
                "host_ra_units",
                "host_dec_units",
                "host_name",
            ],
        }

        groupby_key_for_default_map = {
            "name": None,
            "date_reference": "date_type",
            "coordinate": "coordinate_type",
            "distance": "distance_type",
            "filter_alias": None,
            "schema_version": None,
            "photometry": None,
            "classification": None,
            "host": None,
        }

        # create a blank dictionary since we don't want to overwrite this object
        out = {}

        # find the keys that are
        merge_keys = list(
            self.keys() & other.keys()
        )  # in both t1 and t2 so we need to merge these keys
        only_in_t1 = list(self.keys() - other.keys())  # only in t1
        only_in_t2 = list(other.keys() - self.keys())  # only in t2

        # now let's handle the merge keys
        for key in merge_keys:
            # reference_alias is special
            # we ALWAYS should combine these two
            if key == "reference_alias":
                out[key] = self[key]
                if self[key] != other[key]:
                    # only add t2 values if they aren't already in it
                    bibcodes = {ref["name"] for ref in self[key]}
                    for val in other[key]:
                        if val["name"] not in bibcodes:
                            out[key].append(val)
                continue

            # we can skip this merge process and just add the values from t1
            # if they are equal. We should still add the new reference though!
            if self[key] == other[key]:
                # set the value
                # we don't need to worry about references because this will
                # only be true if the reference is also equal!
                out[key] = deepcopy(self[key])
                continue

            # There are some special keys that we are expecting
            if key in allowed_keywords:
                Transient._merge_arbitrary(
                    key,
                    self,
                    other,
                    out,
                    merge_subkeys=merge_subkeys_map[key],
                    groupby_key=groupby_key_for_default_map[key],
                )
            else:
                # this is an unexpected key!
                if strict_merge:
                    # since this is a strict merge we don't want unexpected data!
                    raise TransientMergeError(f"{key} was not expected! Can not merge!")
                else:
                    # Throw a warning and only keep the old stuff
                    logger.warning(
                        f"{key} was not expected! Only keeping the old information!"
                    )
                    out[key] = deepcopy(self[key])

        # and now combining out with the stuff only in t1 and t2
        out = out | dict(self[only_in_t1]) | dict(other[only_in_t2])

        # now return out as a Transient Object
        return Transient(out)

    def get_meta(self, keys=None) -> Self:
        """
        Get the metadata (no photometry or spectra)

        This essentially just wraps on __getitem__ but with some checks

        Args:
            keys (list[str]) : list of keys to get the metadata for from the transient

        Returns:
            A Transient object of just the meta data
        """
        if keys is None:
            keys = list(self.keys())

            # note: using the remove method is safe here because dict keys are unique
            if "photometry" in keys:
                keys.remove("photometry")
            if "spectra" in keys:
                keys.remove("spectra")
        else:
            # run some checks
            if "photometry" in keys:
                logger.warning("Not returing the photometry!")
                _ = keys.pop("photometry")
            if "spectra" in keys:
                logger.warning("Not returning the spectra!")
                _ = keys.pop("spectra")

            curr_keys = self.keys()
            for key in keys:
                if key not in curr_keys:
                    keys.remove(key)
                    logger.warning(
                        f"Not returning {key} because it is not in this transient!"
                    )

        return self[keys]

    def get_skycoord(self, coord_format="icrs") -> SkyCoord:
        """
        Convert the coordinates to an astropy SkyCoord

        Args:
            coord_format (str): Astropy coordinate format to convert the SkyCoord to
                                defaults to icrs.

        Returns:
            Astropy.coordinates.SkyCoord of the default coordinate for the transient
        """

        # now we can generate the SkyCoord
        f = "df['coordinate_type'] == 'equatorial'"
        coord_dict = self._get_default("coordinate", filt=f)
        coordin = self._reformat_coordinate(coord_dict)
        coord = SkyCoord(**coordin).transform_to(coord_format)

        return coord

    def get_discovery_date(self) -> Time:
        """
        Get the default discovery date for this Transient

        Returns:
            astropy.time.Time of the default discovery date
        """
        key = "date_reference"
        try:
            date = self._get_default(key, filt='df["date_type"] == "discovery"')
        except KeyError:
            return None

        if date is None:
            return date

        if "date_format" in date:
            f = date["date_format"]
        else:
            f = "mjd"

        return Time(str(date["value"]).strip(), format=f)

    def get_redshift(self) -> float:
        """
        Get the default redshift of this Transient

        Returns:
            Float value of the default redshift
        """
        f = "df['distance_type']=='redshift'"
        default = self._get_default("distance", filt=f)
        if default is None:
            return default
        else:
            return float(default["value"])  # cast the redshift to a float

    def get_classification(self) -> tuple(str, float, list):
        """
        Get the default classification of this Transient.
        This normally corresponds to the highest confidence classification that we have
        stored for the transient.

        Returns:
            The default object class as a string, the confidence level in that class,
            and a list of the bibcodes corresponding to that classification. Or, None
            if there is no classification.
        """
        default = self._get_default("classification/value")
        if default is None:
            return default
        return default.object_class, default.confidence, default.reference

    def get_host(self, max_hosts=3, search=False, **kwargs) -> list[Host]:
        """
        Gets the default host information of this Transient. This returns an otter.Host
        object. If search=True, it will also check the BLAST host association database
        for the best match and return it as well. Note that if search is True then
        this has the potential to return max_hosts + 1, if BLAST also returns a result.
        The BLAST result will always be the last value in the returned list.

        Args:
            max_hosts [int] : The maximum number of hosts to return, default is 3
            **kwargs : keyword arguments to be passed to getGHOST

        Returns:
            A list of otter.Host objects. This is useful becuase the Host objects have
            useful methods for querying public catalogs for data of the host.
        """
        # first try to get the host information from our local database
        host = []
        if "host" in self:
            max_hosts = min([max_hosts, len(self["host"])])
            for h in self["host"][:max_hosts]:
                # only return hosts with their ra and dec stored
                if (
                    "host_ra" not in h
                    or "host_dec" not in h
                    or "host_ra_units" not in h
                    or "host_dec_units" not in h
                ):
                    continue

                # now we can construct a host object from this
                host.append(Host(transient_name=self.default_name, **dict(h)))

        # then try BLAST
        if search:
            logger.warning(
                "Trying to find a host with BLAST/astro-ghost. Note\
                 that this won't work for older targets! See https://blast.scimma.org"
            )

            # default_name should always be the TNS name if we have one
            print(self.default_name)
            blast_host = Host.query_blast(self.default_name.replace(" ", ""))
            print(blast_host)
            if blast_host is not None:
                host.append(blast_host)

        return host

    def _get_default(self, key, filt=None):
        """
        Get the default of key

        Args:
            key [str]: key in self to look for the default of
            filt [str]: a valid pandas dataframe filter to index a pandas dataframe
                        called df.
        """
        if key not in self:
            raise KeyError(f"This transient does not have {key} associated with it!")

        df = pd.DataFrame(self[key])
        if len(df) == 0:
            raise KeyError(f"This transient does not have {key} associated with it!")

        if filt is not None:
            df = df[eval(filt)]  # apply the filters

        if "default" in df:
            # first try to get the default
            df_filtered = df[df.default == True]
            if len(df_filtered) == 0:
                df_filtered = df
        else:
            df_filtered = df

        if len(df_filtered) == 0:
            return None

        return df_filtered.iloc[0]

    def _reformat_coordinate(self, item):
        """
        Reformat the coordinate information in item
        """
        coordin = None
        if "ra" in item and "dec" in item:
            # this is an equatorial coordinate
            coordin = {
                "ra": item["ra"],
                "dec": item["dec"],
                "unit": (item["ra_units"], item["dec_units"]),
            }
        elif "l" in item and "b" in item:
            coordin = {
                "l": item["l"],
                "b": item["b"],
                "unit": (item["l_units"], item["b_units"]),
                "frame": "galactic",
            }

        return coordin

    def clean_photometry(
        self,
        flux_unit: u.Unit = "mag(AB)",
        date_unit: u.Unit = "MJD",
        freq_unit: u.Unit = "GHz",
        wave_unit: u.Unit = "nm",
        obs_type: str = None,
        deduplicate: Callable | None = None,
        correct_for_mw_dust: bool = True,
    ) -> pd.DataFrame:
        """
        Ensure the photometry associated with this transient is all in the same
        units/system/etc

        Args:
            flux_unit (astropy.unit.Unit): Either a valid string to convert
                                           or an astropy.unit.Unit, this can be either
                                           flux, flux density, or magnitude unit. This
                                           supports any base units supported by
                                           synphot
                                           (https://synphot.readthedocs.io/en/latest/synphot/units.html#flux-units).
            date_unit (str): Valid astropy date format string. See
                             https://docs.astropy.org/en/stable/time/index.html#time-format
            freq_unit (astropy.unit.Unit): The astropy unit or string representation of
                                           an astropy unit to convert and return the
                                           frequency as. Must have a base unit of
                                           1/time (Hz).
            wave_unit (astropy.unit.Unit): The astropy unit or string representation of
                                           an astropy unit to convert and return the
                                           wavelength as. Must have a base unit of
                                           length.
            obs_type (str): "radio", "xray", or "uvoir". If provided, it only returns
                            data taken within that range of wavelengths/frequencies.
                            Default is None which will return all of the data.
            deduplicate (Callable|None): A function to be used to remove duplicate
                                         reductions of the same data that produces
                                         different flux values. The default is the
                                         otter.deduplicate_photometry  method,
                                         but you can pass
                                         any callable that takes the output pandas
                                         dataframe as input. Set this to False if you
                                         don't want deduplication to occur.
            correct_for_mw_dust (bool): If True we will automatically correct photometry
                                        for MW dust extinction using the SFD dustmaps
                                        and the Gordon+23 extinction curve assuming
                                        R_V=3.1. Note that this will only correct
                                        photometry in the range 0.0912 - 32 um!
        Returns:
            A pandas DataFrame of the cleaned up photometry in the requested units
        """
        if deduplicate is None:
            deduplicate = self.deduplicate_photometry

        warn_filt = _DuplicateFilter()
        logger.addFilter(warn_filt)

        # these imports need to be here for some reason
        # otherwise the code breaks
        from synphot.units import VEGAMAG, convert_flux
        from synphot.spectrum import SourceSpectrum

        # variable so this warning only displays a single time each time this
        # function is called
        source_map_warning = True

        # turn the photometry key into a pandas dataframe
        if "photometry" not in self:
            raise FailedQueryError("No photometry for this object!")

        dfs = []
        for item in self["photometry"]:
            max_len = 0
            for key, val in item.items():
                if isinstance(val, list) and key != "reference":
                    max_len = max(max_len, len(val))

            for key, val in item.items():
                if not isinstance(val, list) or (
                    isinstance(val, list) and len(val) != max_len
                ):
                    item[key] = [val] * max_len

            df = pd.DataFrame(item)
            dfs.append(df)

        if len(dfs) == 0:
            raise FailedQueryError("No photometry for this object!")
        c = pd.concat(dfs)

        # extract the filter information and substitute in any missing columns
        # because of how we handle this later, we just need to make sure the effective
        # wavelengths are never nan
        def fill_wave(row):
            if "wave_eff" not in row or (
                pd.isna(row.wave_eff) and not pd.isna(row.freq_eff)
            ):
                freq_eff = row.freq_eff * u.Unit(row.freq_units)
                wave_eff = freq_eff.to(u.Unit(wave_unit), equivalencies=u.spectral())
                return wave_eff.value, wave_unit
            elif not pd.isna(row.wave_eff):
                return row.wave_eff, row.wave_units
            else:
                raise ValueError("Missing frequency or wavelength information!")

        filters = pd.DataFrame(self["filter_alias"])
        res = filters.apply(fill_wave, axis=1)
        filters["wave_eff"], filters["wave_units"] = zip(*res)
        # merge the photometry with the filter information
        df = c.merge(filters, on="filter_key")

        # drop irrelevant obs_types before continuing
        if obs_type is not None:
            valid_obs_types = {"radio", "uvoir", "xray"}
            if obs_type not in valid_obs_types:
                raise IOError("Please provide a valid obs_type")
            df = df[df.obs_type == obs_type]

        # add some mockup columns if they don't exist
        if "value" not in df:
            df["value"] = np.nan
            df["value_err"] = np.nan
            df["value_units"] = "NaN"

        # fix some bad units that are old and no longer recognized by astropy
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            df.raw_units = df.raw_units.str.replace("ergs", "erg")
            df.raw_units = ["mag(AB)" if uu == "AB" else uu for uu in df.raw_units]
            df.value_units = df.value_units.str.replace("ergs", "erg")
            df.value_units = ["mag(AB)" if uu == "AB" else uu for uu in df.value_units]

        # merge the raw and value keywords based on the requested flux_units
        # first take everything that just has `raw` and not `value`
        df_raw_only = df[df.value.isna()]
        remaining = df[df.value.notna()]
        if len(remaining) == 0:
            df_raw = df_raw_only
            df_value = []  # this tricks the code later
        else:
            # then take the remaining rows and figure out if we want the raw or value
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                flux_unit_astropy = u.Unit(flux_unit)

                val_unit_filt = np.array(
                    [
                        u.Unit(uu).is_equivalent(flux_unit_astropy)
                        for uu in remaining.value_units
                    ]
                )

            df_value = remaining[val_unit_filt]
            df_raw_and_value = remaining[~val_unit_filt]

            # then merge the raw dataframes
            df_raw = pd.concat([df_raw_only, df_raw_and_value], axis=0)

        # then add columns to these dataframes to convert stuff later
        df_raw = df_raw.assign(
            _flux=df_raw["raw"].values,
            _flux_units=df_raw["raw_units"].values,
            _flux_err=(
                df_raw["raw_err"].values
                if "raw_err" in df_raw
                else [np.nan] * len(df_raw)
            ),
        )

        if len(df_value) == 0:
            df = df_raw
        else:
            df_value = df_value.assign(
                _flux=df_value["value"].values,
                _flux_units=df_value["value_units"].values,
                _flux_err=(
                    df_value["value_err"].values
                    if "value_err" in df_value
                    else [np.nan] * len(df_value)
                ),
            )

            # then merge df_value and df_raw back into one df
            df = pd.concat([df_raw, df_value], axis=0)

        # then, for the rest of the code to work, set the "by" variables to _flux
        by = "_flux"

        # skip rows where 'by' is nan
        df = df[df[by].notna()]

        # filter out anything that has _flux_units == "ct" because we can't convert that
        try:
            # this is a test case to see if we can convert ct -> flux_unit
            convert_flux(
                [1 * u.nm, 2 * u.nm], 1 * u.ct, u.Unit(flux_unit), area=1 * u.m**2
            )
        except u.UnitsError:
            bad_units = df[df._flux_units == "ct"]
            if len(bad_units) > 0:
                logger.warning(
                    f"""Removing {len(bad_units)} photometry points from
                    {self.default_name} because we can't convert them from ct ->
                    {flux_unit}"""
                )
            df = df[df._flux_units != "ct"]

        # convert the ads bibcodes to a string of human readable sources here
        def mappedrefs(row):
            if isinstance(row.reference, list):
                return "<br>".join([self.srcmap[bibcode] for bibcode in row.reference])
            else:
                return self.srcmap[row.reference]

        try:
            df["human_readable_refs"] = df.apply(mappedrefs, axis=1)
        except Exception as exc:
            if source_map_warning:
                source_map_warning = False
                logger.warning(f"Unable to apply the source mapping because {exc}")

            df["human_readable_refs"] = df.reference

        # Figure out what columns are good to groupby in the photometry
        outdata = []

        if "telescope" in df:
            tele = True
            to_grp_by = ["obs_type", by + "_units", "telescope"]
        else:
            tele = False
            to_grp_by = ["obs_type", by + "_units"]

        # Do the conversion based on what we decided to group by
        for groupedby, data in df.groupby(to_grp_by, dropna=False):
            if tele:
                obstype, unit, telescope = groupedby
            else:
                obstype, unit = groupedby
                telescope = None

            # get the photometry in the right type
            unit = data[by + "_units"].unique()
            if len(unit) > 1:
                raise OtterLimitationError(
                    "Can not apply multiple units for different obs_types"
                )

            unit = unit[0]
            isvegamag = "vega" in unit.lower()
            try:
                if isvegamag:
                    astropy_units = VEGAMAG
                elif unit == "AB":
                    # In astropy "AB" is a magnitude SYSTEM not unit and while
                    # u.Unit("AB") will succeed without error, it will not produce
                    # the expected result!
                    # We can assume here that this unit really means astropy's "mag(AB)"
                    astropy_units = u.Unit("mag(AB)")
                else:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        astropy_units = u.Unit(unit)

            except ValueError:
                # this means there is something likely slightly off in the input unit
                # string. Let's try to fix it!
                # here are some common mistakes
                unit = unit.replace("ergs", "erg")
                unit = unit.replace("AB", "mag(AB)")

                astropy_units = u.Unit(unit)

            except ValueError:
                raise ValueError(
                    "Could not coerce your string into astropy unit format!"
                )

            # get the flux data and find the type
            indata = np.array(data[by].astype(float))
            err_key = by + "_err"
            if err_key in data:
                indata_err = np.array(data[by + "_err"].astype(float))
            else:
                indata_err = np.zeros(len(data))

            # convert to an astropy quantity
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                q = indata * u.Unit(astropy_units)
                q_err = indata_err * u.Unit(
                    astropy_units
                )  # assume error and values have the same unit

            # get and save the effective wavelength
            # because of cleaning we did to the filter dataframe above wave_eff
            # should NEVER be nan!
            if np.any(pd.isna(data["wave_eff"])):
                raise ValueError("Flushing out the effective wavelength array failed!")

            zz = zip(data["wave_eff"], data["wave_units"])
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                wave_eff = u.Quantity([vv * u.Unit(uu) for vv, uu in zz], wave_unit)
                freq_eff = wave_eff.to(freq_unit, equivalencies=u.spectral())

            data["converted_wave"] = wave_eff.value
            data["converted_wave_unit"] = wave_unit
            data["converted_freq"] = freq_eff.value
            data["converted_freq_unit"] = freq_unit

            # convert using synphot
            # stuff has to be done slightly differently for xray than for the others
            if obstype == "xray":
                if telescope is not None:
                    try:
                        area = XRAY_AREAS[telescope.lower()]
                    except KeyError:
                        raise OtterLimitationError(
                            "Did not find an area corresponding to "
                            + "this telescope, please add to util!"
                        )
                else:
                    raise OtterLimitationError(
                        "Can not convert x-ray data without a telescope"
                    )

                # we also need to make this wave_min and wave_max
                # instead of just the effective wavelength like for radio and uvoir
                zz = zip(data["wave_min"], data["wave_max"], data["wave_units"])
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore")
                    wave_eff = u.Quantity(
                        [np.array([m, M]) * u.Unit(uu) for m, M, uu in zz],
                        u.Unit(wave_unit),
                    )

            else:
                area = None

            if obstype == "xray" or isvegamag:
                # we unfortunately have to loop over the points here because
                # syncphot does not work with a 2D array of min max wavelengths
                # for converting counts to other flux units. It also can't convert
                # vega mags with a wavelength array because it interprets that as the
                # wavelengths corresponding to the SourceSpectrum.from_vega()

                flux, flux_err = [], []
                for wave, xray_point, xray_point_err in zip(wave_eff, q, q_err):
                    with warnings.catch_warnings():
                        warnings.filterwarnings("ignore")
                        f_val = convert_flux(
                            wave,
                            xray_point,
                            u.Unit(flux_unit),
                            vegaspec=SourceSpectrum.from_vega(),
                            area=area,
                        ).value

                    # approximate the uncertainty as dX = dY/Y * X
                    f_err = np.multiply(
                        f_val, np.divide(xray_point_err.value, xray_point.value)
                    )

                    # then we take the average of the minimum and maximum values
                    # computed by syncphot
                    flux.append(np.mean(f_val))
                    flux_err.append(np.mean(f_err))

            else:
                # this will be faster and cover most cases
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore")
                    flux = convert_flux(wave_eff, q, u.Unit(flux_unit)).value

                # since the error propagation is different between logarithmic units
                # and linear units, unfortunately
                if isinstance(u.Unit(flux_unit), u.LogUnit):
                    # approximate the uncertainty as dX = dY/Y * |ln(10)/2.5|
                    prefactor = np.abs(np.log(10) / 2.5)  # this is basically 1
                else:
                    # approximate the uncertainty as dX = dY/Y * X
                    prefactor = flux

                flux_err = np.multiply(prefactor, np.divide(q_err.value, q.value))

            flux = np.array(flux) * u.Unit(flux_unit)
            flux_err = np.array(flux_err) * u.Unit(flux_unit)

            data["converted_flux"] = flux.value
            data["converted_flux_err"] = flux_err.value
            outdata.append(data)

        if len(outdata) == 0:
            raise FailedQueryError()
        outdata = pd.concat(outdata)

        # copy over the flux units
        outdata["converted_flux_unit"] = flux_unit

        # make sure all the datetimes are in the same format here too!!
        times = [
            Time(d, format=f).to_value(date_unit.lower())
            for d, f in zip(outdata.date, outdata.date_format.str.lower())
        ]
        outdata["converted_date"] = times
        outdata["converted_date_unit"] = date_unit

        # compute the upperlimit value based on a 3 sigma detection
        # this is just for rows where we don't already know if it is an upperlimit
        if isinstance(u.Unit(flux_unit), u.LogUnit):
            # this uses the following formula (which is surprising because it means
            # magnitude upperlimits are independent of the actual measurement!)
            # sigma_m > (1/3) * (ln(10)/2.5)
            def is_upperlimit(row):
                if "upperlimit" in row and pd.isna(row.upperlimit):
                    return row.converted_flux_err > np.log(10) / (3 * 2.5)
                else:
                    return row.upperlimit
        else:

            def is_upperlimit(row):
                if "upperlimit" in row and pd.isna(row.upperlimit):
                    return row.converted_flux < 3 * row.converted_flux_err
                else:
                    return row.upperlimit

        outdata["upperlimit"] = outdata.apply(is_upperlimit, axis=1)

        # clean up filter names
        outdata.loc[outdata.obs_type == "uvoir", "filter_name"] = outdata.loc[
            outdata.obs_type == "uvoir", "filter_name"
        ].apply(self._standardize_filter_names)

        # perform some more complex deduplication of the dataset
        if deduplicate:
            outdata = deduplicate(outdata)

        # perform MW dust extinction correction
        if correct_for_mw_dust:
            outdata = self._correct_for_mw_dust(outdata)

        # throw a warning if the output dataframe has UV/Optical/IR or Radio data
        # where we don't know if the dataset has been host corrected or not
        if ("corr_host" not in outdata) or (
            len(outdata[pd.isna(outdata.corr_host) * (outdata.obs_type != "xray")]) >= 0
        ):
            logger.warning(
                f"{self.default_name} has at least one photometry point where it is "
                + "unclear if a host subtraction was performed. This can be especially "
                + "detrimental for UV data. Please consider filtering out UV/Optical/IR"
                + " or radio rows where the corr_host column is null/None/NaN."
            )

        # rearrange the columns so it's more clear
        firstcols = [
            "converted_flux",
            "converted_flux_err",
            "converted_date",
            "converted_wave",
            "converted_freq",
            "converted_flux_unit",
            "converted_date_unit",
            "converted_wave_unit",
            "converted_freq_unit",
            "reference",
            "human_readable_refs",
            "filter_name",
            "obs_type",
            "upperlimit",
        ]

        if "upperlimit" not in outdata:
            outdata["upperlimit"] = False

        othercols = set(outdata.columns) - set(firstcols)
        cols_in_order = firstcols + list(othercols)

        logger.removeFilter(warn_filt)
        return outdata[cols_in_order]

    def get_ebv(self):
        """
        Get the E(B-V) for this transient from SFD.

        Returns:
            The SFD E(B-V) for the MW along this line of sight to the transient.
        """
        skycoord = self.get_skycoord()
        sfd = SFDQuery()
        return sfd(skycoord)

    def _correct_for_mw_dust(
        self, outdata: pd.DataFrame, dust_model: BaseExtRvModel = G23, rv: float = 3.1
    ) -> pd.DataFrame:
        extmod = dust_model(Rv=rv)
        ebv = self.get_ebv()
        wave_unit = u.Unit(outdata.converted_wave_unit.values[0])
        waves = outdata.converted_wave.values * wave_unit
        is_log_flux_unit = isinstance(
            u.Unit(outdata.converted_flux_unit.values[0]), u.LogUnit
        )

        minwav = 1 / max(dust_model.x_range) * u.um
        maxwav = 1 / min(dust_model.x_range) * u.um
        where_wav = np.where((waves > minwav) * (waves < maxwav))[0]
        df_idx = outdata.iloc[where_wav].index

        # first we need to redden any previously corrected values
        # to make sure things are done consistently
        subset = outdata.loc[df_idx]
        if "val_av" in subset and "corr_av" in subset:
            # if it isn't we can assume that corrections are needed
            for val_av, grp in subset[outdata.corr_av == True].groupby("val_av"):
                corr = extmod.extinguish(
                    grp.converted_wave.values * wave_unit, Av=val_av
                )
                if is_log_flux_unit:
                    outdata.loc[grp.index, "converted_flux"] = grp.converted_flux + corr
                else:
                    outdata.loc[grp.index, "converted_flux"] = grp.converted_flux * corr

        # then we need to de-redden the converted flux column
        corr = extmod.extinguish(waves[where_wav], Ebv=ebv)
        if is_log_flux_unit:
            outdata.loc[df_idx, "converted_flux"] = (
                outdata.loc[df_idx, "converted_flux"] - corr
            )
        else:
            outdata.loc[df_idx, "converted_flux"] = (
                outdata.loc[df_idx, "converted_flux"] / corr
            )

        outdata.loc[df_idx, "corr_av"] = True

        return outdata

    def _standardize_filter_names(
        self, filt: str, delimiters: list[str] = [".", "-", " "]
    ) -> list[str]:
        """
        This (private) method is used to clean up the filter names. As an example,
        we want "r.ztf" = "r.ZTF" = "r" but NOT EQUAL to "R" (since capital
        filter names mean something different!). But here's another fun one,
        we want "UVM2.uvot" = "uvm2.uvot" = "uvm2.UVOT" = "UVM2". That one is tricky
        since the sdss and johnson-counsins filters *are* different capitalizations
        but these aren't.
        """
        uppercase_filters = {
            "uvw1",
            "uvw2",
            "uvm2",
            "w1",
            "w2",
            "w3",
            "w4",
        }  # these lowercase filters should be converted to upper case

        newfilt = filt
        for delim in delimiters:
            newfilt = newfilt.split(delim)[0]

        # some additional cleaning
        newfilt = newfilt.strip()
        if newfilt in uppercase_filters:
            newfilt = newfilt.upper()

        return newfilt

    @classmethod
    def deduplicate_photometry(cls, phot: pd.DataFrame, date_tol: int | float = 1):
        """
        This deduplicates a pandas dataframe of photometry that could potentially
        have rows/datasets that are the result of different reductions of the same
        data. This is especially relevant for X-ray and UV observations where different
        reductions can produce different flux values from the same observation.

        The algorithm used here first finds duplicates by normalizing the telescope
        names, then  grouping the dataframe by transient name, norm telescope name,
        filter_key, and the obs_type. It then assumes that data from the same
        reference will not produce duplicated data. Finally, it finds the overlapping
        regions within date +/- date_tol (or between date_min and date_max for binned
        data), and uses any data within that region as duplicated. From there, it
        first tries to choose the reduction that is host subtracted (if only one is
        host subtracted), then if neither or more than one of the datasets are host
        subtracted then it just takes the most recent reduction.

        Args:
            phot (pd.DataFrame): A pandas dataframe of the photometry with keys defined
                                 by the OTTER schema
            date_tol (int|float): The default tolerance (or "uncertainty") to use on the
                                  dates in the "date" column of phot. In days. Defaults
                                  to 1 day.
        """
        # we need to reset the index to keep track of things appropriately
        phot = phot.reset_index(drop=True)

        if "telescope" not in phot:
            phot["telescope"] = np.nan

        # we first have to standardize some columns given some basic assumptions
        phot["_ref_str"] = phot.reference.astype(str)

        # normalize the telescope name so we can group by it
        phot["_norm_tele_name"] = phot.telescope.apply(cls._normalize_tele_name)

        # now find the duplicated data
        dups = []
        phot_grpby = phot.groupby(
            ["_norm_tele_name", "filter_key", "obs_type"], dropna=False
        )
        for (tele, filter_key, obs_type), grp in phot_grpby:
            # by definition, there can only be dups if the name, telescope, and filter
            # are the same

            # if there is only one reference in this group of data, there's no way
            # there are duplicate reductions of the same dataset
            if len(grp._ref_str.unique()) <= 1:
                continue

            # the next trick is that the dates don't need to be the same, but need to
            # fall inside the same range
            grp["_mean_dates"] = grp.apply(cls._convert_dates, axis=1)

            if "date_min" in grp and not np.all(pd.isna(grp.date_min)):
                grp["min_dates"] = grp.apply(
                    lambda row: cls._convert_dates(row, date_key="date_min"), axis=1
                ).astype(float)
                grp["max_dates"] = grp.apply(
                    lambda row: cls._convert_dates(row, date_key="date_max"), axis=1
                ).astype(float)

                # in case any of the min_date and max_date in the grp are nan
                grp.fillna(
                    {
                        "min_dates": grp._mean_dates - date_tol,
                        "max_dates": grp._mean_dates + date_tol,
                    },
                    inplace=True,
                )

            elif "date_err" in grp and not np.any(pd.isna(grp.date_err)):
                grp["min_dates"] = (grp._mean_dates - grp.date_err).astype(float)
                grp["max_dates"] = (grp._mean_dates + grp.date_err).astype(float)
            else:
                # then assume some uncertainty on the date
                grp["min_dates"] = (grp._mean_dates - date_tol).astype(float)
                grp["max_dates"] = (grp._mean_dates + date_tol).astype(float)

            ref_ranges = [
                (subgrp.min_dates.min(), subgrp.max_dates.max())
                for _, subgrp in grp.groupby("_ref_str")
            ]

            overlaps = cls._find_overlapping_regions(ref_ranges)

            if len(overlaps) == 0:
                continue  # then there are no dups

            for min_overlap, max_overlap in overlaps:
                dup_data = grp[
                    (grp.min_dates >= min_overlap) * (grp.max_dates <= max_overlap)
                ]

                if len(dup_data) == 0:
                    continue  # no data falls in this range!

                dups.append(dup_data)

        # now that we've found the duplicated datasets, we can iterate through them
        # and choose the "default"
        phot_res = deepcopy(phot)
        undupd = []
        for dup in dups:
            try:
                phot_res = phot_res.drop(dup.index)  # we'll append back in the non dup
            except KeyError:
                continue  # we already deleted these ones

            # first, check if only one of the dup reductions host subtracted
            if "corr_host" in dup:
                dup_host_corr = dup[dup.corr_host.astype(bool)]
                host_corr_refs = dup_host_corr.human_readable_refs.unique()
                if len(host_corr_refs) == 1:
                    # then one of the reductions is host corrected and the other isn't!
                    undupd.append(dup[dup.human_readable_refs == host_corr_refs[0]])
                    continue

            bibcodes_sorted_by_year = sorted(dup._ref_str.unique(), key=cls._find_year)
            dataset_to_use = dup[dup._ref_str == bibcodes_sorted_by_year[0]]
            undupd.append(dataset_to_use)

        # then return the full photometry dataset but with the dups removed!
        return pd.concat([phot_res] + undupd).reset_index()

    @staticmethod
    def _normalize_tele_name(tele_name):
        if pd.isna(tele_name):
            return tele_name

        common_delims = ["-", "/", " ", "."]
        for delim in common_delims:
            tele_name = tele_name.replace(delim, ":*:")

        # this assumes that the telescope name will almost always be first,
        # before other delimiters
        return tele_name.split(":*:")[0].lower()

    @staticmethod
    def _convert_dates(row, date_key="date"):
        """Make sure the dates are in MJD"""
        if pd.isna(row[date_key]):
            return row[date_key]

        return Time(row[date_key], format=row.date_format.lower()).mjd

    @staticmethod
    def _find_overlapping_regions(intervals):
        """Find the overlaps in a list of tuples of mins and maxs. This is relatively
        inefficient but the len(intervals) should be < 10 so it should be fine"""
        overlap_ranges = []
        for ii, (start_ii, end_ii) in enumerate(intervals):
            for jj, (start_jj, end_jj) in enumerate(intervals):
                if ii <= jj:
                    continue

                if start_ii > start_jj:
                    start = start_ii
                else:
                    start = start_jj

                if end_ii > end_jj:
                    end = end_jj
                else:
                    end = end_ii

                if start < end:
                    # then there is an overlap!
                    overlap_ranges.append((start, end))

        return overlap_ranges

    @staticmethod
    def _find_year(s):
        match = re.search(r"\d{4}", s)
        return int(match.group()) if match else float("inf")

    def _merge_names(t1, t2, out):  # noqa: N805
        """
        Private method to merge the name data in t1 and t2 and put it in out
        """
        key = "name"
        out[key] = {}

        # first deal with the default_name key
        # we are gonna need to use some regex magic to choose a preferred default_name
        if t1[key]["default_name"] == t2[key]["default_name"]:
            out[key]["default_name"] = t1[key]["default_name"]
        else:
            # we need to decide which default_name is better
            # it should be the one that matches the TNS style
            # let's use regex
            n1 = t1[key]["default_name"]
            n2 = t2[key]["default_name"]

            # write some discriminating regex expressions
            # exp1: starts with a number, this is preferred because it is TNS style
            exp1 = "^[0-9]"
            # exp2: starts with any character, also preferred because it is TNS style
            exp2 = ".$"
            # exp3: checks if first four characters are a number, like a year :),
            # this is pretty strict though
            exp3 = "^[0-9]{3}"
            # exp4: # checks if it starts with AT like TNS names
            exp4 = "^AT"

            # combine all the regex expressions, this makes it easier to add more later
            exps = [exp1, exp2, exp3, exp4]

            # score each default_name based on this
            score1 = 0
            score2 = 0
            for e in exps:
                re1 = re.findall(e, n1)
                re2 = re.findall(e, n2)
                if re1:
                    score1 += 1
                if re2:
                    score2 += 1

            # assign a default_name based on the score
            if score1 > score2:
                out[key]["default_name"] = t1[key]["default_name"]
            elif score2 > score1:
                out[key]["default_name"] = t2[key]["default_name"]
            else:
                logger.warning(
                    "Names have the same score! Just using the existing default_name"
                )
                out[key]["default_name"] = t1[key]["default_name"]

        # now deal with aliases
        # create a reference mapping for each
        t1map = {}
        for val in t1[key]["alias"]:
            ref = val["reference"]
            if isinstance(ref, str):
                t1map[val["value"]] = [ref] if isinstance(ref, str) else list(ref)
            else:
                t1map[val["value"]] = [ref] if isinstance(ref, str) else list(ref)

        t2map = {}
        for val in t2[key]["alias"]:
            ref = val["reference"]
            if isinstance(ref, str):
                t2map[val["value"]] = [ref] if isinstance(ref, str) else list(ref)
            else:
                t2map[val["value"]] = [ref] if isinstance(ref, str) else list(ref)

        # figure out which ones we need to be careful with references in
        inboth = list(
            t1map.keys() & t2map.keys()
        )  # in both so we'll have to merge the reference key
        int1 = list(t1map.keys() - t2map.keys())  # only in t1
        int2 = list(t2map.keys() - t1map.keys())  # only in t2

        # add ones that are not in both first, these are easy
        line1 = [{"value": k, "reference": t1map[k]} for k in int1]
        line2 = [{"value": k, "reference": t2map[k]} for k in int2]
        bothlines = [{"value": k, "reference": t1map[k] + t2map[k]} for k in inboth]
        out[key]["alias"] = line2 + line1 + bothlines

    def _merge_filter_alias(t1, t2, out):  # noqa: N805
        """
        Combine the filter alias lists across the transient objects
        """

        key = "filter_alias"

        out[key] = deepcopy(t1[key])
        keys1 = {filt["filter_key"] for filt in t1[key]}
        for filt in t2[key]:
            if filt["filter_key"] not in keys1:
                out[key].append(filt)

    def _merge_schema_version(t1, t2, out):  # noqa: N805
        """
        Just keep whichever schema version is greater
        """
        key = "schema_version/value"
        if "comment" not in t1["schema_version"]:
            t1["schema_version/comment"] = ""

        if "comment" not in t2["schema_version"]:
            t2["schema_version/comment"] = ""

        if key in t1 and key in t2 and int(t1[key]) > int(t2[key]):
            out["schema_version"] = deepcopy(t1["schema_version"])
        else:
            out["schema_version"] = deepcopy(t2["schema_version"])

        out["schema_version"]["comment"] = (
            t1["schema_version/comment"] + ";" + t2["schema_version/comment"]
        )

    def _merge_photometry(t1, t2, out):  # noqa: N805
        """
        Combine photometry sources
        """

        key = "photometry"

        out[key] = deepcopy(t1[key])
        refs = []  # np.array([d["reference"] for d in out[key]])
        # merge_dups = lambda val: np.sum(val) if np.any(val.isna()) else val.iloc[0]
        for val in out[key]:
            if isinstance(val, list):
                refs += val
            elif isinstance(val, np.ndarray):
                refs += list(val)
            else:
                refs.append(val)

        for val in t2[key]:
            # first check if t2's reference is in out
            if val["reference"] not in refs:
                # it's not here so we can just append the new photometry!
                out[key].append(val)
            else:
                # we need to merge it with other photometry
                i1 = np.where(val["reference"] == refs)[0][0]
                df1 = pd.DataFrame(out[key][i1])
                df2 = pd.DataFrame(val)

                # only substitute in values that are nan in df1 or new
                # the combined keys of the two
                mergeon = list(set(df1.keys()) & set(df2.keys()))
                df = df1.merge(df2, on=mergeon, how="outer")
                # convert to a dictionary
                newdict = df.reset_index().to_dict(orient="list")
                del newdict["index"]

                newdict["reference"] = newdict["reference"][0]

                out[key][i1] = newdict  # replace the dictionary at i1 with the new dict

    def _merge_class(t1, t2, out):  # noqa: N805
        """
        Combine the classification attribute
        """
        key = "classification"
        subkey = "value"
        out[key] = deepcopy(t1[key])
        classes = np.array([item["object_class"] for item in out[key][subkey]])

        for item in t2[key][subkey]:
            if item["object_class"] in classes:
                i = np.where(item["object_class"] == classes)[0][0]
                if int(item["confidence"]) > int(out[key][subkey][i]["confidence"]):
                    out[key][subkey][i]["confidence"] = item[
                        "confidence"
                    ]  # we are now more confident

                if not isinstance(out[key][subkey][i]["reference"], list):
                    out[key][subkey][i]["reference"] = [
                        out[key][subkey][i]["reference"]
                    ]

                if not isinstance(item["reference"], list):
                    item["reference"] = [item["reference"]]

                newdata = list(
                    np.unique(out[key][subkey][i]["reference"] + item["reference"])
                )
                out[key][subkey][i]["reference"] = newdata

            else:
                out[key][subkey].append(item)

        # now that we have all of them we need to figure out which one is the default
        maxconf = max(out[key][subkey], key=lambda d: d["confidence"])
        for item in out[key][subkey]:
            if item == maxconf:
                item["default"] = True
            else:
                item["default"] = False

        # then rederive the classification flags
        out = Transient._derive_classification_flags(out)

    @classmethod
    def _derive_classification_flags(cls, out):
        """
        Derive the classification flags based on the confidence flags. This will find
        - spec_classed
        - unambiguous

        See the paper for a detailed description of how this algorithm makes its
        choices
        """

        if "classification" not in out or "value" not in out["classification"]:
            # this means that the transient doesn't have any classifications
            # just return itself without any changes
            return out

        # get the confidences of all of the classifications of this transient
        confs = np.array(
            [item["confidence"] for item in out["classification"]["value"]]
        ).astype(float)

        all_class_roots = np.array(
            [
                _fuzzy_class_root(item["object_class"])
                for item in out["classification"]["value"]
            ]
        )

        if np.any(confs >= 3):
            unambiguous = len(np.unique(all_class_roots)) == 1
            if np.any(confs == 3) or np.any(confs == 3.3):
                # this is a "gold spectrum"
                spec_classed = 3
            elif np.any(confs == 3.2):
                # this is a silver spectrum
                spec_classed = 2
            elif np.any(confs == 3.1):
                # this is a bronze spectrum
                spec_classed = 1
            else:
                raise ValueError("Not prepared for this confidence flag!")

        elif np.any(confs == 2):
            # these always have spec_classed = True, by definition
            # They also have unambiguous = False by definition because they don't
            # have a peer reviewed citation for their classification
            spec_classed = 1
            unambiguous = False

        elif np.any(confs == 1):
            spec_classed = 0  # by definition
            unambiguous = len(np.unique(all_class_roots)) == 1

        else:
            spec_classed = 0
            unambiguous = False

        # finally, set these keys in the classification dict
        out["classification"]["spec_classed"] = spec_classed
        out["classification"]["unambiguous"] = unambiguous

        return out

    @staticmethod
    def _merge_arbitrary(key, t1, t2, out, merge_subkeys=None, groupby_key=None):
        """
        Merge two arbitrary datasets inside the json file using pandas

        The datasets in t1 and t2 in "key" must be able to be forced into
        a NxM pandas dataframe!
        """

        if key == "name":
            t1._merge_names(t2, out)
        elif key == "filter_alias":
            t1._merge_filter_alias(t2, out)
        elif key == "schema_version":
            t1._merge_schema_version(t2, out)
        elif key == "photometry":
            t1._merge_photometry(t2, out)
        elif key == "classification":
            t1._merge_class(t2, out)
        else:
            # this is where we can standardize some of the merging
            df1 = pd.DataFrame(t1[key])
            df2 = pd.DataFrame(t2[key])

            merged_with_dups = pd.concat([df1, df2]).reset_index(drop=True)

            # have to get the indexes to drop using a string rep of the df
            # this is cause we have lists in some cells
            # We also need to deal with merging the lists of references across rows
            # that we deem to be duplicates. This solution to do this quickly is from
            # https://stackoverflow.com/questions/36271413/ \
            # pandas-merge-nearly-duplicate-rows-based-on-column-value
            if merge_subkeys is None:
                merge_subkeys = merged_with_dups.columns.tolist()
                merge_subkeys.remove("reference")
            else:
                for k in merge_subkeys:
                    if k not in merged_with_dups:
                        merge_subkeys.remove(k)

            merged = (
                merged_with_dups.astype(str)
                .groupby(merge_subkeys)["reference"]
                .apply(lambda x: x.sum())
                .reset_index()
            )

            # then we have to turn the merged reference strings into a string list
            merged["reference"] = merged.reference.str.replace("][", ",")

            # then eval the string of a list to get back an actual list of sources
            merged["reference"] = merged.reference.apply(
                lambda v: np.unique(eval(v)).tolist()
            )

            # decide on default values
            if groupby_key is None:
                iterate_through = [(0, merged)]
            else:
                iterate_through = merged.groupby(groupby_key)

            # we will make whichever value has more references the default
            outdict = []
            for data_type, df in iterate_through:
                lengths = df.reference.map(len)
                max_idx_arr = np.argmax(lengths)

                if isinstance(max_idx_arr, np.int64):
                    max_idx = max_idx_arr
                elif len(max_idx_arr) == 0:
                    raise ValueError("Something went wrong with deciding the default")
                else:
                    max_idx = max_idx_arr[0]  # arbitrarily choose the first

                defaults = np.full(len(df), False, dtype=bool)
                defaults[max_idx] = True

                df["default"] = defaults
                outdict.append(df)
            outdict = pd.concat(outdict)

            # from https://stackoverflow.com/questions/52504972/ \
            # converting-a-pandas-df-to-json-without-nan
            outdict = outdict.replace("nan", np.nan)
            outdict_cleaned = [{**x[i]} for i, x in outdict.stack().groupby(level=0)]

            out[key] = outdict_cleaned


def _fuzzy_class_root(s):
    """
    Extract the fuzzy classification root name from the string s
    """
    s = s.upper()
    # first split the class s using regex
    for root in _KNOWN_CLASS_ROOTS:
        if s.startswith(root):
            remaining = s[len(root) :]
            if remaining and root == "SN":
                # we want to be able to distinguish between SN Ia and SN II
                # we will use SN Ia to indicate thoes and SN to indicate CCSN
                if "IA" in remaining or "1A" in remaining:
                    return "SN Ia"
            return root
    return s
