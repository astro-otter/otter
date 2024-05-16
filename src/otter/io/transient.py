"""
Class for a transient,
basically just inherits the dict properties with some overwriting
"""

from __future__ import annotations
import warnings
from copy import deepcopy
import re
from collections.abc import MutableMapping
from typing_extensions import Self

import numpy as np
import pandas as pd

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord

from synphot.units import VEGAMAG, convert_flux
from synphot.spectrum import SourceSpectrum

from ..exceptions import (
    FailedQueryError,
    IOError,
    OtterLimitationError,
    TransientMergeError,
)
from ..util import XRAY_AREAS

warnings.simplefilter("once", RuntimeWarning)
warnings.simplefilter("once", UserWarning)
np.seterr(divide="ignore")


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
            return Transient({key: (self[key] if key in self else []) for key in keys})
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

    def __repr__(self, html=False):
        if not html:
            return f"Transient(\n\tName: {self.default_name},\n\tKeys: {self.keys()}\n)"
        else:
            html = ""

            coord = self.get_skycoord()

            # add the ra and dec
            # These are required so no need to check if they are there
            html += f"""
            <tr>
            <td style="text-align:left">RA [hrs]:</td>
            <td style="text-align:left">{coord.ra}
            </tr>
            <tr>
            <td style="text-align:left">DEC [deg]:</td>
            <td style="text-align:left">{coord.dec}
            </tr>
            """

            if "date_reference" in self:
                discovery = self.getDiscoveryDate().to_value("datetime")
                if discovery is not None:
                    # add the discovery date
                    html += f"""
                    <tr>
                    <td style="text-align:left">Discovery Date [MJD]:</td>
                    <td style="text-align:left">{discovery}
                    </tr>
                    """

            if "distance" in self:
                # add the redshift
                html += f"""
                <tr>
                <td style="text-align:left">Redshift:</td>
                <td style="text-align:left">{self['distance']['redshift'][0]['value']}
                </tr>
                """

            if "reference_alias" in self:
                srcs = ""
                for bibcode, src in self.srcmap.items():
                    srcs += f"<a href='https://ui.adsabs.harvard.edu/abs/{bibcode}'"
                    srcs += f"target='_blank'>{src}</a><br>"

                html += f"""
                <tr>
                <td style="text-align:left">Sources:</td>
                <td style="text-align:left">{srcs}
                </tr>
                """

            return html

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
            if key == "name":
                self._merge_names(other, out)
            elif key == "coordinate":
                self._merge_coords(other, out)
            elif key == "date_reference":
                self._merge_date(other, out)
            elif key == "distance":
                self._merge_distance(other, out)
            elif key == "filter_alias":
                self._merge_filter_alias(other, out)
            elif key == "schema_version":
                self._merge_schema_version(other, out)
            elif key == "photometry":
                self._merge_photometry(other, out)
            elif key == "spectra":
                self._merge_spectra(other, out)
            elif key == "classification":
                self._merge_class(other, out)
            else:
                # this is an unexpected key!
                if strict_merge:
                    # since this is a strict merge we don't want unexpected data!
                    raise TransientMergeError(
                        f"{key} was not expected! Only keeping the old information!"
                    )
                else:
                    # Throw a warning and only keep the old stuff
                    warnings.warn(
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
                warnings.warn("Not returing the photometry!")
                _ = keys.pop("photometry")
            if "spectra" in keys:
                warnings.warn("Not returning the spectra!")
                _ = keys.pop("spectra")

            curr_keys = self.keys()
            for key in keys:
                if key not in curr_keys:
                    keys.remove(key)
                    warnings.warn(
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
        f = "df['coordinate_type'] == 'equitorial'"
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
        date = self._get_default(key, filt='df["date_type"] == "discovery"')
        if "date_format" in date:
            f = date["date_format"]
        else:
            f = "mjd"

        return Time(date["value"], format=f)

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
            return default["value"]

    def _get_default(self, key, filt=""):
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
            # this is an equitorial coordinate
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
        by: str = "raw",
        obs_type: str = None,
    ) -> pd.DataFrame:
        """
        Ensure the photometry associated with this transient is all in the same
        units/system/etc

        Args:
            flux_unit (astropy.unit.Unit): The astropy unit or string representation of
                                           an astropy unit to convert and return the
                                           flux as.
            date_unit (str): Valid astropy date format string.
            freq_unit (astropy.unit.Unit): The astropy unit or string representation of
                                           an astropy unit to convert and return the
                                           frequency as.
            wave_unit (astropy.unit.Unit): The astropy unit or string representation of
                                           an astropy unit to convert and return the
                                           wavelength as.
            by (str): Either 'raw' or 'value'. 'raw' is the default and is highly
                      recommended! If 'value' is used it may skip some photometry.
                      See the schema definition to understand this keyword completely
                      before using it.
            obs_type (str): "radio", "xray", or "uvoir". If provided, it only returns
                            data taken within that range of wavelengths/frequencies.
                            Default is None which will return all of the data.

        Returns:
            A pandas DataFrame of the cleaned up photometry in the requested units
        """

        # check inputs
        if by not in {"value", "raw"}:
            raise IOError("Please choose either value or raw!")

        # turn the photometry key into a pandas dataframe
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

        c = pd.concat(dfs)

        filters = pd.DataFrame(self["filter_alias"])
        df = c.merge(filters, on="filter_key")

        # make sure 'by' is in df
        if by not in df:
            if by == "value":
                by = "raw"
            else:
                by = "value"

        # skip rows where 'by' is nan
        df = df[df[by].notna()]

        # drop irrelevant obs_types before continuing
        if obs_type is not None:
            valid_obs_types = {"radio", "uvoir", "xray"}
            if obs_type not in valid_obs_types:
                raise IOError("Please provide a valid obs_type")
            df = df[df.obs_type == obs_type]

        # convert the ads bibcodes to a string of human readable sources here
        def mappedrefs(row):
            if isinstance(row.reference, list):
                return "<br>".join([self.srcmap[bibcode] for bibcode in row.reference])
            else:
                return self.srcmap[row.reference]

        try:
            df["human_readable_refs"] = df.apply(mappedrefs, axis=1)
        except Exception as exc:
            warnings.warn(f"Unable to apply the source mapping because {exc}")
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
            try:
                if "vega" in unit.lower():
                    astropy_units = VEGAMAG
                else:
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
            q = indata * u.Unit(astropy_units)
            q_err = indata_err * u.Unit(
                astropy_units
            )  # assume error and values have the same unit

            # get the effective wavelength
            if "freq_eff" in data and not np.isnan(data["freq_eff"].iloc[0]):
                freq_units = data["freq_units"]
                if len(np.unique(freq_units)) > 1:
                    raise OtterLimitationError(
                        "Can not convert different units to the same unit!"
                    )

                freq_eff = np.array(data["freq_eff"]) * u.Unit(freq_units.iloc[0])
                wave_eff = freq_eff.to(u.AA, equivalencies=u.spectral())

            elif "wave_eff" in data and not np.isnan(data["wave_eff"].iloc[0]):
                wave_units = data["wave_units"]
                if len(np.unique(wave_units)) > 1:
                    raise OtterLimitationError(
                        "Can not convert different units to the same unit!"
                    )

                wave_eff = np.array(data["wave_eff"]) * u.Unit(wave_units.iloc[0])

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
                        "Can not convert x-ray data without a " + "telescope"
                    )

                # we also need to make this wave_min and wave_max
                # instead of just the effective wavelength like for radio and uvoir
                wave_eff = np.array(
                    list(zip(data["wave_min"], data["wave_max"]))
                ) * u.Unit(wave_units.iloc[0])

            else:
                area = None

            # we unfortunately have to loop over the points here because
            # syncphot does not work with a 2D array of min max wavelengths
            # for converting counts to other flux units. It also can't convert
            # vega mags with a wavelength array because it then interprets that as the
            # wavelengths corresponding to the SourceSpectrum.from_vega()
            flux, flux_err = [], []
            for wave, xray_point, xray_point_err in zip(wave_eff, q, q_err):
                f_val = convert_flux(
                    wave,
                    xray_point,
                    u.Unit(flux_unit),
                    vegaspec=SourceSpectrum.from_vega(),
                    area=area,
                )
                f_err = convert_flux(
                    wave,
                    xray_point_err,
                    u.Unit(flux_unit),
                    vegaspec=SourceSpectrum.from_vega(),
                    area=area,
                )

                # then we take the average of the minimum and maximum values
                # computed by syncphot
                flux.append(np.mean(f_val).value)
                flux_err.append(np.mean(f_err).value)

            flux = np.array(flux) * u.Unit(flux_unit)
            flux_err = np.array(flux_err) * u.Unit(flux_unit)

            data["converted_flux"] = flux.value
            data["converted_flux_err"] = flux_err.value
            outdata.append(data)

        if len(outdata) == 0:
            raise FailedQueryError()
        outdata = pd.concat(outdata)

        # copy over the flux units
        outdata["converted_flux_unit"] = [flux_unit] * len(outdata)

        # make sure all the datetimes are in the same format here too!!
        times = [
            Time(d, format=f).to_value(date_unit.lower())
            for d, f in zip(outdata.date, outdata.date_format.str.lower())
        ]
        outdata["converted_date"] = times
        outdata["converted_date_unit"] = [date_unit] * len(outdata)

        # same with frequencies and wavelengths
        freqs = []
        waves = []

        for _, row in df.iterrows():
            if "freq_eff" in row and not np.isnan(row["freq_eff"]):
                val = row["freq_eff"] * u.Unit(row["freq_units"])
            elif "wave_eff" in df and not np.isnan(row["wave_eff"]):
                val = row["wave_eff"] * u.Unit(row["wave_units"])
            else:
                raise ValueError("No known frequency or wavelength, please fix!")

            freqs.append(val.to(freq_unit, equivalencies=u.spectral()).value)
            waves.append(val.to(wave_unit, equivalencies=u.spectral()).value)

        outdata["converted_freq"] = freqs
        outdata["converted_wave"] = waves
        outdata["converted_wave_unit"] = [wave_unit] * len(outdata)
        outdata["converted_freq_unit"] = [freq_unit] * len(outdata)

        return outdata

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
                warnings.warn(
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

    def _merge_coords(t1, t2, out):  # noqa: N805
        """
        Merge the coordinates subdictionaries for t1 and t2 and put it in out

        Use pandas to drop any duplicates
        """
        key = "coordinate"

        Transient._merge_arbitrary(key, t1, t2, out)

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
        if int(t1[key]) > int(t2[key]):
            out["schema_version"] = deepcopy(t1["schema_version"])
        else:
            out["schema_version"] = deepcopy(t2["schema_version"])

    def _merge_photometry(t1, t2, out):  # noqa: N805
        """
        Combine photometry sources
        """

        key = "photometry"

        out[key] = deepcopy(t1[key])
        refs = np.array([d["reference"] for d in out[key]])
        # merge_dups = lambda val: np.sum(val) if np.any(val.isna()) else val.iloc[0]

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

    def _merge_spectra(t1, t2, out):  # noqa: N805
        """
        Combine spectra sources
        """
        pass

    def _merge_class(t1, t2, out):  # noqa: N805
        """
        Combine the classification attribute
        """
        key = "classification"
        out[key] = deepcopy(t1[key])
        classes = np.array([item["object_class"] for item in out[key]])
        for item in t2[key]:
            if item["object_class"] in classes:
                i = np.where(item["object_class"] == classes)[0][0]
                if int(item["confidence"]) > int(out[key][i]["confidence"]):
                    out[key][i]["confidence"] = item[
                        "confidence"
                    ]  # we are now more confident

                if not isinstance(out[key][i]["reference"], list):
                    out[key][i]["reference"] = [out[key][i]["reference"]]

                if not isinstance(item["reference"], list):
                    item["reference"] = [item["reference"]]

                newdata = list(np.unique(out[key][i]["reference"] + item["reference"]))
                out[key][i]["reference"] = newdata

            else:
                out[key].append(item)

        # now that we have all of them we need to figure out which one is the default
        maxconf = max(out[key], key=lambda d: d["confidence"])
        for item in out[key]:
            if item == maxconf:
                item["default"] = True
            else:
                item["default"] = False

    def _merge_date(t1, t2, out):  # noqa: N805
        """
        Combine epoch data across two transients and write it to "out"
        """
        key = "date_reference"

        Transient._merge_arbitrary(key, t1, t2, out)

    def _merge_distance(t1, t2, out):  # noqa: N805
        """
        Combine distance information for these two transients
        """
        key = "distance"

        Transient._merge_arbitrary(key, t1, t2, out)

    @staticmethod
    def _merge_arbitrary(key, t1, t2, out):
        """
        Merge two arbitrary datasets inside the json file using pandas

        The datasets in t1 and t2 in "key" must be able to be forced into
        a NxM pandas dataframe!
        """

        df1 = pd.DataFrame(t1[key])
        df2 = pd.DataFrame(t2[key])

        merged_with_dups = pd.concat([df1, df2]).reset_index(drop=True)

        # have to get the indexes to drop using a string rep of the df
        # this is cause we have lists in some cells
        to_drop = merged_with_dups.astype(str).drop_duplicates().index

        merged = merged_with_dups.iloc[to_drop].reset_index(drop=True)

        outdict = merged.to_dict(orient="records")

        outdict_cleaned = Transient._remove_nans(
            outdict
        )  # clear out the nans from pandas conversion

        out[key] = outdict_cleaned

    @staticmethod
    def _remove_nans(d):
        """
        Remove nans from a record dictionary

        THIS IS SLOW: O(n^2)!!! WILL NEED TO BE SPED UP LATER
        """

        outd = []
        for item in d:
            outsubd = {}
            for key, val in item.items():
                if not isinstance(val, float):
                    # this definitely is not NaN
                    outsubd[key] = val

                else:
                    if not np.isnan(val):
                        outsubd[key] = val
            outd.append(outsubd)

        return outd
