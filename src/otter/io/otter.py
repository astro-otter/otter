"""
This is the primary class for user interaction with the catalog
"""

from __future__ import annotations
import os
import json
import glob
from warnings import warn

import pandas as pd

from astropy.coordinates import SkyCoord, search_around_sky
from astropy.table import Table
from astropy import units as u

from .transient import Transient
from ..exceptions import FailedQueryError, OtterLimitationError

import warnings

warnings.simplefilter("once", RuntimeWarning)
warnings.simplefilter("once", UserWarning)
warnings.simplefilter("once", u.UnitsWarning)


class Otter(object):
    """
    This is the primary class for users to access the otter backend database

    Args:
        datadir (str): Path to the data directory with the otter data. If not provided
                       will default to a ".otter" directory in the CWD where you call
                       this class from.
        debug (bool): If we should just debug and not do anything serious.

    """

    def __init__(self, datadir: str = None, debug: bool = False) -> None:
        # save inputs
        if datadir is None:
            self.CWD = os.path.dirname(os.path.abspath("__FILE__"))
            self.DATADIR = os.path.join(self.CWD, ".otter")
        else:
            self.CWD = os.path.dirname(datadir)
            self.DATADIR = datadir

        self.debug = debug

        # make sure the data directory exists
        if not os.path.exists(self.DATADIR):
            try:
                os.makedirs(self.DATADIR)
            except FileExistsError:
                warn(
                    "Directory was created between the if statement and trying "
                    + "to create the directory!"
                )
                pass

    def get_meta(self, **kwargs) -> Table:
        """
        Get the metadata of the objects matching the arguments

        Args:
            **kwargs : Arguments to pass to Otter.query()
        Return:
           The metadata for the transients that match the arguments. Will be an astropy
           Table by default, if raw=True will be a dictionary.
        """
        metakeys = [
            "name",
            "coordinate",
            "date_reference",
            "distance",
            "classification",
        ]

        return [t[metakeys] for t in self.query(**kwargs)]

    def cone_search(
        self, coords: SkyCoord, radius: float = 5, raw: bool = False
    ) -> Table:
        """
        Performs a cone search of the catalog over the given coords and radius.

        Args:
            coords (SkyCoord): An astropy SkyCoord object with coordinates to match to
            radius (float): The radius of the cone in arcseconds, default is 0.05"
            raw (bool): If False (the default) return an astropy table of the metadata
                        for matching objects. Otherwise, return the raw json dicts

        Return:
            The metadata for the transients in coords+radius. Will return an astropy
            Table if raw is False, otherwise a dict.
        """

        transients = self.query(coords=coords, radius=radius, raw=raw)

        return transients

    def get_phot(
        self,
        flux_unit="mag(AB)",
        date_unit="MJD",
        return_type="astropy",
        obs_type=None,
        keep_raw=False,
        wave_unit="nm",
        freq_unit="GHz",
        **kwargs,
    ) -> Table:
        """
        Get the photometry of the objects matching the arguments. This will do the
        unit conversion for you!

        Args:
            flux_units (astropy.unit.Unit): Either a valid string to convert
                                            or an astropy.unit.Unit
            date_units (astropy.unit.Unit): Either a valid string to convert to a date
                                            or an astropy.unit.Unit
            return_type (str): Either 'astropy' or 'pandas'. If astropy, returns an
                               astropy Table. If pandas, returns a pandas DataFrame.
                               Default is 'astropy'.
            obs_type (str): Either 'radio', 'uvoir', or 'xray'. Will only return that
                            type of photometry if not None. Default is None and will
                            return any type of photometry.
            keep_raw (bool): If True, keep the raw flux/date/freq/wave associated with
                             the dataset. Else, just keep the converted data. Default
                             is False.
            **kwargs : Arguments to pass to Otter.query(). Can be::

                       names (list[str]): A list of names to get the metadata for
                       coords (SkyCoord): An astropy SkyCoord object with coordinates
                                          to match to
                       radius (float): The radius in arcseconds for a cone search,
                                       default is 0.05"
                       minZ (float): The minimum redshift to search for
                       maxZ (float): The maximum redshift to search for
                       refs (list[str]): A list of ads bibcodes to match to. Will only
                                         return metadata for transients that have this
                                         as a reference.
                       hasSpec (bool): if True, only return events that have spectra.

        Return:
            The photometry for the requested transients that match the arguments.
            Will be an astropy Table sorted by transient default name.

        Raises:
            FailedQueryError: When the query returns no results
            IOError: if one of your inputs is incorrect
        """
        queryres = self.query(hasphot=True, **kwargs)

        dicts = []
        for transient in queryres:
            # clean the photometry
            default_name = transient["name/default_name"]

            try:
                phot = transient.clean_photometry(
                    flux_unit=flux_unit,
                    date_unit=date_unit,
                    wave_unit=wave_unit,
                    freq_unit=freq_unit,
                    obs_type=obs_type,
                )

                phot["name"] = [default_name] * len(phot)

                dicts.append(phot)

            except FailedQueryError:
                # This is fine, it just means that there is no data associated
                # with this one transient. We'll check and make sure there is data
                # associated with at least one of the transients later!
                pass

        if len(dicts) == 0:
            raise FailedQueryError()
        fullphot = pd.concat(dicts)

        # remove some possibly confusing keys
        keys_to_keep = [
            "name",
            "converted_flux",
            "converted_flux_err",
            "converted_date",
            "converted_wave",
            "converted_freq",
            "converted_flux_unit",
            "converted_date_unit",
            "converted_wave_unit",
            "converted_freq_unit",
            "obs_type",
            "upperlimit",
            "reference",
        ]

        if not keep_raw:
            if "telescope" in fullphot:
                fullphot = fullphot[keys_to_keep + ["telescope"]]
            else:
                fullphot = fullphot[keys_to_keep]

        if return_type == "astropy":
            return Table.from_pandas(fullphot)
        elif return_type == "pandas":
            return fullphot
        else:
            raise IOError("return_type can only be pandas or astropy")

    def load_file(self, filename: str) -> dict:
        """
        Loads an otter JSON file

        Args:
            filename (str): The path to the OTTER JSON file to load
        """

        # read in files from summary
        with open(filename, "r") as f:
            to_ret = Transient(json.load(f))

        return to_ret

    def query(
        self,
        names: list[str] = None,
        coords: SkyCoord = None,
        radius: float = 5,
        minz: float = None,
        maxz: float = None,
        refs: list[str] = None,
        hasphot: bool = False,
        hasspec: bool = False,
        raw: bool = False,
    ) -> dict:
        """
        Searches the summary.csv table and reads relevant JSON files

        WARNING! This does not do any conversions for you!
        This is how it differs from the `get_meta` method. Users should prefer to use
        `get_meta`, `getPhot`, and `getSpec` independently because it is a better
        workflow and can return the data in an astropy table with everything in the
        same units.

        Args:
            names (list[str]): A list of names to get the metadata for
            coords (SkyCoord): An astropy SkyCoord object with coordinates to match to
            radius (float): The radius in arcseconds for a cone search, default is 0.05"
            minz (float): The minimum redshift to search for
            maxz (float): The maximum redshift to search for
            refs (list[str]): A list of ads bibcodes to match to. Will only return
                              metadata for transients that have this as a reference.
            hasphot (bool): if True, only returns transients which have photometry.
            hasspec (bool): if True, only return transients that have spectra.

        Return:
           Get all of the raw (unconverted!) data for objects that match the criteria.
        """
        if (
            all(arg is None for arg in [names, coords, maxz, minz, refs])
            and not hasphot
            and not hasspec
        ):
            # there's nothing to query!
            # read in the metdata from all json files
            # this could be dangerous later on!!
            allfiles = glob.glob(os.path.join(self.DATADIR, "*.json"))
            jsondata = []

            # read the data from all the json files and convert to Transients
            for jsonfile in allfiles:
                with open(jsonfile, "r") as j:
                    t = Transient(json.load(j))
                    jsondata.append(t.get_meta())

            return jsondata

        # check if the summary table exists, if it doen't create it
        summary_table = os.path.join(self.DATADIR, "summary.csv")
        if not os.path.exists(summary_table):
            self.generate_summary_table(save=True)

        # then read and query the summary table
        summary = pd.read_csv(summary_table)

        # coordinate search first
        if coords is not None:
            if not isinstance(coords, SkyCoord):
                raise ValueError("Input coordinate must be an astropy SkyCoord!")
            summary_coords = SkyCoord(
                summary.ra.tolist(), summary.dec.tolist(), unit=(u.deg, u.deg)
            )

            try:
                summary_idx, _, _, _ = search_around_sky(
                    summary_coords, coords, seplimit=radius * u.arcsec
                )
            except ValueError:
                summary_idx, _, _, _ = search_around_sky(
                    summary_coords,
                    SkyCoord([coords]),
                    seplimit=radius * u.arcsec,
                )

            summary = summary.iloc[summary_idx]

        # redshift
        if minz is not None:
            summary = summary[summary.z.astype(float) >= minz]

        if maxz is not None:
            summary = summary[summary.z.astype(float) <= maxz]

        # check photometry and spectra
        if hasphot:
            summary = summary[summary.hasPhot == True]

        if hasspec:
            summary = summary[summary.hasSpec == True]

        # check names
        if names is not None:
            if isinstance(names, str):
                n = {names}
            else:
                n = set(names)

            checknames = []
            for alias_row in summary.alias:
                rs = set(eval(alias_row))
                intersection = list(n & rs)
                checknames.append(len(intersection) > 0)

            summary = summary[checknames]

        # check references
        if refs is not None:
            checkrefs = []

            if isinstance(refs, str):
                n = {refs}
            else:
                n = set(refs)

            for ref_row in summary.refs:
                rs = set(eval(ref_row))
                intersection = list(n & rs)
                checkrefs.append(len(intersection) > 0)

            summary = summary[checkrefs]

        outdata = [self.load_file(path) for path in summary.json_path]

        return outdata

    def save(self, schema: list[dict], testing=False) -> None:
        """
        Upload all the data in the given list of schemas.

        Args:
            schema (list[dict]): A list of json dictionaries
            testing (bool): Should we just enter test mode? Default is False

        Raises:
            OtterLimitationError: If some objects in OTTER are within 5" we can't figure
                                  out which ones to merge with which ones.
        """

        if not isinstance(schema, list):
            schema = [schema]

        for transient in schema:
            # convert the json to a Transient
            if not isinstance(transient, Transient):
                transient = Transient(transient)

            print(transient["name/default_name"])

            coord = transient.get_skycoord()
            res = self.cone_search(coords=coord)

            if len(res) == 0:
                # This is a new object to upload
                print("Adding this as a new object...")
                self._save_document(dict(transient), test_mode=testing)

            else:
                # We must merge this with existing data
                print("Found this object in the database already, merging the data...")
                if len(res) == 1:
                    # we can just add these to merge them!
                    combined = res[0] + transient
                    self._save_document(combined, test_mode=testing)
                else:
                    # for now throw an error
                    # this is a limitation we can come back to fix if it is causing
                    # problems though!
                    raise OtterLimitationError("Some objects in Otter are too close!")

        # update the summary table appropriately
        self.generate_summary_table(save=True)

    def _save_document(self, schema, test_mode=False):
        """
        Save a json file in the correct format to the OTTER data directory
        """
        # check if this documents key is in the database already
        # and if so remove it!
        jsonpath = os.path.join(self.DATADIR, "*.json")
        aliases = {item["value"].replace(" ", "-") for item in schema["name"]["alias"]}
        filenames = {
            os.path.basename(fname).split(".")[0] for fname in glob.glob(jsonpath)
        }
        todel = list(aliases & filenames)

        # now save this data
        # create a new file in self.DATADIR with this
        if len(todel) > 0:
            outfilepath = os.path.join(self.DATADIR, todel[0] + ".json")
            if test_mode:
                print("Renaming the following file for backups: ", outfilepath)
        else:
            if test_mode:
                print("Don't need to mess with the files at all!")
            fname = schema["name"]["default_name"] + ".json"
            fname = fname.replace(" ", "-")  # replace spaces in the filename
            outfilepath = os.path.join(self.DATADIR, fname)

        # format as a json
        if isinstance(schema, Transient):
            schema = dict(schema)

        out = json.dumps(schema, indent=4)
        # out = '[' + out
        # out += ']'

        if not test_mode:
            with open(outfilepath, "w") as f:
                f.write(out)
        else:
            print(f"Would write to {outfilepath}")
            # print(out)

    def generate_summary_table(self, save=False) -> pd.DataFrame:
        """
        Generate a summary table for the JSON files in self.DATADIR

        args:
            save (bool): if True, save the summary file to "summary.csv"
                         in self.DATADIR. Default is False and is just returned.

        returns:
            pandas.DataFrame of the summary meta information of the transients
        """
        allfiles = glob.glob(os.path.join(self.DATADIR, "*.json"))

        # read the data from all the json files and convert to Transients
        rows = []
        for jsonfile in allfiles:
            with open(jsonfile, "r") as j:
                t = Transient(json.load(j))
                skycoord = t.get_skycoord()

                row = {
                    "name": t.default_name,
                    "alias": [alias["value"] for alias in t["name"]["alias"]],
                    "ra": skycoord.ra,
                    "dec": skycoord.dec,
                    "refs": [ref["name"] for ref in t["reference_alias"]],
                }

                if "date_reference" in t:
                    row["discovery_date"] = t.get_discovery_date()

                if "distance" in t:
                    row["z"] = t.get_redshift()

                row["hasPhot"] = "photometry" in t
                row["hasSpec"] = "spectra" in t

                row["json_path"] = os.path.abspath(jsonfile)

                rows.append(row)

        alljsons = pd.DataFrame(rows)
        if save:
            alljsons.to_csv(os.path.join(self.DATADIR, "summary.csv"))

        return alljsons
