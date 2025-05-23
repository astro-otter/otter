"""
This is the primary class for user interaction with the catalog
"""

from __future__ import annotations
from typing import Optional
import os
import json
import glob
from copy import deepcopy
import logging

from pyArango.connection import Connection
from pyArango.database import Database
from pyArango.document import Document

import pandas as pd
import numpy as np

from astropy.coordinates import SkyCoord, search_around_sky
from astropy.table import Table
from astropy import units as u

from .transient import Transient
from ..exceptions import FailedQueryError, OtterLimitationError, TransientMergeError
from ..util import bibcode_to_hrn, freq_to_obstype, freq_to_band, _DuplicateFilter

import warnings

warnings.simplefilter("once", RuntimeWarning)
warnings.simplefilter("once", UserWarning)
warnings.simplefilter("once", u.UnitsWarning)

logger = logging.getLogger(__name__)


def _np_encoder(object):
    """
    Numpy data type encoder for json.dump
    """
    if isinstance(object, (np.generic, np.ndarray)):
        return object.item()


class Otter(Database):
    """
    This is the primary class for users to access the otter backend database

    Args:
        url (str): The url where the database api endpoints are located
        username (str): The username to log into the database with
        password (str): The password to log into the database with
        gen_summary (bool): Generate a local summary table, this should generally be
                            left as False!
        datadir (str): Path to the data directory with the otter data. If not provided
                       will default to a ".otter" directory in the CWD where you call
                       this class from.
        debug (bool): If we should just debug and not do anything serious.

    Returns:
        An Otter object that is connected to the otter database
    """

    def __init__(
        self,
        url: str = "http://127.0.0.1:8529",
        username: str = os.environ.get("ARANGO_USER_USERNAME", "user-guest"),
        password: str = os.environ.get("ARANGO_USER_PASSWORD", ""),
        gen_summary: bool = False,
        datadir: str = None,
        debug: bool = False,
        **kwargs,
    ) -> None:
        print("Attempting to login with the following credentials:")
        print(f"username: {username}")
        print(f"password: {password}")

        # save inputs
        if datadir is None:
            self.CWD = os.path.dirname(os.path.abspath("__FILE__"))
            self.DATADIR = os.path.join(self.CWD, ".otter")
        else:
            self.CWD = os.path.dirname(datadir)
            self.DATADIR = datadir

        self.debug = debug

        if gen_summary:
            self.generate_summary_table(save=True)

        # make sure the data directory exists
        if not os.path.exists(self.DATADIR):
            try:
                os.makedirs(self.DATADIR)
            except FileExistsError:
                logger.warning(
                    "Directory was created between the if statement and trying "
                    + "to create the directory!"
                )
                pass

        connection = Connection(username=username, password=password, arangoURL=url)
        super().__init__(connection, "otter", **kwargs)

    def get_meta(self, **kwargs) -> Table:
        """
        Get the metadata of the objects matching the arguments

        Args:
            **kwargs : Arguments to pass to Otter.query(). See that documentation with
                       `help(otter.Otter.query)`.
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
            "reference_alias",
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
            flux_unit (astropy.unit.Unit): Either a valid string to convert
                                            or an astropy.unit.Unit
            date_unit (astropy.unit.Unit): Either a valid string to convert to a date
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
            wave_unit (str): The astropy wavelength unit to return with
            freq_unit (str): The astropy frequency unit to return with`
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
        warn_filt = _DuplicateFilter()
        logger.addFilter(warn_filt)

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
            "filter_name",
            "obs_type",
            "upperlimit",
            "reference",
            "human_readable_refs",
        ]

        if "upperlimit" not in fullphot:
            fullphot["upperlimit"] = False

        if not keep_raw:
            if "telescope" in fullphot:
                fullphot = fullphot[keys_to_keep + ["telescope"]]
            else:
                fullphot = fullphot[keys_to_keep]

        logger.removeFilter(warn_filt)
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

        Returns:
            dictionary with the otter JSON file contents
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
        spec_classed: bool = False,
        unambiguous: bool = False,
        classification: str = None,
        class_confidence_threshold: float = 0,
        query_private=False,
        **kwargs,
    ) -> dict:
        """
        Searches the arango database table and reads relevant JSON files

        WARNING! This does not do any conversions for you!
        This is how it differs from the `get_meta` method. Users should prefer to use
        `get_meta`, and `get_phot` independently because it is a better
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
            hasspec (bool): NOT IMPLEMENTED! Will return False for all targets!
            spec_classed (bool): If True, only returns transients that have been
                                 specotroscopically classified/confirmed
            unambiguous (bool): If True, only returns transients that only have a single
                                published classification in OTTER. If classifications
                                disagree for a transient, it will be filtered out.
            classification (str): A classification string to search for
            class_confidence_threshold (float): classification confidence cutoff for
                                                query, between 0 and 1. Default is 0.
            query_private (bool): Set to True if you would like to also query the
                                  dataset located at whatever you set datadir to

        Return:
           Get all of the raw (unconverted!) data for objects that match the criteria.
        """
        # write some AQL filters based on the inputs
        query_filters = ""

        if hasphot is True:
            query_filters += "FILTER 'photometry' IN ATTRIBUTES(transient)\n"

        if hasspec is True:
            query_filters += "FILTER 'spectra' IN ATTRIBUTES(transient)\n"

        if spec_classed:
            query_filters += "FILTER transient.classification.spec_classed >= 1"

        if unambiguous:
            query_filters += "FILTER transient.classification.unambiguous"

        if classification is not None:
            query_filters += f"""
            FOR subdoc IN transient.classification.value
                FILTER subdoc.confidence > TO_NUMBER({class_confidence_threshold})
                FILTER subdoc.object_class LIKE '%{classification}%'
            """

        if minz is not None:
            sfilt = f"""
            FILTER 'redshift' IN transient.distance[*].distance_type
            LET redshifts1 = (
                FOR val IN transient.distance
                FILTER val.distance_type == 'redshift'
                FILTER TO_NUMBER(val.value) >= {minz}
                RETURN val
            )
            FILTER COUNT(redshifts1) > 0
            """
            query_filters += sfilt
        if maxz is not None:
            sfilt = f"""
            FILTER 'redshift' IN transient.distance[*].distance_type
            LET redshifts2 = (
                FOR val IN transient.distance
                FILTER val.distance_type == 'redshift'
                FILTER TO_NUMBER(val.value) <= {maxz}
                RETURN val
            )
            FILTER COUNT(redshifts2) > 0
            """
            query_filters += sfilt

        if names is not None:
            if isinstance(names, str):
                query_filters += f"""
            FILTER UPPER(transient.name) LIKE UPPER('%{names}%')\n
                """
            elif isinstance(names, list):
                namefilt = f"""
            FOR name IN {names}
                FILTER name IN transient.name.alias[*].value\n
                """
                query_filters += namefilt
            else:
                raise Exception("Names must be either a string or list")

        if refs is not None:
            if isinstance(refs, str):  # this is just a single bibcode
                query_filters += f"FILTER {refs} IN transient.reference_alias[*].name"
            elif isinstance(refs, list):
                query_filters += f"""
                FOR ref IN {refs}
                    FILTER ref IN transient.reference_alias[*].name
                """
            else:
                raise Exception("reference list must be either a string or a list")

        # define the query
        query = f"""
        FOR transient IN transients
            {query_filters}
            RETURN transient
        """

        # set batch size to 100 million (for now at least)
        result = self.AQLQuery(query, rawResults=True, batchSize=100_000_000)

        # now that we have the query results do the RA and Dec queries if they exist
        if coords is not None:
            # get the catalog RAs and Decs to compare against
            query_coords = coords
            good_tdes = []

            for tde in result:
                for coordinfo in tde["coordinate"]:
                    if "ra" in coordinfo and "dec" in coordinfo:
                        coord = SkyCoord(
                            coordinfo["ra"],
                            coordinfo["dec"],
                            unit=(coordinfo["ra_units"], coordinfo["dec_units"]),
                        )
                    elif "l" in coordinfo and "b" in coordinfo:
                        # this is galactic
                        coord = SkyCoord(
                            coordinfo["l"],
                            coordinfo["b"],
                            unit=(coordinfo["l_units"], coordinfo["b_units"]),
                            frame="galactic",
                        )
                    else:
                        raise ValueError(
                            "Either needs to have ra and dec or l and b as keys!"
                        )
                    if query_coords.separation(coord) < radius * u.arcsec:
                        good_tdes.append(tde)
                        break  # we've confirmed this tde is in the cone!

            arango_query_results = [Transient(t) for t in good_tdes]

        else:
            arango_query_results = [Transient(res) for res in result.result]

        if not query_private:
            return arango_query_results

        private_results = self._query_datadir(
            names=names,
            coords=coords,
            radius=radius,
            minz=minz,
            maxz=maxz,
            refs=refs,
            hasphot=hasphot,
            hasspec=hasspec,
        )

        partially_merged = deepcopy(arango_query_results)
        new_transients = []
        for jj, t_private in enumerate(private_results):
            for ii, t_public in enumerate(arango_query_results):
                try:
                    partially_merged[ii] += t_private
                    break
                except TransientMergeError:
                    continue
            else:
                new_transients.append(t_private)

        return partially_merged + new_transients

    def _query_datadir(
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
        This is a private method and is here just for the pipeline!!!
        This should not be used by end users!

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
            jsondata = [self.load_file(jsonfile) for jsonfile in allfiles]

            return jsondata

        # check if the summary table exists, if it doen't create it
        summary_table = os.path.join(self.DATADIR, "summary.csv")
        if not os.path.exists(summary_table):
            self.generate_summary_table(save=True)

        # then read and query the summary table
        summary = pd.read_csv(summary_table)
        if len(summary) == 0:
            return []

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

    def upload(self, json_data, collection="vetting", testing=False) -> Document:
        """
        Upload json_data to collection WITHOUT deduplication!

        Args:
            json_data [dict] : A dictionary of the json data to upload to Otter
            collection [str] : The collection to upload to
            testing [bool] : Default is False

        Returns:
            The pyArango document that was uplaoded
        """

        # now add the document
        doc = self[collection].createDocument(json_data)
        if not testing:
            doc.save()
        return doc

    def upload_private(self, collection="vetting", testing=False) -> None:
        """
        Upload the local/private data stored in self.DATADIR to the vetting collection
        (like a SQL table) in the central arangodb document database.

        WARNING! This will make any data in self.DATADIR public! Please double check
        before using this method!

        Args:
            collection (str) : The collection to add the documents to. Default is
                               "vetting" where the documents will then be vetted by
                               our team.
            testing (bool) : True if don't actually upload, False is default

        Returns:
            If testing is false (the default), returns the arangodb upload result. If
            testing is true, returns the list of merged dictionaries that would get
            uploaded.

        Raises:
            OtterLimitationError: If some objects in OTTER are within 5" we can't figure
                                  out which ones to merge with which ones.

        """

        if not self.hasCollection(collection):
            raise ValueError(f"{collection} not in {self}!")

        local_data = self._query_datadir()
        docs = []
        for t in local_data:
            res = self.query(coords=t.get_skycoord())

            if len(res) > 1:
                raise OtterLimitationError("Some objects in Otter are too close!")

            elif len(res) == 1:
                # this object exists in otter already, let's grab the transient data and
                # merge the files
                merged = t + res[0]

                # copy over the special arangodb keys
                merged["_key"] = res[0]["_key"]
                merged["_id"] = res[0]["_id"]

                # we also have to delete the document from the OTTER database
                doc = self.fetchDocument(merged["_id"])
                if not testing:
                    doc.delete()
                else:
                    print(f"Would delete\n{doc}")

            else:
                # this means the object doesn't exist in otter already
                merged = t

            docs.append(self.upload(merged, collection=collection, testing=testing))

        return docs

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
            res = self._query_datadir(coords=coord)

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

        out = json.dumps(schema, indent=4, default=_np_encoder)
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
                    date_types = {d["date_type"] for d in t["date_reference"]}
                    if "discovery" in date_types:
                        row["discovery_date"] = t.get_discovery_date()

                if "distance" in t:
                    dist_types = {d["distance_type"] for d in t["distance"]}
                    if "redshift" in dist_types:
                        row["z"] = t.get_redshift()

                row["hasPhot"] = "photometry" in t
                row["hasSpec"] = "spectra" in t

                row["json_path"] = os.path.abspath(jsonfile)

                rows.append(row)

        alljsons = pd.DataFrame(rows)
        if save:
            alljsons.to_csv(os.path.join(self.DATADIR, "summary.csv"))

        return alljsons

    @staticmethod
    def from_csvs(
        metafile: str,
        photfile: str = None,
        local_outpath: Optional[str] = None,
        db: Otter = None,
    ) -> Otter:
        """
        Converts private metadata and photometry csvs to an Otter object stored
        *locally* so you don't need to worry about accidentally uploading them to the
        real Otter database.

        Args:
            metafile (str) : String filepath or string io csv object of the csv metadata
            photfile (str) : String filepath or string io csv object of the csv
                                          photometry
            local_outpath (str) : The outpath to write the OTTER json files to
            db (Otter) : An Otter instance to add the local_outpath to for querying.
                         This keyword can be useful if you have special permission for
                         the otter database and want to upload your private data

        Returns:
            An Otter object where the json files are stored locally
        """
        # read in the metadata and photometry file
        meta = pd.read_csv(metafile)
        meta.columns = meta.columns.str.strip()  # clean up the col names
        phot = None

        required_phot_cols = [
            "name",
            "date",
            "date_format",
            "filter",
            "filter_eff",
            "filter_eff_units",
            "flux",
            "flux_err",
            "flux_unit",
        ]

        if photfile is not None:
            phot_unclean = pd.read_csv(photfile)
            phot_unclean.columns = phot_unclean.columns.str.strip()  # cleanup colnames

            phot = phot_unclean.dropna(subset=required_phot_cols)
            if len(phot) != len(phot_unclean):
                logger.warning("""
                Filtered out rows with nan in the photometry file! Make sure you
                expect this behaviour!
                """)

            if "bibcode" not in phot:
                phot["bibcode"] = "private"
                logger.warning("""
                Setting the bibcode column to the special keyword 'private'!
                """)

            # we need to generate columns of wave_eff and freq_eff
            wave_eff = []
            freq_eff = []
            wave_eff_unit = u.nm
            freq_eff_unit = u.GHz
            for val, unit in zip(phot.filter_eff, phot.filter_eff_units):
                wave_eff.append(
                    (val * u.Unit(unit))
                    .to(wave_eff_unit, equivalencies=u.spectral())
                    .value
                )
                freq_eff.append(
                    (val * u.Unit(unit))
                    .to(freq_eff_unit, equivalencies=u.spectral())
                    .value
                )

            phot["band_eff_wave"] = wave_eff
            phot["band_eff_wave_unit"] = str(wave_eff_unit)
            phot["band_eff_freq"] = freq_eff
            phot["band_eff_freq_unit"] = str(freq_eff_unit)

        if not os.path.exists(local_outpath):
            os.makedirs(local_outpath)

        # drop duplicated names in meta and keep the first
        meta = meta.drop_duplicates(subset="name", keep="first")

        # merge the meta and phot data
        if phot is not None:
            data = pd.merge(phot, meta, on="name", how="inner")
        else:
            data = meta

        # perform some data checks
        assert (
            len(data[pd.isna(data.ra)].name.unique()) == 0
        ), "Missing some RA and Decs, please check the input files!"
        if phot is not None:
            for name in meta.name:
                assert len(data[data.name == name]) == len(
                    phot[phot.name == name]
                ), f"failed on {name}"

        # actually do the data conversion to OTTER
        all_jsons = []
        for name, tde in data.groupby("name"):
            json = {}
            tde = tde.reset_index()

            # name first
            json["name"] = dict(
                default_name=name,
                alias=[dict(value=name, reference=[tde.coord_bibcode[0]])],
            )

            # coordinates
            json["coordinate"] = [
                dict(
                    ra=tde.ra[0],
                    dec=tde.dec[0],
                    ra_units=tde.ra_unit[0],
                    dec_units=tde.dec_unit[0],
                    reference=[tde.coord_bibcode[0]],
                    coordinate_type="equatorial",
                )
            ]

            ### distance info
            json["distance"] = []

            # redshift
            if "redshift" in tde and not np.any(pd.isna(tde["redshift"])):
                json["distance"].append(
                    dict(
                        value=tde.redshift[0],
                        reference=[tde.redshift_bibcode[0]],
                        computed=False,
                        distance_type="redshift",
                    )
                )

            # luminosity distance
            if "luminosity_distance" in tde and not np.any(
                pd.isna(tde["luminosity_distance"])
            ):
                json["distance"].append(
                    value=tde.luminosity_distance[0],
                    reference=[tde.luminosity_distance_bibcode[0]],
                    unit=tde.luminosity_distance_unit[0],
                    computed=False,
                    distance_type="luminosity",
                )

            # comoving distance
            if "comoving_distance" in tde and not np.any(
                pd.isna(tde["comoving_distance"])
            ):
                json["distance"].append(
                    value=tde.comoving_distance[0],
                    reference=[tde.comoving_distance_bibcode[0]],
                    unit=tde.comoving_distance_unit[0],
                    computed=False,
                    distance_type="comoving",
                )

            # remove the distance list if it is empty still
            if len(json["distance"]) == 0:
                del json["distance"]

            ### Classification information that is in the csvs
            # classification
            if "classification" in tde:
                class_flag = 0
                if "classification_flag" in tde:
                    class_flag = tde.classification_flag[0]
                json["classification"] = dict(
                    value=[
                        dict(
                            object_class=tde.classification[0],
                            confidence=class_flag,
                            reference=[tde.classification_bibcode[0]],
                        )
                    ]
                )

            # discovery date
            # print(tde)
            if "discovery_date" in tde and not np.any(pd.isna(tde.discovery_date)):
                json["date_reference"] = [
                    dict(
                        value=str(tde.discovery_date.tolist()[0]).strip(),
                        date_format=tde.discovery_date_format.tolist()[0].lower(),
                        reference=tde.discovery_date_bibcode.tolist(),
                        computed=False,
                        date_type="discovery",
                    )
                ]

            # host information
            if "host_ref" in tde and not np.any(pd.isna(tde.host_ref)):
                host_info = dict(
                    host_name=tde.host_name.tolist()[0].strip(),
                    host_ra=tde.host_ra.tolist()[0],
                    host_dec=tde.host_dec.tolist()[0],
                    host_ra_units=tde.host_ra_unit.tolist()[0],
                    host_dec_units=tde.host_dec_unit.tolist()[0],
                    reference=[tde.host_ref.tolist()[0]],
                )

                if not pd.isna(tde.host_redshift.tolist()[0]):
                    host_info["host_z"] = tde.host_redshift.tolist()[0]

                if "host" in json:
                    json["host"].append(host_info)
                else:
                    json["host"] = [host_info]

            # comments
            if "comment" in tde and not np.any(pd.isna(tde.comment)):
                if "schema_version" not in json:
                    json["schema_version"] = {}
                json["schema_version"]["comment"] = tde.comment.tolist()[0]

            # skip the photometry code if there is no photometry file
            # if there is a photometry file then we want to convert it below
            phot_sources = []
            if phot is not None:
                tde["obs_type"] = [
                    freq_to_obstype(vv * u.Unit(uu))
                    for vv, uu in zip(
                        tde.band_eff_freq.values,
                        tde.band_eff_freq_unit.values,
                    )
                ]

                unique_filter_keys = []
                index_for_match = []
                json["photometry"] = []

                if "telescope" in tde:
                    to_grpby = ["bibcode", "telescope", "obs_type"]
                else:
                    to_grpby = ["bibcode", "obs_type"]

                for grp_keys, p in tde.groupby(to_grpby, dropna=False):
                    if len(grp_keys) == 3:
                        src, tele, obstype = grp_keys
                    else:
                        src, obstype = grp_keys
                        tele = None

                    if src not in phot_sources:
                        phot_sources.append(src)

                    if len(np.unique(p.flux_unit)) == 1:
                        raw_units = p.flux_unit.tolist()[0]
                    else:
                        raw_units = p.flux_unit.tolist()

                    # add a column to phot with the unique filter key
                    if obstype == "radio":
                        filter_uq_key = (
                            p.band_eff_freq.astype(str)
                            + p.band_eff_freq_unit.astype(str)
                        ).tolist()

                    elif obstype in ("uvoir", "xray"):
                        filter_uq_key = p["filter"].astype(str).tolist()

                    else:
                        raise ValueError("not prepared for this obstype!")

                    unique_filter_keys += filter_uq_key
                    index_for_match += p.index.tolist()

                    if "upperlimit" not in p:
                        p["upperlimit"] = False

                    json_phot = dict(
                        reference=src,
                        raw=p.flux.astype(float).tolist(),
                        raw_err=p.flux_err.astype(float).tolist(),
                        raw_units=raw_units,
                        date=p.date.tolist(),
                        date_format=p.date_format.tolist(),
                        upperlimit=p.upperlimit.tolist(),
                        filter_key=filter_uq_key,
                        obs_type=obstype,
                    )

                    if not pd.isna(tele):
                        json_phot["telescope"] = tele

                    if pd.isna(tele) and obstype == "xray":
                        raise ValueError("The telescope is required for X-ray data!")

                    # check the minimum and maximum filter values
                    if obstype == "xray" and (
                        "filter_min" not in p or "filter_max" not in p
                    ):
                        raise ValueError(
                            "Minimum and maximum filters required for X-ray data!"
                        )

                    # check optional keys
                    optional_keys = [
                        "date_err",
                        "sigma",
                        "instrument",
                        "phot_type",
                        "exptime",
                        "aperature",
                        "observer",
                        "reducer",
                        "pipeline",
                    ]
                    for k in optional_keys:
                        if k in p and not np.all(pd.isna(p[k])):
                            # fill the nan values
                            # this is to match with the official json format
                            # and works with arangodb document structure
                            p[k].fillna("null", inplace=True)

                            json_phot[k] = p[k].tolist()

                    # handle more detailed uncertainty information
                    raw_err_detail = {}
                    for key in ["statistical_err", "systematic_err", "iss_err"]:
                        if key in p and not np.all(pd.isna(p[key])):
                            k = key.split("_")[0]

                            # fill the nan values
                            # this is to match with the official json format
                            # and works with arangodb document structure
                            p[key].fillna(0, inplace=True)

                            raw_err_detail[k] = p[key].tolist()

                    if len(raw_err_detail) > 0:
                        json_phot["raw_err_detail"] = raw_err_detail

                    # check the possible corrections
                    corrs = ["val_k", "val_s", "val_host", "val_av", "val_hostav"]
                    for c in corrs:
                        bool_v_key = c.replace("val", "corr")
                        json_phot[c] = False

                        if c in p:
                            # fill the nan values
                            # this is to match with the official json format
                            # and works with arangodb document structure
                            p[c].fillna("null", inplace=True)

                            json_phot[c] = p[c].tolist()
                            json_phot[bool_v_key] = [v != "null" for v in json_phot[c]]

                    json["photometry"].append(json_phot)

                tde["filter_uq_key"] = pd.Series(
                    unique_filter_keys, index=index_for_match
                )

                # filter alias
                # radio filters first
                filter_keys1 = ["filter_uq_key", "band_eff_wave", "band_eff_wave_unit"]
                if "filter_min" in tde:
                    filter_keys1.append("filter_min")
                if "filter_max" in tde:
                    filter_keys1.append("filter_max")

                filter_map = (
                    tde[filter_keys1].drop_duplicates().set_index("filter_uq_key")
                )  # .to_dict(orient='index')
                try:
                    filter_map_radio = filter_map.to_dict(orient="index")
                except Exception:
                    print(filter_map)
                    print(name)
                    raise Exception

                json["filter_alias"] = []
                for filt, val in filter_map_radio.items():
                    obs_type = freq_to_obstype(
                        float(val["band_eff_wave"]) * u.Unit(val["band_eff_wave_unit"])
                    )
                    if obs_type == "radio":
                        filter_name = freq_to_band(
                            (
                                float(val["band_eff_wave"])
                                * u.Unit(val["band_eff_wave_unit"])
                            ).to(u.GHz, equivalencies=u.spectral())
                        )
                    else:
                        filter_name = filt

                    filter_alias_dict = dict(
                        filter_key=filt,
                        filter_name=filter_name,
                        wave_eff=float(val["band_eff_wave"]),
                        wave_units=val["band_eff_wave_unit"],
                    )

                    if "filter_min" in val:
                        filter_alias_dict["wave_min"] = (
                            val["filter_min"] * u.Unit(phot.filter_eff_units)
                        ).to(
                            u.Unit(
                                val["band_eff_wave_unit"], equivalencies=u.spectral()
                            )
                        )

                    if "filter_max" in val:
                        filter_alias_dict["wave_max"] = (
                            val["filter_max"] * u.Unit(phot.filter_eff_units)
                        ).to(
                            u.Unit(
                                val["band_eff_wave_unit"], equivalencies=u.spectral()
                            )
                        )

                    json["filter_alias"].append(filter_alias_dict)

            # reference alias
            # gather all the bibcodes
            all_bibcodes = [tde.coord_bibcode[0]] + phot_sources
            if (
                "redshift_bibcode" in tde
                and tde.redshift_bibcode[0] not in all_bibcodes
                and not np.any(pd.isna(tde.redshift))
            ):
                all_bibcodes.append(tde.redshift_bibcode[0])

            if (
                "luminosity_distance_bibcode" in tde
                and tde.luminosity_distance_bibcode[0] not in all_bibcodes
                and not np.any(pd.isna(tde.luminosity_distance))
            ):
                all_bibcodes.append(tde.luminosity_distance_bibcode[0])

            if (
                "comoving_distance_bibcode" in tde
                and tde.comoving_distance_bibcode[0] not in all_bibcodes
                and not np.any(pd.isna(tde.comoving_distance))
            ):
                all_bibcodes.append(tde.comoving_distance_bibcode[0])

            if (
                "discovery_date_bibcode" in tde
                and tde.discovery_date_bibcode[0] not in all_bibcodes
                and not np.any(pd.isna(tde.discovery_date))
            ):
                all_bibcodes.append(tde.discovery_date_bibcode[0])

            if (
                "classification_bibcode" in tde
                and tde.classification_bibcode[0] not in all_bibcodes
                and not np.any(pd.isna(tde.classification))
            ):
                all_bibcodes.append(tde.classification_bibcode[0])

            if (
                "host_bibcode" in tde
                and tde.host_bibcode not in all_bibcodes
                and not np.any(pd.isna(tde.host_bibcode))
            ):
                all_bibcodes.append(tde.host_bibcode[0])

            # find the hrn's for all of these bibcodes
            uq_bibcodes, all_hrns = bibcode_to_hrn(all_bibcodes)

            # package these into the reference alias
            json["reference_alias"] = [
                dict(name=name, human_readable_name=hrn)
                for name, hrn in zip(uq_bibcodes, all_hrns)
            ]

            all_jsons.append(Transient(json))

        if db is None:
            db = Otter(datadir=local_outpath)
        else:
            db.DATADIR = local_outpath

        # always save this document as a new one
        db.save(all_jsons)
        db.generate_summary_table(save=True)
        return db
