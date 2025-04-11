"""
Host object that stores information on the Transient DataFinder and provides utility
methods for pulling in data corresponding to that host
"""

from __future__ import annotations
import os
import csv
import io
import re
import time
import math
from urllib.request import urlopen
import requests

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.table import Table
from astropy.io.votable import parse_single_table
from astropy.io import ascii

import numpy as np
import pandas as pd
import logging

from fundamentals.stats import rolling_window_sigma_clip
from operator import itemgetter

from ..util import VIZIER_LARGE_CATALOGS
from ..exceptions import MissingEnvVarError

logger = logging.getLogger(__name__)


class DataFinder(object):
    def __init__(
        self,
        ra: str | float,
        dec: str | float,
        ra_units: str | u.Unit,
        dec_units: str | u.Unit,
        name: str = None,
        redshift: float = None,
        reference: list[str] = None,
        **kwargs,
    ) -> None:
        """
        Object to store DataFinder info to query public data sources of host galaxies

        Args:
            ra (str|float) : The RA of the host to be passed to an astropy SkyCoord
            dec (str|float) : The declination of the host to be passed to an
                                   astropy SkyCoord
            ra_units (str|astropy.units.Unit) : units of the RA, to be passed to
                                                     the unit keyword of SkyCoord
            dec_units (str|astropy.units.Unit) : units of the declination, to be
                                                      passed to the unit keyword of
                                                      SkyCoord
            name (str) : The name of the host galaxy
            redshift (float) : The redshift of the host galaxy
            reference (list[str]) : a list of bibcodes that found this to be the host
            kwargs : Just here so we can pass **Transient['host'] into this constructor
                     and any extraneous properties will be ignored.
        """
        self.coord = SkyCoord(ra, dec, unit=(ra_units, dec_units))
        self.name = name
        self.z = redshift
        self.redshift = redshift  # just here for ease of use
        self.bibcodes = reference

    def __repr__(self) -> str:
        """
        String representation of the DataFinder for printing
        """

        if self.name is None:
            print_name = "No Name DataFinder"
        else:
            print_name = self.name

        return f"{print_name} @ (RA, Dec)=({self.coord.ra},{self.coord.dec})"

    def __iter__(self) -> dict:
        """
        Provides an iterator for the properties of this DataFinder. Yields (key, value)
        """
        out = dict(
            host_ra=self.coord.ra.value,
            host_dec=self.coord.dec.value,
            host_ra_units="deg",
            host_dec_units="deg",
        )

        if self.name is not None:
            out["host_name"] = self.name

        if self.z is not None:
            out["host_redshift"] = self.z

        if self.bibcodes is not None:
            out["reference"] = self.bibcodes

        for k, v in out.items():
            yield (k, v)

    ###################################################################################
    ################### CONVENIENCE METHODS FOR QUERYING HOST METADATA ################
    ###################################################################################

    @staticmethod
    def _wrap_astroquery(module, *args, **kwargs):
        """
        Private convenience method that just standardizes how we call the query_region
        method in astroquery
        """
        return module.query_region(*args, **kwargs)

    def query_simbad(self, radius="5 arcsec", **kwargs):
        """
        Query SIMBAD through astroquery to provide any other "meta" information on this
        host that may not be stored in the OTTER

        Args:
            radius (str|astropy.quantity.Quantity) : search radius for astroquery
            **kwargs : any other arguments for astroquery.vizier.Vizier.query_region

        Returns:
            astropy Table of the simbad results.
        """
        from astroquery.simbad import Simbad

        return DataFinder._wrap_astroquery(Simbad, self.coord, radius=radius, **kwargs)

    def query_vizier(self, radius="5 arcsec", **kwargs):
        """
        Query the ViZier catalog for TIME-AVERAGED data from their major/large catalogs.

        ViZier Catalogs Queried:
           - 2MASS-PSC
           - 2MASX
           - AC2000.2
           - AKARI
           - ALLWISE
           - ASCC-2.5
           - B/DENIS
           - CMC14
           - Gaia-DR1
           - GALEX
           - GLIMPSE
           - GSC-ACT
           - GSC1.2
           - GSC2.2
           - GSC2.3
           - HIP
           - HIP2
           - IRAS
           - NOMAD1
           - NVSS
           - PanSTARRS-DR1
           - PGC
           - Planck-DR1
           - PPMX
           - PPMXL
           - SDSS-DR12
           - SDSS-DR7
           - SDSS-DR9
           - Tycho-2
           - UCAC2
           - UCAC3
           - UCAC4
           - UKIDSS
           - USNO-A2
           - USNO-B1
           - WISE

        Args:
            radius (str|astropy.quantity.Quantity) : search radius for astroquery
            **kwargs : any other arguments for astroquery.vizier.Vizier.query_region

        Returns:
            astropy TableList of the time-averaged photometry associated with this host.
        """
        from astroquery.vizier import Vizier

        return DataFinder._wrap_astroquery(
            Vizier, self.coord, radius=radius, catalog=VIZIER_LARGE_CATALOGS
        )

    ###################################################################################
    ######### CONVENIENCE METHODS FOR QUERYING HOST TIME SERIES PHOTOMETRY  ###########
    ###################################################################################

    def query_atlas(
        self, days_ago: int = 365, disc_date: float = None, clip_sigma: float = 2.0
    ) -> pd.DataFrame:
        """
        Query ATLAS forced photometry for photometry for this host

        Args:
            days_ago (int) : Number of days before the transients discovery date
                             (or today if no disc_date is given) to get ATLAS
                             forced photometry for.
            disc_date (float) : The discovery date of the transient in MJD.
            clip_sigma (float) : amount to sigma clip the ATLAS data by

        Return:
            pandas DataFrame of the ATLAS forced photometry for this host
        """
        base_url = "https://fallingstar-data.com/forcedphot"

        token = os.environ.get("ATLAS_API_TOKEN", None)
        if token is None:
            logger.warn(
                "Getting your token from ATLAS. Please add ATLAS_API_TOKEN to your \
                environment variables to avoid this!"
            )

            uname = os.environ.get("ATLAS_UNAME", default=None)
            pword = os.environ.get("ATLAS_PWORD", default=None)

            if uname is None and pword is None:
                raise MissingEnvVarError(["ATLAS_UNAME", "ATLAS_PWORD"], base_url)
            elif uname is None and pword is not None:
                raise MissingEnvVarError(["ATLAS_UNAME"], base_url)
            elif uname is not None and pword is None:
                raise MissingEnvVarError(["ATLAS_PWORD"], base_url)

            resp = requests.post(
                url=f"{base_url}/api-token-auth/",
                data={"username": uname, "password": pword},
            )

            token = resp.json()["token"]

        headers = {"Authorization": f"Token {token}", "Accept": "application/json"}

        # compute the query start
        if disc_date is None:
            t_queryend = Time.now().mjd
            logger.warn(
                "Since no transient name is given we are using today \
                as the query end!"
            )
        else:
            t_queryend = Time(disc_date, format="mjd").mjd

        t_querystart = t_queryend - days_ago

        # submit the query to the ATLAS forced photometry server
        task_url = None
        while not task_url:
            with requests.Session() as s:
                resp = s.post(
                    f"{base_url}/queue/",
                    headers=headers,
                    data={
                        "ra": self.coord.ra.value,
                        "dec": self.coord.ra.value,
                        "send_email": False,
                        "mjd_min": t_querystart,
                        "mjd_max": t_queryend,
                        "use_reduced": False,
                    },
                )
                if resp.status_code == 201:  # success
                    task_url = resp.json()["url"]
                    logger.info(f"The task URL is {task_url}")
                elif resp.status_code == 429:  # throttled
                    message = resp.json()["detail"]
                    logger.info(f"{resp.status_code} {message}")
                    t_sec = re.findall(r"available in (\d+) seconds", message)
                    t_min = re.findall(r"available in (\d+) minutes", message)
                    if t_sec:
                        waittime = int(t_sec[0])
                    elif t_min:
                        waittime = int(t_min[0]) * 60
                    else:
                        waittime = 10
                    logger.info(f"Waiting {waittime} seconds")
                    time.sleep(waittime)
                else:
                    raise Exception(f"ERROR {resp.status_code}\n{resp.text}")

        # Now wait for the result
        result_url = None
        taskstarted_printed = False
        while not result_url:
            with requests.Session() as s:
                resp = s.get(task_url, headers=headers)

                if resp.status_code == 200:  # HTTP OK
                    if resp.json()["finishtimestamp"]:
                        result_url = resp.json()["result_url"]
                        logger.info(
                            f"Task is complete with results available at {result_url}"
                        )
                    elif resp.json()["starttimestamp"]:
                        if not taskstarted_printed:
                            print(
                                f"Task is running (started at\
                                {resp.json()['starttimestamp']})"
                            )
                            taskstarted_printed = True
                        time.sleep(2)
                    else:
                        # print(f"Waiting for job to start (queued at {timestamp})")
                        time.sleep(4)
                else:
                    raise Exception(f"ERROR {resp.status_code}\n{resp.text}")

        # get and clean up the result
        with requests.Session() as s:
            textdata = s.get(result_url, headers=headers).text

        atlas_phot = DataFinder._atlas_stack(textdata, clipping_sigma=clip_sigma)

        return pd.DataFrame(atlas_phot)

    def query_ptf(self, radius: str | u.Quantity = "5 arcsec", **kwargs) -> Table:
        """
        Query the palomer transient facility's light curve catalog for this host

        Args:
            radius (str|astropy.quantity.Quantity) : search radius
            **kwargs : other optional arguments for astroquery's query_region

        Returns:
            An astropy Table of the resulting light curve
        """
        from astroquery.ipac.irsa import Irsa

        ptf_lc_catalog = "ptf_lightcurves"
        return DataFinder._wrap_astroquery(
            Irsa, self.coord, radius=radius, catalog=ptf_lc_catalog
        )

    def query_ztf(self, radius: float = 5):
        """
        Query ZTF photometry/forced photometry for photometry for this host

        Args:
            radius (float) : The search radius in arcseconds

        Returns:
            An astropy table of the time series data from the cone search in ZTF
        """

        base_url = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?"

        ra, dec = self.coord.ra.value, self.coord.dec.value
        search_radius_arcseconds = radius  # in arcseconds
        search_radius_degree = search_radius_arcseconds / 3600

        query_url = f"{base_url}POS=CIRCLE%20{ra}%20{dec}%20{search_radius_degree}"

        resp = urlopen(query_url)

        votab = parse_single_table(io.BytesIO(resp.read()))

        return Table(votab.array)

    def query_asassn(self, radius: float = 5.0, nthreads: int = 2) -> pd.DataFrame:
        """
        Query ASASSN photometry/forced photometry for photometry for this host

        Args:
            radius (float) : search radius in arcseconds
            nthreads (int) : number of threads to utilize during download, default is 2

        Returns:
            A pandas dataframe with the ASASSN lightcurve for this object
        """
        from pyasassn.client import SkyPatrolClient

        client = SkyPatrolClient()
        light_curve = client.cone_search(
            self.coord.ra.value,
            self.coord.dec.value,
            radius=radius,
            units="arcsec",
            download=True,
            threads=nthreads,
        )
        return light_curve.data

    def query_wise(
        self,
        radius: float = 5,
        datadir: str = "ipac/",
        overwrite: bool = False,
        verbose=False,
        **kwargs,
    ) -> pd.DataFrame:
        """
        Query NEOWISE for their multiepoch photometry

        The method used to query wise here was taken from this github repo:
        https://github.com/HC-Hwang/wise_light_curves/tree/master
        and you should cite this other paper that the authors of this code developed
        it for: https://ui.adsabs.harvard.edu/abs/2020MNRAS.493.2271H/abstract

        This will download the ipac data files to the "datadir" argument. by default,
        these will go into os.getcwd()/ipac

        Args:
            radius (float) : The cone search radius in arcseconds
            overwrite (bool) : Overwrite the existing datasets downloaded from wise
            **kwargs : Other optional arguments for the astroquery query_region
        Returns:
            An astropy Table of the multiepoch wise data for this host
        """
        # from https://www.cambridge.org/core/journals/
        # publications-of-the-astronomical-society-of-australia/article/
        # recalibrating-the-widefield-infrared-survey-explorer-wise-w4-filter/
        # B238BFFE19A533A2D2638FE88CCC2E89
        band_vals = {"w1": 3.4, "w2": 4.6, "w3": 12, "w4": 22}  # in um

        ra, dec = self.coord.ra.value, self.coord.dec.value

        fbasename = f"wise_{self.name}"
        allwise_name = f"{fbasename}_allwise.ipac"
        neowise_name = f"{fbasename}_neowise.ipac"

        if not os.path.exists(datadir):
            os.makedirs(datadir)

        self._download_single_data(
            name=fbasename,
            ra=ra,
            dec=dec,
            root_path=datadir,
            radius=radius,
            overwrite=overwrite,
        )

        allwise = ascii.read(f"ipac/{allwise_name}", format="ipac")
        neowise = ascii.read(f"ipac/{neowise_name}", format="ipac")

        allwise, neowise = self._only_good_data(allwise, neowise, verbose=verbose)
        if verbose and (allwise is None or neowise is None):
            print(f"Limited good infrared data for {self.name}, skipping!")

        mjd, mag, mag_err, filts = self._make_full_lightcurve_multibands(
            allwise, neowise, bands=["w1", "w2", "w3", "w4"]
        )

        df = pd.DataFrame(
            dict(
                name=[self.name] * len(mjd),
                date_mjd=mjd,
                filter=filts,
                filter_eff=[band_vals[f] for f in filts],
                filter_eff_unit=["um"] * len(mjd),
                flux=mag,
                flux_err=mag_err,
                flux_unit=["mag(AB)"] * len(mjd),
                upperlimit=[False] * len(mjd),
            )
        )

        # clean up the wise data by filtering out negative flux
        wise = df[df.flux > 0].reset_index(drop=True)
        return wise

    def query_alma(self, radius: float = 5, **kwargs) -> Table:
        """
        Query ALMA to see if there are observations of this host.

        NOTE: Since this is radio/mm data, it is unlikely that the output table will
        simply have fluxes in it. Instead you will need to use the access_url column
        to download and reduce this data.

        Args:
            radius (float) : The cone search radius in arcseconds
            **kwargs : Other optional arguments for the astroquery query_region
        Returns:
            An astropy Table of the multiepoch wise data for this host
        """

        logger.warn(
            "This method may not work if you are using a conda environment!\
            This is a known issue in setuptools that is not resolved!"
        )

        from astroquery.alma import Alma

        res = DataFinder._wrap_astroquery(
            Alma, self.coord, radius=5 * u.arcsec, **kwargs
        )
        return res

    def query_first(
        self, radius: u.Quantity = 5 * u.arcmin, get_image: bool = False, **kwargs
    ) -> list:
        """
        Query the FIRST radio survey and return an astropy table of the flux density

        This queries Table 2 from Ofek & Frail (2011); 2011ApJ...737...45O

        Args:
            radius (u.Quantity) : An astropy Quantity with the image height/width
            get_image (bool) : If True, download and return a list of the associated
                               images too.
            **kwargs : any other arguments to pass to the astroquery.image_cutouts
                       get_images method

        Returns:
            Astropy table of the flux densities. If get_image is True, it also returns
            a list of FIRST radio survey images
        """
        from astroquery.vizier import Vizier

        res = DataFinder._wrap_astroquery(
            Vizier, self.coord, radius=radius, catalog="J/ApJ/737/45/table2"
        )

        if get_image:
            from astroquery.image_cutouts.first import First

            res_img = First.get_images(self.coord, image_size=radius, **kwargs)
            return res, res_img

        return res

    def query_nvss(self, radius: u.Quantity = 5 * u.arcsec, **kwargs) -> Table:
        """
        Query the NRAO VLA Sky Survey (NVSS) and return a table list of the
        result

        This queries Table 1 from Ofek & Frail (2011); 2011ApJ...737...45O

        Args:
            radius (u.Quantity) : An astropy Quantity with the radius
            **kwargs : Any other arguments to pass to query_region
        """
        from astroquery.vizier import Vizier

        res = DataFinder._wrap_astroquery(
            Vizier, self.coord, radius=radius, catalog="J/ApJ/737/45/table1"
        )
        return res

    def query_heasarc(self, radius: u.Quantity = 5 * u.arcsec, **kwargs) -> Table:
        """
        Query Heasarc by the argument "heasarc_key" for the ra/dec associated with this
        DataLoader object.

        Args:
            radius (u.Quantity) : An astropy Quantity with the radius
            heasarc_table (str) : String with name of heasarc table to query. Default is
                                'xray' which queries the heasarc master x-ray catalog,
                                'radio' will query the heasarc master radio catalog. See
                                https://heasarc.gsfc.nasa.gov/cgi-bin/W3Browse/w3catindex.pl
                                for a complete list.
            **kwargs : Any other arguments to pass to query_region

        Returns:
            Astropy table of the rows in `heasarc_table` that match self.coord.
        """
        from astroquery.heasarc import Heasarc

        res = DataFinder._wrap_astroquery(Heasarc, self.coord, radius=radius, **kwargs)

        return res

    ###################################################################################
    ######### CONVENIENCE METHODS FOR QUERYING HOST SPECTR  ###########################
    ###################################################################################

    def query_sparcl(
        self, radius: u.Quantity = 5 * u.arcsec, include: str | list = "DEFAULT"
    ) -> Table:
        """
        Query the NOIRLab DataLabs Sparcl database for spectra for this host

        Args:
            radius (Quantity) : search radius as an Astropy.unit.Quantity
            include [list|str] : list or string of columns to include in the result. See
                                 the sparcl client documentation for more info. The
                                 default returns specid, ra, dec, sparcl_id, flux,
                                 wavelength, and the spectroscopic surveyu (_dr)

        Returns:
            astropy Table of the results, one row per spectrum
        """

        from sparcl.client import SparclClient
        from dl import queryClient as qc  # noqa: N813

        client = SparclClient()

        # first do a cone search on sparcl.main
        ra, dec = self.coord.ra.value, self.coord.dec.value
        radius_deg = radius.to(u.deg).value

        adql = f"""
        SELECT *
        FROM sparcl.main
        WHERE 't'=Q3C_RADIAL_QUERY(ra,dec,{ra},{dec},{radius_deg})
        """
        cone_search_res = qc.query(adql=adql, fmt="pandas")

        # then retrieve all of the spectra corresponding to those sparcl_ids
        sparcl_ids = cone_search_res.sparcl_id.tolist()
        res = client.retrieve(uuid_list=sparcl_ids, include=include)
        all_spec = pd.concat([pd.DataFrame([record]) for record in res.records])
        return Table.from_pandas(all_spec)

    ###################################################################################
    ######### PRIVATE HELPER METHODS FOR THE QUERYING #################################
    ###################################################################################
    @staticmethod
    def _atlas_stack(filecontent, clipping_sigma, log=logger):
        """
        Function adapted from David Young's :func:`plotter.plot_single_result`
        https://github.com/thespacedoctor/plot-results-from-atlas-force-photometry-service/blob/main/plot_atlas_fp.py

        And again adapted from https://github.com/SAGUARO-MMA/kne-cand-vetting/blob/master/kne_cand_vetting/survey_phot.py
        """
        epochs = DataFinder._atlas_read_and_sigma_clip_data(
            filecontent, log=log, clipping_sigma=clipping_sigma
        )

        # c = cyan, o = arange
        magnitudes = {
            "c": {"mjds": [], "mags": [], "magErrs": [], "lim5sig": []},
            "o": {"mjds": [], "mags": [], "magErrs": [], "lim5sig": []},
            "I": {"mjds": [], "mags": [], "magErrs": [], "lim5sig": []},
        }

        # SPLIT BY FILTER
        for epoch in epochs:
            if epoch["F"] in ["c", "o", "I"]:
                magnitudes[epoch["F"]]["mjds"].append(epoch["MJD"])
                magnitudes[epoch["F"]]["mags"].append(epoch["uJy"])
                magnitudes[epoch["F"]]["magErrs"].append(epoch["duJy"])
                magnitudes[epoch["F"]]["lim5sig"].append(epoch["mag5sig"])

        # STACK PHOTOMETRY IF REQUIRED
        stacked_magnitudes = DataFinder._stack_photometry(magnitudes, binningdays=1)

        return stacked_magnitudes

    @staticmethod
    def _atlas_read_and_sigma_clip_data(filecontent, log, clipping_sigma=2.2):
        """
        Function adapted from David Young's :func:`plotter.read_and_sigma_clip_data`
        https://github.com/thespacedoctor/plot-results-from-atlas-force-photometry-service/blob/main/plot_atlas_fp.py

        And again adapted from
        https://github.com/SAGUARO-MMA/kne-cand-vetting/blob/master/kne_cand_vetting/survey_phot.py

        *clean up rouge data from the files by performing some basic clipping*
        **Key Arguments:**
        - `fpFile` -- path to single force photometry file
        - `clippingSigma` -- the level at which to clip flux data
        **Return:**
        - `epochs` -- sigma clipped and cleaned epoch data
        """

        # CLEAN UP FILE FOR EASIER READING
        fpdata = (
            filecontent.replace("###", "")
            .replace(" ", ",")
            .replace(",,", ",")
            .replace(",,", ",")
            .replace(",,", ",")
            .replace(",,", ",")
            .splitlines()
        )

        # PARSE DATA WITH SOME FIXED CLIPPING
        oepochs = []
        cepochs = []
        csvreader = csv.DictReader(
            fpdata, dialect="excel", delimiter=",", quotechar='"'
        )

        for row in csvreader:
            for k, v in row.items():
                try:
                    row[k] = float(v)
                except Exception:
                    pass
            # REMOVE VERY HIGH ERROR DATA POINTS, POOR CHI SQUARED, OR POOR EPOCHS
            if row["duJy"] > 4000 or row["chi/N"] > 100 or row["mag5sig"] < 17.0:
                continue
            if row["F"] == "c":
                cepochs.append(row)
            if row["F"] == "o":
                oepochs.append(row)

        # SORT BY MJD
        cepochs = sorted(cepochs, key=itemgetter("MJD"), reverse=False)
        oepochs = sorted(oepochs, key=itemgetter("MJD"), reverse=False)

        # SIGMA-CLIP THE DATA WITH A ROLLING WINDOW
        cdataflux = []
        cdataflux[:] = [row["uJy"] for row in cepochs]
        odataflux = []
        odataflux[:] = [row["uJy"] for row in oepochs]

        masklist = []
        for flux in [cdataflux, odataflux]:
            fullmask = rolling_window_sigma_clip(
                log=log, array=flux, clippingSigma=clipping_sigma, windowSize=11
            )
            masklist.append(fullmask)

        try:
            cepochs = [e for e, m in zip(cepochs, masklist[0]) if m == False]
        except Exception:
            cepochs = []

        try:
            oepochs = [e for e, m in zip(oepochs, masklist[1]) if m == False]
        except Exception:
            oepochs = []

        logger.info("Completed the ``read_and_sigma_clip_data`` function")
        # Returns ordered dictionary of all parameters
        return cepochs + oepochs

    @staticmethod
    def _stack_photometry(magnitudes, binningdays=1.0):
        """
        Function adapted from David Young's :func:`plotter.stack_photometry`
        https://github.com/thespacedoctor/plot-results-from-atlas-force-photometry-service/blob/main/plot_atlas_fp.py

        And again adapted from
        https://github.com/SAGUARO-MMA/kne-cand-vetting/blob/master/kne_cand_vetting/survey_phot.py

        *stack the photometry for the given temporal range*
        **Key Arguments:**
            - `magnitudes` -- dictionary of photometry divided into filter sets
            - `binningDays` -- the binning to use (in days)
        **Return:**
            - `summedMagnitudes` -- the stacked photometry
        """

        # IF WE WANT TO 'STACK' THE PHOTOMETRY
        summed_magnitudes = {
            "c": {"mjds": [], "mags": [], "magErrs": [], "n": [], "lim5sig": []},
            "o": {"mjds": [], "mags": [], "magErrs": [], "n": [], "lim5sig": []},
            "I": {"mjds": [], "mags": [], "magErrs": [], "n": [], "lim5sig": []},
        }

        # MAGNITUDES/FLUXES ARE DIVIDED IN UNIQUE FILTER SETS - SO ITERATE OVER
        # FILTERS
        alldata = []
        for fil, data in list(magnitudes.items()):
            # WE'RE GOING TO CREATE FURTHER SUBSETS FOR EACH UNQIUE MJD
            # (FLOORED TO AN INTEGER)
            # MAG VARIABLE == FLUX (JUST TO CONFUSE YOU)
            distinctmjds = {}
            for mjd, flx, err, lim in zip(
                data["mjds"], data["mags"], data["magErrs"], data["lim5sig"]
            ):
                # DICT KEY IS THE UNIQUE INTEGER MJD
                key = str(int(math.floor(mjd / float(binningdays))))
                # FIRST DATA POINT OF THE NIGHTS? CREATE NEW DATA SET
                if key not in distinctmjds:
                    distinctmjds[key] = {
                        "mjds": [mjd],
                        "mags": [flx],
                        "magErrs": [err],
                        "lim5sig": [lim],
                    }
                # OR NOT THE FIRST? APPEND TO ALREADY CREATED LIST
                else:
                    distinctmjds[key]["mjds"].append(mjd)
                    distinctmjds[key]["mags"].append(flx)
                    distinctmjds[key]["magErrs"].append(err)
                    distinctmjds[key]["lim5sig"].append(lim)

            # ALL DATA NOW IN MJD SUBSETS. SO FOR EACH SUBSET (I.E. INDIVIDUAL
            # NIGHTS) ...
            for k, v in list(distinctmjds.items()):
                # GIVE ME THE MEAN MJD
                meanmjd = sum(v["mjds"]) / len(v["mjds"])
                summed_magnitudes[fil]["mjds"].append(meanmjd)
                # GIVE ME THE MEAN FLUX
                meanflux = sum(v["mags"]) / len(v["mags"])
                summed_magnitudes[fil]["mags"].append(meanflux)
                # GIVE ME THE COMBINED ERROR
                sum_of_squares = sum(x**2 for x in v["magErrs"])
                comberror = math.sqrt(sum_of_squares) / len(v["magErrs"])
                summed_magnitudes[fil]["magErrs"].append(comberror)
                # 5-sigma limits
                comb5siglimit = 23.9 - 2.5 * math.log10(5.0 * comberror)
                summed_magnitudes[fil]["lim5sig"].append(comb5siglimit)
                # GIVE ME NUMBER OF DATA POINTS COMBINED
                n = len(v["mjds"])
                summed_magnitudes[fil]["n"].append(n)
                alldata.append(
                    {
                        "mjd": meanmjd,
                        "uJy": meanflux,
                        "duJy": comberror,
                        "F": fil,
                        "n": n,
                        "mag5sig": comb5siglimit,
                    }
                )
        print("completed the ``stack_photometry`` method")

        return alldata

    """
    The following code was taken and modified for the purposes of this package from
    https://github.com/HC-Hwang/wise_light_curves/blob/master/wise_light_curves.py

    Original Authors:
    - Matthew Hill
    - Hsiang-Chih Hwang

    Update Author:
    - Noah Franz
    """

    @staticmethod
    def _get_by_position(ra, dec, radius=2.5):
        allwise_cat = "allwise_p3as_mep"
        neowise_cat = "neowiser_p1bs_psd"
        query_url = "http://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query"
        payload = {
            "catalog": allwise_cat,
            "spatial": "cone",
            "objstr": " ".join([str(ra), str(dec)]),
            "radius": str(radius),
            "radunits": "arcsec",
            "outfmt": "1",
        }
        r = requests.get(query_url, params=payload)
        allwise = ascii.read(r.text)
        payload = {
            "catalog": neowise_cat,
            "spatial": "cone",
            "objstr": " ".join([str(ra), str(dec)]),
            "radius": str(radius),
            "radunits": "arcsec",
            "outfmt": "1",
            "selcols": "ra,dec,sigra,sigdec,sigradec,glon,glat,elon,elat,w1mpro,w1sigmpro,w1snr,w1rchi2,w2mpro,w2sigmpro,w2snr,w2rchi2,rchi2,nb,na,w1sat,w2sat,satnum,cc_flags,det_bit,ph_qual,sso_flg,qual_frame,qi_fact,saa_sep,moon_masked,w1frtr,w2frtr,mjd,allwise_cntr,r_allwise,pa_allwise,n_allwise,w1mpro_allwise,w1sigmpro_allwise,w2mpro_allwise,w2sigmpro_allwise,w3mpro_allwise,w3sigmpro_allwise,w4mpro_allwise,w4sigmpro_allwise",  # noqa: E501
        }
        r = requests.get(query_url, params=payload)

        neowise = ascii.read(r.text, guess=False, format="ipac")

        return allwise, neowise

    @staticmethod
    def _download_single_data(
        name, ra, dec, root_path="ipac/", radius=2.5, overwrite=False
    ):
        # ra, dec: in degree
        # name, ra, dec = row['Name'], row['RAJ2000'], row['DEJ2000']
        # name = 'J' + ra + dec
        if root_path[-1] != "/":
            root_path += "/"
        if (
            not overwrite
            and os.path.isfile(root_path + name + "_allwise.ipac")
            and os.path.isfile(root_path + name + "_neowise.ipac")
        ):
            pass
        else:
            allwise, neowise = DataFinder._get_by_position(ra, dec, radius=radius)
            allwise.write(
                root_path + name + "_allwise.ipac", format="ascii.ipac", overwrite=True
            )
            neowise.write(
                root_path + name + "_neowise.ipac", format="ascii.ipac", overwrite=True
            )

    @staticmethod
    def _get_data_arrays(table, t, mag, magerr):
        """Get the time series from a potentially masked astropy table"""
        if table.masked:
            full_mask = table[t].mask | table[mag].mask | table[magerr].mask
            t = table[t].data
            mag = table[mag].data
            magerr = table[magerr].data

            t.mask = full_mask
            mag.mask = full_mask
            magerr.mask = full_mask

            return t.compressed(), mag.compressed(), magerr.compressed()

        else:
            return table[t].data, table[mag].data, table[magerr].data

    @staticmethod
    def _make_full_lightcurve(allwise, neowise, band):
        """band = 'w1', 'w2', 'w3', or 'w4'"""
        """Get a combined AllWISE and NEOWISE lightcurve from their Astropy tables"""

        if band not in ["w1", "w2", "w3", "w4"]:
            raise ValueError("band can only be w1, w2, w3, or w4")

        use_neowise = band in {"w1", "w2"}
        use_allwise = allwise is not None

        if use_neowise and use_allwise:
            t, m, e = DataFinder._get_data_arrays(
                allwise, "mjd", band + "mpro_ep", band + "sigmpro_ep"
            )
            t_n, m_n, e_n = DataFinder._get_data_arrays(
                neowise, "mjd", band + "mpro", band + "sigmpro"
            )
            t, m, e = (
                np.concatenate((t, t_n)),
                np.concatenate((m, m_n)),
                np.concatenate((e, e_n)),
            )

        elif use_neowise and not use_allwise:
            t, m, e = DataFinder._get_data_arrays(
                neowise, "mjd", band + "mpro", band + "sigmpro"
            )

        elif not use_neowise and use_allwise:
            t, m, e = DataFinder._get_data_arrays(
                allwise, "mjd", band + "mpro_ep", band + "sigmpro_ep"
            )

        else:
            raise Exception("No good allwise or neowise data!")

        t_index = t.argsort()
        t, m, e = map(lambda e: e[t_index], [t, m, e])

        return t, m, e

    @staticmethod
    def _make_full_lightcurve_multibands(allwise, neowise, bands=["w1", "w2"]):
        t, m, e = DataFinder._make_full_lightcurve(allwise, neowise, bands[0])
        filts = [bands[0] for i in range(len(t))]
        for band in bands[1:]:
            try:
                t_tmp, m_tmp, e_tmp = DataFinder._make_full_lightcurve(
                    allwise, neowise, band
                )
            except Exception:
                continue
            t = np.concatenate((t, t_tmp))
            m = np.concatenate((m, m_tmp))
            e = np.concatenate((e, e_tmp))
            filts += [band for i in range(len(t_tmp))]
        return t, m, e, np.array(filts)

    @staticmethod
    def _cntr_to_source_id(cntr):
        cntr = str(cntr)

        # fill leanding 0s
        if len(cntr) < 19:
            num_leading_zeros = 19 - len(cntr)
            cntr = "0" * num_leading_zeros + cntr

        pm = "p"
        if cntr[4] == "0":
            pm = "m"

        t = chr(96 + int(cntr[8:10]))

        return "%s%s%s_%cc%s-%s" % (
            cntr[0:4],
            pm,
            cntr[5:8],
            t,
            cntr[11:13],
            cntr[13:19],
        )

    @staticmethod
    def _only_good_data(allwise, neowise, verbose=False):
        """
        Select good-quality data. The criteria include:
        - matching the all-wise ID

        To be done:
        - deal with multiple cntr

        This filtering is described here:
        https://wise2.ipac.caltech.edu/docs/release/neowise/expsup/sec2_3.html
        """

        neowise_prefilter_n = len(neowise)
        neowise = neowise[
            (neowise["qual_frame"] > 0.0)
            * (neowise["qi_fact"] > 0.9)
            * (neowise["saa_sep"] > 0)
            * (neowise["moon_masked"] == "00")
        ]
        neowise_postfilter_n = len(neowise)
        if verbose:
            print(
                f"Filtered out {neowise_prefilter_n-neowise_postfilter_n} neowise \
                points, leaving {neowise_postfilter_n}"
            )

        cntr_list = []
        for data in neowise:
            if data["allwise_cntr"] not in cntr_list and data["allwise_cntr"] > 10.0:
                cntr_list.append(data["allwise_cntr"])

        if len(cntr_list) >= 2:
            print("multiple cntr:")
            print(cntr_list)
            return None, neowise

        if len(cntr_list) == 0:
            # import pdb; pdb.set_trace()
            # raise Exception('No center!')
            return None, neowise

        cntr = cntr_list[0]

        source_id = DataFinder._cntr_to_source_id(cntr)

        allwise_prefilter_n = len(allwise)
        allwise = allwise[
            (allwise["source_id_mf"] == source_id)
            * (allwise["saa_sep"] > 0.0)
            * (allwise["moon_masked"] == "0000")
            * (allwise["qi_fact"] > 0.9)
        ]
        allwise_postfilter_n = len(neowise)
        if verbose:
            print(
                f"Filtered out {allwise_prefilter_n-allwise_postfilter_n} allwise \
                points, leaving {allwise_postfilter_n}"
            )

        return allwise, neowise
