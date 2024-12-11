"""
Host object that stores information on the Transient Host and provides utility methods
for pulling in data corresponding to that host
"""

from __future__ import annotations

from urllib.request import urlopen
import json

import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u

from .data_finder import DataFinder
from ..exceptions import OtterLimitationError

import logging

logger = logging.getLogger(__name__)


class Host(DataFinder):
    def __init__(
        self,
        host_ra: str | float,
        host_dec: str | float,
        host_ra_units: str | u.Unit,
        host_dec_units: str | u.Unit,
        host_name: str = None,
        host_redshift: float = None,
        redshift_type: str = None,
        reference: list[str] = None,
        transient_name: str = None,
        **kwargs,
    ) -> None:
        """
        Object to store host information and query public data sources of host galaxies

        Subclass of the data scraper class to allow for these queries to happen

        Args:
            host_ra (str|float) : The RA of the host to be passed to an astropy SkyCoord
            host_dec (str|float) : The declination of the host to be passed to an
                                   astropy SkyCoord
            host_ra_units (str|astropy.units.Unit) : units of the RA, to be passed to
                                                     the unit keyword of SkyCoord
            host_dec_units (str|astropy.units.Unit) : units of the declination, to be
                                                      passed to the unit keyword of
                                                      SkyCoord
            host_name (str) : The name of the host galaxy
            host_redshift (float) : The redshift of the host galaxy
            redshift_type (str) : Either "phot" or "spec", tells you if the redshift
                                  is a phot-z or spec-z.
            reference (list[str]) : a list of bibcodes that found this to be the host
            transient_name (str) : the name of the transient associated with this host
            kwargs : Just here so we can pass **Transient['host'] into this constructor
                     and any extraneous properties will be ignored.
        """
        self.coord = SkyCoord(host_ra, host_dec, unit=(host_ra_units, host_dec_units))
        self.name = host_name
        self.z = host_redshift
        self.redshift = host_redshift  # just here for ease of use
        self.bibcodes = reference
        self.transient_name = transient_name
        self.redshift_type = redshift_type

    def pcc(self, transient_coord: SkyCoord, mag: float = None):
        """
        Compute the Probability of Chance Coincindence as described in
        Bloom et al. (2002) "Offset Distribution of Gamma-Ray Bursts.

        This computes the probability that this galaxy is by chance nearby to the
        transient on the sky. Or, in simpler terms this essentially computes the
        probability that we are wrong about this being the transient host. So, a
        smaller probability is better!

        Note: This probability was initially defined for GRB afterglows, which tend to
        be redder transients (supernova too). So, be cautious when using this algorithm
        for TDEs!

        Args:
            transient_coord (astropy.coordinates.SkyCoord) : The coordinates of the
                                                             transient object.
            mag (float) : An r-band magnitude to compute from. Default is None which
                          will prompt us to check SDSS for one within 10".
        Returns:
            A float probability in the range [0,1]
        """

        # first get the separation r, in arcseconds
        r = self.coord.separation(transient_coord).arcsec

        # next get the host r magnitude
        if mag is None:
            res = self.query_vizier(radius=10 * u.arcsec)

            if len(res) == 0:
                raise OtterLimitationError(
                    "No magnitude found in SDSS! Please provide a magnitude via the \
                    `mag` keyword to make this calculation!"
                )

            sdss = [k for k in res.keys() if "sdss" in k]
            use = max(sdss, key=lambda k: int(k.split("sdss")[-1]))
            print(f"Using the r magnitude from the {use} table")
            mag = res[use]["rmag"][0]

        # then compute the probability
        sigma_prefactor = 1 / (3600**2 * 0.334 * np.log(10))
        sigma_pow = 0.334 * (mag - 22.963) + 4.320
        sigma = sigma_prefactor * 10**sigma_pow

        eta = np.pi * r**2 * sigma

        prob = 1 - np.exp(-eta)

        return prob

    ###################################################################################
    ######### METHODS FOR FINDING HOSTS  ###########################
    ###################################################################################
    @staticmethod
    def query_blast(tns_name: str) -> dict:
        """
        Query the BLAST host galaxy service

        Args:
            tns_name (str) : The TNS target name to grab the data from BLAST for

        Returns:
            best match BLAST ID, host ra, host dec, host redshift. Redshift will be
            spectroscopic if available, otherwise photometric.
        """

        # clean up the input name a little
        if tns_name[:2] == "AT":
            tns_name = tns_name.replace("AT", "")

        failed_query_res = None

        # do the query
        blast_base_url = "https://blast.scimma.org"
        blast_query_url = f"{blast_base_url}/api/transient/?name={tns_name}&format=json"
        with urlopen(blast_query_url) as response:
            if response.status != 200:
                logger.warn(f"BLAST query failed with response code {response.status}!")
                return failed_query_res
            else:
                blast_data = json.loads(response.read())
        if len(blast_data) == 0:
            logger.warn("BLAST query returned no results!")
            return failed_query_res

        blast_data = blast_data[0]
        if "host" not in blast_data:
            logger.warn("BLAST query found the object but it has no host associated!")
            return failed_query_res

        blast_host = blast_data["host"]

        if blast_host["redshift"] is not None:
            # prefer spec-z over phot-z
            z = blast_host["redshift"]
            z_type = "spec"
        else:
            # well I guess we need to use phot-z
            z = blast_host["photometric_redshift"]
            z_type = "phot"

        refs = [
            "2021ApJ...908..170G",  # GHOST citation
            "2024arXiv241017322J",  # BLAST citation
        ]

        return Host(
            host_ra=blast_host["ra_deg"],
            host_dec=blast_host["dec_deg"],
            host_ra_units="deg",
            host_dec_units="deg",
            host_name=blast_host["id"],
            host_redshift=z,
            redshift_type=z_type,
            reference=refs,
            transient_name=tns_name,
        )
