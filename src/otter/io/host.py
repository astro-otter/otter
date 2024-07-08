"""
Host object that stores information on the Transient Host and provides utility methods
for pulling in data corresponding to that host
"""

from __future__ import annotations
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u

from .data_finder import DataFinder
from ..exceptions import OtterLimitationError


class Host(DataFinder):
    def __init__(
        self,
        host_ra: str | float,
        host_dec: str | float,
        host_ra_units: str | u.Unit,
        host_dec_units: str | u.Unit,
        host_name: str = None,
        host_redshift: float = None,
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
