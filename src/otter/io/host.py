"""
Host object that stores information on the Transient Host and provides utility methods
for pulling in data corresponding to that host
"""

from .scraper import Scraper
from astropy.coordinates import SkyCoord
from astropy import units as u


class Host(Scraper):
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
