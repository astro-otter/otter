'''
This is the primary class for user interaction with the catalog
'''

from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy import units as u

from pyArango.database import Database
from pyArango.connection import Connection

from .transient import Transient

class Otter(Database):
    '''
    This is the primary class for users to access the otter backend database

    Args:
        username [str]: Your connection username to the database, default is the user
                        login which only has read permission.
        password [str]: Your password corresponding to your username.
        db [str]: The database name to connect to. This is default to 'otter' which is
                  the only database so far.
        collection [str]: The collection to read data from. Right now the only 
                          collection is 'tdes'.
        debug [bool]: debug mode, set to true to limit reading from database.
    '''
    
    def __init__(self,
                 username:str='user@otter',
                 password:str='insecure',
                 db:str='otter',
                 collection:str='tdes',
                 debug:bool=False
                 ) -> None:

        # save inputs
        self.dbName = db
        self.collectionName = collection
        
        c = Connection(username=username, password=password)
        
        # initiate the tdes database
        super().__init__(c, db)

    def getMeta(self, **kwargs) -> Table:
        '''
        Get the metadata of the objects matching the arguments

        Args:
            **kwargs : Arguments to pass to Otter.query()
        Return:
           The metadata for the transients that match the arguments. Will be an astropy
           Table by default, if raw=True will be a dictionary.
        '''
        metakeys = ['name', 'coordinate', 'epoch', 'distance', 'classification']
        
        return [t[metakeys] for t in self.query(**kwargs)]     
        
        
    def coneSearch(self,
                   coords:SkyCoord,
                   radius:float=0.05,
                   raw:bool=False
                   ) -> Table:
        '''
        Performs a cone search of the catalog over the given coords and radius.

        Args:
            coords [SkyCoord]: An astropy SkyCoord object with coordinates to match to
            radius [float]: The radius of the cone in arcseconds, default is 0.05"
            raw [bool]: If False (the default) return an astropy table of the metadata
                        for matching objects. Otherwise, return the raw json dicts
        
        Return:
            The metadata for the transients in coords+radius. Will return an astropy 
            Table if raw is False, otherwise a dict.
        '''

        transients = self.query(coords=coords,
                                radius=radius,
                                raw=False)
        
        
    def query(self,
              names:list[str]=None,
              coords:SkyCoord=None,
              radius:float=0.05,
              minZ:float=0,
              maxZ:float=None,
              refs:list[str]=None,
              hasPhot:bool=False,
              hasSpec:bool=False,
              raw:bool=False
              ) -> dict:
        '''
        Wraps on the super.AQLQuery

        This is how it differs from the `getMeta` method. Users should prefer to use
        `getMeta`, `getPhot`, and `getSpec` independently because it is a better 
        workflow and can return the data in an astropy table.

        Args:
            names [list[str]]: A list of names to get the metadata for
            coords [SkyCoord]: An astropy SkyCoord object with coordinates to match to
            radius [float]: The radius in arcseconds for a cone search, default is 0.05"
            minZ [float]: The minimum redshift to search for
            maxZ [float]: The maximum redshift to search for
            refs [list[str]]: A list of ads bibcodes to match to. Will only return 
                              metadata for transients that have this as a reference.
            hasPhot [bool]: if True, only returns transients which have photometry.
            hasSpec [bool]: if True, only return transients that have spectra.

        Return:
           Get all of the raw json data for objects that match the criteria. 
        '''

        # write some AQL filters based on the inputs
        queryFilters = ''

        if hasPhot is True:
            queryFilters += "FILTER 'photometry' IN ATTRIBUTES(tde)\n"

        if hasSpec is True:
            queryFilters += "FILTER 'spectra' IN ATTRIBUTES(tde)\n"
        
        if minZ > 0:
            sfilt = f'''
            FILTER TO_NUMBER(tde.distance.redshift[*].value) >= {minZ} \n
            '''
            queryFilters += sfilt
        if maxZ is not None:
            sfilt = f'''
            FILTER TO_NUMBER(tde.distance.redshift[*].value) <= {maxZ} \n
            '''
            queryFilters += sfilt 
            
        if names is not None:
            if isinstance(names, str):
                queryFilters += f"FILTER tde.name LIKE '%{names}%'\n"
            elif isinstance(names, list):
                namefilt = f'''
            FOR name IN {names} 
                FILTER name IN tde.name.alias[*].value\n
                '''
                queryFilters += namefilt
            else:
                raise Exception('Names must be either a string or list')

        if refs is not None:
            if isinstance(refs, str): # this is just a single bibcode
                queryFilters += f"FILTER {refs} IN tde.reference_alias[*].name"
            elif isinstance(refs, list):
                queryFilters += f'''
                FOR ref IN {refs}
                    FILTER ref IN tde.reference_alias[*].name
                '''
            else:
                raise Exception('reference list must be either a string or a list')
            
        # define the query
        query = f'''
        FOR tde IN tdes
            {queryFilters}
            RETURN tde
        '''
        
        result = self.AQLQuery(query, rawResults=True)
        
        # now that we have the query results do the RA and Dec queries if they exist
        if coords is not None:
            # get the catalog RAs and Decs to compare against
            queryCoords = coords
            goodTDEs = []

            for tde in result:
                for coordinfo in tde['coordinate']['equitorial']:
                    coord = SkyCoord(coordinfo['ra'], coordinfo['dec'],
                                     unit=(coordinfo['ra_units'], coordinfo['dec_units']))
                    if queryCoords.separation(coord) < radius*u.arcsec:
                        goodTDEs.append(tde)
                        break # we've confirmed this tde is in the cone!

            if not raw:
                return [Transient(t) for t in goodTDEs]
            else:
                return goodTDEs

        else:
            if not raw:
                return [Transient(res) for res in result.result]
            else:
                return result.result
            
    def _close(self) -> None:
        '''
        Close the database connection. We just need this for flask!
        '''
        del self
