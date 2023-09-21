'''
Class for to make queries in the catalog
'''
import warnings
from pyArango.database import Database
from pyArango.connection import Connection
from .tde import TDE

import astropy.units as u
from astropy.coordinates import SkyCoord

class TDECatalog(Database):

    def __init__(self,
                 username='user@tide',
                 password='password',
                 db='tide',
                 collection='tdes'):

        c = Connection(username=username, password=password)

        # initiate the tdes database
        super().__init__(c, db)

        # for now, just get all data
        # THIS IS NOT SCALABLE, FIX LATER
        self.rawData = self[collection].fetchAll(rawResults=True)

        # clean up the raw data for easier parsing and return a list of TDEs
        self.tdes = self._clean()
        
    def _clean(self, dataDict=None) -> list[TDE]:
        '''
        Get the data from the database and clean it 
    
        Args:
            tdes [list[dict]]: output of connect
        
        Returns: names, sources, ra, dec, z
        '''

        if dataDict is not None:
            data = dataDict
        else:
            data = self.rawData
        
        names = [tde['name'] for tde in data]        
            
        tdes = {}
        sourceWarningThrown = False
        coordWarningThrown = False
        zWarningThrown = False
        for i in range(len(data)):

            tde = data[i]
            
            # try to get the name of the TDE, this will be the key
            try:
                name = tde['name']
            except KeyError as err:
                raise Exception('Name must be provided for every TDE!!!') from err

            # try to get any sources for the TDE
            try:
                sources = tde['sources']
            except KeyError:
                if not sourceWarningThrown:
                    warnings.warn('No sources found at least one TDE!')
                    sourceWarningThrown = True
                sources = None

            # try to get the ra and dec, these aren't necessary since they
            # aren't measured well for some
            try:
                ra = tde['ra'][0]['value']
                dec = tde['dec'][0]['value']
            except KeyError as err:
                if not coordWarningThrown:
                    warnings.warn(f"No RA and DEC found for {tde['name']}!!")
                    coordWarningThrown = True
                ra = None
                dec = None

            # try to get the redshift
            try:
                z = tde['redshift'][0]['value']
            except KeyError as err:
                try:
                    z = tde['z'][0]['value']
                except KeyError as err:
                    if not zWarningThrown:
                        warnings.warn(f"No redshift found for {tde['name']}!!")
                        zWarningThrown = True
                    z = None

            # try to pull out the photometry
            try:
                photo = tde['photometry']
            except KeyError as err:
                photo = None

            # try pulling out the spectra
            try:
                spec = tde['spectra']
            except KeyError as err:
                spec = None            
                    
            thisTDE = TDE(name, ra, dec, z, sources, photometry=photo, spectra=spec)
            tdes[thisTDE.name] = thisTDE
        
        return tdes
    
    def upload(self, tdes:list[TDE]) -> None:
        '''
        Will be used to import new JSON files
        '''

    def query(self,
              names:list[str]=None,
              z:list[float]=None,
              minZ:float=None,
              maxZ:float=None,
              ra:list[str]=None,
              dec:list[str]=None,
              photometryType:str=None,
              spectraType:str=None
              ) -> list[TDE]:
        '''
        Wrapper on AQLQuery.

        If no arguments are provided it returns everything. 

        Args:
            names [list]: A list of the TDE names to grab, default is None. If just a 
                          string is provided it is interpreted as the only names to 
                          get data for.
            z [list] : A list of redshifts as floats. If only one float is provided then
                       it gets all TDEs with that redshift.
            minZ [float]: minimum redshift 
            maxZ [float]: max redshift
            ra [list] : A astring with the exact RA
            dec [list]: A string with the exact Dec
        '''

        # clean up inputs
        queryFilters = ''
        
        if minZ is not None:
            sfilt = f'''
            FOR z1 in tde.z
                FILTER TO_NUMBER(z1.value) >= {minZ}\n
            '''
            queryFilters += sfilt
        if maxZ is not None:
            sfilt = f'''
            FOR z2 in tde.z
                FILTER TO_NUMBER(z2.value) <= {maxZ}\n
            '''
            queryFilters += sfilt
                    
        if names is not None:
            if isinstance(names, str):
                queryFilters += f"FILTER tde.name LIKE '%{names}%'\n"
            elif isinstance(names, list):
                queryFilters += f'FILTER tde.name IN {names}\n'
            else:
                raise Exception('Names must be either a string or list')
            
        if z is not None:
            if isinstance(z, float):
                sfilt = f'''
                FOR z3 in tde.z
                    FILTER TO_NUMBER(z3.value) == {z} \n
                '''
                queryFilters += sfilt
            elif isinstance(names, list):
                sfilt = f'''
                FOR z3 in tde.z
                    FILTER TO_NUMBER(z3.value) IN {z} \n
                '''
                queryFilters += sfilt
            else:
                raise Exception('Redshifts must be either a float or list')

        if ra is not None:
            filt = f'''
            FOR ra in tde.ra
                FILTER ra.value == '{ra}'\n
            '''
            queryFilters += filt
            
        if dec is not None:
            filt = f'''
            FOR d in tde.dec
                FILTER d.value == '{dec}'\n
            '''
            queryFilters += filt

        if photometryType is not None:
            queryFilters += f"FILTER '{photometryType}' IN ATTRIBUTES(tde.photometry)"

        if spectraType is not None:
            queryFilters += f"FILTER '{spectraType}' IN ATTRIBUTES(tde.spectra)"

        # define the query
        query = f'''
        FOR tde IN tdes
            {queryFilters}
            RETURN tde
        '''

        print(query)
        result = self.AQLQuery(query, rawResults=True)
        return self._clean(result)
            
    def close(self) -> None:
        '''
        closes the database connection
        '''
        del self
    
    def __str__(self):
        ret = ''
        for tde in self.tdes:
            ret += str(tde)
            ret += '\n'
        return ret

