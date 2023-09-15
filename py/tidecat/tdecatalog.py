'''
Class for to make queries in the catalog
'''
import warnings
from pyArango.database import Database
from pyArango.connection import Connection
from .tde import TDE

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
        
    def _clean(self) -> list[TDE]:
        '''
        Get the data from the database and clean it 
    
        Args:
            tdes [list[dict]]: output of connect
        
        Returns: names, sources, ra, dec, z
        '''
        names = [tde['name'] for tde in self.rawData]
        
            
        tdes = []
        sourceWarningThrown = False
        coordWarningThrown = False
        zWarningThrown = False
        for i in range(len(self.rawData)):

            tde = self.rawData[i]
            
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

            thisTDE = TDE(name, ra, dec, z, sources)
            tdes.append(thisTDE)
        
        return tdes
    
    def upload(self, tdes:list[TDE]) -> None:
        '''
        Will be used to import new JSON files
        '''

    def query(self) -> list[TDE]:
        '''
        wrapper on AQLQuery
        '''

    def __str__(self):
        ret = ''
        for tde in self.tdes:
            ret += str(tde)
            ret += '\n'
        return ret
