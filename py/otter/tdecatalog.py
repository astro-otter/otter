'''
Class for to make queries in the catalog
'''
from copy import deepcopy
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
                 collection='tdes',
                 debug=False):

        # save inputs
        self.dbName = db
        self.collectionName = collection
        
        c = Connection(username=username, password=password)
        
        # initiate the tdes database
        super().__init__(c, db)

        if not debug:
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
            
        tdes = {}
        for tde in data:
            thisTDE = TDE(tde)
            tdes[thisTDE.name.name] = thisTDE
        
        return tdes
    
    def upload(self, tdes:list[TDE], test=False) -> None:
        '''
        Will be used to import new JSON files
        '''

        c = self[self.collectionName] # collection to add data to

        for tde in tdes:
            inCatalog = tde.name.name in self.tdes
            if inCatalog:
                query = f"FOR tde IN tdes FILTER tde.name == '{tde.name.name}' RETURN tde"
                catdata = self.AQLQuery(query, rawResults=True)
                if len(catdata) != 0:
                    catjson = TDE(catdata[0]).tojson()
                    out = self._append(tde, catjson)
                    out['_key'] = catdata[0]['_key']
                    replaceQuery = f"REPLACE {out} IN tdes"
                    if test:
                        print(replaceQuery)
                        print()
                    else:
                        self.AQLQuery(replaceQuery)
                else:
                    inCatalog = False
            if not inCatalog:
                json = tde.tojson()
                doc = c.createDocument(json)
            
                if not test:
                    doc.save() # update the database
                else:
                    print(doc['name'])
                    print()
            
    def _append(self, tde, doc):
        '''
        Appends data to an existing entry in the catalog and returns the new 
        collection for upload
        '''
        aliases = [int(s['alias']) for s in doc['sources']]
        bibcodes = {s['bibcode']:s['alias'] for s in doc['sources']}
        if len(aliases) == 0:
            newAlias = 1
        else:
            newAlias = max(aliases)+1    
            
        json = tde.tojson()
        # get the correct alias to handle sources well
        sourcemap = {} # key is old source number and value is new source number
        for val in json['sources']:
            if val['bibcode'] in bibcodes:
                # this source is already in the database!!!
                alias = bibcodes[val['bibcode']]
            else:
                alias = newAlias
                newAlias += 1

            sourcemap[str(val['alias'])] = str(alias)

            val['alias'] = str(alias)
            
        # append or add keys now
        for key in json:
            if key == 'name': continue            
            if key in doc:                
                for val in json[key]:

                    if key == 'sources' and val['bibcode'] in bibcodes: continue
                    
                    # update the source alias for specific content
                    if 'source' in val:
                        src = val['source']
                        if ',' in str(src):
                            srcs = src.split(',')
                            val['source'] = ','.join([sourcemap[s] for s in srcs])
                        else:
                            val['source'] = sourcemap[src]

                        
                    # append the value to the corresponding key!
                    if key == 'photometry' or key == 'spectra':
                        for filt in json[key]: # loop over the filters
                            if filt in doc[key]: # already have data for this band
                                for item in json[key][filt]:
                                    doc[key][filt].append(item)
                            else: # no data for this filter
                                doc[key][filt] = deepcopy(json[key][filt])
                    else:
                        if isinstance(val, dict):
                            doc[key].append(val)
                        else:
                            for item in val:
                                doc[key].append(item)
            else:
                doc[key] = json[key]
        print(doc['photometry'])
        return doc

    def query(self,
              names:list[str]=None,
              z:list[float]=None,
              minZ:float=None,
              maxZ:float=None,
              ra:list[str]=None,
              dec:list[str]=None,
              searchRadius:int=5,
              photometryType:str=None,
              spectraType:str=None,
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
        
        result = self.AQLQuery(query, rawResults=True)

        # now that we have the query results do the RA and Dec queries if they exist
        cleanResult = self._clean(result) 
        if ra is not None and dec is not None:
            # get the catalog RAs and Decs to compare against
            queryCoords = SkyCoord(ra, dec, unit=u.deg)
            goodTDEs = {}
            for tde in cleanResult.values():
                for ra, dec in zip(tde.ra, tde.dec):
                    coord = SkyCoord(ra.valueString, dec.valueString)
                    if queryCoords.separation(coord) < searchRadius*u.arcsec:
                        goodTDEs[tde.name.name] = tde
            return goodTDEs

        elif ra is not None or dec is not None:
            raise Exception('You must provide both an RA and Dec')

        else:
            return cleanResult
            
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
