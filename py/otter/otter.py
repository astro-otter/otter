'''
This is the primary class for user interaction with the catalog
'''
import os
import json
import glob
import uuid
from zipfile import ZipFile
import pandas as pd

from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy import units as u

import ads

from pyArango.database import Database
from pyArango.connection import Connection

from .transient import Transient
from .constants import *
from .helpers import * 

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
                   radius:float=5,
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
                                raw=raw)
        return transients

    def getPhot(self, flux_unit='mag(AB)', date_unit='MJD',
                return_type='astropy', **kwargs
                ) -> Table:
        '''
        Get the photometry of the objects matching the arguments. This will do the
        unit conversion for you!

        Args:
            flux_units [astropy.unit.Unit]: Either a valid string to convert 
                                            or an astropy.unit.Unit
            date_units [astropy.unit.Unit]: Either a valid string to convert to a date 
                                            or an astropy.unit.Unit
            return_type [str]: Either 'astropy' or 'pandas'. If astropy, returns an
                               astropy Table. If pandas, returns a pandas DataFrame.
                               Default is 'astropy'.
            
            **kwargs : Arguments to pass to Otter.query(). Can be:
                       names [list[str]]: A list of names to get the metadata for
                       coords [SkyCoord]: An astropy SkyCoord object with coordinates to match to
                       radius [float]: The radius in arcseconds for a cone search, default is 0.05"
                       minZ [float]: The minimum redshift to search for
                       maxZ [float]: The maximum redshift to search for
                       refs [list[str]]: A list of ads bibcodes to match to. Will only return 
                              metadata for transients that have this as a reference.
                       hasSpec [bool]: if True, only return transients that have spectra.
        
        Return:
           The photometry for the requested transients that match the arguments. 
           Will be an astropy Table sorted by transient default name.
        '''
        queryres = self.query(hasPhot=True, **kwargs)

        dicts = []
        for transient in queryres:

            # clean the photometry
            default_name = transient['name/default_name']
            phot = transient.cleanPhotometry(flux_unit=flux_unit, date_unit=date_unit)
            phot['name'] = [default_name]*len(phot)
            
            dicts.append(phot)
            
        fullphot = pd.concat(dicts)

        # remove some possibly confusing keys
        del fullphot['raw']
        del fullphot['raw_units']
        del fullphot['date_format']
        
        if return_type == 'astropy':
            return Table.from_pandas(fullphot)
        elif return_type == 'pandas':
            return fullphot
        else:
            raise ValueError('return_type can only be pandas or astropy')

    
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
        Wraps on the super.AQLQuery and queries the OTTER database more intuitively.

        WARNING! This does not do any conversions for you! 
        This is how it differs from the `getMeta` method. Users should prefer to use
        `getMeta`, `getPhot`, and `getSpec` independently because it is a better 
        workflow and can return the data in an astropy table with everything in the
        same units.

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
           Get all of the raw (unconverted!) json data for objects that match the criteria. 
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
                for coordinfo in tde['coordinate']:
                    if 'ra' in coordinfo and 'dec' in coordinfo:
                        coord = SkyCoord(coordinfo['ra'], coordinfo['dec'],
                                         unit=(coordinfo['ra_units'], coordinfo['dec_units']))
                    elif 'l' in coordinfo and 'b' in coordinfo:
                        # this is galactic
                        coord = SkyCoord(coordinfo['l'], coordinfo['b'],
                                         unit=(coordinfo['l_units'], coordinfo['b_units']),
                                         frame='galactic'
                                         )
                    else:
                        raise ValueError('Either needs to have ra and dec or l and b as keys!')
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


    def upload_zip(self, zipfile:str, testing:bool=False) -> None:
        '''
        Upload a zipfile of information about transients. See the README for info on 
        formatting this zipfile.

        Wraps on self.upload which uploads a json-style file

        Args:
            zipfile [str]: path to the zipfile you want to upload
        '''

        datadir = os.path.dirname(zipfile)
        datapath = os.path.join(datadir, os.path.basename(zipfile).replace('.zip', ''))
        metapath = os.path.join(datapath, 'meta.csv')
        
        # read in the zipfile and extract all of the files into datapath
        with ZipFile(zipfile) as z:
            z.extractall(datadir)

        # read in the metadata
        # now read in the meta file
        if not os.path.exists(metapath):
            raise ValueError('meta.csv not found, can not upload this data!')
        df = pd.read_csv(metapath)
            
        # convert the rows to a json file
        schema = [self._row_to_json(row, datapath, testing=testing) for _, row in df.iterrows()]
        
        # now upload this json file
        self.upload(schema, outpath=datapath, test_mode=testing)
            
    def upload(self, schema:list[dict], outpath:str=os.getcwd(), **kwargs) -> None:
        '''
        Upload all the data in the given list of schemas.

        Args:
            schema [list[dict]]: A list of json dictionaries
            outpath [str]: The path to the directory where to write the json files.
                           Default is current working directory.
        '''

        if not isinstance(schema, list):
            schema = [schema]
        
        for json in schema:

            # convert the json to a Transient
            if not isinstance(json, Transient):
                json = Transient(json)
                
            coord = json.getSkyCoord()
            res = self.coneSearch(coords=coord)

            if len(res) == 0:
                # This is a new object to upload
                print('Adding this as a new object...')    
                self._upload_document(dict(json), **kwargs)
                
            else:
                # We must merge this with existing data
                print('Found this object in the database already, merging the data...')
                if len(res) == 1:
                    # we can just add these to merge them!
                    combined = res[0] + json
                    self._upload_document(combined, **kwargs)
                else:
                    # for now throw an error
                    # this is a limitation we can come back to fix if it is causing
                    # problems though!
                    raise Exception('Current Limitation Found: Some objects in Otter are too close!')

    def _upload_document(self, schema, test_mode=False):
        '''
        Upload a json file in the correct format to the OTTER database
        '''
        # check if this documents key is in the database already
        # and if so remove it!
        jsonpath = os.path.join(DATADIR, '*.json')
        aliases = {item['value'] for item in schema['name']['alias']}
        filenames = {os.path.basename(fname).split('.')[0] for fname in glob.glob(jsonpath)}
        todel = list(aliases & filenames)
        if len(todel) > 0 and not test_mode:
            os.remove(os.path.join(DATADIR, todel[0]+'.json'))
        else:
            # for testing
            print('Deleting the following file: ', todel)
        
        # now do two things to save this data
        # 1) create a new file in "base" with this
        outfilepath = os.path.join(DATADIR,schema['name']['default_name']+'.json')
        
        # format as a json
        if isinstance(schema, Transient):
            schema = dict(schema)

        out = json.dumps(schema, indent=4)
        out = '[' + out
        out += ']'

        if not test_mode:
            with open(outfilepath, 'w') as f:
                f.write(out)

            # 2) upload to the database
            c = self[self.collectionName]
            newdoc = c.createDocument(schema)
            newdoc.save()
        else:
            print(out)
        
    def _row_to_json(self, dfrow:pd.Series, datapath:str, testing=False) -> None:
        '''
        Add a new transient to otter because the input doesn't match any existing transients

        Args:
            json [dict]: a dictionary in the correct json format with the correct keys
        '''
        dfrow.dropna(inplace=True)
        uu = uuid.uuid4()
        
        schema = {
            'schema_version': {'value': '0'},
            'name': {'default_name': dfrow['name'],
                     'alias': [{'value': dfrow['name'], 'reference': dfrow['reference']}]
                     },
        }

        # add coordinates
        # these were required so we don't have to check anything
        galactic = SkyCoord(dfrow.ra, dfrow.dec, unit=(dfrow.ra_units, dfrow.dec_units),
                            frame='icrs').galactic
        
        schema['coordinate'] = [{
            'ra': dfrow.ra,
            'dec': dfrow.dec,
            'ra_units': dfrow.ra_units,
            'dec_units': dfrow.dec_units,
            'computed': False,
            'uuid': str(uu),
            'default': True,
            'reference': dfrow.reference,
            'coord_type': 'equitorial'
        },
        {
            'l': float(galactic.l.value),
            'b': float(galactic.b.value),
            'l_units': 'deg',
            'b_units': 'deg',
            'reference': str(uu),
            'computed': True,
            'coord_type': 'galactic'
        }]
        
        # now add distance measurements if they have any
        dist_keys = ['redshift', 'luminosity_distance', 'dispersion_measure']
        dist_dict = {}
        for key in dist_keys:
            if key in dfrow:
                dist_dict[key] = [{'value':dfrow[key],
                                   'reference': dfrow.reference,
                                   'computed': False,
                                   }]

        if len(dist_dict) > 0:
            schema['distance'] = dist_dict

        # do the same thing with the classification
        if 'class' in dfrow:
            schema['classification'] = [{'object_class': dfrow['class'],
                                         'confidence': 1.0, # THIS IS DANGEROUS! FIX LATER
                                         'reference': dfrow.reference,
                                         'defualt': True
                                         }]
            
            
        # do the same with epoch info
        epoch_keys = ['date_discovery', 'date_peak', 'date_explosion']
        epoch_dict = []
        for key in epoch_keys:
            if key in dfrow:
                epoch_dict.append(
                    {
                        'value':dfrow[key],
                        'date_format': dfrow['date_format'],
                        'reference': dfrow.reference,
                        'computed': False,
                        'measurement_type': key
                    }
                )
        if len(epoch_dict) > 0:
            schema['date_reference'] = epoch_dict

        # create the reference_alias        
        if not testing: # we don't want to use up all our queries
            adsquery = list(ads.SearchQuery(bibcode=dfrow.reference))[0]
            authors = adsquery.author
            year = adsquery.year

            if len(authors) == 0:
                raise ValueError('This ADS bibcode does not exist!')
            elif len(authors) == 1:
                author = authors[0]
            elif len(authors) == 2:
                author = authors[0] + ' & ' + authors [1]
            else: # longer than 2
                author = authors[0] + ' et al.'

            # generate the human readable name
            hrn = author + ' (' + year + ')'
            schema['reference_alias'] = [{'name': dfrow.reference,
                                          'human_readable_name':hrn
                                          }]
        else:
            print(f'We would be querying for bibcode={dfrow.reference}')
        
        # check if there is a photometry file path
        if 'phot_path' in dfrow:
            phot = pd.read_csv(os.path.join(datapath, dfrow.phot_path))

            # remove any annoying columns from phot
            for key in phot:
                if 'Unnamed' in key:
                    del phot[key]

            # replace all NaNs with null
            phot.fillna('null', inplace=True)

            # rename the filter key to filter_key
            phot['filter_key'] = phot['filter']
            del phot['filter']
            
            # get the observation types
            phot['obs_type'] = phot['filter_key'].apply(filter_to_obstype)
            
            # convert to a dictionary
            phot_dict = phot.to_dict(orient='list')
            phot_dict['reference'] = dfrow.reference
            
            # put this in the schema
            schema['photometry'] = [phot_dict]

            # create a filter alias
            filteralias = []
            filteralias_keys = []
            for idx in range(len(phot_dict['filter_key'])):
                if phot_dict['filter_key'][idx] in filteralias_keys: continue

                filtername = phot_dict['filter_key'][idx]
                filteralias_keys.append(filtername) # to make sure we don't duplicate
                
                indict = {
                    'filter_key': filtername
                }
                
                if phot_dict['obs_type'][idx] == 'radio':
                    if 'filter_eff' not in phot_dict.keys():
                        eff = FILTER_MAP_FREQ[filtername]
                        eff_units = 'THz'
                    else:
                        eff = phot_dict['filter_eff'][idx]
                        eff_units = phot_dict['filter_eff_units'][idx]

                else:
                    if 'filter_eff' not in phot_dict.keys():
                        eff = FILTER_MAP_WAVE[filtername]
                        eff_units = 'nm'
                    else:
                        eff = phot_dict['filter_eff'][idx]
                        eff_units = phot_dict['filter_eff_units'][idx]

                if 'hz' in eff_units.lower():
                    indict['freq_eff'] = eff
                    indict['freq_units'] = eff_units
                else:
                    indict['wave_eff'] = eff
                    indict['wave_units'] = eff_units    

                filteralias.append(indict)

            # add these filteraliases to the schema
            schema['filter_alias'] = filteralias
                    
        # check if there is a spectra file path
        if 'spec_path' in dfrow:
            pass # ADD THIS CODE ONCE WE HANDLE SPECTRA
        
        return schema
