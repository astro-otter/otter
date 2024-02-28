'''
This is the primary class for user interaction with the catalog
'''
import os
import json
import glob
from warnings import warn
import uuid
from zipfile import ZipFile
from collections.abc import Iterable

import pandas as pd
import numpy as np

from astropy.coordinates import SkyCoord, search_around_sky
from astropy.table import Table
from astropy import units as u

import ads

from .transient import Transient
from ..constants import *
from .helpers import * 

class Otter(object):
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
                 datadir:str=None,
                 debug:bool=False
                 ) -> None:
        
        # save inputs
        if datadir is None:
            self.CWD = os.path.dirname(os.path.abspath('__FILE__'))
            self.DATADIR = os.path.join(self.CWD, '.otter')
        else:
            self.CWD = os.path.dirname(datadir)
            self.DATADIR = datadir

        self.debug = debug

        # make sure the data directory exists
        if not os.path.exists(self.DATADIR):
            try:
                os.makedirs(self.DATADIR)
            except FileExistsError:
                warn('Directory was created between the if statement and trying to create the directory!')
                pass
        
    def getMeta(self, **kwargs) -> Table:
        '''
        Get the metadata of the objects matching the arguments

        Args:
            **kwargs : Arguments to pass to Otter.query()
        Return:
           The metadata for the transients that match the arguments. Will be an astropy
           Table by default, if raw=True will be a dictionary.
        '''
        metakeys = ['name', 'coordinate', 'date_reference', 'distance', 'classification']
        
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


    def load_file(self, filename:str) -> dict:
        '''
        Loads an otter JSON file
        
        Args:
            filename [str]: The path to the file to load
        '''

        # read in files from summary
        with open(filename, 'r') as f:
            to_ret = Transient(json.load(f))

        return to_ret

        
    def query(self,
              names:list[str]=None,
              coords:SkyCoord=None,
              radius:float=5,
              minZ:float=0,
              maxZ:float=None,
              refs:list[str]=None,
              hasPhot:bool=False,
              hasSpec:bool=False,
              raw:bool=False
              ) -> dict:
        '''
        Searches the summary.csv table and reads relevant JSON files

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
        if all(arg is None for arg in [names, coords, maxZ, refs]) and not hasPhot and not hasSpec:
            # there's nothing to query!
            # read in the metdata from all json files
            # this could be dangerous later on!!
            allfiles = glob.glob(os.path.join(self.DATADIR, '*.json'))
            jsondata = []

            # read the data from all the json files and convert to Transients
            for jsonfile in allfiles:
                with open(jsonfile, 'r') as j:
                    t = Transient(json.load(j))
                    jsondata.append(t.getMeta())

            return jsondata
        
        # check if the summary table exists, if it doen't create it
        summary_table = os.path.join(self.DATADIR, 'summary.csv')
        if not os.path.exists(summary_table):
            self.generate_summary_table(save=True)

        # then read and query the summary table
        summary = pd.read_csv(summary_table)

        # coordinate search first
        if coords is not None:
            if not isinstance(coords, SkyCoord):
                raise ValueError('Input coordinate must be an astropy SkyCoord!')
            summary_coords = SkyCoord(summary.ra.tolist(), summary.dec.tolist(), unit=(u.hourangle, u.deg))

            try:
                summary_idx, _, _, _ = search_around_sky(summary_coords,
                                                         coords,
                                                         seplimit=radius*u.arcsec)
            except ValueError:
                summary_idx, _, _, _ = search_around_sky(summary_coords,
                                                         SkyCoord([coords]),
                                                         seplimit=radius*u.arcsec)
                
            summary = summary.iloc[summary_idx]
        
        # redshift
        summary = summary[summary.z.astype(float) >= minZ]
        if maxZ is not None:
            summary = summary[summary.z.astype(float) <= maxZ]
            
        # check photometry and spectra
        if hasPhot:
            summary = summary[summary.hasPhot == True]
            
        if hasSpec:
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
        self.save(schema, outpath=datapath, test_mode=testing)
            
    def save(self, schema:list[dict], **kwargs) -> None:
        '''
        Upload all the data in the given list of schemas.

        Args:
            schema [list[dict]]: A list of json dictionaries
        '''

        if not isinstance(schema, list):
            schema = [schema]
        
        for json in schema:

            print(json['name/default_name'])
            
            # convert the json to a Transient
            if not isinstance(json, Transient):
                json = Transient(json)
                
            coord = json.getSkyCoord()
            res = self.coneSearch(coords=coord)

            if len(res) == 0:
                # This is a new object to upload
                print('Adding this as a new object...')    
                self._save_document(dict(json), **kwargs)
                
            else:
                # We must merge this with existing data
                print('Found this object in the database already, merging the data...')
                if len(res) == 1:
                    # we can just add these to merge them!
                    combined = res[0] + json
                    self._save_document(combined, **kwargs)
                else:
                    # for now throw an error
                    # this is a limitation we can come back to fix if it is causing
                    # problems though!
                    raise Exception('Current Limitation Found: Some objects in Otter are too close!')

        # update the summary table appropriately
        self.generate_summary_table(save=True)
                
    def _save_document(self, schema, test_mode=False):
        '''
        Save a json file in the correct format to the OTTER data directory
        '''
        # check if this documents key is in the database already
        # and if so remove it!
        jsonpath = os.path.join(self.DATADIR, '*.json')
        aliases = {item['value'].replace(' ', '-') for item in schema['name']['alias']}
        filenames = {os.path.basename(fname).split('.')[0] for fname in glob.glob(jsonpath)}
        todel = list(aliases & filenames)

        # now save this data
        # create a new file in self.DATADIR with this
        if len(todel) > 0:
            outfilepath = os.path.join(self.DATADIR, todel[0]+'.json')
            if test_mode:
                print('Renaming the following file for backups: ', outfilepath)
            else:
                os.rename(outfilepath, outfilepath+'.backup')
        else:
            if test_mode:
                print("Don't need to mess with the files at all!")
            fname = schema['name']['default_name']+'.json'
            fname = fname.replace(' ', '-') # replace spaces in the filename
            outfilepath = os.path.join(self.DATADIR, fname)
            
        
        # format as a json
        if isinstance(schema, Transient):
            schema = dict(schema)

        out = json.dumps(schema, indent=4)
        #out = '[' + out
        #out += ']'

        if not test_mode:
            with open(outfilepath, 'w') as f:
                f.write(out)
        else:
            print(f'Would write to {outfilepath}')
            print(out)

    def generate_summary_table(self, save=False):
        '''
        Generate a summary table for the JSON files in self.DATADIR

        args:
            save [bool]: if True, save the summary file to "summary.csv"
                         in self.DATADIR. Default is False.
        '''
        allfiles = glob.glob(os.path.join(self.DATADIR, '*.json'))
        jsondata = []
        
        # read the data from all the json files and convert to Transients
        rows = []
        for jsonfile in allfiles:
            with open(jsonfile, 'r') as j:
                t = Transient(json.load(j))
                skycoord = t.getSkyCoord()

                row = {'name': t.default_name,
                       'alias': [alias['value'] for alias in t['name']['alias']], 
                       'ra': skycoord.ra,
                       'dec': skycoord.dec,
                       'refs': [ref['name'] for ref in t['reference_alias']]
                       }
                
                if 'date_reference' in t:
                    row['discovery_date'] = t.getDiscoveryDate()

                if 'distance' in t:
                    row['z'] = t.getRedshift()

                row['hasPhot'] = 'photometry' in t
                row['hasSpec'] = 'spectra' in t
                    
                row['json_path'] = os.path.abspath(jsonfile)

                rows.append(row)
                
        alljsons = pd.DataFrame(rows)
        if save:
            alljsons.to_csv(os.path.join(self.DATADIR, 'summary.csv'))

        return alljsons
        
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
