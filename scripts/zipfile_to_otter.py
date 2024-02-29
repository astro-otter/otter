'''
Convert a specifically formatted zipfile directory to merge with otter
'''
import os
import uuid
from zipfile import ZipFile

import pandas as pd
import ads

from astropy.coordinates import SkyCoord

from otter import Otter, Transient
from otter.io.helpers import filter_to_obstype
from otter.constants import FILTER_MAP_WAVE

def upload_zip(outdir:str, zipfile:str, testing:bool=False) -> None:
    '''
    Upload a zipfile of information about transients. See the README for info on 
    formatting this zipfile.

    Wraps on self.upload which uploads a json-style file

    Args:
        zipfile [str]: path to the zipfile you want to upload
    '''

    otter = Otter(outdir)
    
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
    schema = [_row_to_json(row, datapath, testing=testing) for _, row in df.iterrows()]

    # now upload this json file
    otter.save(schema, outpath=datapath, test_mode=testing)

def _row_to_json(dfrow:pd.Series, datapath:str, testing=False) -> None:
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
        'coordinate_type': 'equitorial'
    },
    {
        'l': float(galactic.l.value),
        'b': float(galactic.b.value),
        'l_units': 'deg',
        'b_units': 'deg',
        'reference': str(uu),
        'computed': True,
        'coordinate_type': 'galactic'
    }]

    # now add distance measurements if they have any
    dist_keys = ['redshift', 'luminosity_distance', 'dispersion_measure']
    dist_dict = []
    for key in dist_keys:
        if key in dfrow:
            dist_dict.append({'value':dfrow[key],
                               'reference': dfrow.reference,
                               'computed': False,
                               'distance_type': key.replace('_distance', '')
                               })

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
                    'date_type': key.replace('date_', '')
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

    return Transient(schema)

def main():
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument('--otterdir', help='The directory where otter data is stored')
    p.add_argument('--zipfile', help='The path to the zipfile')
    p.add_argument('--testing', dest='testing', action='store_true', help='run in testing mode')
    p.set_defaults(testing=False)
    args = p.parse_args()

    upload_zip(args.otterdir, args.zipfile, testing=args.testing)

if __name__ == '__main__':
    main()
