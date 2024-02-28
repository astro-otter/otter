'''
Reformat some original raw json files of radio data for OTTER
'''

import os, glob
import json
import uuid
from copy import deepcopy
import numpy as np

from astropy.coordinates import SkyCoord
import astropy.units as u

from otter import Otter

# some important meta params
t_disc_map = {
    'AT2022cmc': 59621.45,
    'AT2020vwl': 59132.73,
    'ASASSN-15oi' : 57248.2,
    'XMMSL1 J0740-85': 56748,
    'AT2019azh' : 58536.02,
    'Sw J1644+57': 55648.54,
    'AT2018hyz':  58428.64,
    'AT2019dsg': 58582.46,
    'IGR J12580+0134' : 55572,
    'ASASSN-14li': 56983,
    'AT2020opy' : 59038.23
}

thisfile_dir = os.path.dirname(os.path.realpath(__file__)) #"scripts" directory
otter_dir = os.path.join(thisfile_dir, '.otter') # set this as the default but it probably shouldn't be used

def main():

    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('--indir', help='Directory for the basic radio data')
    p.add_argument('--outdir', help='Directory where the otter json files will go',
                   default=otter_dir)
    args = p.parse_args()

    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)
    
    # convert the input arguments
    jsonpath = args.indir
    jsons = glob.glob(os.path.join(jsonpath, '*.json'))
    
    # reformat these appropriately
    for j in jsons:
        with open(j, 'r') as s:
            data = json.load(s)[0]

        # get names organized
        names = [data['name']]
        nameRefs = []
        for name in names:
            if 'AT' in name:
                nameRefs.append('TNS')
            elif 'Sw' in name:
                nameRefs.append('Swift')
            elif 'XMM' in name:
                nameRefs.append('XMM-Newton')
            elif 'CNSS' in name:
                nameRefs.append('CNSS')
            elif 'ASASSN' in name:
                nameRefs.append('ASASSN')
            elif 'IGR' in name:
                nameRefs.append('IGR')

        # get coordinates organized
        coords = [(ra['value'], dec['value']) for ra, dec in zip(data['ra'], data['dec'])]
        units = [('hourangle', 'deg') for i in range(len(data['ra']))]
        sources = [source['bibcode'] for source in data['sources']]
        default = [False]*len(data['ra'])
        default[0] = True
        uuids = [uuid.uuid4() for i in range(len(data['ra']))]
        galactic = SkyCoord(coords, unit=(u.hourangle, u.deg), frame='icrs').galactic

        # get redshifts
        zs = [z['value'] for z in data['z']]
        defaultz = [False]*len(zs)
        defaultz[0] = True

        # define the general info schema
        schema = {
             "schema_version": {"value":'0',
                        'comment': 'Original Dataset'
                        },
            "name": {'default_name': data['name'],
                     'alias': [{'value': name, 'reference':ref} for name, ref in zip(names, nameRefs)]
                    },
            "coordinate": [
                {'ra': coord[0], 
                 'dec': coord[1], 
                 'epoch': 'J2000',
                 'system': 'ICRS',
                 'ra_units': unit[0],
                 'dec_units': unit[1],
                 'reference': sources,
                 'computed': False,
                 'default': defa,
                 'uuid': str(uu),
                 'coordinate_type': 'equitorial'
                } for coord, unit, defa, uu in zip(coords, units, default, uuids)
            ] + [
                {'l': float(row.l.value),
                 'b': float(row.b.value),
                 'l_units': 'deg',
                 'b_units': 'deg',
                 'reference': str(uu),
                 'computed': True,
                 'coordinate_type': 'galactic'
                } for row, uu in zip(galactic, uuids)],

            "distance": [{'value': z,
                          'reference': sources,
                          'computed': False,
                          'default': d,
                          'distance_type': 'redshift'
                          } for z, d in zip(zs, defaultz)],

            "classification": [{
                'object_class': 'TDE',
                'confidence': 1.0,
                'reference': sources,
                'default': True
            }],

            'reference_alias': [{"name":src['bibcode'],
                                 "human_readable_name":src['name']
                                } for src in data['sources'] 
                               ],      

             }

        # epoch info if available
        if 't_disc' in data:
            schema['date_reference'] = [{"value":date['value'],
                                "date_format":"mjd", 
                                "reference": sources,
                                "computed": False,
                                "date_type": "discovery"
                                } for date in data['t_disc']]
        else:
            schema['date_reference'] = [{"value":t_disc_map[names[0]],
                                "date_format":"mjd", 
                                "reference": sources,
                                "computed": False,
                                "date_type": "discovery"
                               }]

        # add in the photometry
        schema['photometry'] = {}

        ## first find the different sources
        photsources = {}
        for photo in data['photometry']['radio']:
            if photo['source'] not in photsources:

                alias = str(photo['source'])
                if alias == 'unknown':
                    photsources[alias] = sources
                # get the bibcode
                for src in data['sources']:
                    if src['alias'] == alias:
                        # put in the source mapping
                        photsources[alias] = src['bibcode']
                        break     


        ## then iterate over the possible sources and find all photometry values with that
        i = 0
        samplephot = {}
        for alias, bibcode in photsources.items():
            for photo in data['photometry']['radio']:
                if photo['source'] == alias:
                    name = f'phot_{i}'
                    if name not in schema['photometry']:
                        schema['photometry'][name] = {'reference':bibcode}
                        samplephot[name] = []

                    samplephot[name].append(photo)
            i += 1

        filteraliases = []
        schema['photometry'] = []
        for phot in samplephot:
            source = photsources[samplephot[phot][0]['source']]
            if not isinstance(source, list):
                source = [source] # FIX ME LATER!!!! SHOULD NOT JUST TAKE THE FIRST!!!!!
            photdict = {'reference': source,
                        'raw': [],
                        'raw_units': [],
                        'date': [],
                        'date_format': [],
                        'filter_key': [],
                        'computed': [],
                        'obs_type': [],
                        'upperlimit': []
                        }
            for point in samplephot[phot]: # all these points belong together
                if 'nu_GHz' in data:
                    filteralias = data['nu_GHz'][0]['value'] + 'GHz'
                else: # it's in the individual photometry
                    filteralias = str(point['nu_GHz']) + 'GHz'

                if filteralias not in filteraliases:
                    filteraliases.append(filteralias)

                photdict['raw'].append(point['F_mJy'])
                photdict['raw_units'].append('mJy')
                photdict['date_format'].append('mjd')
                photdict['filter_key'].append(filteralias)
                photdict['computed'].append(False)
                photdict['obs_type'].append('radio')

                if 'upperlimit' in data:
                    photdict['upperlimit'].append(eval(data['upperlimit'][0]['value']))
                else:
                    photdict['upperlimit'].append(point['upperlimit'])

                if point['t_MJD'] is not None:
                    photdict['date'].append(point['t_MJD'])
                else:
                    mjd = t_disc_map[names[0]] + point['dt_days']
                    photdict['date'].append(mjd)

            schema['photometry'].append(photdict)

        # add a filter_alias property
        schema['filter_alias'] = []
        for filtername in filteraliases:
            forSchema = {'filter_key': filtername,
                         'freq_eff': float(filtername.replace('GHz', '')),
                         'freq_units': 'GHz'
                        }
            schema['filter_alias'].append(forSchema)


        out = json.dumps(schema, indent=4)

        outpath = os.path.join(args.outdir, os.path.basename(j))
        with open(outpath, 'w+') as outfile:
            outfile.write(out)

    otter = Otter(args.outdir)
    otter.generate_summary_table(save=True)
            
if __name__ == '__main__':
    main()
    
