'''
Converts the data from tde.space to the otter json format while also incorporating it
with existing files
'''

# imports
import os, glob, warnings
from copy import deepcopy
import numpy as np
import pandas as pd
import json
from otter import Otter, Transient
from otter import util as otter_const
from otter import util as otter_helper
from astroquery.simbad import Simbad
import astropy.units as u
import astropy.constants as const

# some useful mappings
bandwavelengths = {
    "u"      : 354.,
    "g"      : 475.,
    "r"      : 622.,
    "i"      : 763.,
    "z"      : 905.,
    "u'"     : 354.,
    "g'"     : 475.,
    "r'"     : 622.,
    "i'"     : 763.,
    "z'"     : 905.,
    "u_SDSS" : 354.3,
    "g_SDSS" : 477.0,
    "r_SDSS" : 623.1,
    "i_SDSS" : 762.5,
    "z_SDSS" : 913.4,
    "U"      : 365.,
    "B"      : 445.,
    "V"      : 551.,
    "R"      : 658.,
    "I"      : 806.,
    "Y"      : 1020.,
    "J"      : 1220.,
    "H"      : 1630.,
    "K"      : 2190.,
    "M2"     : 260.,
    "W1"     : 224.6,
    "W2"     : 192.8,
    "w"      : 622.,
}

bandreps = {
    'Ks': ['K_s'],
    'M2': ['uvm2', 'UVM2', 'UVm2', 'Um2', 'm2', 'um2'],
    'W1': ['uvw1', 'UVW1', 'UVw1', 'Uw1', 'w1', 'uw1'],
    'W2': ['uvw2', 'UVW2', 'UVw2', 'Uw2', 'w2', 'uw2']
}

# this isn't correct and can be updated late
# But, just use central energy to compute a central wavelength
xraycodes = {
    "0.3 - 10" : (const.h*const.c / (10-0.3)/u.keV).to(u.nm).value,
    "0.5 - 8": (const.h*const.c / (8-0.5)/u.keV).to(u.nm).value,
    '0.3 - 2.0': (const.h*const.c / (2-0.3)/u.keV).to(u.nm).value,
    '0.2 - 2.0': (const.h*const.c / (2-0.2)/u.keV).to(u.nm).value
}

# from https://imagine.gsfc.nasa.gov/science/toolbox/spectrum_chart.html
bandOptions = {'radio':{'min':1e-3, 'max':1e-100}, 
               'x-ray':{'min':1e-11, 'max':1e-8}, 
               'optical':{'min':4e-7, 'max':7e-7}, 
               'infared':{'min':7e-7, 'max':1e-3}, 
               'uv':{'min':1e-8, 'max':4e-7}
              }

# now the functions

def clean_schema(schema):
    '''
    Clean out Nones and empty lists from the given subschema
    '''
    for key, val in list(schema.items()):
        if val is None or (isinstance(val, (list, dict)) and len(val) == 0):
            del schema[key]
    return schema

def main():

    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('--indir', help='Directory for the basic radio data')
    p.add_argument('--outdir', help='Directory where the otter json files will go')
    args = p.parse_args()

    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)

    exc = {'ASAS-SN Supernovae', 'MAST'}

    mappedsrc = lambda source_map, src: [source_map[key] for key in src if key in source_map]

    allschemas = []
    allsrcs = 0
    badsrcs = 0
    allphot = 0
    badphot = 0

    for file in glob.glob(os.path.join(args.indir, '*.json')):
        print(f'Reformating {file}')
        with open(file, 'r') as f:
            j = json.load(f)[0]

        # copy over the references
        schema = Transient(deepcopy(otter_const.schema))
        source_map = {}

        for src in j['sources']:
            allsrcs += 1 # another source, double are okay I think...
            if 'bibcode' not in src and src['name'] in exc:
                bib = src['name']
            elif 'bibcode' not in src:
                print(f'No bibcode for {src}, skipping!')
                badsrcs += 1
                continue
            else:
                bib = src['bibcode']
            sub = deepcopy(otter_const.subschema['reference_alias'])
            sub['name'] = bib
            if 'reference' in src:
                sub['human_readable_name'] = src['reference']
            else:
                sub['human_readable_name'] = src['name']
            schema['reference_alias'].append(sub)

            source_map[src['alias']] = bib

        # copy over names and aliases
        schema['name/default_name'] = j['name']
        for val in j['alias']:
            sub = deepcopy(otter_const.subschema['name/alias'])
            sub['value'] = val['value']
            sub['reference'] = mappedsrc(source_map, val['source'].split(','))
            schema['name/alias'].append(sub)

        # copy over coordinates
        if 'ra' in j and 'dec' in j:
            key1, key2 = 'ra', 'dec'
        elif 'hostra' in j and 'hostdec' in j:
            key1, key2 = 'hostra', 'hostdec'
        else:
            print(f'Skipping {file} because no ra and dec associated with it!')
            continue

        for ra, dec in zip(j[key1], j[key2]):
            sub = deepcopy(otter_const.subschema['coordinate'])
            if ra['source'] != dec['source']:
                src = list(set(ra['source'].split(',')) or set(ra['source'].split(',')))
            else:
                src = ra['source'].split(',')

            if not any(s in source_map for s in src):
                print(f'Skipping {(ra, dec)} because it does not have a reliable reference!')
                continue

        if ra['u_value'] == 'hours' or ':' in ra['value']:
            ra_u = 'hour'
        elif ra['u_value'] == 'degrees' or ra['u_value'] == 'floatdegrees':
            ra_u = 'deg'
        else:
            ra_u = ra['u_value']

        if dec['u_value'] == 'degrees' or dec['u_value'] == 'floatdegrees':
            dec_u = 'deg'
        else:
            dec_u = dec['u_value']

        sub['ra'] = ra['value']
        sub['dec'] = dec['value']
        sub['ra_units'] = ra_u
        sub['dec_units'] = dec_u
        sub['reference'] = mappedsrc(source_map, src)
        sub['coordinate_type'] = 'equitorial'
        schema['coordinate'].append(clean_schema(sub))

        # copy over distance measurements
        # first redshift
        if 'redshift' in j or 'lumdist' in j or 'comovingdist' in j:    
            if 'redshift' in j:
                for ii, z in enumerate(j['redshift']):
                    sub = deepcopy(otter_const.subschema['distance'])
                    src = z['source'].split(',')
                    if not any(s in source_map for s in src):
                        print(f'Skipping {z} because it does not have a reliable reference!')
                        continue

                    sub['value'] = z['value']
                    sub['reference'] = mappedsrc(source_map, src)
                    sub['computed'] = False
                    sub['distance_type'] = 'redshift'
                    if ii == 0:
                        sub['default'] = True
                    schema['distance'].append(clean_schema(sub))
                del j['redshift']

            if 'lumdist' in j:
                for ii, d in enumerate(j['lumdist']):
                    sub = deepcopy(otter_const.subschema['distance'])
                    src = d['source'].split(',')
                    if not any(s in source_map for s in src):
                        print(f'Skipping {d} because it does not have a reliable reference!')
                        continue

                    sub['value'] = d['value']
                    sub['reference'] = mappedsrc(source_map, src)
                    sub['computed'] = False
                    sub['unit'] = d['u_value']
                    sub['distance_type'] = 'luminosity'
                    if ii == 0:
                        sub['default'] = True
                    schema['distance'].append(clean_schema(sub))
                del j['lumdist']

            if 'comovingdist' in j:
                for ii, d in enumerate(j['comovingdist']):
                    sub = deepcopy(otter_const.subschema['distance'])
                    src = d['source'].split(',')
                    if not any(s in source_map for s in src):
                        print(f'Skipping {d} because it does not have a reliable reference!')
                        continue

                    sub['value'] = d['value']
                    sub['reference'] = mappedsrc(source_map, src)
                    sub['computed'] = False
                    sub['unit'] = d['u_value']
                    sub['distance_type'] = 'comoving'
                    if ii == 0:
                        sub['default'] = True
                    schema['distance'].append(clean_schema(sub))
                del j['comovingdist']

        # dates
        if 'discoverdate' in j:
            for ii, d in enumerate(j['discoverdate']):
                sub = deepcopy(otter_const.subschema['date_reference'])
                src = d['source'].split(',')
                if not any(s in source_map for s in src):
                    print(f'Skipping {d} because it does not have a reliable reference!')
                    continue

                if '/' in d['value']:
                    indate = d['value'].replace('/', '-')
                    sub['date_format'] = 'iso'
                    day = d['value'].split('/')[-1]
                    if '.' in day:
                        # shorten the day so it doesn't have a float
                        indate = indate.replace(day, day.split('.')[0])

                        # then this has a time of day associated with it
                        intime = float(day)%1 # just get the decimal point 
                        hrs =  int(np.floor(intime * 24))
                        min_leftover = (intime * 24)%1
                        minutes = int(np.floor(min_leftover * 60))
                        sec_leftover = (min_leftover * 60)%1 
                        sec = sec_leftover * 60 # astropy allows this to be a float!



                        indate += f' {hrs}:{minutes}:{sec}'

                elif isinstance(d['value'], str) and len(d['value']) == 4: # this is just a year
                    indate = f'{d["value"]}-01-01'
                    sub['date_format'] = 'iso'
                else:
                    print(j['discoverdate'])
                    raise ValueError('New Type of date found, please fix!')

                sub['value'] = indate
                sub['reference'] = mappedsrc(source_map, src)
                sub['measurement_type'] = 'discovery'
                if 'derived' in d:
                    sub['computed'] = d['derived']
                else:
                    sub['computed'] = False

                sub['computed'] = False
                if ii == 0:
                    sub['default'] = True
                sub['date_type'] = 'discovery'
                schema['date_reference'].append(clean_schema(sub))
            del j['discoverdate']

        # copy classification
        if 'claimedtype' in j:
            #print(json.dumps(j, indent=4))
            for ii, d in enumerate(j['claimedtype']):
                src = d['source'].split(',')
                if not any(s in source_map for s in src):
                    print(f'Skipping {d} because it does not have a reliable reference!')
                    continue

                c = d['value']
                if c == 'TDE':
                    conf = 1.0 # we can trust this
                elif c == 'TDE?':
                    conf = 0.5 # it's probably a TDE
                else:
                    conf = 0.1 # it's might be a TDE


                sub = deepcopy(otter_const.subschema['classification'])
                sub['object_class'] = 'TDE'
                sub['reference'] = mappedsrc(source_map, src)
                sub['confidence'] = conf
                schema['classification'].append(clean_schema(sub))
            del j['claimedtype']

        # copy photometry
        if 'photometry' in j:
            phot = pd.DataFrame(j['photometry'])
            allphot += len(phot)

            if 'telescope' in phot and 'u_time' in phot:
                gby = ['source', 'telescope', 'u_time']
            elif 'u_time' in phot:
                gby = ['source', 'u_time']
            else:
                gby = ['source']

            for grouped_by, group in phot.groupby(gby):

                if len(grouped_by) == 3:
                    ref, telescope, tfmt = grouped_by
                elif len(grouped_by) == 2:
                    ref, tfmt = grouped_by
                    telescope = None
                else:
                    ref = grouped_by
                    telescope = None
                    raise ValueError('Time format uncertain') # for now, hopefully we never come to this anyways...

                # clean up group
                group = group.dropna(axis=1, how='any')
                #print(group)
                sub = deepcopy(otter_const.subschema['photometry'])
                src = ref.split(',')
                if not any(s in source_map for s in src):
                    print(f'Skipping {d} because it does not have a reliable reference!')
                    badphot += len(group)
                    continue

                sub['reference'] = mappedsrc(source_map, src)
                sub['telescope'] = telescope

                usedFluxForRaw = False
                if 'magnitude' in group:
                    sub['raw'] = list(group.magnitude)
                    if 'e_magnitude' in group:
                        sub['raw_err'] = list(group.e_magnitude)
                    if 'system' in group:
                        sub['raw_units'] = list(group.system)
                    else:
                        sub['raw_units'] = 'AB' # hopefully this assumption is okay!
                    if 'band' in group:
                        sub['filter_key'] = list(group.band)
                    else:
                        print(f'Skipping {group} because no filter given so unreliable!')
                        badphot += len(group)
                        continue # don't add this one
                elif 'countrate' in group:
                    sub['raw'] = list(group.countrate)
                    if 'e_countrate' in group:
                        sub['raw_err'] = list(group.e_countrate)
                    if 'u_countrate' in group:
                        sub['raw_units'] = list(group.u_countrate)
                    else:
                        print(f'Skipping {group} because no units given so unreliable!')
                        badphot += len(group)
                        continue # don't add this one
                    if 'energy' in group:
                        sub['filter_key'] = [' - '.join(e) if isinstance(e, list) else e for e in group.energy]
                    else:
                        raise ValueError()
                        print(f'Skipping {group} because no filter given so unreliable!')
                        continue # don't add this one

                elif 'flux' in group and 'u_flux' in group: # we will just use the flux for the raw data
                    usedFluxForRaw = True
                    sub['raw'] = list(group.flux)
                    sub['raw_units'] = list(group.u_flux)
                    if 'e_flux' in group:
                        sub['raw_err'] = list(group.e_flux)
                    if 'band' in group:
                        sub['filter_key'] = list(group.band)
                    elif 'energy' in group:
                        sub['filter_key'] = [' - '.join(e) if isinstance(e, list) else e for e in group.energy]
                    else:
                        print(f'Skipping {group} because no filter given so unreliable!')
                        raise ValueError()
                        continue # don't add this one
                else:
                    # for radio data, throw an error for now
                    #print(group)
                    #raise ValueError('Unknown type of photometry! Please Fix!')
                    print('Skipping this photometry point because it is an unknown type! Please fix!')
                    badphot += len(group)

                if 'flux' in group and 'u_flux' in group and not usedFluxForRaw:
                    sub['value'] = list(group.flux)
                    sub['value_units'] = list(group.u_flux)
                    if 'e_flux' in group:
                        sub['value_err'] = list(group.e_flux)

                if 'upperlimit' in group:
                    sub['upperlimit'] = list(group['upperlimit'])
                else:
                    sub['upperlimit'] = [False]*len(group) # james tended to only store this if it is True

                sub['date'] = [t[0] if isinstance(t, list) else t for t in group['time']]
                sub['date_format'] = tfmt

                # add the observation type

                if sub['filter_key'] is None: continue
                sub['obs_type'] = []
                for filt in sub['filter_key']:
                    if filt in bandwavelengths:
                        sub['obs_type'].append('uvoir')
                    elif filt in xraycodes:
                        sub['obs_type'].append('xray')
                    else:
                        sub['obs_type'].append(otter_helper.filter_to_obstype(filt))

                schema['photometry'].append(clean_schema(sub))

            del j['photometry']

        # create the filter alias
        schema['filter_alias'] = []
        for val in schema['photometry']:
            filter_keys = np.unique(val['filter_key'])
            for key in filter_keys:

                curr_filts = [filt['filter_key'] for filt in schema['filter_alias']]
                if key in curr_filts: continue

                sub = deepcopy(otter_const.subschema['filter_alias'])
                sub['filter_key'] = key
                sub['wave_units'] = 'nm'
                if key in bandwavelengths:
                    sub['wave_eff'] = bandwavelengths[key]
                elif key in otter_const.FILTER_MAP_WAVE:
                    sub['wave_eff'] = otter_const.FILTER_MAP_WAVE[key]
                elif key in xraycodes:
                    sub['wave_eff'] = xraycodes[key]
                else:
                    raise ValueError('Can not add filter {key} because we dont know wave_eff')

                for key, val in list(sub.items()):
                    if val is None:
                        del sub[key]

                schema['filter_alias'].append(sub)

        # ADD SPECTRA ONCE WE HAVE A BETTER IDEA OF FORMATTING

        # ADD HOST INFO ONCE WE HAVE A BETTER IDEA OF FORMATTING

        # remove everything we've added from j to help keep track
        del j['name']
        del j['alias']
        del j['sources']
        del j[key1]
        del j[key2]

        #print(j.keys())
        #print(json.dumps(j, indent=4))

        # some checks of the outputs
        assert len(schema['coordinate']) > 0
        assert len(schema['name/alias']) > 0

        # get rid of dispersion measure cause none of these have it
        if len(schema['distance']) == 0:
            del schema['distance']

        del schema['spectra']

        if len(schema['classification']) == 0:
            # give it a low classification score
            schema["classification"] = [{
                    'object_class': 'TDE',
                    'confidence': 0.1, # we don't trust this classification very much
                    'reference': '2017ApJ...835...64G',
                    'default': True
                }]

        if len(schema['photometry']) == 0:
            del schema['photometry']
            del schema['filter_alias'] # this will be empty too

        if 'date_reference' in schema and len(schema['date_reference']) == 0:
            del schema['date_reference']

        if 'distance' in schema and len(schema['distance']) == 0:
            del schema['distance']

        allschemas.append(schema)

        #json_schema = json.dumps(dict(schema), indent=4) # clean_schema(dict(schema))
        print(schema)
        print()

    # print some useful stats
    print()
    print('#######################################################################')
    print(f'Skipped {badsrcs/allsrcs * 100 : .2f}% of sources because we did not recognize them')
    print(f'Skipped {badphot/allphot * 100 : .2f}% of photometry points because we did not trust them')
    print('#######################################################################')
    print()

    # add the data to the output directory
    p = args.outdir
    otter = Otter(p)
    
    otter.save(allschemas, test_mode=False)

if __name__ == '__main__':
    main()
