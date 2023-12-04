'''
Class for a transient, 
basically just inherits the dict properties with some overwriting
'''
import warnings
from copy import deepcopy
import re
from collections.abc import MutableMapping

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.io import to_html

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord

from .unit_types import *

class Transient(MutableMapping):

    def __init__(self, d={}, name=None):
        '''
        Overwrite the dictionary init
        
        Args:
            d [dict]: A transient dictionary
        '''
        self.data = d

        if 'reference_alias' in self:
            self.srcmap = {ref['name']:ref['human_readable_name'] for ref in self['reference_alias']}
        else:
            self.srcmap = {}
            
        if 'name' in self:
            if 'default_name' in self['name']:
                self.default_name = self['name']['default_name']
            else:
                raise AttributeError('Missing the default name!!')
        elif name is not None:
            self.default_name = name
        else:
            self.default_name = 'Missing Default Name'

        # Make it so all coordinates are astropy skycoords
            
    def __getitem__(self, keys):
        if isinstance(keys, (list, tuple)):
            return Transient({key:self[key] for key in keys})
        elif isinstance(keys, str) and '/' in keys: # this is for a path
            s = "']['".join(keys.split('/'))
            s = "['" + s
            s += "']"
            return eval(f"self{s}") 
        else:
            return self.data[keys]
            
    def __setitem__(self, key, value):
        if isinstance(key, str) and '/' in key: # this is for a path
            s = "']['".join(key.split('/'))
            s = "['" + s
            s += "']"
            exec(f"self{s} = value")
        else:
            self.data[key] = value

    def __delitem__(self, keys):
        del self.data[keys]

    def __iter__(self):
        return iter(self.data)

    def __len__(self):
        return len(self.data)

    def __iter__(self):
        for key, item in self.data.items():
            yield key, item
    
    def __repr__(self, html=False):
        if not html:
            return str(self.data)
        else:
            html = ''
            
            # add the ra and dec
            # These are required so no need to check if they are there
            html += f'''
            <tr>
            <td style="text-align:left">RA [hrs]:</td>
            <td style="text-align:left">{self['coordinate']['equitorial'][0]['ra']}
            </tr>
            <tr>
            <td style="text-align:left">DEC [deg]:</td>
            <td style="text-align:left">{self['coordinate']['equitorial'][0]['dec']}
            </tr>
            '''

            if 'epoch' in self and 'date_discovery' in self['epoch']:
                # add the discovery date
                html += f'''
                <tr>
                <td style="text-align:left">Discovery Date [MJD]:</td>
                <td style="text-align:left">{self['epoch']['date_discovery'][0]['value']}
                </tr>
                '''

            if 'distance' in self and 'redshift' in self['distance']:
                # add the redshift
                html += f'''
                <tr>
                <td style="text-align:left">Redshift:</td>
                <td style="text-align:left">{self['distance']['redshift'][0]['value']}
                </tr>
                '''

            if 'reference_alias' in self:
                srcs = ''
                for bibcode, src in self.srcmap.items():
                    srcs += f"<a href='https://ui.adsabs.harvard.edu/abs/{bibcode}' target='_blank'>{src}</a><br>"

                html += f'''
                <tr>
                <td style="text-align:left">Sources:</td>
                <td style="text-align:left">{srcs}
                </tr>
                '''
                
            return html

    def keys(self):
        return self.data.keys()

    def __add__(self, other, strict_merge=True):
        '''
        Merge this transient object with another transient object

        Args:
            other [Transient]: A Transient object to merge with
            strict_merge [bool]: If True it won't let you merge objects that
                                 intuitively shouldn't be merged (ie. different
                                 transient events).
        '''

        # first check that this object is within a good distance of the other object
        if strict_merge and self.getSkyCoord().separation(other.getSkyCoord()) > 5*u.arcsec: 
            raise ValueError('These two transients are not within 5 arcseconds!' +
                             ' They probably do not belong together! If they do' +
                             ' You can set strict_merge=False to override the check')
            
        # create a blank dictionary since we don't want to overwrite this object
        out = {}

        # find the keys that are 
        merge_keys = list(self.keys() & other.keys()) # in both t1 and t2 so we need to merge these keys
        only_in_t1 = list(self.keys() - other.keys()) # only in t1
        only_in_t2 = list(other.keys() - self.keys()) # only in t2

        # now let's handle the merge keys
        for key in merge_keys:

            # reference_alias is special
            # we ALWAYS should combine these two
            if key == 'reference_alias':
                out[key] = self[key]
                if self[key] != other[key]:
                    # only add t2 values if they aren't already in it
                    bibcodes = {ref['name'] for ref in self[key]}
                    for val in other[key]:
                        if val['name'] not in bibcodes:
                            out[key].append(val)
                continue

            # we can skip this merge process and just add the values from t1 
            # if they are equal. We should still add the new reference though!
            if self[key] == other[key]:
                # set the value
                # we don't need to worry about references because this will
                # only be true if the reference is also equal!
                out[key] = deepcopy(self[key])
                continue

            # There are some special keys that we are expecting
            if key == 'name':
                self._merge_names(other, out)               
            elif key == 'coordinate':
                self._merge_coords(other, out)
            elif key == 'epoch':
                self._merge_epoch(other, out)
            elif key == 'distance':
                self._merge_distance(other, out)
            elif key == 'filter_alias':
                self._merge_filter_alias(other, out)
            elif key == 'schema_version':
                self._merge_schema_version(other, out)
            elif key == 'photometry':
                self._merge_photometry(other, out)
            elif key == 'spectra':
                self._merge_spectra(other, out)
            elif key == 'classification':
                self._merge_class(other, out)
            else:
                # this is an unexpected key!
                if strict_merge:
                    # since this is a strict merge we don't want unexpected data!
                    raise Exception(f'{key} was not expected! Only keeping the old information!')
                else:
                    # Throw a warning and only keep the old stuff
                    warnings.warn(f'{key} was not expected! Only keeping the old information!')
                    out[key] = deepcopy(self[key])

        # and now combining out with the stuff only in t1 and t2
        out = out | dict(self[only_in_t1]) | dict(other[only_in_t2])
        
        # now return out as a Transient Object
        return Transient(out)
        
    def getMeta(self, keys=None):
        '''
        Get the metadata (no photometry or spectra)
            
        This essentially just wraps on __getitem__ but with some checks
        
        Args:
            keys [list[str]] : list of keys
        '''
        if keys is None:
            keys = list(self.keys())
            _=keys.pop('photometry')
            _=keys.pop('spectra')
        else:
            # run some checks
            if 'photometry' in keys:
                warnings.warn('Not returing the photometry!')
                _=keys.pop('photometry')
            if 'spectra' in keys:
                warnings.warn('Not returning the spectra!')
                _=keys.pop('spectra')

        return self[keys]        

    def getSkyCoord(self, coord_type='equitorial', idx=0):
        '''
        Convert the coordinates to an astropy SkyCoord
        '''

        # first run some checks that the data we need is there
        req_keys = {'ra', 'dec', 'ra_units', 'dec_units'}
        key1 = f'coordinate/{coord_type}'
        if key1 not in self:
            raise KeyError(f'This transient does not have {key1} associated with it!')
        if len(self[key1])-1 < idx:
            raise KeyError(f'This transient does not have a coordinate with index {idx}')
        diff = list(req_keys-self[key1][idx].keys())
        if len(diff) > 0:
            raise KeyError(f'The following keys are missing, cant convert! {diff}')

        # now we can generate the SkyCoord
        if coord_type == 'equitorial':
            coordin = {'ra': self[key1][idx]['ra'],
                       'dec': self[key1][idx]['dec'],
                       'unit': (self[key1][idx]['ra_units'],
                              self[key1][idx]['dec_units'])
                        }
        elif coord_type == 'galactic':
            coordin = {'l': self[key1][idx]['l'],
                       'b': self[key1][idx]['b'],
                       'unit': (self[key1][idx]['l_units'],
                                self[key1][idx]['b_units'])
                       }
        else:
            raise ValueError('coord_type must be either equitorial or galactic')
            
        coord = SkyCoord(**coordin)
        return coord
        
    def cleanPhotometry(self, flux_unit='mag(AB)', date_unit='MJD'):
        '''
        Ensure the photometry associated with this transient is all in the same units/system/etc

        THIS IS HELL
        '''

        # turn the photometry key into a pandas dataframe
        dictlist = []
        for phot in self['photometry']:
            meta = deepcopy(self[f'photometry/{phot}'])
            del meta['flux']
            for ddict in self[f'photometry/{phot}/flux']:
                for key in meta:
                    ddict[key] = meta[key]
                ddict['phot_num'] = phot
                dictlist.append(ddict)

        df = pd.DataFrame(dictlist)

        # combine with the filter_alias
        filters = pd.DataFrame(self['filter_alias'])
        df = df.merge(filters, on='filter_key')

        # convert the ads bibcodes to a string of human readable sources here
        def mappedrefs(row):
            if isinstance(row.reference, list):
                return '<br>'.join([self.srcmap[bibcode] for bibcode in row.reference])
            else:
                return self.srcmap[row.reference] 
        try:
            df['human_readable_refs'] = df.apply(mappedrefs, axis=1)
        except Exception as exc:
            warnings.warn(f'Unable to apply the source mapping because {exc}')
            df['human_readable_refs'] = df.reference
        
        # figure out the units of the photometry
        outdata = []
        for obstype, data in df.groupby('obs_type'):

            # get the photometry in the right type
            unit = data.raw_units.unique()
            if len(unit) > 1:
                raise ValueError('Can not apply multiple units for different obs_types')
            astropy_units = u.Unit(unit[0])
            phot = get_type(np.asarray(data.raw) * astropy_units)

            # convert this to a flux
            if 'freq_eff' in data:
                system = 'freq_eff'
                system_units = 'freq_units'
            elif 'wave_eff' in data:
                system = 'wave_eff'
                system_units = 'wave_units'

            uniteff = data[system_units].unique()
            if len(uniteff) > 1:
                raise ValueError('Applying multiple units to the effective frequencies is currently not supported!')
            astropy_uniteff = u.Unit(uniteff[0])
            conversion = {system:np.asarray(data[system])*astropy_uniteff,
                          'out_units':flux_unit}

            if Flux.isflux(1*u.Unit(flux_unit)):
                flux = phot.toflux(**conversion)
            elif FluxDensity.isfluxdensity(1*u.Unit(flux_unit)):
                flux = phot.tofluxdensity(**conversion)
            else:
                raise ValueError('The y-axis units must be either flux or fluxdensity!')

            # add a new column called converted_flux 
            data['converted_flux'] = flux.value
            outdata.append(data)

        outdata = pd.concat(outdata)
        
        # make sure all the datetimes are in the same format here too!!
        times = [Time(d,format=f).to_value(date_unit.lower()) for d, f in zip(outdata.date, outdata.date_format.str.lower())]
        outdata['converted_date'] = times
        
        return outdata


    def plotPhotometry(self, flux_unit='mag(AB)', date_unit='datetime', **kwargs):
        '''
        Plot the photometry associated with this transient (if any)
        
        Args:
            flux_unit [str]: Valid astropy unit string for the flux (y-axis) units.
                             Default: 'ABmag'
            date_unit [str]: Valid astropy unit string for the date (x-axis) units.
                             Default: 'MJD'
        '''
        if 'photometry' not in self:
            raise AttributeError('This transient does not have any photometry associated with it')
        if 'filter_alias' not in self:
            raise AttributeError('filter_alias is required to plot photometry!')
        if 'reference_alias' not in self:
            warnings.warn('Warning! References will be given as bibcodes, not in human readable format!')

        # clean the photometry to prepare it for plotting
        cleanPhot = self.cleanPhotometry(flux_unit=flux_unit, date_unit=date_unit)

        # add symbol for upperlimit or not
        m = lambda row: 'triangle-down' if row.upperlimit else 'circle'
        cleanPhot['symbols'] = cleanPhot.apply(m, axis=1)

        # assign a numerical value to each unique filter
        filters = cleanPhot.filter_key.unique()
        
        filtercolors = ['#FF0000',  # Red
                        '#FF7F00',  # Orange
                        '#FFFF00',  # Yellow
                        '#00FF00',  # Green
                        '#0000FF',  # Blue
                        '#4B0082',  # Indigo
                        '#9400D3'  # Violet
                        ]*len(filters)
        colormap = dict(zip(filters, filtercolors))
        cm = lambda row: colormap[row.filter_key]
        cleanPhot['color'] = cleanPhot.apply(cm, axis=1)
        
        # now plot using plotly
        fig = go.Figure()
        visible = True
        buttons = []
        ii = 0
        for obstype, phot in cleanPhot.groupby('obs_type'):
            fig.add_trace(go.Scatter(x=cleanPhot.converted_date,
                                     y=cleanPhot.converted_flux,
                                     marker_symbol=cleanPhot.symbols,
                                     marker_color=cleanPhot.color,
                                     mode='markers',
                                     customdata=np.stack((cleanPhot.human_readable_refs,
                                                          cleanPhot.filter_key
                                                          ), axis=-1), 
                                     hovertemplate=
                                     '<b>Sources</b>: %{customdata[0]}<br>'+
                                     '<b>Filter</b>: %{customdata[1]}'+
                                     '<extra></extra>',
                                     visible=visible,
                                     **kwargs
                                     ))
        
            visible = False

            vis = [False] * len(cleanPhot.obs_type.unique())
            vis[ii] = True
            d = dict(label=obstype,
                     method='update',
                     args=[{'visible':vis},
                            {'title':'',
                             'annotations':[]}]
                     )
            buttons.append(d)

            ii += 1


        if 'mag' in flux_unit:
            fig.update_yaxes(autorange="reversed")
            ylabel = 'Magnitude'
        else:
            ylabel = f'Flux [{flux_unit}]'

        xlabel = f'Time [{date_unit}]'
            
        fig.update_layout(
            xaxis_title=xlabel,
            yaxis_title=ylabel,
            updatemenus=[go.layout.Updatemenu(
                x = 1.25,
                y = 1,
                active = 0,
                buttons = buttons,
                pad = dict(l=10)
            )    
            ]
        )

        fig.update_yaxes(exponentformat='none')
        fig.update_xaxes(exponentformat='none')

        return to_html(fig, full_html=False, default_width='500px',
                       default_height='500px')
        

    def _merge_names(t1, t2, out):
        '''
        Private method to merge the name data in t1 and t2 and put it in out
        '''
        key = 'name'
        out[key] = {}
        
        # first deal with the default_name key
        # we are gonna need to use some regex magic to choose a preferred default_name
        if t1[key]['default_name'] == t2[key]['default_name']:
            out[key]['default_name'] = t1[key]['default_name']
        else:
            # we need to decide which default_name is better
            # it should be the one that matches the TNS style
            # let's use regex
            n1 = t1[key]['default_name']
            n2 = t2[key]['default_name']

            # write some discriminating regex expressions
            exp1 = '^[0-9]' # starts with a number, this is preferred because it is TNS style
            exp2 = '.$' # ends with any character, this is also preferred because it is TNS style
            exp3 = '^[0-9]{3}' # checks if first four characters are a number, like a year :), this is pretty strict though
            exp4 = '^AT' # checks if it starts with AT like TNS names
            exps = [exp1, exp2, exp3, exp4]

            # score each default_name based on this
            score1 = 0
            score2 = 0
            for e in exps:
                re1 = re.findall(e, n1)
                re2 = re.findall(e, n2)
                if re1:
                    score1 += 1
                if re2:
                    score2 += 1

            # assign a default_name based on the score
            if score1 > score2: 
                out[key]['default_name'] = t1[key]['default_name']
            elif score2 > score1:
                out[key]['default_name'] = t2[key]['default_name']
            else:
                warnings.warn('Names have the same score! Just using the existing default_name')
                out[key]['default_name'] = t1[key]['default_name']

        # now deal with aliases
        # create a reference mapping for each
        t1map = {}
        for val in t1[key]['alias']:
            ref = val['reference']
            if isinstance(ref, str):
                t1map[val['value']] = [ref]
            else:
                t1map[val['value']] = [ref]

        t2map = {}
        for val in t2[key]['alias']:
            ref = val['reference']
            if isinstance(ref, str):
                t2map[val['value']] = [ref]
            else:
                t2map[val['value']] = [ref]

        # figure out which ones we need to be careful with references in        
        inboth = list(t1map.keys() & t2map.keys()) # in both so we'll have to merge the reference key
        int1 = list(t1map.keys() - t2map.keys()) # only in t1
        int2 = list(t2map.keys() - t1map.keys()) # only in t2

        # add ones that are not in both first, these are easy
        L1 = [{'value':k, 'reference':t1map[k]} for k in int1]
        L2 = [{'value':k, 'reference':t2map[k]} for k in int2]
        Lboth = [{'value':k, 'reference':t1map[k]+t2map[k]} for k in inboth]
        out[key]['alias'] =  L1+L2+Lboth

    def _merge_coords(t1, t2, out):
        '''
        Merge the coordinates subdictionaries for t1 and t2 and put it in out
        '''
        key = 'coordinate'
        out[key] = {}

        # first deal with equitorial and then galactic
        subkeys = ['equitorial', 'galactic']
        cnames = [('ra', 'dec', 'icrs'), ('l', 'b', 'galactic')]
        for subkey, c in zip(subkeys, cnames):

            c1, c2, frame = c
            c1_units, c2_units = f'{c1}_units', f'{c2}_units'

            if subkey in t1[key] and subkey in t2[key]:
                out[key][subkey] = t1[key][subkey]
                curr_coords = np.array([SkyCoord(val[c1], val[c2], unit=(val[c1_units], val[c2_units]), frame=frame) for val in t1[key][subkey]])
                for coord in t2[key][subkey]:
                    coorddict = {c1:coord[c1],
                                 c2:coord[c2],
                                 'unit':(coord[c1_units], coord[c2_units]),
                                 'frame': frame
                                }
                    skycoord = SkyCoord(**coorddict)
                    if skycoord not in curr_coords:
                        out[key][subkey].append(coord)
                    else:
                        idx = np.where(skycoord == curr_coords)[0][0] # we only need the first value
                        ref = out[key][subkey][idx]['reference']
                        if not isinstance(ref, list):
                            out[key][subkey][idx]['reference'] = [ref]

                        if not isinstance(coord['reference'], list):
                            coord['reference'] = [coord['reference']]

                        newdata = list(np.unique(out[key][subkey][idx]['reference']+coord['reference']))
                        out[key][subkey][idx]['reference'] = newdata

            elif subkey in t1[key]:
                out[key][subkey] = t1[key][subkey]

            elif subkey in t2[key]:
                out[key][subkey] = t2[key][subkey]

    def _merge_filter_alias(t1, t2, out):
        '''
        Combine the filter alias lists across the transient objects
        '''

        key = 'filter_alias'

        out[key] = deepcopy(t1[key])
        keys1 = {filt['filter_key'] for filt in t1[key]}
        for filt in t2[key]:
            if filt['filter_key'] not in keys1:
                out[key].append(filt)

    def _merge_schema_version(t1, t2, out):
        '''
        Just keep whichever schema version is greater
        '''
        key = 'schema_version/value'
        if int(t1[key]) > int(t2[key]):
            out['schema_version'] = deepcopy(t1['schema_version'])
        else:
            out['schema_version'] = deepcopy(t2['schema_version'])

    def _merge_photometry(t1, t2, out):
        '''
        Combine photometry sources
        '''

        key = 'photometry'

        out[key] = deepcopy(t1[key])

        idx = int(list(out[key].keys())[-1][-1])+1
        telescopes = np.array([phot['telescope'] for phot in out[key].values() if 'telescope' in phot])
        refs = np.concatenate([phot['reference'] for phot in out[key].values() if 'reference' in phot], axis=None)
        for phot in t2[key].values():

            if len(telescopes) > 0 and 'telescope' in phot and phot['telescope'] in telescopes:
                i = np.where(phot['telescope'] == telescopes)[0][0]
                toappend = out[key][f'phot_{i}']
            elif len(refs) > 0 and 'reference' in phot and phot['reference'] in refs:
                i = np.where(phot['reference'] == refs)[0][0]
                toappend = out[key][f'phot_{i}']
            else:
                # nothing with this telescope has been added
                out[key][f'phot_{idx}'] = phot
                idx += 1
                continue

            # if the code has gotten here we need to append to an existing list of photometry
            for point in phot['flux']:
                if point not in toappend['flux']:
                    toappend['flux'].append(point)
                else:
                    if not isinstance(toappend['reference'], list):
                        toappend['reference'] = [toappend['reference']]

                    if not isinstance(phot['reference'], list):
                        phot['reference'] = [phot['reference']]

                    if phot['reference'] not in toappend['reference']:    
                        newdata = list(np.unique(toappend['reference']+phot['reference']))
                        toappend['reference'] = newdata

    def _merge_spectra(t1, t2, out):
        '''
        Combine spectra sources
        '''
        pass

    def _merge_class(t1, t2, out):
        '''
        Combine the classification attribute
        '''
        key = 'classification'
        out[key] = deepcopy(t1[key])
        classes = np.array([item['object_class'] for item in out[key]])
        for item in t2[key]:
            if item['object_class'] in classes:
                i = np.where(item['object_class'] == classes)[0][0]
                if int(item['confidence']) > int(out[key][i]['confidence']):
                    out[key][i]['confidence'] = item['confidence'] # we are now more confident

                if not isinstance(out[key][i]['reference'], list):
                    out[key][i]['reference'] = [out[key][i]['reference']]

                if not isinstance(item['reference'], list):
                    item['reference'] = [item['reference']]

                newdata = list(np.unique(out[key][i]['reference']+item['reference']))
                out[key][i]['reference'] = newdata

            else:
                out[key].append(item)

        # now that we have all of them we need to figure out which one is the default
        maxconf = max(out[key], key=lambda d: d['confidence'])  
        for item in out[key]:
            if item == maxconf:
                item['default'] = True
            else:
                item['default'] = False

    def _merge_epoch(t1, t2, out):
        '''
        Combine epoch data across two transients and write it to "out"
        '''
        key = 'epoch'
        subkeys = ['date_explosion', 'date_peak', 'date_discovery']

        out[key] = {}

        for subkey in subkeys:
            if subkey in t1[key] and subkey in t2[key]:
                out[key][subkey] = t1[key][subkey]
                values = np.array([val['value'] for val in out[key][subkey]])
                for item in t2[key][subkey]:
                    if item['value'] in values:
                        i = np.where(item['value'] == values)[0][0]
                        if not isinstance(out[key][subkey][i]['reference'], list):
                            out[key][subkey][i]['reference'] = [out[key][subkey][i]['reference']]
                        if not isinstance(item['reference'], list):
                            item['reference'] = [item['reference']]

                        out[key][subkey][i]['reference'] = list(np.unique(out[key][subkey][i]['reference']+item['reference']))
                    else:
                        out[key][subkey].append(item)

            elif subkey in t1[key]:
                out[key][subkey] = t1[key][subkey]

            elif subkey in t2[key]:
                out[key][subkey] = t2[key][subkey]

    def _merge_distance(t1, t2, out):
        '''
        Combine distance information for these two transients
        '''
        key = 'distance'
        subkeys = ['redshift', 'luminosity_distance', 'dispersion_measure']
        out[key] = {}
        for subkey in subkeys:
            if subkey in t1[key] and subkey in t2[key]:
                out[key][subkey] = t1[key][subkey]
                values = np.array([val['value'] for val in out[key][subkey]])
                for item in t2[key][subkey]:
                    if item['value'] in values:
                        i = np.where(item['value'] == values)[0][0]
                        if not isinstance(out[key][subkey][i]['reference'], list):
                            out[key][subkey][i]['reference'] = [out[key][subkey][i]['reference']]
                        if not isinstance(item['reference'], list):
                            item['reference'] = [item['reference']]

                        out[key][subkey][i]['reference'] = list(np.unique(out[key][subkey][i]['reference']+item['reference']))
                    else:
                        out[key][subkey].append(item)

            elif subkey in t1[key]:
                out[key][subkey] = t1[key][subkey]

            elif subkey in t2[key]:
                out[key][subkey] = t2[key][subkey]
