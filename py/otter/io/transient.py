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

from ..unit_types import *

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
        '''
        Override getitem to recursively access Transient elements
        '''
        
        if isinstance(keys, (list, tuple)):
            return Transient({key:self[key] for key in keys})
        elif isinstance(keys, str) and '/' in keys: # this is for a path
            s = "']['".join(keys.split('/'))
            s = "['" + s
            s += "']"
            return eval(f"self{s}")
        elif isinstance(keys, int) or keys.isdigit() or (keys[0] == '-' and keys[1:].isdigit()):
            # this is for indexing a sublist
            return self[int(keys)]
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

            coord = self.getSkyCoord()
            
            # add the ra and dec
            # These are required so no need to check if they are there
            html += f'''
            <tr>
            <td style="text-align:left">RA [hrs]:</td>
            <td style="text-align:left">{coord.ra}
            </tr>
            <tr>
            <td style="text-align:left">DEC [deg]:</td>
            <td style="text-align:left">{coord.dec}
            </tr>
            '''

            if 'date_reference' in self:
                discovery = self.getDiscoveryDate().to_value('datetime')
                if discovery is not None:
                    # add the discovery date
                    html += f'''
                    <tr>
                    <td style="text-align:left">Discovery Date [MJD]:</td>
                    <td style="text-align:left">{discovery}
                    </tr>
                    '''

            if 'distance' in self:
                
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
            elif key == 'date_reference':
                self._merge_date(other, out)
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

            # note: using the remove method is safe here because dict keys are unique
            if 'photometry' in keys:
                keys.remove('photometry')
            if 'spectra' in keys:
                keys.remove('spectra')
        else:
            # run some checks
            if 'photometry' in keys:
                warnings.warn('Not returing the photometry!')
                _=keys.pop('photometry')
            if 'spectra' in keys:
                warnings.warn('Not returning the spectra!')
                _=keys.pop('spectra')

            curr_keys = self.keys()
            for key in keys:
                if key not in curr_keys:
                    keys.remove(key)
                    warnings.warn(f'Not returning {key} because it is not in this transient!')
                    
        return self[keys]        

    def getSkyCoord(self, coord_format='icrs'):
        '''
        Convert the coordinates to an astropy SkyCoord
        '''

        # now we can generate the SkyCoord
        f = "df['coordinate_type'] == 'equitorial'"
        coord_dict = self.get_default("coordinate", filt=f)
        coordin = self._reformat_coordinate(coord_dict)
        coord = SkyCoord(**coordin).transform_to(coord_format)

        return coord

    def getDiscoveryDate(self):
        '''
        Get the default discovery date
        '''
        key = 'date_reference'
        date = self.get_default(key, filt='df["date_type"] == "discovery"')
        if 'date_format' in date:
            f = date['date_format']
        else:
            f = 'mjd'

        return Time(date['value'], format=f)

    def getRedshift(self):
        '''
        Get the default redshift
        '''
        f = "df['distance_type']=='redshift'"
        default = self.get_default('distance', filt=f)
        if default is None:
            return default
        else:
            return default['value']
    
    def get_default(self, key, filt=''):
        '''
        Get the default of key

        Args:
            key [str]: key in self to look for the default of
            filt [str]: a valid pandas dataframe filter to index a pandas dataframe called df.
        '''
        if key not in self:
            raise KeyError(f'This transient does not have {key} associated with it!')

        df = pd.DataFrame(self[key])
        df = df[eval(filt)] # apply the filters
        
        if 'default' in df:
            # first try to get the default
            df_filtered = df[df.default == True]
            if len(df_filtered) == 0:
                df_filtered = df
        else:
            df_filtered = df

        if len(df_filtered) == 0:
            return None
        return df_filtered.iloc[0]
            
    def _reformat_coordinate(self, item):
        '''
        Reformat the coordinate information in item
        '''
        coordin = None
        if 'ra' in item and 'dec' in item:
            # this is an equitorial coordinate
            coordin = {'ra': item['ra'],
                       'dec': item['dec'],
                       'unit': (item['ra_units'],
                                item['dec_units'])
                       }
        elif 'l' in item and 'b' in item:
            coordin = {'l': item['l'],
                       'b': item['b'],
                       'unit': (item['l_units'],
                                item['b_units']),
                       'frame': 'galactic'
                       }
            
        return coordin

    
    def cleanPhotometry(self, flux_unit='mag(AB)', date_unit='MJD', by='raw'):
        '''
        Ensure the photometry associated with this transient is all in the same units/system/etc
        '''

        # check inputs
        if by not in {'value', 'raw'}:
            raise ValueError('Please choose either value or raw!')
        
        # turn the photometry key into a pandas dataframe
        dfs = []
        for item in self['photometry']:
            max_len = 0
            for val in item.values():
                if isinstance(val, list):
                    max_len = max(max_len, len(val))
                    
            for key, val in item.items():
                if not isinstance(val, list) or isinstance(val, list) and len(val) < max_len:
                    item[key] = [val]*max_len
                    
            df = pd.DataFrame(item)
            dfs.append(df)
            
        c = pd.concat(dfs)

        filters = pd.DataFrame(self['filter_alias'])
        df = c.merge(filters, on='filter_key')

        # make sure 'by' is in df
        if by not in df:
            if by == 'value':
                by = 'raw'
            else:
                by = 'value'
        
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
        for groupedby, data in df.groupby(['obs_type', by+'_units']):

            obstype, unit = groupedby
            
            # get the photometry in the right type
            unit = data[by+'_units'].unique()
            if len(unit) > 1:
                raise ValueError('Can not apply multiple units for different obs_types')

            unit = unit[0]
            try:
                astropy_units = u.Unit(unit)
            except ValueError:
                # this means there is something likely slightly off in the input unit
                # string. Let's try to fix it!
                # here are some common mistakes
                unit = unit.replace('ergs', 'erg')
                unit = unit.replace('AB', 'mag(AB)')
                
                warnings.warn('Attempting to coerce vega mag to AB mag, this is potentially dangerous!')
                unit = unit.lower().replace('vega', 'mag(AB)')
                
                
                astropy_units = u.Unit(unit)
            except ValueError:
                raise ValueError('Could not coerce your string into astropy unit format!')

            indata = np.array(data[by].astype(float))
            q = indata*u.Unit(astropy_units)
            phot = get_type(q)

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
            elif CountRate.iscountrate(1*u.Unit(flux_unit)):
                flux = phot.tocountrate(**conversion)
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

        Use pandas to drop any duplicates
        '''
        key = 'coordinate'

        Transient._merge_arbitrary(key, t1, t2, out)
        
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
        refs = np.array([d['reference'] for d in out[key]])
        merge_dups = lambda val: np.sum(val) if np.any(val.isna()) else val.iloc[0]
        
        for val in t2[key]:

            # first check if t2's reference is in out
            if val['reference'] not in refs:
                # it's not here so we can just append the new photometry!
                out[key].append(val)
            else:
                # we need to merge it with other photometry
                i1 = np.where(val['reference'] == refs)[0][0]
                df1 = pd.DataFrame(out[key][i1])
                df2 = pd.DataFrame(val)
                
                # only substitute in values that are nan in df1 or new
                mergeon = list(df1.keys() & df2.keys()) # the combined keys of the two
                df = df1.merge(df2, on=mergeon, how='outer')
                
                # convert to a dictionary
                newdict = df.reset_index().to_dict(orient='list')
                del newdict['index']
                
                newdict['reference'] = newdict['reference'][0]
                
                out[key][i1] = newdict # replace the dictionary at i1 with the new dict
                
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

    def _merge_date(t1, t2, out):
        '''
        Combine epoch data across two transients and write it to "out"
        '''
        key = 'date_reference'

        Transient._merge_arbitrary(key, t1, t2, out)
        
    def _merge_distance(t1, t2, out):
        '''
        Combine distance information for these two transients
        '''
        key = 'distance'

        Transient._merge_arbitrary(key, t1, t2, out)

    @staticmethod
    def _merge_arbitrary(key, t1, t2, out):
        '''
        Merge two arbitrary datasets inside the json file using pandas
        
        The datasets in t1 and t2 in "key" must be able to be forced into 
        a NxM pandas dataframe!
        '''

        df1 = pd.DataFrame(t1[key])
        df2 = pd.DataFrame(t2[key])

        merged_with_dups = pd.concat([df1, df2]).reset_index(drop=True)

        # have to get the indexes to drop using a string rep of the df
        # this is cause we have lists in some cells
        to_drop = merged_with_dups.astype(str).drop_duplicates().index

        merged = merged_with_dups.iloc[to_drop].reset_index(drop=True)
        
        outdict = merged.to_dict(orient='records')

        outdict_cleaned = Transient._remove_nans(outdict) # clear out the nans from pandas conversion

        out[key] = outdict_cleaned

    @staticmethod
    def _remove_nans(d):
        '''
        Remove nans from a record dictionary

        THIS IS SLOW: O(n^2)!!! WILL NEED TO BE SPED UP LATER
        '''

        outd = []
        for item in d:
            outsubd = {}
            for key, val in item.items():
                if not isinstance(val, float):
                    # this definitely is not NaN
                    outsubd[key] = val
                
                else:
                    if not np.isnan(val):
                        outsubd[key] = val
            outd.append(outsubd)

        return outd
