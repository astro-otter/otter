'''
Class for a transient, 
basically just inherits the dict properties with some overwriting
'''
import warnings
from copy import deepcopy
from collections.abc import MutableMapping

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.io import to_html

import astropy.units as u
from astropy.time import Time

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
            return Transient({key:self.data[key] for key in keys})
        elif isinstance(keys, str) and '/' in keys: # this is for a path
            s = "']['".join(keys.split('/'))
            s = "['" + s
            s += "']"
            return eval(f"self{s}") 
        else:
            return self.data[keys]
            
    def __setitem__(self, key, value):
        self.data[key] = value

    def __delitem__(self, keys):
        del self.data[keys]

    def __iter__(self):
        return iter(self.data)

    def __len__(self):
        return len(self.data)
        
    def __repr__(self, html=False):
        if not html:
            return f'''
            Transient({self.default_name})\n
            '''
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
