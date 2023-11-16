'''
Class for a transient, 
basically just inherits the dict properties with some overwriting
'''
import warnings
from collections.abc import MutableMapping

import pandas as pd
import plotly.graph_objects as go
from plotly.io import to_html

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
            
    def __getitem__(self, keys):
        if isinstance(keys, list):
            return Transient({key:self.data[key] for key in keys})
        elif isinstance(keys, tuple):
            s = "']['".join(keys)
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

            if 'date_discovery' in self['epoch']:
                # add the discovery date
                html += f'''
                <tr>
                <td style="text-align:left">Discovery Date [MJD]:</td>
                <td style="text-align:left">{self['epoch']['date_discovery'][0]['value']}
                </tr>
                '''

            if 'redshift' in self['distance']:
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

    def plotPhotometry(self, **kwargs):
        '''
        Plot the photometry associated with this transient (if any)
        '''
        if 'photometry' not in self:
            raise AttributeError('This transient does not have any photometry associated with it')
        if 'filter_alias' not in self:
            raise AttributeError('filter_alias is required to plot photometry!')
        if 'reference_alias' not in self:
            warnings.warn('Warning! References will be given as bibcodes, not in human readable format!')

        # clean the photometry to prepare it for plotting
        cleanPhot = self.cleanPhotometry()

        print(cleanPhot)
        # add symbol for upperlimit or not
        m = lambda row: 'triangle-down' if row.upperlimit else 'circle'
        cleanPhot['symbols'] = cleanPhot.apply(m, axis=1)
        
        # now plot using plotly
        fig = go.Figure()
        visible = True
        buttons = []
        ii = 0
        for obstype, phot in cleanPhot.groupby('obs_type'):
            fig.add_trace(go.Scatter(x=cleanPhot.date,
                                     y=cleanPhot.raw,
                                     marker_symbol=cleanPhot.symbols,
                                     mode='markers',
                                     customdata=cleanPhot.reference, 
                                     hovertemplate=
                                     'Sources: %customdata'+
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

        fig.update_layout(
            updatemenus=[go.layout.Updatemenu(
                x = 1.25,
                y = 1,
                active = 0,
                buttons = buttons,
                pad = dict(l=10)
            )    
            ]
        )

        return to_html(fig, full_html=False, default_width='500px',
                       default_height='500px')
        
            
    def cleanPhotometry(self):
        '''
        Ensure the photometry associated with this transient is all in the same units/system/etc

        THIS IS HELL
        '''

        # turn the photometry key into a pandas dataframe
        dictlist = []
        for phot in self['photometry']:
            ref = self['photometry', phot, 'reference']
            for ddict in self['photometry', phot, 'flux']:
                ddict['reference'] = ref
                ddict['phot_num'] = phot
                dictlist.append(ddict)

        df = pd.DataFrame(dictlist)

        # combine with the filter_alias
        filters = pd.DataFrame(self['filter_alias'])
        df = df.merge(filters, on='filter_key')

        # SKIP THE NEXT STUFF FOR NOW
        
        # convert the ads bibcodes to a string of human readable sources here
        
        # figure out the units

        # make sure all the datetimes are in the same format here too!!
        
        return df
        
