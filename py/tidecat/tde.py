'''
Simple TDE class with information about an individual TDE
'''
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.io import to_html

class TDE:

    def __init__(self, name, ra, dec, z, sources, photometry=None, spectra=None):

        self.name = name
        self.ra = ra
        self.dec = dec
        self.z = z
        self.sources = sources
        if self.sources is not None:
            self.fancySources = self._fancySources()
        else:
            self.fancySources = ''

        self.path = self._getpath()

        if photometry is not None:
            self.photometry = photometry

        if spectra is not None:
            self.spectra = spectra

        # create a mapping from the source aliases to the source names
        self._sourcemap = {s['alias']:s['name'] for s in sources}
        self._sourcemap['unknown'] = 'unknown' # THIS IS FOR PHOTOMETRY WE ARENT SURE ABOUT
        
    def _fancySources(self):

        html = ''
        for s in self.sources:
            html += f"<a href='https://ui.adsabs.harvard.edu/abs/{s['bibcode']}' target='_blank'>{s['name']}</a><br>"

        return html

    def _getpath(self):
        return '' # for now

    def plotPhotometry(self,  **kwargs):
        '''
        Plots the photometry for this TDE using plotly

        The TDE MUST have photometry associated with it!
        '''

        keywords = ['telescope', 
                'observatory',
                'instrument']
        
        def create_layout_button(key):
            return dict(label = key,
                        method = 'update',
                        args = [{'visible': key in self.photometry.keys(),
                                 'title': key,
                                 'showlegend': True}])
        
        if self.photometry is None:
            raise Exception("There is no photometry associated with this object!")

        fig = go.Figure()
        
        for key in self.photometry:

            data = {'Time [MJD]': [],
                    'Luminosity [erg/s]': [],
                    'symbols': [],
                    'telescope': [],
                    'observatory': [],
                    'instrument': [],
                    'source': []}
            
            for d in self.photometry[key]:
            
                if 'time' not in d: continue
                
                time = d['time']
                if isinstance(time, list):
                    time = float(d['time'][0])
                else:
                    time = float(d['time'])
                    
                if 'magnitude' in d:
                    lum = float(d['magnitude'])
                    ylabel = 'Apparent Magnitude'
                elif 'luminosity' in d:
                    lum = float(d['luminosity'])
                    ylabel = 'Luminosity [erg/s]'
                    #ax.set_yscale('log')
                else:
                    continue
            
                for key in keywords:
                    if key in d:
                        data[key].append(d[key])
                    else:
                        data[key].append('')
            
                data['Time [MJD]'].append(time)
                data['Luminosity [erg/s]'].append(lum)
            
                if 'source' in d:
                    print(self.sources)
                    data['source'].append(self._sourcemap[d['source']])
                else:
                    data['source'].append('')
                
                if 'upperlimit' in d and d['upperlimit'] == True:
                    data['symbols'].append('triangle-down')
                else:
                    data['symbols'].append('circle')

            df = pd.DataFrame(data)
            fig.add_trace(go.Scatter(x=df['Time [MJD]'],
                                     y=df['Luminosity [erg/s]'],
                                     marker_symbol=df['symbols'],
                                     customdata=np.stack((df['telescope'],
                                                          df['observatory'],
                                                          df['instrument'],
                                                          df['source']
                                                          ), axis=-1),
                                     hovertemplate=
                                         'Telescope: %{customdata[0]}<br>'+
                                         'Observatory: %{customdata[1]}<br>'+
                                         'Instrument: %{customdata[2]} <br>'+
                                         'Sources: %{customdata[3]}'+
                                         '<extra></extra>',
                                     **kwargs))
        #fig.update_traces(mode="markers+lines")

        fig.update_layout(
            updatemenus=[go.layout.Updatemenu(
                x = 1.25,
                y = 1,
                active = 0,
                buttons = [create_layout_button(key) for key in self.photometry],
                pad = dict(l=10)
            )
                         ])
        
        return to_html(fig, full_html=False,
                       default_width='500px',
                       default_height='500px')

        
    def plotSpectra(self):
        '''
        Plots the spectra for this TDE using plotly
        '''
    
    def __str__(self, html=False):

        if not html:
            return f'''
            TDE: {self.name}
            --------------------------------
            RA      : {self.ra}
            DEC     : {self.dec}
            Z       : {self.z}
            Sources : {self.fancySources}
            '''

        else:
            html = ''

            infoKeys = {'ra': 'RA [hrs]',
                        'dec': 'Dec [deg]',
                        'z': 'Redshift',
                        'fancySources': 'Sources'
                        }
            
            for key in infoKeys:

                info = ''
                if isinstance(getattr(self,key), list):
                    for val in getattr(self,key):
                        info += f"{val}<br>\n"
                else:
                    info = getattr(self,key)    

                html += f'''
                <tr>
                <td style="text-align=left">{infoKeys[key]}</td>
                <td style="text-align=right">{info}</td>
                <tr>
                '''
            
            return html
            
        
