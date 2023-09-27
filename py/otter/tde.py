'''
Simple TDE class with information about an individual TDE
'''
import json
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.io import to_html

from .dbattr import *

class TDE:

    def __init__(self, attrs):
        '''        
        `attrs` should be a dictionary with the following keys:
        
        Required Keys:
        - name
        - ra
        - dec
        - sources
        - discovery_date
        
        Optional Keys:
        - redshift
        - spectra
        - photometry

        NOTE: if the TDE is already in the database, only the name is required! Other 
              information will be appended to the existing database object.
        
        All other data will be ignored because it is not stored in the
        database. If photometry is provided it must have the following properties. All
        numerical values should be in CGS!
        - source (corresponding to a source in "sources")
        - flux
        - luminosity
        - time (MJD since discovery)
        - telescope (this can be the telescope, observatory, instrument, or some combination)
        - wavelength (the central wavelength of the band)

        If spectra is provided it must have the following keys
        - source (corresponding to a source in "sources")
        - time (MJD since discovery date)
        - telescope (instrument/observatory used)
        - wavelength (Array of wavelengths in REST FRAME corresponding to flux)
        - flux (array of fluxes)

        Here is a list of acceptable units:
        - Time: MJD
        - Flux: ergs/s/cm^2
        - Luminosity: ergs/s
        '''

        self.requiredInputs = {Name, RA, Dec, Source}
        self.optionalInputs = {Redshift, Spectra, Photometry, DiscoveryDate}

        # define the sources
        self.sources = [Source(s) for s in attrs['sources']]
        
        # create a mapping from the source aliases to the source names
        self._sourcemap = {s.alias:s for s in self.sources}
        self._sourcemap['unknown'] = 'unknown' # THIS IS FOR PHOTOMETRY WE ARENT SURE ABOUT

        # unpack the rest of the inputs
        self._unpackInput(attrs)
        
        if self.sources is not None:
            self.fancySources = self._fancySources()
        else:
            self.fancySources = ''

        self.path = self._getpath()

        # now clean up the photometry and make sure it complies
        
    def _unpackInput(self, attrs):
        for a in self.requiredInputs:
            if a.strname == 'sources': continue
            try:
                if a.strname == 'name':
                    attr = attrs[a.strname]
                    if 'alias' in attrs:
                        attrToSet = a(attr, aliases=attrs['alias'])
                        del attrs['alias']
                    else:
                        attrToSet = a(attr)
                else:
                    # because there can be multiple
                    attrToSet = [a(val, sourcemap=self._sourcemap) for val in attrs[a.strname]] 

                setattr(self, a.strname, attrToSet)

                del attrs[a.strname] # remove from the dictionary
                
            except KeyError:
                # verify input
                print(attrs)
                raise Exception(f'Missing "{a.strname}" from {attrs["name"]}, which is a required input key')

        # now unpack optional inputs
        for a in self.optionalInputs:
            if a.strname in attrs:
                if a.strname == 'photometry' or a.strname == 'spectra':
                    setattr(self, a.strname, [a(val, key, sourcemap=self._sourcemap) for key, val in attrs[a.strname].items()])
                  
                else:
                    setattr(self, a.strname, [a(val, sourcemap=self._sourcemap) for val in attrs[a.strname]])

                del attrs[a.strname] # remove from the dictionary
        
        # save the rest of the data in the "other" attribute
        if len(attrs) > 0:
            self.other = [attrs]
        else:
            self.other = None
                
    def _cleanPhotometry(self):

        # possible instrument types
        instrumentTypes = {'radio', 'optical', 'xray', 'infared', 'uv'}

    def _cleanSpectra(self):
        # I'm not sure how to do these things yet
        return
    
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

        bands = [p.band for p in self.photometry]
        
        if self.photometry is None:
            raise Exception("There is no photometry associated with this object!")

        fig = go.Figure()
        visible = True
        for photo in self.photometry:

            data = {'Time Since Discovery [MJD]': photo.time,
                    'Luminosity [erg/s]': photo.luminosity,
                    'symbols': [],
                    'source': photo.source}
            
            for b in photo.upperlimit:
                if b:
                    data['symbols'].append('triangle-down')
                else:
                    data['symbols'].append('circle')

            df = pd.DataFrame(data)
            fig.add_trace(go.Scatter(x=df['Time Since Discovery [MJD]'],
                                     y=df['Luminosity [erg/s]'],
                                     marker_symbol=df['symbols'],
                                     mode='markers',
                                     customdata=df['source'],
                                     hovertemplate=
                                         'Sources: %{customdata}'+
                                         '<extra></extra>',
                                     visible = visible,
                                     **kwargs))
            visible = False


        buttons = []
        for ii, p in enumerate(self.photometry):
            vis = [False]*len(self.photometry)
            vis[ii] = True
            d = dict(label=p.band,
                     method='update',
                     args=[{"visible": vis},
                           {'title':'',
                            'annotations':[]}]
                     )
            buttons.append(d)
            
        fig.update_layout(
            updatemenus=[go.layout.Updatemenu(
                x = 1.25,
                y = 1,
                active = 0,
                buttons = buttons,
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

    def tojson(self):
        '''
        Writes the attributes of the TDE to a JSON formatted string
        '''
        
        json = {'name':self.name.name,
                'ra': [s.tojson() for s in self.ra],
                'dec': [s.tojson() for s in self.dec],
                'sources':[s.tojson() for s in self.sources]
                }
    
        for opt in self.optionalInputs:
            if opt.strname in dir(self):
                if opt.strname == 'spectra': continue # JUST FOR NOW
                if opt.strname == 'photometry':
                    attr = getattr(self, opt.strname)
                    json[opt.strname] = {a.band:a.json for a in attr} 
                else:
                    json[opt.strname] = [s.tojson() for s in getattr(self, opt.strname)]

        if self.other is not None:
            json['other'] = self.other

        return json
    
    def __str__(self, html=False):

        if not html:
            return f'''
            TDE: {self.name}
            --------------------------------
            RA      : {self.ra.value}
            DEC     : {self.dec.value}
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
                if not hasattr(self, key): continue
                if isinstance(getattr(self,key), list):
                    for val in getattr(self,key):
                        info += f"{val['value']}<br>\n"
                else:
                    info = getattr(self,key)    

                html += f'''
                <tr>
                <td style="text-align=left">{infoKeys[key]}</td>
                <td style="text-align=right">{info}</td>
                <tr>
                '''
            
            return html
            
        
