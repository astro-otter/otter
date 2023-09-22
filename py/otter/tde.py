'''
Simple TDE class with information about an individual TDE
'''
import json
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.io import to_html

class Source:
    def __init__(self, sourceDict):
        '''
        Args:
            source [dict]: Contains at least keys with "name", "bibcode", "alias"
        '''

        self.name = sourceDict['name']
        self.bibcode = sourceDict['bibcode']
        self.alias = sourceDict['alias']
        self.strname = 'source'
        
class RAs:
    def __init__(self, ra):
        '''
        Args:
            ra [dict]: RA dict that can be in any of the following formats and must contain
                       the keys "ra" and "source"
                      - HH MM SS.SSSS
                      - HH:MM:SS.SSSS
        '''
        self.strname = 'ra'
        self.ra = []
        self.sources = []
        for r in ra:
            self.ra.append(r[self.strname])
            if 'source' in r:
                src = Source(r['source'])
                self.sources.append(src)
            else:
                self.sources.append(None)
            
class Decs:
    def __init__(self, dec):
        '''
        Args:
            dec [dict]: Declination dict that can be in any of the following formats and must contain
                       the keys "dec" and "source"
                      - HH MM SS.SSSS
                      - HH:MM:SS.SSSS
        '''
        self.strname = 'dec'
        self.dec = []
        self.sources = []
        for r in dec:
            self.dec.append(r[self.strname])
            if 'source' in r:
                src = Source(r['source'])
                self.sources.append(src)
            else:
                self.sources.append(None)
        
class Redshifts:
    def __init__(self, zs):
        '''
        Args:
            z [list[dict]]: List of redshift dictionaries with keys 'z' and 'source'
        '''
        self.strname = 'z'
        self.z = []
        self.sources = []
        for z in zs:
            self.z.append(float(z[self.strname]))
            if 'source' in z:
                src = Source(z['source'])
                self.source.append(src)
            else:
                self.source.append(None)

class Name:
    def __init__(self, name, aliases):
        '''
        Args:
            name [str]: Primary key of the TDE 
            aliases [list[str]]: other names the TDE goes by
        '''
        self.strname = 'name'
        self.name = name
        self.aliases = aliases

class DiscoveryDates:
    def __init__(self, discoveryDates):
        '''
        Args:
            discoveryDate[list[dict]]: list of dictionaries of discover dates
        '''
        self.strname = 'discover_date'
        self.discoveryDates = []
        self.sources = []
        for dd in discoveryDates:
            self.discoveryDates.append(dd[self.strname])
            if 'source' in dd:
                src = Source(dd['source'])
                self.source.append(src)
            else:
                self.source.append(None)

        
class Spectra:
    def __init__(self, spectra):
        '''
        Args:
            spectra [list[dict]]: list of dictioanries of spectra
        '''
        pass

class Photometry:
    def __init__(self, photometry):
        '''
        Args:
            photometry [list[dict]]: list of photometry dictionaries
        '''
        attrs = ['flux', 'luminosity', 'time', 'telescope', 'wavelength']
        
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

        self.requiredInputs = {Name, RAs, Decs, Source}
        self.optionalInputs = {Z, Spectra, Photometry, DiscoveryDate}
        
        self._unpackInput(attrs)
        
        if self.sources is not None:
            self.fancySources = self._fancySources()
        else:
            self.fancySources = ''

        self.path = self._getpath()

        # create a mapping from the source aliases to the source names
        self._sourcemap = {s['alias']:s['name'] for s in self.sources if 'alias' in s}
        self._sourcemap['unknown'] = 'unknown' # THIS IS FOR PHOTOMETRY WE ARENT SURE ABOUT

        # now clean up the photometry and make sure it complies

    def _unpackInput(self, attrs):

        # verify input
        if not self.requiredInputs.issubset(attrs.keys()):
            raise Exception(f'The required input keys are {requiredInputs}')

        # now unpack
        for a in attrs:
            if a in self.requiredInputs or a in self.optionalInputs:

                # make sure sources is formatted correctly
                if a == 'sources':
                    for val in attrs[a]:
                        # these things are needed for later analysis
                        # even if they are meaningless
                        assert 'bibcode' in val
                        assert 'name' in val
                        assert 'alias' in val
                
                setattr(self, a, attrs[a])
        
    def _unpackInput_new(self, attrs):

        # verify input
        if not self.requiredInputs.issubset(attrs.keys()):
            raise Exception(f'The required input keys are {requiredInputs}')

        # first unpack the required inputs
        self.sources = [Source(s) for s in attrs['sources']]

        for a in self.requiredInputs:
            if a.strname == 'source': continue
            setattr(self, a.strname, a(attrs[a.strname]))
        
        # now unpack optional inputs
        for a in self.optionalInputs:
            if a.strname in attrs:
                setattr(self, a.strname, a(attrs[a.strname]))
        
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

    def tojson(self):
        '''
        Writes the attributes of the TDE to a JSON formatted string
        '''

        json = {'name':self.name,
                'sources':self.sources,
                'ra': self.ra,
                'dec': self.dec
                }
        
        for opt in self.optionalInputs:
            if opt in dir(self):
                json[opt] = getattr(self, opt)

        return json
    
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
            
        
