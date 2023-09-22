'''
Functions to generate different plotly express plots with summary values for the catalog
'''

# imports
import numpy as np
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from plotly.io import to_html
import astropy.coordinates as coord
import astropy.units as u

def skymap(fig, tdes):
    '''
    Generates a skymap of the locations of the tdes
    '''

    fig.update_geos(projection_type="mollweide", showland=False,
                    showcoastlines=False, showframe=True,
                    lonaxis=dict(showgrid=True, dtick=30),
                    lataxis=dict(showgrid=True, dtick=30),
                    row=1, col=1
                    )
    
    info = {'RA [deg]':[],
            'Dec [deg]':[],
            'TDE Name':[],
            'Z': []
            }

    for t in tdes.values():
        if t.ra is not None and t.dec is not None:
            info['RA [deg]'].append(t.ra[0]['value'])
            info['Dec [deg]'].append(t.dec[0]['value'])
            info['TDE Name'].append(t.name)
            info['Z'].append(float(t.z[0]['value']))
            
    info['RA [deg]'] = coord.Angle(info['RA [deg]'], unit=u.hourangle)
    info['Dec [deg]'] = coord.Angle(info['Dec [deg]'], unit=u.deg)
    
    info['RA [deg]'] = info['RA [deg]'].deg
    info['Dec [deg]'] = info['Dec [deg]'].deg
    
    for lon in np.arange(0, 360, 30):
        fig.add_trace(go.Scattergeo(lon=[lon, lon, lon],
                                    lat=[90, 0, -90],
                                    mode='lines+text',
                                    showlegend=False,
                                    text=["", f"{lon}", ""],
                                    textposition='bottom center',
                                    line=go.scattergeo.Line(width=1, dash='dot', color='grey')))
    for lat in [-90, -60, -30, 30, 60, 90]:
        fig.add_trace(go.Scattergeo(lon=[0, 180, 360], lat=[lat, lat, lat],
                                    mode='lines+text',
                                    text=[f"{lat}"],
                                    textposition="middle right",
                                    showlegend=False,
                                    line=go.scattergeo.Line(width=1, dash='dot', color='grey')))
        
    fig.add_trace(go.Scattergeo(lon=info['RA [deg]'],
                                lat=info['Dec [deg]'],
                                hovertext=info['TDE Name'],
                                hovertemplate="TDE Name: %{hovertext}\nRA: %{lon}\nDec: %{lat}", 
                                name='Coordinates',
                                marker=dict(color=info['Z'],
                                            colorbar=dict(thickness=20),
                                            colorscale="Thermal"
                                            )
                                ),
                  row=1, col=1)
        
    return fig
    
def redshifts(fig, tdes):
    '''
    Generates a plot of the redshifts
    '''

    z = np.array([t.z[0]['value'] for t in tdes.values() if t.z is not None]).astype(float)
    
    fig.add_trace(go.Histogram(x=z, nbinsx=5, name='Redshifts',
                               marker={"color": "grey"}),
                  row=2, col=1)
    return fig    

def plotAll(tdes):
    '''
    makes a subplot with all of them and returns the html
    '''
    fig = make_subplots(rows=2, cols=1, specs=[[{"type":"scattergeo"}],
                                              [{"type":"histogram"}]
                                              ]
                        )
    fig.update_xaxes(title_text="Right Ascension [deg]", row=1, col=1)
    fig.update_xaxes(title_text="Redshift", row=2, col=1)
    fig.update_yaxes(title_text="Declination [deg]", row=1, col=1)
    fig.update_yaxes(title_text="Number of TDEs", row=2, col=1)
    
    #fig.update_layout(title = {'text':f'Number of TDEs: {len(tdes)}'})
    
    fig = skymap(fig, tdes)
    fig = redshifts(fig, tdes)

    return to_html(fig, full_html=False,
                       default_width='750px',
                       default_height='1000px')
